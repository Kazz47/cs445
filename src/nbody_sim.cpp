#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "packet.hpp"
#include "node.hpp"
#include "event.hpp"

#include <glog/logging.h>

// Stats vectors
std::vector<double> simulation_times;

// Output File
std::ofstream outfile;

// Print packet
std::ostream& operator<< ( std::ostream& out, Packet& packet) {
    out << "[x:" << packet.getX() << " y:" << packet.getY() << " z:" << packet.getZ();
    out << std::right << "]";
    return out;
}

// Print node
std::ostream& operator<< ( std::ostream& out, Node& node) {
    out << "[x:" << node.getX() << " y:" << node.getY() << " z:" << node.getZ();
    out << std::right << "]";
    return out;
}

class CompareEvent {
    public:
        bool operator()(Event* &e1, Event* &e2) {
            return e1->time > e2->time;
        }
};

// Print event
std::ostream& operator<< ( std::ostream& out, Event& event) {
    out << "[sim_time: " << std::setw(10) << std::setprecision(3) << std::fixed << event.time << ", type: " << std::setw(4) << static_cast<int>(event.type);
    if (event.type < EventType::SIZE) {
        out << " - " << std::left << std::setw(50) << static_cast<int>(event.type);
    } else {
        out << std::left << std::setw(50) << " - UNKNOWN";
    }
    out << std::right << "]";
    return out;
}

double runSimulation(size_t nodes, size_t data, double compute_time, double latency, Topology topology, size_t iterations, size_t pass_num) {
    size_t nodes_x = 0;
    size_t nodes_y = 0;
    size_t nodes_z = 0;
    size_t total_nodes = nodes;

    if (topology == Topology::RING) {
        nodes_x = nodes;
        nodes_y = 1;
        nodes_z = 1;
        total_nodes = nodes;
    } else if (topology == Topology::GRID) {
        nodes_x = nodes;
        nodes_y = nodes;
        nodes_z = 1;
        total_nodes = nodes * nodes;
    } else if (topology == Topology::CUBE) {
        nodes_x = nodes;
        nodes_y = nodes;
        nodes_z = nodes;
        total_nodes = nodes * nodes * nodes;
    } else {
        LOG(FATAL) << "Unknown topology: '" << static_cast<size_t>(topology) << "'";
    }

    double complete_nodes = 0;
    double simulation_time_s = 0;
    double previous_time_s = 0;

    // Init RNG
    std::random_device rd;
    std::mt19937 generator(rd());
    std::exponential_distribution<double> compute_distribution(1/compute_time);
    std::exponential_distribution<double> send_distribution(1/latency);

    std::priority_queue<Event*, std::vector<Event*>, CompareEvent> heap;

    size_t packets_per_node = ceil(data/static_cast<float>(total_nodes));

    std::vector<std::vector<std::vector<Node*>*>*> node_mat;

    // Put initial events in the heap
    for (size_t x = 0; x < nodes_x; x++) {
        std::vector<std::vector<Node*>*> *x_vec = new std::vector<std::vector<Node*>*>();
        node_mat.push_back(x_vec);
        for (size_t y = 0; y < nodes_y; y++) {
            std::vector<Node*> *y_vec = new std::vector<Node*>();
            x_vec->push_back(y_vec);
            for (size_t z = 0; z < nodes_z; z++) {
                Node *node = new Node(x, y, z, 0, packets_per_node, packets_per_node*(total_nodes - 1));
                y_vec->push_back(node);
                double compute_duration = compute_distribution(generator) * data;
                node->updateComputeTime(simulation_time_s, compute_duration);
                LOG_IF(FATAL, node == nullptr);
                heap.push(new Event(simulation_time_s + compute_duration, EventType::COMPUTE, node));
            }
        }
    }

    VLOG(2) << "Start Simulation...";

    while (complete_nodes != total_nodes) {
        LOG_IF(FATAL, heap.empty()) << "Simulation heap was empty.";
        Event *current_event = heap.top();
        heap.pop();

        if (current_event == nullptr) {
            LOG(ERROR) << "Simulation not complete and there are no events in the min heap.";
            LOG(ERROR) << "\tsimulation time: " << simulation_time_s;
            exit(0);
        }

        previous_time_s = simulation_time_s;
        simulation_time_s = current_event->time;

        // Calculate average time stuff
        double time_since_last = simulation_time_s - previous_time_s;

        Node *node = current_event->node;
        Stream* stream = current_event->stream;
        Packet* packet = nullptr;
        Node *next_node = nullptr;

        switch (current_event->type) {
            case EventType::COMPUTE:
                VLOG(2) << "Compute Begin: Node" << *node;
                VLOG(2) << "Iteration: " << node->getIteration();
                for (size_t x = 0; x < node_mat.size(); x++) {
                    std::vector<std::vector<Node*>*> *y_vec = node_mat.at(x);
                    for (size_t y = 0; y < y_vec->size(); y++) {
                        std::vector<Node*> *z_vec = y_vec->at(y);
                        for (size_t z = 0; z < z_vec->size(); z++) {
                            if (*node != *(z_vec->at(z))) {
                                Node *dest = z_vec->at(z);
                                LOG_IF(FATAL, dest == nullptr);
                                for (size_t i = 0; i < node->getNumPacketsToSend(); i++) {
                                    Packet* new_packet = new Packet(dest->getX(), dest->getY(), dest->getZ(), current_event->iteration);
                                    Stream* stream = node->enqueuePacket(new_packet, topology, nodes);
                                    if (stream->needsInit()) {
                                        stream->transferPacket();
                                        heap.push(new Event(simulation_time_s + send_distribution(generator), EventType::PASS_PACKET, node, stream));
                                    }
                                }
                            }
                        }
                    }
                }
                VLOG(2) << "Compute End";
                break;
            case EventType::PASS_PACKET:
                // This needs work. Should continue creating PASS_PACKET events
                // unless the current stream is empty.
                VLOG(4) << "Pass Packet Begin";
                LOG_IF(FATAL, stream == nullptr) << "Stream was null.";
                packet = stream->servicePacket();
                LOG_IF(FATAL, packet == nullptr) << "Stream was empty: " << *node << static_cast<size_t>(stream->getDirection());

                next_node = node->getNextNode(stream, topology, nodes, node_mat);
                LOG_IF(FATAL, next_node == nullptr);
                if (packet->atDest(next_node->getX(), next_node->getY(), next_node->getZ())) {
                    next_node->receivePacket(packet);
                    delete packet;
                    packet = nullptr;
                    if (next_node->readyToCompute()) {
                        double compute_duration = compute_distribution(generator) * data;
                        next_node->updateComputeTime(simulation_time_s, compute_duration);
                        LOG_IF(FATAL, next_node == nullptr);

                        next_node->startNextIteration();
                        next_node->writeMeanCompute(pass_num);
                        if (next_node->getIteration() == iterations) {
                            complete_nodes++;
                        }

                        heap.push(new Event(simulation_time_s + compute_duration, EventType::COMPUTE, next_node));
                        VLOG(1) << "Iteration: " << next_node->getIteration();
                    }
                } else { // Packet needs to be passed on to next node
                    Stream* next_stream = next_node->enqueuePacket(packet, topology, nodes);
                    if (next_stream->needsInit()) {
                        next_stream->transferPacket();
                        heap.push(new Event(simulation_time_s + send_distribution(generator), EventType::PASS_PACKET, next_node, next_stream));
                    }
                }

                if (stream->needsInit()) {
                    stream->transferPacket();
                    heap.push(new Event(simulation_time_s + send_distribution(generator), EventType::PASS_PACKET, node, stream));
                } else {
                    VLOG(3) << "Stream is empty: " << *node << static_cast<size_t>(stream->getDirection());
                }

                VLOG(4) << "Pass Packet End";
                break;
            default:
                LOG(ERROR) << "Simulation had an event with an unknown type: " << static_cast<size_t>(current_event->type);
                LOG(ERROR) << "\tsimulation time: " << simulation_time_s;
                exit(0);
        }

        delete current_event; // Event's are created with new, so we need to delete them when we're done with them
        current_event = nullptr;
    }

    LOG(INFO) << "Ending simulation...";

    LOG_IF(ERROR, !heap.empty());
    while(!node_mat.empty()) {
        std::vector<std::vector<Node*>*> *y_vec = node_mat.back();
        while(!y_vec->empty()) {
            std::vector<Node*> *z_vec = y_vec->back();
            while(!z_vec->empty()) {
                Node *node = z_vec->back();
                z_vec->pop_back();
                LOG(INFO) << "Compute/Transfer Ratio: " << node->getComputeRatio() * 100 << "%";
                LOG(INFO) << "Wait Time: : " << node->getWaitTime(); 
                LOG(INFO) << "Compute Time: : " << node->getComputeTime();
                delete node;
            }
            delete y_vec->back();
            y_vec->pop_back();
        }
        delete node_mat.back();
        node_mat.pop_back();
    }

    VLOG(1) << "The simulation ended at time: " << simulation_time_s;
    return simulation_time_s;
}

int main(int argc, char **argv) {
    // Initialize Google Logging
    google::InitGoogleLogging(argv[0]);
    // Log to Stderr
    FLAGS_logtostderr = 1;

    // Intiate
    size_t iterations = 500;
    size_t nodes = 5;
    size_t data = 200;
    double compute_time= 0.005;
    double latency = 0.00001;
    Topology topology = Topology::RING;
    //Topology topology = Topology::GRID;
    //Topology topology = Topology::CUBE;

    if (argc >= 2) {
        nodes = atoi(argv[1]);
    }
    if (argc >= 3) {
        data = atoi(argv[2]);
    }
    if (argc >= 4) {
        compute_time = atof(argv[3]);
    }
    if (argc >= 5) {
        latency = atof(argv[4]);
    }

    LOG(INFO) << "Iterations: " << iterations;
    LOG(INFO) << "Nodes: " << nodes;
    LOG(INFO) << "Data: " << data;
    LOG(INFO) << "Compute Time: " << compute_time;
    LOG(INFO) << "Latency: " << latency;

    // Open files
    std::stringstream filename;
    filename << "nbody_stats" << ".dat";
    outfile.open(filename.str());

    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 10; j++) {
            double simulation_time = runSimulation(nodes, data, compute_time, latency, topology, iterations, j);
            simulation_times.push_back(simulation_time);
        }

        // Collect and print statistics.
        double sum_time = std::accumulate(simulation_times.begin(), simulation_times.end(), 0.0);
        double min_time = *std::min_element(simulation_times.begin(), simulation_times.end());
        double max_time = *std::max_element(simulation_times.begin(), simulation_times.end());
        double avg_time = sum_time / simulation_times.size();
        double sqr_time = std::inner_product(simulation_times.begin(), simulation_times.end(), simulation_times.begin(), 0.0);
        double std_time = std::sqrt(sqr_time/ simulation_times.size() - avg_time* avg_time);

        outfile << "------------------------------------------" << std::endl;
        outfile << "-------------------- " << i << " -------------------" << std::endl;
        outfile << "------------------------------------------" << std::endl;
        outfile << std::endl;

        outfile << "------------------------------------------" << std::endl;
        outfile << "---------------- RUN TIMES ---------------" << std::endl;
        outfile << "------------------------------------------" << std::endl;
        outfile << "Min: " << std::fixed <<  min_time << std::endl;
        outfile << "Max: " << std::fixed << max_time << std::endl;
        outfile << "Mean: " << std::fixed << avg_time << std::endl;
        outfile << "Stdev: " << std::fixed << std_time << std::endl;
        outfile << std::endl;

    }
    outfile.close();
}
