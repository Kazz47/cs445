#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "node.hpp"
#include "event.hpp"

#include <glog/logging.h>

enum class Topology {RING, GRID, CUBE, SIZE};

// Stats vectors
std::vector<double> simulation_times;

// Output File
std::ofstream outfile;

// Print packet
std::ostream& operator<< ( std::ostream& out, Packet& packet) {
    out << "[x:" << packet.x << " y:" << packet.y << " z:" << packet.z;
    out << std::right << "]";
    return out;
}

// Print node
std::ostream& operator<< ( std::ostream& out, Node& node) {
    out << "[x:" << node.x << " y:" << node.y << " z:" << node.z;
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

double runSimulation(size_t nodes, size_t data, double compute_time, double latency, Topology topology, size_t iterations) {
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
                Node *node = new Node(x, y, z, 0, packets_per_node);
                y_vec->push_back(node);
                double compute_duration = (compute_distribution(generator) * data);
                node->updateComputeTime(simulation_time_s, compute_duration);
                heap.push(new Event(simulation_time_s + compute_duration, EventType::COMPUTE, node));
            }
        }
    }

    VLOG(2) << "Start Simulation...";

    while (complete_nodes != total_nodes) {
        Event *current_event = heap.top();
        heap.pop();

        if (current_event == NULL) {
            LOG(ERROR) << "Simulation not complete and there are no events in the min heap.";
            LOG(ERROR) << "\tsimulation time: " << simulation_time_s;
            exit(0);
        }

        previous_time_s = simulation_time_s;
        simulation_time_s = current_event->time;

        // Calculate average time stuff
        double time_since_last = simulation_time_s - previous_time_s;

        Node *node = current_event->node;
        std::queue<Packet*> *stream= current_event->stream;
        Packet *packet = nullptr;
        Node *next_node = nullptr;

        switch (current_event->type) {
            case 0: // COMPUTE
                VLOG(3) << "Compute Begin";
                LOG_IF(ERROR, stream != nullptr) << "Stream was not freed!";
                for (size_t x = 0; x < node_mat.size(); x++) {
                    std::vector<std::vector<Node*>*> *y_vec = node_mat.at(x);
                    for (size_t y = 0; y < y_vec->size(); y++) {
                        std::vector<Node*> *z_vec = y_vec->at(y);
                        for (size_t z = 0; z < z_vec->size(); z++) {
                            if (*node != *(z_vec->at(z))) {
                                Node *dest = z_vec->at(z);
                                for (size_t i = 0; i < node->getNumPacketsToSend(); i++) {
                                    Packet *new_packet = new Packet(dest->x, dest->y, dest->z, node->iteration);
                                    bool need_init = false;
                                    std::queue<Packet*> *queue = enqueuePacket(node, new_packet, topology, nodes, need_init);
                                    if (need_init) {
                                        heap.push(new Event(simulation_time_s + send_distribution(generator), EventType::PASS_PACKET, node, queue));
                                    }
                                }
                            }
                        }
                    }
                }
                break;
            case 1: // PASS_PACKET
                // This needs work. Should continue creating PASS_PACKET events
                // unless the current stream is empty.
                VLOG(4) << "Pass Packet Begin";
                LOG_IF(FATAL, stream->empty()) << "Stream is empty!";
                packet = stream->front();
                stream->pop();

                next_node = getNextNode(node, stream, topology, nodes, node_mat);
                if (packet->x == next_node->x && packet->y == next_node->y && packet->z == next_node->z) { // Packet is at its destination.
                    VLOG(2) << "Packet is at its destination";
                    if (packet->iteration == node->iteration) { // Received a packet for this iteration
                        node->current_iteration_packets++;
                        VLOG(2) << "Received packet for this iteration: " << node->getNumPacketsToSend()*(total_nodes-1) << " : " << node->current_iteration_packets;
                        if (node->getNumPacketsToSend()*(total_nodes-1) == node->current_iteration_packets) {
                            VLOG(1) << "Node ready to compute: " << node->getNumPacketsToSend()*(total_nodes-1) << " : " << node->current_iteration_packets;
                            double compute_duration = (compute_distribution(generator) * data);
                            void updateComputeTime(double sim_time, double compute_duration) {
                            heap.push(new Event(simulation_time_s + compute_duration, COMPUTE, node));

                            node->startNextIteration();
                            if (node->iteration == iterations) {
                                complete_nodes++;
                            }
                            VLOG(1) << "Iteration: " << node->iteration;
                        }
                    } else if (packet->iteration == node->iteration + 1) { // Received a packet for the next interation
                        VLOG(2) << "Received packet for the next iteration";
                        node->next_iteration_packets++;
                    } else {
                        LOG(ERROR) << "Incorrect packet iteration for node (Packet: " << packet->iteration << ", Node: " << node->iteration << ")";
                    }
                    delete packet;
                    packet = nullptr;
                } else { // Packet needs to be passed on to next node
                    bool need_init = false;
                    std::queue<Packet*> *queue = enqueuePacket(next_node, packet, topology, nodes, need_init);
                    LOG_IF(FATAL, queue == stream);
                    if (need_init) {
                        heap.push(new Event(simulation_time_s + send_distribution(generator), PASS_PACKET, next_node, queue));
                    }
                }

                // Check the current stream and process the next packet.
                if (!stream->empty()) {
                    heap.push(new Event(simulation_time_s + send_distribution(generator), PASS_PACKET, node, stream));
                }

                VLOG(4) << "Pass Packet End";
                break;
            default:
                LOG(ERROR) << "Simulation had an event with an unknown type: " << current_event->type;
                LOG(ERROR) << "\tsimulation time: " << simulation_time_s;
                exit(0);
        }

        delete current_event; // Event's are created with new, so we need to delete them when we're done with them
    }

    LOG_IF(ERROR, !heap.empty());
    while(!node_mat.empty()) {
        std::vector<std::vector<Node*>*> *y_vec = node_mat.back();
        while(!y_vec->empty()) {
            std::vector<Node*> *z_vec = y_vec->back();
            while(!z_vec->empty()) {
                Node *node = z_vec->back();
                z_vec->pop_back();
                LOG(INFO) << "Compute/Transfer Ratio: " << node->total_wait_time/(node->total_wait_time + node->total_compute_time) * 100 << "%";
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
    size_t iterations = 100;
    size_t nodes = 5;
    size_t data = 200;
    double compute_time= 0.005;
    double latency = 0.00001;
    //unsigned int topology = RING;
    unsigned int topology = GRID;
    //unsigned int topology = CUBE;

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
            double simulation_time = runSimulation(nodes, data, compute_time, latency, topology, iterations);
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
