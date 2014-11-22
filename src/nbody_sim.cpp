#include <iostream>
#include <fstream>
#include <queue>
#include <string>
#include <limits>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include <glog/logging.h>

#define ONLY_EVENT_TYPE     0
#define NUMBER_EVENT_TYPES  2
#define NUMBER_TOPOLOGY_TYPES  3

const static unsigned int COMPUTE = 0;
const static unsigned int PASS_PACKET = 1;
const std::string event_names[NUMBER_EVENT_TYPES] = {
    "COMPUTE",
    "PASS_PACKET"
};

const static unsigned int RING = 0;
const static unsigned int GRID = 1;
const static unsigned int CUBE = 2;
const std::string topology_names[NUMBER_TOPOLOGY_TYPES] = {
    "RING",
    "GRID",
    "HYPERCUBE"
};

// Stats vectors
std::vector<double> simulation_times;

// Output File
std::ofstream outfile;

class Node {
    public:
        size_t x;
        size_t y;
        size_t z;
        size_t iteration;
        size_t total_packets;

        size_t current_iteration_packets = 0;
        size_t next_iteration_packets = 0;

        Node(size_t x, size_t y, size_t z, size_t iteration, size_t total_packets) : x(x), y(y), z(z), iteration(iteration), total_packets(total_packets) {}

        bool operator==(const Node& j) const {
            return (x == j.x && y == j.y && z == j.z);
        }

        bool operator!=(const Node& j) const {
            return (x != j.x || y != j.y || z != j.z);
        }

        friend std::ostream& operator<< (std::ostream& out, Node& node);
};

// Print event
std::ostream& operator<< ( std::ostream& out, Node& node) {
    out << "[x:" << node.x << " y:" << node.y << " z:" << node.z;
    out << std::right << "]";
    return out;
}

class Packet {
    public:
        size_t x;
        size_t y;
        size_t z;
        size_t size;
        size_t iteration;

        Packet(size_t x, size_t y, size_t z, size_t iteration) : x(x), y(y), z(z), iteration(iteration) {}

        bool operator==(const Packet& j) const {
            return (x == j.x && y == j.y && z == j.z);
        }

        bool operator!=(const Packet& j) const {
            return (x != j.x || y != j.y || z != j.z);
        }

        friend std::ostream& operator<< (std::ostream& out, Packet& packet);
};

// Print event
std::ostream& operator<< ( std::ostream& out, Packet& packet) {
    out << "[x:" << packet.x << " y:" << packet.y << " z:" << packet.z;
    out << std::right << "]";
    return out;
}

/**
 *  A class for an event, which holds an int for the event type, a double for simulation time and possibly other data.
 *  We may want to subclass this with our own events
 */
class Event {
    public:
        const double time;
        const int type;
        Node *node;
        Packet *packet = nullptr;

        Event(double time, int type) : time(time), type(type) {
            node = nullptr;
            packet = nullptr;
            //cout << "created an event with simulation time: " << this->time << endl;
            //cout << "created an event with event type: " << this->type << endl;
        }

        Event(double time, int type, Node *node) : time(time), type(type), node(node) {
            packet = nullptr;
            //cout << "created an event with simulation time: " << this->time << endl;
            //cout << "created an event with event type: " << this->type << endl;
        }

        Event(double time, int type, Node *node, Packet *packet) : time(time), type(type), node(node), packet(packet) {
            //cout << "created an event with simulation time: " << this->time << endl;
            //cout << "created an event with event type: " << this->type << endl;
        }

        bool operator<( const Event& e ) const {
            return time < e.time;
        }

        bool operator>( const Event& e ) const {
            return time > e.time;
        }

        bool less( const Event& e) const {
            return time < e.time;
        }

        friend std::ostream& operator<< (std::ostream& out, Event& event);
};

class CompareEvent {
    public:
        bool operator()(Event* &e1, Event* &e2) {
            return e1->time > e2->time;
        }
};

// Print event
std::ostream& operator<< ( std::ostream& out, Event& event) {
    out << "[sim_time: " << std::setw(10) << std::setprecision(3) << std::fixed << event.time << ", type: " << std::setw(4) << event.type;
    if (event.type < NUMBER_EVENT_TYPES) {
        out << " - " << std::left << std::setw(50) << event_names[event.type];
    } else {
        out << std::left << std::setw(50) << " - UNKNOWN";
    }
    out << std::right << "]";
    return out;
}

Node* getNextNode(Node *node, Packet *packet, unsigned int topology, size_t num_nodes, std::vector<std::vector<std::vector<Node*>*>*> node_mat) {
    Node *next_node = node;
    VLOG(2) << "NODE: " << *node;
    VLOG(2) << "PCKT: " << *packet;
    switch (topology) {
        case 0: // RING
            if (packet->x > node->x) {
                if (packet->x - node->x > num_nodes - packet->x + node->x) {
                    // Get Node with x - 1
                    if (node->x == 0) {
                        next_node = node_mat.at(num_nodes-1)->at(node->y)->at(node->z);
                    } else {
                        next_node = node_mat.at(node->x-1)->at(node->y)->at(node->z);
                    }
                } else {
                    // Get Node with x + 1
                    if (node->x == num_nodes - 1) {
                        next_node = node_mat.at(0)->at(node->y)->at(node->z);
                    } else {
                        next_node = node_mat.at(node->x+1)->at(node->y)->at(node->z);
                    }
                }
            } else {
                if (node->x - packet->x > num_nodes - node->x + packet->x) {
                    // Get Node with x + 1
                    if (node->x == num_nodes - 1) {
                        next_node = node_mat.at(0)->at(node->y)->at(node->z);
                    } else {
                        next_node = node_mat.at(node->x+1)->at(node->y)->at(node->z);
                    }
                } else {
                    // Get Node with x - 1
                    if (node->x == 0) {
                        next_node = node_mat.at(num_nodes-1)->at(node->y)->at(node->z);
                    } else {
                        next_node = node_mat.at(node->x-1)->at(node->y)->at(node->z);
                    }
                }
            }
            break;
        case 1: // GRID
            if (node->x != packet->x) {
                // Move packet in x direction
                if (packet->x > node->x) {
                    // Down
                    next_node = node_mat.at(node->x+1)->at(node->y)->at(node->z);
                } else {
                    // Up
                    next_node = node_mat.at(node->x-1)->at(node->y)->at(node->z);
                }
            } else if (node->y != packet->y) {
                // Move packet in y direction
                if (packet->y > node->y) {
                    // Right
                    next_node = node_mat.at(node->x)->at(node->y+1)->at(node->z);
                } else {
                    // Left
                    next_node = node_mat.at(node->x)->at(node->y-1)->at(node->z);
                }
            } else {
                LOG(FATAL) << "Packet already at destination!";
            }
            break;
        case 2: // CUBE
            if (node->x != packet->x) {
                // Move packet in x direction
                if (packet->x > node->x) {
                    // Down
                    next_node = node_mat.at(node->x+1)->at(node->y)->at(node->z);
                } else {
                    // Up
                    next_node = node_mat.at(node->x-1)->at(node->y)->at(node->z);
                }
            } else if (node->y != packet->y) {
                // Move packet in y direction
                if (packet->y > node->y) {
                    // Right
                    next_node = node_mat.at(node->x)->at(node->y+1)->at(node->z);
                } else {
                    // Left
                    next_node = node_mat.at(node->x)->at(node->y-1)->at(node->z);
                }
            } else if (node->z != packet->z) {
                // Move packet in z direction
                if (packet->z > node->z) {
                    // Foreward
                    next_node = node_mat.at(node->x)->at(node->y)->at(node->z+1);
                } else {
                    // Back
                    next_node = node_mat.at(node->x)->at(node->y)->at(node->z-1);
                }
            } else {
                LOG(FATAL) << "Packet already at destination!";
            }
            break;
        default:
            LOG(FATAL) << "No topology with id '" << topology  << "' exists.";
    }
    VLOG(2) << "NEXT: " << *next_node;
    return next_node;
}

double run_simulation(size_t nodes, size_t data, double compute_time, double latency, unsigned int topology, size_t iterations) {
    size_t nodes_x;
    size_t nodes_y;
    size_t nodes_z;

    if (topology == RING) {
        nodes_x = nodes;
        nodes_y = 1;
        nodes_z = 1;
    } else if (topology == GRID) {
        nodes_x = nodes;
        nodes_y = nodes;
        nodes_z = 1;
    } else if (topology == CUBE) {
        nodes_x = nodes;
        nodes_y = nodes;
        nodes_z = nodes;
    }


    double iteration = 0;
    double simulation_time_s = 0;
    double previous_time_s = 0;

    // Init RNG
    std::random_device rd;
    std::mt19937 generator(rd());
    std::exponential_distribution<double> compute_distribution(1/compute_time);
    std::exponential_distribution<double> send_distribution(1/latency);

    std::priority_queue<Event*, std::vector<Event*>, CompareEvent> heap;

    size_t packets_per_node = ceil(data/nodes);

    std::vector<std::vector<std::vector<Node*>*>*> node_mat;

    // Put initial events in the heap and set worker error delays
    for (size_t x = 0; x < nodes_x; x++) {
        std::vector<std::vector<Node*>*> *x_vec = new std::vector<std::vector<Node*>*>();
        node_mat.push_back(x_vec);
        for (size_t y = 0; y < nodes_y; y++) {
            std::vector<Node*> *y_vec = new std::vector<Node*>();
            x_vec->push_back(y_vec);
            for (size_t z = 0; z < nodes_z; z++) {
                Node *node = new Node(x, y, z, 0, packets_per_node);
                y_vec->push_back(node);
                heap.push(new Event(simulation_time_s + (compute_distribution(generator) * data), COMPUTE, node));
            }
        }
    }

    VLOG(2) << "Start Simulation...";

    while (!heap.empty() && iteration < iterations) {
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
        Packet *packet = current_event->packet;

        switch (current_event->type) {
            case 0: // COMPUTE
                VLOG(2) << "Compute Begin";
                LOG_IF(ERROR, packet != nullptr) << "Packet was not freed!";
                for (size_t x = 0; x < node_mat.size(); x++) {
                    std::vector<std::vector<Node*>*> *y_vec = node_mat.at(x);
                    for (size_t y = 0; y < y_vec->size(); y++) {
                        std::vector<Node*> *z_vec = y_vec->at(y);
                        for (size_t z = 0; z < z_vec->size(); z++) {
                            if (*node != *(z_vec->at(z))) {
                                Node *dest = z_vec->at(z);
                                for (size_t i = 0; i < node->total_packets; i++) {
                                    Packet *new_packet = new Packet(dest->x, dest->y, dest->z, node->iteration + 1);
                                    heap.push(new Event(simulation_time_s + send_distribution(generator), PASS_PACKET, node, new_packet));
                                }
                            }
                        }
                    }
                }
                VLOG(2) << "Compute iteration " << iteration << " End";
                break;
            case 1: // PASS_PACKET
                VLOG(2) << "Pass Packet Begin";
                if (packet->x == node->x && packet->y == node->y && packet->z == node->z) { // Packet is at its destination.
                    VLOG(2) << "Packet is at its destination";
                    if (packet->iteration == node->iteration) { // Received a packet for this iteration
                        node->current_iteration_packets++;
                        if (node->total_packets*(nodes-1) == node->current_iteration_packets - node->total_packets) {
                            VLOG(2) << "Node ready to compute, start compute.";
                            heap.push(new Event(simulation_time_s + (compute_distribution(generator) * data), COMPUTE));
                        }
                    } else if (packet->iteration == node->iteration + 1) { // Received a packet for the next interation
                        node->next_iteration_packets++;
                    } else {
                        LOG(ERROR) << "Incorrect packet iteration for current node!";
                    }
                    delete packet;
                    packet = nullptr;
                } else { // Packet needs to be passed on to next node
                    Node *next_node = getNextNode(node, packet, topology, nodes, node_mat);
                    VLOG(2) << *next_node;
                    heap.push(new Event(simulation_time_s + send_distribution(generator), PASS_PACKET, next_node, packet));
                }
                VLOG(2) << "Pass Packet End";
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
                delete z_vec->back();
                z_vec->pop_back();
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
    size_t iterations = 100000;
    size_t nodes = 5;
    size_t data = 200;
    double compute_time= 0.005;
    double latency = 0.00001;
    //unsigned int topology = RING;
    //unsigned int topology = GRID;
    unsigned int topology = CUBE;

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
            double simulation_time = run_simulation(nodes, data, compute_time, latency, topology, iterations);
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
