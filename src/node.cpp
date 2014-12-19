#include "node.hpp"

#include <sstream>
#include <fstream>

#include <glog/logging.h>

Node::Node(size_t x, size_t y, size_t z, size_t iteration, size_t num_packets_to_send, size_t num_packets_to_receive) : x(x), y(y), z(z), iteration(iteration), num_packets_to_send(num_packets_to_send), num_packets_to_receive(num_packets_to_receive) {
    for (size_t d = 0; d < static_cast<size_t>(StreamDirection::SIZE); d++) {
        streams.push_back(new Stream(static_cast<StreamDirection>(d)));
    }
}

Node::~Node() {
    while (!streams.empty()) {
        delete streams.back();
        streams.pop_back();
    }
}

size_t Node::getX() {
    return this->x;
}

size_t Node::getY() {
    return this->y;
}

size_t Node::getZ() {
    return this->z;
}

size_t Node::getIteration() {
    return this->iteration;
}

size_t Node::getNumPacketsToSend() {
    return this->num_packets_to_send;
}

bool Node::readyToCompute() {
    return this->ready_to_compute;
}

double Node::getWaitTime() {
    return this->total_wait_time;
}

double Node::getComputeTime() {
    return this->total_compute_time;
}

double Node::getComputeRatio() {
    return this->total_wait_time/(this->total_wait_time + this->total_compute_time);
}

void Node::startNextIteration() {
    this->iteration++;
    this->current_iteration_packets = this->next_iteration_packets;
    this->next_iteration_packets = 0;
    this->ready_to_compute = false;
}

void Node::updateComputeTime(double sim_time, double compute_duration) {
    this->total_compute_time += compute_duration;
    this->total_wait_time = sim_time + compute_duration - total_compute_time;

    this->mean_iteration_time = ((sim_time - previous_compute_time) + (this->iteration * this->mean_iteration_time)) / (this->iteration + 1);
    this->previous_iteration_duration = sim_time - previous_compute_time;
    this->previous_compute_time = sim_time;
}

void Node::receivePacket(Packet* packet) {
    if (packet->getIteration() == this->getIteration()) { // Received a packet for this iteration
        this->current_iteration_packets++;
        VLOG(2) << "Received packet for this iteration: " << *this << " " << this->num_packets_to_receive << " : " << this->current_iteration_packets;
        if (this->num_packets_to_receive == this->current_iteration_packets) {
            VLOG(1) << "Node" << *this << " ready to compute: " << this->num_packets_to_receive << " : " << this->current_iteration_packets;
            this->ready_to_compute = true;
        }
    } else if (packet->getIteration() == this->getIteration() + 1) { // Received a packet for the next interation
        VLOG(2) << "Received packet for the next iteration: " << *this << " " << this->num_packets_to_receive << " : " << this->current_iteration_packets;
        this->next_iteration_packets++;
    } else {
        LOG(ERROR) << "Incorrect packet iteration for node (Packet: " << packet->getIteration() << ", Node: " << this->getIteration() << ")";
    }
}

Stream* Node::enqueuePacket(Packet* packet, Topology topology, size_t num_nodes) {
    LOG_IF(FATAL, packet == nullptr) << "Packet is null!";
    VLOG(3) << "PCKT: " << *packet;
    LOG_IF(FATAL, packet->getX() > 1000);
    LOG_IF(FATAL, packet->getY() > 1000);
    LOG_IF(FATAL, packet->getZ() > 1000);
    std::vector<Stream*> queues;
    switch (topology) {
        case Topology::RING:
            if (packet->getX() > this->x) {
                if (packet->getX() - this->x > num_nodes - packet->getX() + this->x) {
                    // Get Node with x - 1
                    if (this->x == 0) {
                        queues.push_back(this->streams[static_cast<size_t>(StreamDirection::BACK)]);
                    } else {
                        queues.push_back(this->streams[static_cast<size_t>(StreamDirection::BACK)]);
                    }
                } else {
                    // Get Node with x + 1
                    if (this->x == num_nodes - 1) {
                        queues.push_back(this->streams[static_cast<size_t>(StreamDirection::FOREWARD)]);
                    } else {
                        queues.push_back(this->streams[static_cast<size_t>(StreamDirection::FOREWARD)]);
                    }
                }
            } else {
                if (this->x - packet->getX() > num_nodes - this->x + packet->getX()) {
                    // Get Node with x + 1
                    if (this->x == num_nodes - 1) {
                        queues.push_back(this->streams[static_cast<size_t>(StreamDirection::FOREWARD)]);
                    } else {
                        queues.push_back(this->streams[static_cast<size_t>(StreamDirection::FOREWARD)]);
                    }
                } else {
                    // Get Node with x - 1
                    if (this->x == 0) {
                        queues.push_back(this->streams[static_cast<size_t>(StreamDirection::BACK)]);
                    } else {
                        queues.push_back(this->streams[static_cast<size_t>(StreamDirection::BACK)]);
                    }
                }
            }
            break;
        case Topology::GRID:
            if (this->x != packet->getX()) {
                // Move packet in x direction
                if (packet->getX() > this->x) {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::FOREWARD)]); // Foreward
                } else {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::BACK)]); // Back
                }
            }
            if (this->y != packet->getY()) {
                // Move packet in y direction
                if (packet->getY() > this->y) {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::RIGHT)]); // Right
                } else {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::LEFT)]); // Left
                }
            }
            break;
        case Topology::CUBE:
            if (this->x != packet->getX()) {
                // Move packet in x direction
                if (packet->getX() > this->x) {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::FOREWARD)]); // Foreward
                } else {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::BACK)]); // Back
                }
            }
            if (this->y != packet->getY()) {
                // Move packet in y direction
                if (packet->getY() > this->y) {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::RIGHT)]); // Right
                } else {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::LEFT)]); // Left
                }
            }
            if (this->z != packet->getZ()) {
                // Move packet in z direction
                if (packet->getZ() > this->z) {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::DOWN)]); // Down
                } else {
                    queues.push_back(this->streams[static_cast<size_t>(StreamDirection::UP)]); // Up
                }
            }
            break;
        default:
            LOG(FATAL) << "No topology with id '" << static_cast<size_t>(topology)  << "' exists.";
    }
    // Insert packet into viable queue with shortest length
    Stream* stream = nullptr;
    if (!queues.empty()) {
        for (size_t i = 0; i < queues.size(); i++) {
            if (stream == nullptr || queues.at(i)->size() < stream->size()) {
                LOG_IF(FATAL, queues.at(i) == nullptr) << "Queue is null";
                stream = queues.at(i);
            }
        }
        LOG_IF(FATAL, packet->getX() > 1000);
        LOG_IF(FATAL, packet->getY() > 1000);
        LOG_IF(FATAL, packet->getZ() > 1000);
        stream->push(packet);
    } else {
        LOG(FATAL) << "Packet already at destination!";
    }
    return stream;
}

Node *Node::getNextNode(
        Stream *stream,
        Topology topology,
        size_t num_nodes,
        std::vector<std::vector<std::vector<Node*>*>*> node_mat) {
    Node *next_node = nullptr;
    switch (topology) {
        case Topology::RING:
            if (stream->getDirection() == StreamDirection::BACK) {
                // Get Node with x - 1
                if (this->x == 0) {
                    next_node = node_mat.at(num_nodes-1)->at(this->y)->at(this->z);
                } else {
                    next_node = node_mat.at(this->x-1)->at(this->y)->at(this->z);
                }
            } else if (stream->getDirection() == StreamDirection::FOREWARD) {
                // Get Node with x + 1
                if (this->x == num_nodes - 1) {
                    next_node = node_mat.at(0)->at(this->y)->at(this->z);
                } else {
                    next_node = node_mat.at(this->x+1)->at(this->y)->at(this->z);
                }
            } else {
                LOG(ERROR) << "Invalid stream";
            }
            break;
        case Topology::GRID:
            if (stream->getDirection() == StreamDirection::FOREWARD) {
                next_node = node_mat.at(this->x+1)->at(this->y)->at(this->z);
            } else if (stream->getDirection() == StreamDirection::BACK) {
                next_node = node_mat.at(this->x-1)->at(this->y)->at(this->z);
            } else if (stream->getDirection() == StreamDirection::RIGHT) {
                next_node = node_mat.at(this->x)->at(this->y+1)->at(this->z);
            } else if (stream->getDirection() == StreamDirection::LEFT) {
                next_node = node_mat.at(this->x)->at(this->y-1)->at(this->z);
            } else {
                LOG(ERROR) << "Invalid stream";
            }
            break;
        case Topology::CUBE:
            if (stream->getDirection() == StreamDirection::FOREWARD) {
                next_node = node_mat.at(this->x+1)->at(this->y)->at(this->z);
            } else if (stream->getDirection() == StreamDirection::BACK) {
                next_node = node_mat.at(this->x-1)->at(this->y)->at(this->z);
            } else if (stream->getDirection() == StreamDirection::RIGHT) {
                next_node = node_mat.at(this->x)->at(this->y+1)->at(this->z);
            } else if (stream->getDirection() == StreamDirection::LEFT) {
                next_node = node_mat.at(this->x)->at(this->y-1)->at(this->z);
            } else if (stream->getDirection() == StreamDirection::DOWN) {
                next_node = node_mat.at(this->x)->at(this->y)->at(this->z+1);
            } else if (stream->getDirection() == StreamDirection::UP) {
                next_node = node_mat.at(this->x)->at(this->y)->at(this->z-1);
            } else {
                LOG(ERROR) << "Invalid stream";
            }
            break;
        default:
            LOG(FATAL) << "No topology with id '" << static_cast<size_t>(topology)  << "' exists.";
    }
    VLOG(3) << "NEXT: " << *next_node;
    return next_node;
}

void Node::writeStats(size_t pass_num) {
    std::stringstream filename;
    filename << pass_num << ":" << this->getX() << "_" << this->getY() << "_" << getZ() << ".dat";
    std::ofstream outfile(filename.str(), std::ios::app);
    outfile << this->iteration << "\t" << this->previous_iteration_duration << "\t" << this->mean_iteration_time << std::endl;
    outfile.close();
}

bool Node::operator==(const Node& j) const {
    return (x == j.x && y == j.y && z == j.z);
}

bool Node::operator!=(const Node& j) const {
    return (x != j.x || y != j.y || z != j.z);
}
