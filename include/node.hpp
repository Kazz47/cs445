#ifndef NODE_H
#define NODE_H

#include <vector>

#include "stream.hpp"

#include <glog/logging.h>

class Node {
    private:
        size_t x;
        size_t y;
        size_t z;
        size_t iteration;
        size_t num_packets_to_send;

        double previous_compute_time = 0;
        double total_wait_time = 0;
        double total_compute_time = 0;

        size_t current_iteration_packets = 0;
        size_t next_iteration_packets = 0;

        // Out-bound queues
        std::vector<Stream*> streams;

    public:
        Node(size_t x, size_t y, size_t z, size_t iteration, size_t num_packets_to_send) : x(x), y(y), z(z), iteration(iteration), num_packets_to_send(num_packets_to_send) {
            for (size_t d = 0; d < static_cast<size_t>(StreamDirection::SIZE); d++) {
                streams.push_back(new Stream(static_cast<StreamDirection>(d)));
            }
        }

        ~Node() {
            while (!streams.empty()) {
                delete streams.back();
                streams.pop_back();
            }
        }

        size_t getNumPacketsToSend() {
            return this->num_packets_to_send;
        }

        void startNextIteration() {
            this->iteration++;
            this->current_iteration_packets = this->next_iteration_packets;
            this->next_iteration_packets = 0;
        }

        void updateComputeTime(double sim_time, double compute_duration) {
            this->total_compute_time += compute_duration;
            this->total_wait_time += sim_time - this->previous_compute_time;
            this->previous_compute_time = sim_time + compute_duration;
        }

        Stream* enqueuePacket(Packet* packet, unsigned int topology, size_t num_nodes) {
            VLOG(3) << "PCKT: " << *packet;
            std::vector<Stream*> queues;
            switch (topology) {
                case 0: // RING
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
                case 1: // GRID
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
                case 2: // CUBE
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
                    LOG(FATAL) << "No topology with id '" << topology  << "' exists.";
            }
            // Insert packet into viable queue with shortest length
            Stream* stream= nullptr;
            if (!queues.empty()) {
                for (size_t i = 0; i < queues.size(); i++) {
                    if (stream == nullptr || queues.at(i)->size() < stream->size()) {
                        stream= queues.at(i);
                    }
                }
                stream->push(packet);
            } else {
                LOG(FATAL) << "Packet already at destination!";
            }
            return stream;
        }

        Node *getNextNode(
                Stream *stream,
                unsigned int topology,
                size_t num_nodes,
                std::vector<std::vector<std::vector<Node*>*>*> node_mat) {
            Node *next_node = nullptr;
            switch (topology) {
                case 0: // RING
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
                case 1: // GRID
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
                case 2: // CUBE
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
                    LOG(FATAL) << "No topology with id '" << topology  << "' exists.";
            }
            VLOG(3) << "NEXT: " << *next_node;
            return next_node;
        }

        bool operator==(const Node& j) const {
            return (x == j.x && y == j.y && z == j.z);
        }

        bool operator!=(const Node& j) const {
            return (x != j.x || y != j.y || z != j.z);
        }

        friend std::ostream& operator<< (std::ostream& out, Node& node);
};

#endif //NODE_H
