#ifndef NODE_H
#define NODE_H

#include <cstddef>
#include <iostream>
#include <vector>

#include "stream.hpp"

#include <glog/logging.h>

enum class Topology {RING, GRID, CUBE, SIZE};

class Node {
    private:
        size_t x;
        size_t y;
        size_t z;
        size_t iteration;
        size_t num_packets_to_send;
        size_t num_packets_to_receive;

        bool ready_to_compute = false;

        double mean_iteration_time = 0;
        double previous_iteration_duration = 0;
        double previous_compute_time = 0;
        double total_wait_time = 0;
        double total_compute_time = 0;

        size_t current_iteration_packets = 0;
        size_t next_iteration_packets = 0;

        // Out-bound queues
        std::vector<Stream*> streams;

    public:
        Node(size_t x, size_t y, size_t z, size_t iteration, size_t num_packets_to_send, size_t num_packets_to_receive);
        ~Node();
        size_t getX();
        size_t getY();
        size_t getZ();
        size_t getIteration();
        size_t getNumPacketsToSend();
        bool readyToCompute();
        double getWaitTime();
        double getComputeTime();
        double getComputeRatio();
        void startNextIteration();
        void updateComputeTime(double sim_time, double compute_duration);
        void receivePacket(Packet* packet);

        Stream* enqueuePacket(Packet* packet, Topology topology, size_t num_nodes);
        Node* getNextNode(Stream *stream, Topology topology, size_t num_nodes, std::vector<std::vector<std::vector<Node*>*>*> node_mat);

        void writeStats(size_t pass_num);

        bool operator==(const Node& j) const;
        bool operator!=(const Node& j) const;

        friend std::ostream& operator<< (std::ostream& out, Node& node);
};

#endif //NODE_H
