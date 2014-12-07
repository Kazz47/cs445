#ifndef PACKET_H
#define PACKET_H

#include <cstddef>
#include <iostream>

class Packet {
    private:
        size_t x;
        size_t y;
        size_t z;
        size_t iteration;

    public:
        Packet(size_t x, size_t y, size_t z, size_t iteration);

        size_t getX();
        size_t getY();
        size_t getZ();
        size_t getIteration();
        bool atDest(size_t x, size_t y, size_t z);

        friend std::ostream& operator<< (std::ostream& out, Packet& packet);
};

#endif //PACKET_H
