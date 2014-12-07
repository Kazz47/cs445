#include "packet.hpp"
#include "node.hpp"

#include <glog/logging.h>

Packet::Packet(size_t x, size_t y, size_t z, size_t iteration) : x(x), y(y), z(z), iteration(iteration) {
    LOG_IF(FATAL, this->getX() > 1000);
    LOG_IF(FATAL, this->getY() > 1000);
    LOG_IF(FATAL, this->getZ() > 1000);
}

size_t Packet::getX() {
    return this->x;
}

size_t Packet::getY() {
    return this->y;
}

size_t Packet::getZ() {
    return this->z;
}

size_t Packet::getIteration() {
    return this->iteration;
}

bool Packet::atDest(size_t x, size_t y, size_t z) {
    if (x == this->x && y == this->y && z == this->z){
        VLOG(3) << "Packet is at its destination: [" << x << ", " << y << ", " << z << "]";
        return true;
    } else {
        return false;
    }
}

