#include "stream.hpp"

#include <glog/logging.h>

Stream::Stream(StreamDirection dir) {
    this->direction = dir;
    this->stream = new std::queue<Packet*>();
}

Stream::~Stream() {
    delete stream;
}

StreamDirection Stream::getDirection() {
    return this->direction;
}

size_t Stream::size() {
    return this->stream->size();
}

void Stream::push(Packet* packet) {
    this->stream->push(packet);
}

void Stream::transferPacket() {
    LOG_IF(FATAL, this->active_packet != nullptr);
    this->active_packet = this->stream->front();
    this->stream->pop();
}

Packet* Stream::servicePacket() {
    LOG_IF(FATAL, this->active_packet == nullptr);
    Packet* service_packet = this->active_packet;
    this->active_packet = nullptr;
    return service_packet;
}

bool Stream::needsInit() {
    VLOG(4) << "Stream size: " << this->stream->size() << " active: " << (this->active_packet != nullptr);
    if (!this->stream->empty() && this->active_packet == nullptr) {
        return true;
    } else {
        return false;
    }
}
