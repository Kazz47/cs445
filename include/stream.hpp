#ifndef STREAM_H
#define STREAM_H

#include <cstddef>
#include <queue>

#include "packet.hpp"

enum class StreamDirection {BACK, FOREWARD, UP, DOWN, RIGHT, LEFT, SIZE};

class Stream {
    private:
        StreamDirection direction;
        std::queue<Packet*> *stream;
        Packet* active_packet = nullptr;

    public:
        Stream(StreamDirection dir);
        ~Stream();
        StreamDirection getDirection();
        size_t size();
        void push(Packet* packet);

        void transferPacket();
        Packet* servicePacket();

        bool needsInit();
};

#endif //STREAM_H
