#ifndef STREAM_H
#define STREAM_H

#include <queue>

#include "packet.hpp"

enum class StreamDirection {BACK, FOREWARD, UP, DOWN, RIGHT, LEFT, SIZE};

class Stream {
    private:
        StreamDirection direction;
        std::queue<Packet*> *stream;
        bool running = false;

    public:
        Stream(StreamDirection dir) {
            this->direction = dir;
            stream = new std::queue<Packet*>;
        }

        ~Stream() {
            delete stream;
        }

        StreamDirection getDirection() {
            return this->direction;
        }

        size_t size() {
            return stream->size();
        }

        void push(Packet* packet) {
            stream->push(packet);
        }

        Packet* pop() {
            Packet* result = nullptr;
            if (!stream->empty()) {
                result = stream->back();
                stream->pop();
            }
            return result;
        }

        bool needsInit() {
            if (!stream->empty() && !running) {
                return false;
            } else {
                return true;
            }
        }
};

#endif //STREAM_H
