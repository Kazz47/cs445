#ifndef PACKET_H
#define PACKET_H

class Packet {
    private:
        size_t x;
        size_t y;
        size_t z;
        size_t size;
        size_t iteration;

    public:
        Packet(size_t x, size_t y, size_t z, size_t iteration) : x(x), y(y), z(z), iteration(iteration) {}

        size_t getX() {
            return this->x;
        }

        size_t getY() {
            return this->y;
        }

        size_t getZ() {
            return this->z;
        }

        size_t getSize() {
            return this->size;
        }

        size_t getIteration() {
            return this->iteration;
        }

        bool operator==(const Packet& j) const {
            return (x == j.x && y == j.y && z == j.z);
        }

        bool operator!=(const Packet& j) const {
            return (x != j.x || y != j.y || z != j.z);
        }

        friend std::ostream& operator<< (std::ostream& out, Packet& packet);
};

#endif //PACKET_H
