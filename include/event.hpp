#ifndef EVENT_H
#define EVENT_H

enum class EventType {COMPUTE, PASS_PACKET, SIZE};

/**
 *  A class for an event, which holds an int for the event type, a double for simulation time and possibly other data.
 *  We may want to subclass this with our own events
 */
class Event {
    public:
        const double time;
        const EventType type;
        size_t iteration;
        Node *node = nullptr;
        Stream* stream = nullptr;

        Event(double time, EventType type) : time(time), type(type) {
            node = nullptr;
            stream = nullptr;
            //cout << "created an event with simulation time: " << this->time << endl;
            //cout << "created an event with event type: " << this->type << endl;
        }

        Event(double time, EventType type, Node* node) : time(time), type(type), node(node) {
            stream = nullptr;
            iteration = node->getIteration();
            //cout << "created an event with simulation time: " << this->time << endl;
            //cout << "created an event with event type: " << this->type << endl;
        }

        Event(double time, EventType type, Node* node, Stream* stream) : time(time), type(type), node(node), stream(stream) {
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

#endif //NODE_H
