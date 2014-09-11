#include <iostream>
#include <vector>
#include <queue>
#include <string>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>

using boost::variate_generator;
using boost::mt19937;
using boost::exponential_distribution;

#define ONLY_EVENT_TYPE     0
#define NUMBER_EVENT_TYPES  2   //NEED TO MAKE SURE THIS IS 1 MORE THAN THE LAST DEFINED EVENT

const std::string event_names[NUMBER_EVENT_TYPES] = {
    "ARRIVE",
    "DEPART"
};

/**
 *  A class for an event, which holds an int for the event type, a double for simulation time and possibly other data.
 *  We may want to subclass this with our own events
 */
class Event {
    public:
        const double time;
        const int type;

        Event(double time, int type) : time(time), type(type) {
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

class CompareEvent {
    public:
        bool operator()(Event* &e1, Event* &e2) {
            return e1->time < e2->time;
            return false;
        }
};

std::ostream& operator<< ( std::ostream& out, Event& event) {
    out << "[sim_time: " << std::setw(10) << std::setprecision(3) << std::fixed << event.time << ", type: " << std::setw(4) << event.type;
    if (event.type < NUMBER_EVENT_TYPES) {
        out << " - " << std::left << std::setw(50) << event_names[event.type];
    } else {
        out << std::left << std::setw(50) << " - UNKNOWN";
    }
    out << std::right << "]";
    return out;
}

int main(int argc, char **argv) {
    int seed = time(0);
    double distribution_mean = 3.0;
    double simulation_time_s = 0;
    double end_simulation_time_s = 1000;
    double previous_time_s;

    variate_generator< mt19937, exponential_distribution<> > rand_generator(mt19937(seed), exponential_distribution<>(1/distribution_mean));

    std::priority_queue<Event*, std::vector<Event*>, CompareEvent> heap;

    for (int i = 0; i < 100; i++) {
        heap.push(new Event(simulation_time_s + rand_generator(), 0));
    }

    while (simulation_time_s  < end_simulation_time_s) {
        Event *current_event = heap.top();
        heap.pop();

        if (current_event == NULL) {
            std::cerr << "ERROR, simulation not complete and there are no events in the min heap." << std::endl;
            std::cerr << "\tsimulation time: " << simulation_time_s << std::endl;
            exit(0);
        }

        previous_time_s = simulation_time_s;
        simulation_time_s = current_event->time;

        switch (current_event->type) {
            case 0:
                //When we process an event, insert a new one
                heap.push(new Event(simulation_time_s + rand_generator(), 0));
                break;
            case 1:
                break;
            default:
                std::cerr << "ERROR, simulation had an event with an unknown type: " << current_event->type << std::endl;
                std::cerr << "\tsimulation time: " << simulation_time_s << std::endl;
                exit(0);
        }

        std::cout << *current_event
             << ", heap size: " << heap.size()
             << std::endl;

        delete current_event; //Event's are created with new, so we need to delete them when we're done with them
    }

    std::cout << "The simulation ended at time: " << end_simulation_time_s << std::endl;
}
