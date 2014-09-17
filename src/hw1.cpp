#include <iostream>
#include <queue>
#include <string>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>

#include <glog/logging.h>

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
        const double server_start_time;
        const int type;

        Event(double time, int type) : time(time), server_start_time(0), type(type) {
            //cout << "created an event with simulation time: " << this->time << endl;
            //cout << "created an event with event type: " << this->type << endl;
        }

        Event(double time, double server_start_time, int type) : time(time), server_start_time(server_start_time), type(type) {
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
            return e1->time > e2->time;
            return false;
        }
};

// Print event
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
    // Initialize Google Logging
    google::InitGoogleLogging(argv[0]);
    // Log to Stderr
    FLAGS_logtostderr = 1;

    int seed = time(0);
    double distribution_mean = 3.0;
    double simulation_time_s = 0;
    double end_simulation_time_s = atoi(argv[1]);
    double previous_time_s;

    bool server_idle = true;
    double server_work_time = 0;

    double total_queue_time = 0;
    double total_departures = 0;

    double sum_of_time_queue_length = 0;

    variate_generator< mt19937, exponential_distribution<> > rand_generator(mt19937(seed), exponential_distribution<>(1/distribution_mean));

    std::priority_queue<Event*, std::vector<Event*>, CompareEvent> heap;
    std::queue<double> queue;

    // Put initial event in the heap
    heap.push(new Event(simulation_time_s + rand_generator(), 0));

    while (simulation_time_s  < end_simulation_time_s) {
        Event *current_event = heap.top();
        heap.pop();

        if (current_event == NULL) {
            LOG(ERROR) << "Simulation not complete and there are no events in the min heap.";
            LOG(ERROR) << "\tsimulation time: " << simulation_time_s;
            exit(0);
        }

        previous_time_s = simulation_time_s;
        simulation_time_s = current_event->time;

        // Calculate average time stuff
        double time_since_last = simulation_time_s - previous_time_s;
        sum_of_time_queue_length += queue.size() * time_since_last;

        switch (current_event->type) {
            case 0: // ARRIVE
                // If the server is busy add new arrival to waiting queue
                // otherwise set the server as busy and add new departure time.
                if (!server_idle) {
                    queue.push(simulation_time_s);
                } else {
                    server_idle = false;
                    heap.push(new Event(simulation_time_s + rand_generator(), simulation_time_s, 1));
                }
                // Add next arrival event
                heap.push(new Event(simulation_time_s + rand_generator(), 0));
                break;
            case 1: // DEPART
                // Add server busy time.
                server_work_time += simulation_time_s - current_event->server_start_time;
                // If the queue is empty set the server to idle otherwise
                // get the next person from the queue and set their departure
                // time.
                if (queue.empty()) {
                    server_idle = true;
                } else {
                    // Add current simulation time minus time stored in the
                    // queue to the total queue time.
                    total_queue_time += simulation_time_s - queue.front();
                    queue.pop();
                    // Add a new depature to the heap.
                    heap.push(new Event(simulation_time_s + rand_generator(), simulation_time_s, 1));
                }
                total_departures++;
                break;
            default:
                LOG(ERROR) << "Simulation had an event with an unknown type: " << current_event->type;
                LOG(ERROR) << "\tsimulation time: " << simulation_time_s;
                exit(0);
        }

        VLOG(2) << *current_event << ", h: " << heap.size() << ", q: " << queue.size();

        delete current_event; // Event's are created with new, so we need to delete them when we're done with them
    }

    LOG(INFO) << "The simulation ended at time: " << end_simulation_time_s;
    VLOG(1) << "The server was busy for time: " << server_work_time;
    LOG(INFO) << "The server utilization was: " << server_work_time / end_simulation_time_s;
    LOG(INFO) << "The average length of the queue was: " << sum_of_time_queue_length / end_simulation_time_s;
    VLOG(1) << "Total time spent in the queue: " << total_queue_time;
    VLOG(1) << "Total departures: " << total_departures;
    LOG(INFO) << "The average time spent in queue: " << total_queue_time / total_departures;
}

