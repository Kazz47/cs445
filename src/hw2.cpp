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
#define NUMBER_EVENT_TYPES  3   //NEED TO MAKE SURE THIS IS 1 MORE THAN THE LAST DEFINED EVENT

const std::string event_names[NUMBER_EVENT_TYPES] = {
    "ARRIVE",
    "DEPART",
    "CLOSE"
};

bool closed = false;

std::vector<bool> servers_idle;
std::vector< std::queue<double>* > servers_queue;

/**
 *  A class for an event, which holds an int for the event type, a double for simulation time and possibly other data.
 *  We may want to subclass this with our own events
 */
class Event {
    public:
        const double time;
        const double server_start_time;
        const int server_index;
        const int type;

        Event(double time, int type) : time(time), server_start_time(0), server_index(-1), type(type) {
            //cout << "created an event with simulation time: " << this->time << endl;
            //cout << "created an event with event type: " << this->type << endl;
        }

        Event(double time, double server_start_time, int server_index, int type) : time(time), server_start_time(server_start_time), server_index(server_index), type(type) {
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
            //return false;
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


bool server_idle(int &idle_index) {
    for (int i = 0; i < servers_idle.size(); i++) {
        if (servers_idle[i] == true) {
            idle_index = i;
            return true;
        }
    }
    return false;
}

void addToQueue(double enqueue_time) {
    std::queue<double> *shortest = NULL;
    for (int i = 0; i < servers_queue.size(); i++) {
        if (shortest == NULL || servers_queue[i]->size() < shortest->size()) {
            shortest = servers_queue[i];
        }
    }
    shortest->push(enqueue_time);
}

void run_simulation(const double close_time, variate_generator< mt19937, exponential_distribution<> > rand_generator) {
    double simulation_time_s = 0;
    double previous_time_s;

    std::vector<double> servers_work_time;
    std::vector<double> servers_time_queue_length;
    for (int i = 0; i < servers_idle.size(); i++) {
        servers_work_time.push_back(0);
        servers_time_queue_length.push_back(0);
    }

    double total_queue_time = 0;
    double total_departures = 0;


    std::priority_queue<Event*, std::vector<Event*>, CompareEvent> heap;
    //std::queue<double> queue;

    // Put initial event in the heap
    heap.push(new Event(simulation_time_s + rand_generator(), 0));
    heap.push(new Event(close_time, 2));

    VLOG(2) << "Start Simulation...";

    while (!heap.empty()) {
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
        // Calculate sum of time queue length for all servers
        for (int i = 0; i < servers_idle.size(); i++) {
            servers_time_queue_length[i] += servers_queue[i]->size() * time_since_last;
        }

        int idle_index = -1;
        switch (current_event->type) {
            case 0: // ARRIVE
                // If the server is busy add new arrival to waiting queue
                // otherwise set the server as busy and add new departure time.
                VLOG(2) << "Arrival Begin";
                if (!server_idle(idle_index)) {
                    addToQueue(simulation_time_s);
                } else {
                    VLOG(2) << "Idle Index: " << idle_index;
                    servers_idle[idle_index] = false;
                    heap.push(new Event(simulation_time_s + rand_generator(), simulation_time_s, idle_index, 1));
                }
                if (!closed) {
                    // Add next arrival event
                    heap.push(new Event(simulation_time_s + rand_generator(), 0));
                }
                VLOG(2) << "Arrival End";
                break;
            case 1: // DEPART
                // Add server busy time.
                // Which server is busy?
                servers_work_time[current_event->server_index] += simulation_time_s - current_event->server_start_time;
                // If the queue is empty set the server to idle otherwise
                // get the next person from the queue and set their departure
                // time.
                VLOG(2) << "Depart Begin";
                if (servers_queue[current_event->server_index]->empty()) {
                    servers_idle[current_event->server_index] = true;
                } else {
                    // Add current simulation time minus time stored in the
                    // queue to the total queue time.
                    total_queue_time += simulation_time_s - servers_queue[current_event->server_index]->front();
                    servers_queue[current_event->server_index]->pop();
                    // Add a new depature to the heap.
                    heap.push(new Event(simulation_time_s + rand_generator(), simulation_time_s, current_event->server_index, 1));
                }
                total_departures++;
                VLOG(2) << "Depart End";
                break;
            case 2: // CLOSE
                // Set flag to stop arrivals
                VLOG(2) << "Closed";
                closed = true;
                break;
            default:
                LOG(ERROR) << "Simulation had an event with an unknown type: " << current_event->type;
                LOG(ERROR) << "\tsimulation time: " << simulation_time_s;
                exit(0);
        }

        VLOG(3) << *current_event << ", h: " << heap.size();

        delete current_event; // Event's are created with new, so we need to delete them when we're done with them
    }

    LOG(INFO) << "The simulation ended at time: " << simulation_time_s;
    LOG(INFO) << "The simulation ended " << simulation_time_s - close_time << " after close.";
    for (int i = 0; i < servers_idle.size(); i++) {
        VLOG(1) << "The server " << i+1 << " was busy for time: " << servers_work_time[i];
        LOG(INFO) << "The server " << i+1 << " utilization was: " << servers_work_time[i] / simulation_time_s;
        LOG(INFO) << "The average length of server " << i+1 << " queue was: " << servers_time_queue_length[i] / simulation_time_s;
    }
    VLOG(1) << "Total time spent in the queue: " << total_queue_time;
    VLOG(1) << "Total departures: " << total_departures;
    LOG(INFO) << "The average time spent in queue: " << total_queue_time / total_departures;
}

int main(int argc, char **argv) {
    // Initialize Google Logging
    google::InitGoogleLogging(argv[0]);
    // Log to Stderr
    FLAGS_logtostderr = 1;

    int seed = time(0);
    double distribution_mean = 3.0;
    variate_generator< mt19937, exponential_distribution<> > rand_generator(mt19937(seed), exponential_distribution<>(1/distribution_mean));

    // Init number of servers
    int num_servers = 2;
    for (int i = 0; i < num_servers; i++) {
        servers_idle.push_back(true);
        servers_queue.push_back(new std::queue<double>());
    }

    run_simulation(atoi(argv[1]), rand_generator);

    for (int i = 0; i < num_servers; i++) {
        delete servers_queue[i];
    }
}

