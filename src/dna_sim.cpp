#include <iostream>
#include <fstream>
#include <queue>
#include <string>
#include <limits>

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

// Stats vectors
std::vector<double> times_after_close;
std::vector<double> average_times_in_queue;
std::vector< std::vector<double>* > average_server_utilizations;
std::vector< std::vector<double>* > average_length_of_queues;

// Output File
std::ofstream outfile;

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

void parseFile(const std::string &file_name, std::vector<unsigned int> &start_times, std::vector<unsigned int> &end_times, std::vector<float> &cpu_times, std::vector<unsigned int> &result_flags, std::vector<std::string> &wu_names) {
    std::ifstream infile(file_name);
    unsigned int start_time, end_time, result_flag;
    float cpu_time;
    std::string wu_name;

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (!(iss >> start_time >> end_time >> cpu_time >> result_flag >> wu_name)) {
            break;
        }
        start_times.push_back(start_time);
        end_times.push_back(end_time);
        cpu_times.push_back(cpu_time);
        result_flags.push_back(result_flag);
        wu_names.push_back(wu_name);
    }
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

void jockey(const int &departed_queue) {
    for (int i = 1; i < servers_queue.size(); i++) {
        int pos = departed_queue - i;
        if (pos < 0) pos = servers_queue.size() + pos;
        if (servers_queue[departed_queue]->size() < servers_queue[pos]->size()) {
            servers_queue[departed_queue]->push(servers_queue[pos]->front());
            servers_queue[pos]->pop();
        }
    }
}

void run_simulation(
        const double close_time,
        variate_generator< mt19937, exponential_distribution<> > &arrival_generator,
        variate_generator< mt19937, exponential_distribution<> > &departure_generator) {
    double simulation_time_s = 0;
    double previous_time_s = 0;

    std::vector<double> servers_work_time;
    std::vector<double> servers_time_queue_length;
    for (int i = 0; i < servers_idle.size(); i++) {
        servers_work_time.push_back(0);
        servers_time_queue_length.push_back(0);
    }

    double total_queue_time = 0;
    double total_departures = 0;


    std::priority_queue<Event*, std::vector<Event*>, CompareEvent> heap;

    // Put initial event in the heap
    heap.push(new Event(simulation_time_s + arrival_generator(), 0));
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
                    heap.push(new Event(simulation_time_s + departure_generator(), simulation_time_s, idle_index, 1));
                }
                if (!closed) {
                    // Add next arrival event
                    heap.push(new Event(simulation_time_s + arrival_generator(), 0));
                }
                VLOG(2) << "Arrival End";
                break;
            case 1: // DEPART
                // Add server busy time.
                // Which server is busy?
                servers_work_time[current_event->server_index] += simulation_time_s - current_event->server_start_time;
                jockey(current_event->server_index);
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
                    heap.push(new Event(simulation_time_s + departure_generator(), simulation_time_s, current_event->server_index, 1));
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

    assert(heap.empty());

    VLOG(1) << "The simulation ended at time: " << simulation_time_s;
    times_after_close.push_back(simulation_time_s - close_time);
    VLOG(1) << "The simulation ended " << simulation_time_s - close_time << " after close.";
    average_server_utilizations.push_back(new std::vector<double>());
    average_length_of_queues.push_back(new std::vector<double>());
    VLOG(1) << "Server util vector size: " << average_server_utilizations.size();
    for (int i = 0; i < servers_idle.size(); i++) {
        VLOG(1) << "The server " << i+1 << " was busy for time: " << servers_work_time[i];
        average_server_utilizations.back()->push_back(servers_work_time[i] / simulation_time_s);
        VLOG(1) << "The server " << i+1 << " utilization was: " << servers_work_time[i] / simulation_time_s;
        average_length_of_queues.back()->push_back(servers_time_queue_length[i] / simulation_time_s);
        VLOG(1) << "The average length of server " << i+1 << " queue was: " << servers_time_queue_length[i] / simulation_time_s;
        VLOG(1) << "Server util size: " << average_server_utilizations.back()->size();
    }
    VLOG(1) << "Total time spent in the queue: " << total_queue_time;
    VLOG(1) << "Total departures: " << total_departures;
    average_times_in_queue.push_back(total_queue_time / total_departures);
    VLOG(1) << "The average time spent in queue: " << total_queue_time / total_departures;
}

int main(int argc, char **argv) {
    // Initialize Google Logging
    google::InitGoogleLogging(argv[0]);
    // Log to Stderr
    FLAGS_logtostderr = 1;

    int seed = time(0);
    double arrival_mean = 1.0;
    double departure_mean = 4.5;
    variate_generator< mt19937, exponential_distribution<> > arrival_generator(mt19937(seed), exponential_distribution<>(1/arrival_mean));
    variate_generator< mt19937, exponential_distribution<> > departure_generator(mt19937(seed+1), exponential_distribution<>(1/departure_mean));

    // Intiate
    int duration = 480;
    int min_servers = 4;
    int max_servers = 7;
    int num_iterations = 10;

    // Open files
    std::stringstream filename;
    filename << "servers_stats" << ".dat";
    outfile.open(filename.str());

    for (int sev = min_servers; sev <= max_servers; sev++) {
        int num_servers = sev;
        for (int j = 0; j < num_iterations; j++) {
            closed = false;
            servers_idle.clear();
            servers_queue.clear();
            for (int k = 0; k < num_servers; k++) {
                servers_idle.push_back(true);
                servers_queue.push_back(new std::queue<double>());
            }
            run_simulation(duration, arrival_generator, departure_generator);

            for (int j = 0; j < num_servers; j++) {
                delete servers_queue[j];
            }
        }

        // Collect and print statistics.

        double sum_average_queue_time = std::accumulate(average_times_in_queue.begin(), average_times_in_queue.end(), 0.0);
        double min_average_queue_time = *std::min_element(average_times_in_queue.begin(), average_times_in_queue.end());
        double max_average_queue_time = *std::max_element(average_times_in_queue.begin(), average_times_in_queue.end());
        double avg_average_queue_time = sum_average_queue_time / average_times_in_queue.size();
        double sqr_average_queue_time = std::inner_product(average_times_in_queue.begin(), average_times_in_queue.end(), average_times_in_queue.begin(), 0.0);
        double std_average_queue_time = std::sqrt(sqr_average_queue_time / average_times_in_queue.size() - avg_average_queue_time * avg_average_queue_time);

        double sum_times_after_close = std::accumulate(times_after_close.begin(), times_after_close.end(), 0.0);
        double min_times_after_close = *std::min_element(times_after_close.begin(), times_after_close.end());
        double max_times_after_close = *std::max_element(times_after_close.begin(), times_after_close.end());
        double avg_times_after_close = sum_times_after_close / average_times_in_queue.size();
        double sqr_times_after_close = std::inner_product(times_after_close.begin(), times_after_close.end(), times_after_close.begin(), 0.0);
        double std_times_after_close = std::sqrt(sqr_times_after_close / times_after_close.size() - avg_times_after_close * avg_times_after_close);

        double min_server_utilization = std::numeric_limits<double>::max();
        double max_server_utilization = std::numeric_limits<double>::min();
        double sum_server_utilization = 0;
        double sqr_server_utilization = 0;

        double min_average_length = std::numeric_limits<double>::max();
        double max_average_length = std::numeric_limits<double>::min();
        double sum_average_length = 0;
        double sqr_average_length = 0;

        for (int i = 0; i < num_iterations; i++) {
            double temp_avg_server_utilization = 0;
            double temp_avg_average_length = 0;

            // Find the average server utilization of the queues;
            for (int j = 0; j < num_servers; j++) {
                double current_utilization = average_server_utilizations.at(i)->at(j);
                VLOG(2) << "Server " << j+1 << " utilization: " << current_utilization;
                temp_avg_server_utilization += current_utilization;

                double current_average_length = average_length_of_queues.at(i)->at(j);
                VLOG(2) << "Server " << j+1 << " average length: " << current_average_length;
                temp_avg_average_length += current_average_length;
                temp_avg_average_length = temp_avg_average_length / num_servers;
            }
            delete average_server_utilizations.at(i);
            delete average_length_of_queues.at(i);

            temp_avg_server_utilization = temp_avg_server_utilization / num_servers;
            sum_server_utilization += temp_avg_server_utilization;
            if (temp_avg_server_utilization > max_server_utilization) max_server_utilization = temp_avg_server_utilization;
            if (temp_avg_server_utilization < min_server_utilization) min_server_utilization = temp_avg_server_utilization;
            sqr_server_utilization += temp_avg_server_utilization * temp_avg_server_utilization;

            temp_avg_average_length = temp_avg_average_length / num_servers;
            sum_average_length += temp_avg_average_length;
            if (temp_avg_average_length > max_average_length) max_average_length = temp_avg_average_length;
            if (temp_avg_average_length < min_average_length) min_average_length = temp_avg_average_length;
            sqr_average_length += temp_avg_average_length * temp_avg_average_length;
        }
        average_server_utilizations.clear();
        average_length_of_queues.clear();

        double avg_server_utilization = sum_server_utilization / num_iterations;
        double std_server_utilization = std::sqrt(sqr_server_utilization / num_iterations - avg_server_utilization * avg_server_utilization);

        double avg_average_length = sum_average_length / (num_servers * num_iterations);
        double std_average_length = std::sqrt(sqr_average_length / (num_servers * num_iterations) - avg_average_length * avg_average_length);


        outfile << "------------------------------------------" << std::endl;
        outfile << "---------------- " << num_servers << " SERVERS" << " ---------------" << std::endl;
        outfile << "------------------------------------------" << std::endl;
        outfile << std::endl;

        outfile << "-----------------------------------------" << std::endl;
        outfile << "------- AVERAGE SERVER UTILIZATION ------" << std::endl;
        outfile << "-----------------------------------------" << std::endl;
        outfile << "Min server utiliation: " << min_server_utilization << std::endl;
        outfile << "Max server utiliation: " << max_server_utilization << std::endl;
        outfile << "Average server utilization: " << avg_server_utilization << std::endl;
        outfile << "Standard deviation of server utilization: " << std_server_utilization << std::endl;
        outfile << std::endl;

        outfile << "-----------------------------------------" << std::endl;
        outfile << "---------- AVERAGE QUEUE LENGTH ---------" << std::endl;
        outfile << "-----------------------------------------" << std::endl;
        outfile << "Min queue length: " << min_average_length << std::endl;
        outfile << "Max queue length: " << max_average_length << std::endl;
        outfile << "Average queue length: " << avg_average_length << std::endl;
        outfile << "Standard deviation of queue length: " << std_average_length << std::endl;
        outfile << std::endl;

        outfile << "-----------------------------------------" << std::endl;
        outfile << "---------- AVERAGE QUEUE TIMES ----------" << std::endl;
        outfile << "-----------------------------------------" << std::endl;
        outfile << "Min of average queue times: " << min_average_queue_time << std::endl;
        outfile << "Max of average queue times: " << max_average_queue_time << std::endl;
        outfile << "Mean of average queue times: " << avg_average_queue_time << std::endl;
        outfile << "Standard deviation of average queue times: " << std_average_queue_time << std::endl;
        outfile << std::endl;

        outfile << "------------------------------------------" << std::endl;
        outfile << "-------- AVERAGE TIME AFTER CLOSE --------" << std::endl;
        outfile << "------------------------------------------" << std::endl;
        outfile << "Min of times after close: " << min_times_after_close << std::endl;
        outfile << "Max of times after close: " << max_times_after_close << std::endl;
        outfile << "Mean of times after close: " << avg_times_after_close << std::endl;
        outfile << "Standard deviation of times after close: " << std_times_after_close << std::endl;
        outfile << std::endl;
        outfile << std::endl;
        outfile << std::endl;

    }
    outfile.close();
}

