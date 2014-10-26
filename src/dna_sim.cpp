#include <iostream>
#include <fstream>
#include <queue>
#include <string>
#include <limits>
#include <cmath>

#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>

#include <glog/logging.h>

using boost::variate_generator;
using boost::mt19937;
using boost::lognormal_distribution;

#define ONLY_EVENT_TYPE     0
#define NUMBER_EVENT_TYPES  3   //NEED TO MAKE SURE THIS IS 1 MORE THAN THE LAST DEFINED EVENT

const std::string event_names[NUMBER_EVENT_TYPES] = {
    "REQUEST_JOB",
    "RETURN_SUCCESS",
    "RETURN_ERROR"
};

std::vector<bool> job_queue;
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
            return e1->time > e2->time;
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

void parseFile(const std::string &file_name, std::vector<float> &means, std::vector<float> &stdevs, std::vector<float> &errors) {
    LOG(INFO) << "Load file: '" << file_name << "'";
    std::ifstream infile(file_name);

    std::vector<unsigned long> start_times, end_times, result_flags, wu_names;
    std::vector<float> cpu_times;

    std::string start_time, end_time, result_flag, cpu_time, wu_name;

    std::string line;
    std::string temp;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (!(iss >> start_time >> end_time >> cpu_time >> result_flag >> wu_name)) {
            break;
        }
        start_time.pop_back();
        VLOG(1) << "START_TIME: " << start_time;
        end_time.pop_back();
        VLOG(1) << "END_TIME: " << end_time;
        cpu_time.pop_back();
        VLOG(1) << "CPU_TIME: " << cpu_time;
        result_flag.pop_back();
        VLOG(1) << "RESULT_FLAG: " << result_flag;
        start_times.push_back(atoi(start_time.c_str()));
        end_times.push_back(atoi(end_time.c_str()));
        cpu_times.push_back(atof(cpu_time.c_str()));
        result_flags.push_back(atoi(result_flag.c_str()));

        // Order by names with more zeros.
        if (wu_name.find("gibbs_test_hg19_1000") == 0) {
            wu_names.push_back(1000);
        } else if (wu_name.find("gibbs_test_hg19_100") == 0) {
            wu_names.push_back(100);
        } else if (wu_name.find("gibbs_test_hg19_10") == 0) {
            wu_names.push_back(10);
        } else {
            LOG(ERROR) << "No matching name: '" << wu_name << "'";
            wu_names.push_back(0);
        }
    }
    infile.close();

    unsigned long sum[3] = {0, 0, 0};
    unsigned long size[3] = {0, 0, 0};
    unsigned long size_errors[3] = {0, 0, 0};
    for (size_t i = 0; i < wu_names.size(); i++) {
        long duration = end_times[i] - start_times[i];
        if (wu_names[i] == 1000) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                sum[0] += duration;
                size[0]++;
            } else {
                size_errors[0]++;
            }
        } else if(wu_names[i] == 100) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                sum[1] += duration;
                size[1]++;
            } else {
                size_errors[1]++;
            }
        } else if(wu_names[i] == 10) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                sum[2] += duration;
                size[2]++;
            } else {
                size_errors[2]++;
            }
        }
    }

    for (size_t i = 0; i < 3; i++) {
        means.push_back(static_cast<float>(sum[i])/size[i]);
        errors.push_back(static_cast<float>(size_errors[i])/(size_errors[i] + size[i]));
    }

    unsigned long stdev_sum[3] = {0, 0, 0};
    for (size_t i = 0; i < wu_names.size(); i++) {
        if (result_flags[i] == 1) {
            int duration = end_times[i] - start_times[i];
            LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
            if (wu_names[i] == 1000) {
                stdev_sum[0] += (duration - means[0]) * (duration - means[0]);
            } else if(wu_names[i] == 100) {
                stdev_sum[1] += (duration - means[1]) * (duration - means[1]);
            } else if(wu_names[i] == 10) {
                stdev_sum[2] += (duration - means[2]) * (duration - means[2]);
            }
        }
    }

    for (size_t i = 0; i < 3; i++) {
        stdevs.push_back(sqrt(static_cast<float>(stdev_sum[i])/size[i]));
    }
}

void run_simulation(
        variate_generator< mt19937, lognormal_distribution<> > &duration_generator,
        variate_generator< mt19937, std::uniform_real_distribution<> > &error_generator) {
    double simulation_time_s = 0;
    double previous_time_s = 0;

    //std::vector<double> servers_work_time;
    //std::vector<double> servers_time_queue_length;
    /*
    for (int i = 0; i < servers_idle.size(); i++) {
        servers_work_time.push_back(0);
        servers_time_queue_length.push_back(0);
    }
    */

    std::priority_queue<Event*, std::vector<Event*>, CompareEvent> heap;

    // Put initial event in the heap
    heap.push(new Event(simulation_time_s + duration_generator(), 0));

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

        switch (current_event->type) {
            case 0: // REQUESET_JOB
                // If there is a job push a new event for either success or
                // error. If there is no job in the queue then push another job
                // request onto the queue.
                VLOG(2) << "Request Job Begin";
                if (!job_queue.empty()) {
                    // Add another request event
                    heap.push(new Event(simulation_time_s + error_generator(), 0));
                }
                VLOG(2) << "Request Job End";
                break;
            case 1: // RETURN_SUCCESS
                // Server checks if there the quarum for this sample is
                // complete if it is then create jobs for the next sample in
                // the walk. If there are jobs available create a new success
                // or error event otherwise a new request event.
                VLOG(2) << "Success Begin";
                heap.push(new Event(simulation_time_s + duration_generator(), simulation_time_s));
                VLOG(2) << "Success End";
                break;
            case 2: // RETURN_ERROR
                // Server creates a new job for the errored sample. If there
                // are available jobs create a success or error event. Otherwise create a
                // request event.
                VLOG(2) << "Error Begin";
                // Do stuff here
                VLOG(2) << "Error End";
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
}

int main(int argc, char **argv) {
    // Initialize Google Logging
    google::InitGoogleLogging(argv[0]);
    // Log to Stderr
    FLAGS_logtostderr = 1;

    std::vector<float> means;
    std::vector<float> stdevs;
    std::vector<float> errors;
    parseFile("../data/dna_workunit_transit.csv", means, stdevs, errors);

    int seed = time(0);
    variate_generator< mt19937, lognormal_distribution<> > duration_generator(mt19937(seed), lognormal_distribution<>(means[0], stdevs[0]));
    variate_generator< mt19937, std::uniform_real_distribution<> > error_generator(mt19937(seed+1), std::uniform_real_distribution<>(0, 1));

    for (size_t i = 0; i < 3; i++) {
        LOG(INFO) << i << ":Mean:" << means[i];
        LOG(INFO) << i << ":Stdevs:" << stdevs[i];
        LOG(INFO) << i << ":Errors:" << errors[i];
    }

    // Intiate
    size_t num_samples = 1;
    size_t sample_size = 1;
    size_t num_workers = 2;
    size_t num_jobs_per_sample = 2;
    size_t quorum = 2;

    // Open files
    std::stringstream filename;
    filename << "servers_stats" << ".dat";
    outfile.open(filename.str());

    for (int i = 2; i < 3; i++) {
        for (int j = 0; j < 1; j++) {
            // Clear arrays.
            //servers_idle.clear();
            //servers_queue.clear();
            for (int k = 0; k < num_workers; k++) {
                //servers_idle.push_back(true);
                //servers_queue.push_back(new std::queue<double>());
            }
            run_simulation(duration_generator, error_generator);

            for (int j = 0; j < num_workers; j++) {
                //delete servers_queue[j];
            }
        }

        // Collect and print statistics.

        /*
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
        */

    }
    outfile.close();
}

