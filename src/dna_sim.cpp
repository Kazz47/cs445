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
#define NUMBER_EVENT_TYPES  3

const static unsigned int REQUEST_JOB = 0;
const static unsigned int RETURN_SUCCESS = 1;
const static unsigned int RETURN_ERROR = 2;
const std::string event_names[NUMBER_EVENT_TYPES] = {
    "REQUEST_JOB",
    "RETURN_SUCCESS",
    "RETURN_ERROR"
};

std::vector<bool> job_queue;
std::vector< std::queue<double>* > servers_queue;

size_t total_jobs = 0;
size_t jobs_complete = 0;
// Stats vectors
std::vector<double> times_after_close;
std::vector<double> average_times_in_queue;
std::vector< std::vector<double>* > average_server_utilizations;
std::vector< std::vector<double>* > average_length_of_queues;

// Output File
std::ofstream outfile;

class Job {
    public:
        bool done;
        size_t walk;
        size_t sample;
        size_t job_id;
        std::vector<size_t> failed_clients;

        Job() : walk(0), sample(0), job_id(0) {}
        Job(size_t walk, size_t sample, size_t job) : done(false), walk(walk), sample(sample), job_id(job) {}

        void successClient(size_t client) {
            failed_clients.push_back(client);
            done = true;
            jobs_complete++;
        }

        void failClient(size_t client) {
            //failed_clients.push_back(client);
        }

        bool clientOkay(size_t client) {
            for (size_t i = 0; i < failed_clients.size(); i++) {
                if (client == failed_clients.at(i)) {
                    VLOG(1) << "Client not okay: (" << client << ")";
                    return false;
                }
            }
            return true;
        }

        bool operator==(const Job& j) const {
            return (walk == j.walk && sample == j.sample && job_id == j.job_id);
        }

        friend std::ostream& operator<< (std::ostream& out, Job& job);
};

// Print event
std::ostream& operator<< ( std::ostream& out, Job& job) {
    out << "[w:" << job.walk << " s:" << job.sample << " id:" << job.job_id;
    out << std::right << "]";
    return out;
}

/**
 *  A class for an event, which holds an int for the event type, a double for simulation time and possibly other data.
 *  We may want to subclass this with our own events
 */
class Event {
    public:
        const double time;
        const int type;
        const size_t client;
        size_t wait_time;
        Job *job;

        Event(double time, int type, size_t client, size_t prev_wait_time) : time(time), type(type), client(client), job(nullptr) {
            if (prev_wait_time == 1) {
                wait_time = 2;
            } else if (prev_wait_time == 2) {
                wait_time = 4;
            } else if (prev_wait_time == 4) {
                wait_time = 8;
            } else if (prev_wait_time == 8) {
                wait_time = 15;
            } else if (prev_wait_time == 15) {
                wait_time = 15;
            } else {
                wait_time = 1;
            }
            //cout << "created an event with simulation time: " << this->time << endl;
            //cout << "created an event with event type: " << this->type << endl;
        }

        Event(double time, int type, size_t client, Job *job) : time(time), type(type), client(client), job(job) {
            wait_time = 1;
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

void parseFile(const std::string &file_name, std::vector<float> &means, std::vector<float> &stdevs, std::vector<float> &errors, std::vector<float> &err_means, std::vector<float> &err_stdevs) {
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
        VLOG(3) << "START_TIME: " << start_time;
        end_time.pop_back();
        VLOG(3) << "END_TIME: " << end_time;
        cpu_time.pop_back();
        VLOG(3) << "CPU_TIME: " << cpu_time;
        result_flag.pop_back();
        VLOG(3) << "RESULT_FLAG: " << result_flag;
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
    unsigned long err_sum[3] = {0, 0, 0};
    unsigned long err_size[3] = {0, 0, 0};
    for (size_t i = 0; i < wu_names.size(); i++) {
        long duration = end_times[i] - start_times[i];
        if (wu_names[i] == 1000) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                sum[0] += duration;
                size[0]++;
            } else {
                if (duration > 0) {
                    err_sum[0] += duration;
                }
                err_size[0]++;
            }
        } else if(wu_names[i] == 100) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                sum[1] += duration;
                size[1]++;
            } else {
                if (duration > 0) {
                    err_sum[1] += duration;
                }
                err_size[1]++;
            }
        } else if(wu_names[i] == 10) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                sum[2] += duration;
                size[2]++;
            } else {
                if (duration > 0) {
                    err_sum[2] += duration;
                }
                err_size[2]++;
            }
        }
    }

    for (size_t i = 0; i < 3; i++) {
        means.push_back(static_cast<float>(sum[i])/size[i]);
        err_means.push_back(static_cast<float>(err_sum[i])/err_size[i]);
        errors.push_back(static_cast<float>(err_size[i])/(err_size[i] + size[i]));
    }

    unsigned long stdev_sum[3] = {0, 0, 0};
    unsigned long err_stdev_sum[3] = {0, 0, 0};
    for (size_t i = 0; i < wu_names.size(); i++) {
        int duration = end_times[i] - start_times[i];
        if (wu_names[i] == 1000) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                stdev_sum[0] += (duration - means[0]) * (duration - means[0]);
            } else if (duration > 0) {
                err_stdev_sum[0] += (duration - err_means[0]) * (duration - err_means[0]);
            }
        } else if(wu_names[i] == 100) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                stdev_sum[1] += (duration - means[1]) * (duration - means[1]);
            } else if (duration > 0) {
                err_stdev_sum[1] += (duration - err_means[1]) * (duration - err_means[1]);
            }
        } else if(wu_names[i] == 10) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                stdev_sum[2] += (duration - means[2]) * (duration - means[2]);
            } else if (duration > 0) {
                err_stdev_sum[2] += (duration - err_means[2]) * (duration - err_means[2]);
            }
        }
    }

    for (size_t i = 0; i < 3; i++) {
        stdevs.push_back(sqrt(static_cast<float>(stdev_sum[i])/size[i]));
        err_stdevs.push_back(sqrt(static_cast<float>(err_stdev_sum[i])/err_size[i]));
    }
}

bool nextAvailableJob(const size_t client, std::queue<Job*> &available_jobs, Job *&job) {
    for (size_t i = 0; i < available_jobs.size(); i++) {
        job = available_jobs.front();
        if (job->clientOkay(client)) {
            available_jobs.pop();
            return true;
        } else {
            VLOG(2) << "Unfit client request.";
            return false;
        }
    }
    VLOG(3) << "No more available jobs.";
    return false;
}

void checkSample(std::vector<std::vector<std::vector<Job*>*>*> &job_queue, std::queue<Job*> &available_jobs, const size_t &walk_id, const size_t &sample_id) {
    bool addNextSample = true;
    std::vector<std::vector<Job*>*> *walk = job_queue.at(walk_id);
    std::vector<Job*> *sample = walk->at(sample_id);
    for (size_t i = 0; i < sample->size(); i++) {
        if (!sample->at(i)->done) {
            addNextSample = false;
            break;
        }
    }

    if (addNextSample && sample_id+1 < walk->size()) {
        VLOG(2) << "Make next sample available";
        sample = walk->at(sample_id+1);
        for (size_t i = 0; i < sample->size(); i++) {
            available_jobs.push(sample->at(i));
        }
    }
}

bool jobsDone(std::vector<std::vector<std::vector<Job*>*>*> &job_queue) {
    //TODO Fix this!
    if (jobs_complete < total_jobs) {
        return false;
    }
    return true;
    /*
    for (size_t i = 0; i < job_queue.size(); i++) {
        std::vector<std::vector<Job*>*> *walk = job_queue[i];
        for (size_t j = 0; j < walk->size(); j++) {
            std::vector<Job*> *sample = walk->at(j);
            for (size_t k = 0; k < sample->size(); k++) {
                if (sample->at(k)->done == false) {
                    return false;
                }
            }
        }
    }
    return true;
    */
}

void run_simulation(
        size_t num_walks,
        size_t samples_per_walk,
        size_t num_workers,
        size_t num_jobs_per_sample,
        size_t quorum,
        float error_percent,
        variate_generator< mt19937, lognormal_distribution<> > &duration_generator,
        variate_generator< mt19937, lognormal_distribution<> > &err_duration_generator,
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


    size_t complete_samples = 0;
    size_t complete_jobs_in_sample = 0;

    // This needs to be freed
    std::vector<std::vector<std::vector<Job*>*>*> job_queue;
    for (size_t i = 0; i < num_walks; i++) {
        job_queue.push_back(new std::vector<std::vector<Job*>*>());
        std::vector<std::vector<Job*>*> *walk = job_queue[i];
        for (size_t j = 0; j < samples_per_walk; j++) {
            walk->push_back(new std::vector<Job*>());
            std::vector<Job*> *sample = walk->at(j);
            for (size_t k = 0; k < num_jobs_per_sample; k++) {
                sample->push_back(new Job(i, j, k));
            }
        }
    }

    // For each walk make the first sample available.
    std::queue<Job*> available_jobs;
    for (size_t i = 0; i < num_walks; i++) {
        for (size_t k = 0; k < num_jobs_per_sample; k++) {
            available_jobs.push(job_queue.at(i)->at(0)->at(k));
        }
    }

    // Need another array to track which clients failed and on which samples.


    // Put initial events in the heap
    for (int client_id = 0; client_id < num_workers; client_id++) {
        heap.push(new Event(simulation_time_s + client_id, REQUEST_JOB, client_id, 1));
    }

    VLOG(2) << "Start Simulation...";

    while (!heap.empty() && !jobsDone(job_queue)) {
        Event *current_event = heap.top();
        size_t current_client = current_event->client;
        Job *current_job = current_event->job;
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

        Job *next_job = nullptr;

        switch (current_event->type) {
            case 0: // REQUESET_JOB
                // If there is a job push a new event for either success or
                // error. If there is no job in the queue then push another job
                // request onto the queue.
                VLOG(2) << "Request Job Begin (" << current_client << ")";
                if (nextAvailableJob(current_client, available_jobs, next_job)) {
                    assert(next_job != nullptr);
                    // Add another request event
                    if (error_generator() <= error_percent) {
                        heap.push(new Event(simulation_time_s + err_duration_generator(), RETURN_ERROR, current_client, next_job));
                    } else {
                        heap.push(new Event(simulation_time_s + duration_generator(), RETURN_SUCCESS, current_client, next_job));
                    }
                } else {
                    heap.push(new Event(simulation_time_s + current_event->wait_time, REQUEST_JOB, current_client, current_event->wait_time));
                }
                VLOG(2) << "Request Job End";
                break;
            case 1: // RETURN_SUCCESS
                // Server checks if there the quarum for this sample is
                // complete if it is then create jobs for the next sample in
                // the walk. If there are jobs available create a new success
                // or error event otherwise a new request event.
                VLOG(2) << "Success Begin";
                VLOG(2) << *current_job << " was successful.";
                current_job->successClient(current_event->client);
                checkSample(job_queue, available_jobs, current_job->walk, current_job->sample);
                if (nextAvailableJob(current_client, available_jobs, next_job)) {
                    // Add another request event
                    if (error_generator() <= error_percent) {
                        heap.push(new Event(simulation_time_s + err_duration_generator(), RETURN_ERROR, current_client, next_job));
                    } else {
                        heap.push(new Event(simulation_time_s + duration_generator(), RETURN_SUCCESS, current_client, next_job));
                    }
                } else {
                    heap.push(new Event(simulation_time_s + current_event->wait_time, REQUEST_JOB, current_client, current_event->wait_time));
                }
                VLOG(2) << "Success End";
                break;
            case 2: // RETURN_ERROR
                // Server creates a new job for the errored sample. If there
                // are available jobs create a success or error event. Otherwise create a
                // request event.
                VLOG(2) << "Error Begin";
                current_job->failClient(current_event->client);
                available_jobs.push(current_job);
                VLOG(2) << "Get next available job";
                if (nextAvailableJob(current_event->client, available_jobs, next_job)) {
                    // Add another request event
                    if (error_generator() <= error_percent) {
                        heap.push(new Event(simulation_time_s + err_duration_generator(), RETURN_ERROR, current_client, next_job));
                    } else {
                        heap.push(new Event(simulation_time_s + duration_generator(), RETURN_SUCCESS, current_client, next_job));
                    }
                } else {
                    heap.push(new Event(simulation_time_s + current_event->wait_time, REQUEST_JOB, current_client, current_event->wait_time));
                }
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

    for (size_t i = 0; i < num_workers; i++) {
        Event *event = heap.top();
        heap.pop();
        delete event;
    }
    assert(heap.empty());

    // This needs to be freed
    for (size_t i = 0; i < num_walks; i++) {
        std::vector<std::vector<Job*>*> *walk = job_queue[i];
        for (size_t j = 0; j < samples_per_walk; j++) {
            std::vector<Job*> *sample = walk->at(j);
            for (size_t k = 0; k < num_jobs_per_sample; k++) {
                delete sample->at(k);
            }
            delete sample;
        }
        delete walk;
    }

    VLOG(1) << "The simulation ended at time: " << simulation_time_s;
}

int main(int argc, char **argv) {
    // Initialize Google Logging
    google::InitGoogleLogging(argv[0]);
    // Log to Stderr
    FLAGS_logtostderr = 1;

    std::vector<float> means;
    std::vector<float> stdevs;
    std::vector<float> err_means;
    std::vector<float> err_stdevs;
    std::vector<float> errors;
    parseFile("../data/dna_workunit_transit.csv", means, stdevs, errors, err_means, err_stdevs);

    int seed = time(0);
    variate_generator< mt19937, std::uniform_real_distribution<> > error_generator(mt19937(seed), std::uniform_real_distribution<>(0, 1));
    variate_generator< mt19937, lognormal_distribution<> > duration_generator(mt19937(seed+1), lognormal_distribution<>(means[0], stdevs[0]));
    variate_generator< mt19937, lognormal_distribution<> > err_duration_generator(mt19937(seed+2), lognormal_distribution<>(err_means[0], err_stdevs[0]));

    for (size_t i = 0; i < 3; i++) {
        LOG(INFO) << i << ":Mean:" << means[i];
        LOG(INFO) << i << ":Stdevs:" << stdevs[i];
        LOG(INFO) << i << ":ErrMean:" << err_means[i];
        LOG(INFO) << i << ":ErrStdevs:" << err_stdevs[i];
        LOG(INFO) << i << ":Errors:" << errors[i];
    }

    std::ofstream outfile("duration.dat");
    for (size_t i = 0; i < 9317; i++) {
        outfile << duration_generator() << std::endl;
    }
    outfile.close();

    // Intiate
    size_t num_walks = 1;
    size_t samples_per_walk = 10;
    size_t num_workers = 2;
    size_t num_jobs_per_sample = 2;
    size_t quorum = 2;

    if (argc >= 2) {
        samples_per_walk = atoi(argv[1]);
    }
    if (argc >= 3) {
        num_workers = atoi(argv[2]);
    }
    if (argc >= 4) {
        num_walks = atoi(argv[3]);
    }
    total_jobs = num_walks * samples_per_walk * num_jobs_per_sample;

    LOG(INFO) << "Number of Walks: " << num_walks;
    LOG(INFO) << "Walk size: " << samples_per_walk;
    LOG(INFO) << "Number of Workers: " << num_workers;
    LOG(INFO) << "Number of simultaneous jobs: " << num_jobs_per_sample;

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
                //servers_idle.uush_back(true);
                //servers_queue.push_back(new std::queue<double>());
            }
            run_simulation(num_walks, samples_per_walk, num_workers, num_jobs_per_sample, quorum, errors[i], duration_generator, err_duration_generator, error_generator);

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
