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
#define NUMBER_EVENT_TYPES  5

const static unsigned int REQUEST_JOB = 0;
const static unsigned int RETURN_SUCCESS = 1;
const static unsigned int RETURN_ERROR = 2;
const static unsigned int CHECK_VALID = 3;
const static unsigned int CHECK_WORKERS = 4;
const std::string event_names[NUMBER_EVENT_TYPES] = {
    "REQUEST_JOB",
    "RETURN_SUCCESS",
    "RETURN_ERROR",
    "CHECK_VALID",
    "CHECK_WORKERS"
};

std::vector<bool> job_queue;
std::vector< std::queue<double>* > servers_queue;

size_t total_jobs = 0;
size_t jobs_complete = 0;
// Stats vectors
std::vector<double> simulation_times;

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
            //failed_clients.push_back(client);
            done = true;
            jobs_complete++;
        }

        void failClient(size_t client) {
            failed_clients.push_back(client);
        }

        bool clientOkay(size_t client) {
            for (size_t i = 0; i < failed_clients.size(); i++) {
                if (client == failed_clients.at(i)) {
                    VLOG(3) << "Client not okay: (" << client << ")";
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

        Event(double time, int type) : time(time), type(type), client(-1), job(nullptr) {
            wait_time = 0;
            //cout << "created an event with simulation time: " << this->time << endl;
            //cout << "created an event with event type: " << this->type << endl;
        }

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

void parseFile(
        const std::string &file_name,
        std::vector<double> &means, std::vector<double> &stdevs,
        std::vector<double> &errors, std::vector<double> &err_means,
        std::vector<double> &err_stdevs,
        std::vector<std::vector<double>> &data) {
    LOG(INFO) << "Load file: '" << file_name << "'";
    std::ifstream infile(file_name);

    std::vector<unsigned long> diff_times, result_flags, wu_names;
    std::vector<double> cpu_times;

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
        diff_times.push_back(atoi(end_time.c_str()) - atoi(start_time.c_str()));
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

    //Adjust data size to three arrays
    data.resize(3);

    unsigned long sum[3] = {0, 0, 0};
    unsigned long size[3] = {0, 0, 0};
    unsigned long err_sum[3] = {0, 0, 0};
    unsigned long err_size[3] = {0, 0, 0};
    for (size_t i = 0; i < wu_names.size(); i++) {
        long duration = diff_times[i];
        if (wu_names[i] == 1000) {
            if (result_flags[i] == 1) {
                LOG_IF(ERROR, duration <= 0) << "Invalid duration: " << duration;
                sum[0] += duration;
                data[0].push_back(duration);
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
                data[1].push_back(duration);
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
                data[2].push_back(duration);
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
        means.push_back(static_cast<double>(sum[i])/size[i]);
        err_means.push_back(static_cast<double>(err_sum[i])/err_size[i]);
        errors.push_back(static_cast<double>(err_size[i])/(err_size[i] + size[i]));
    }

    unsigned long stdev_sum[3] = {0, 0, 0};
    unsigned long err_stdev_sum[3] = {0, 0, 0};
    for (size_t i = 0; i < wu_names.size(); i++) {
        int duration = diff_times[i];
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
        stdevs.push_back(sqrt(static_cast<double>(stdev_sum[i])/size[i]));
        err_stdevs.push_back(sqrt(static_cast<double>(err_stdev_sum[i])/err_size[i]));
    }

    for (size_t i = 0; i < data.size(); i++) {
        std::sort(data.at(i).begin(), data.at(i).end());
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
}

// Input data must be sorted and distribution must be from 0 to 1.
double generator(variate_generator< mt19937, std::uniform_real_distribution<> > &gen, std::vector<double> &data) {
    double val = gen() * (data.size()-1);
    size_t low = data.at(floor(val));
    size_t high = data.at(ceil(val));
    return gen() * (high - low) + low;
}

double run_simulation(
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

    std::priority_queue<Event*, std::vector<Event*>, CompareEvent> heap;
    std::queue<Job*> success_jobs;
    std::queue<Job*> failure_jobs;

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

    // Put initial events in the heap and set worker error delays
    std::vector<size_t> worker_errors;
    std::vector<bool> worker_bans;
    for (int client_id = 0; client_id < num_workers; client_id++) {
        worker_bans.push_back(false);
        worker_errors.push_back(0);
        heap.push(new Event(simulation_time_s + client_id, REQUEST_JOB, client_id, 1));
    }

    // Put initial validator events in heap
    heap.push(new Event(simulation_time_s + 10, CHECK_VALID));
    heap.push(new Event(simulation_time_s + 86400, CHECK_WORKERS));

    // Reset Globals
    jobs_complete = 0;

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
                if (!worker_bans.at(current_client) && nextAvailableJob(current_client, available_jobs, next_job)) {
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
                success_jobs.push(current_job);
                if (!worker_bans.at(current_client) && nextAvailableJob(current_client, available_jobs, next_job)) {
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
                worker_errors.at(current_event->client)++;
                failure_jobs.push(current_job);
                VLOG(2) << "Get next available job";
                if (!worker_bans.at(current_client) && nextAvailableJob(current_event->client, available_jobs, next_job)) {
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
            case 3: // CHECK_VALID
                // Server runs the validator and create new jobs where
                // necessary.
                VLOG(3) << "Validator Begin";
                while (!success_jobs.empty()) {
                    Job *temp_job = success_jobs.front();
                    success_jobs.pop();
                    checkSample(job_queue, available_jobs, temp_job->walk, temp_job->sample);
                }
                while (!failure_jobs.empty()) {
                    Job *temp_job = failure_jobs.front();
                    failure_jobs.pop();
                    available_jobs.push(temp_job);
                }
                heap.push(new Event(simulation_time_s + 10, CHECK_VALID));
                VLOG(3) << "Validator End";
                break;
            case 4: // CHECK_WORKERS
                // Server checks how many errors each client has made each day.
                // If the client's errors exceed 10 for that day then ban them
                // from work for the next day.
                VLOG(2) << "Worker Check Begin";
                for (size_t worker_id = 0; worker_id < num_workers; worker_id++) {
                    if (worker_errors.at(worker_id) >= 10) {
                        worker_bans.at(worker_id) = true;
                    } else {
                        worker_bans.at(worker_id) = false;
                    }
                    worker_errors.at(worker_id) = 0;
                }
                heap.push(new Event(simulation_time_s + 86400, CHECK_WORKERS));
                VLOG(2) << "Worker Check End";
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
    heap.pop(); // One for the validator
    heap.pop(); // One for the worker check
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
    return simulation_time_s;
}

int main(int argc, char **argv) {
    // Initialize Google Logging google::InitGoogleLogging(argv[0]);
    // Log to Stderr
    FLAGS_logtostderr = 1;

    std::vector<double> means;
    std::vector<double> stdevs;
    std::vector<double> err_means;
    std::vector<double> err_stdevs;
    std::vector<double> errors;
    std::vector<std::vector<double>> data;
    parseFile("../data/dna_workunit_transit.csv", means, stdevs, errors, err_means, err_stdevs, data);

    int seed = time(0);
    variate_generator< mt19937, std::uniform_real_distribution<> > error_generator(mt19937(seed), std::uniform_real_distribution<>(0, 1));
    variate_generator< mt19937, lognormal_distribution<> > duration_generator_1000(mt19937(seed+2), lognormal_distribution<>(means[1], stdevs[1]));
    variate_generator< mt19937, lognormal_distribution<> > duration_generator_100(mt19937(seed+2), lognormal_distribution<>(means[1], stdevs[1]));
    variate_generator< mt19937, lognormal_distribution<> > duration_generator_10(mt19937(seed+3), lognormal_distribution<>(means[2], stdevs[2]));
    std::vector<variate_generator< mt19937, lognormal_distribution<> >> duration_generators;
    duration_generators.push_back(duration_generator_1000);
    duration_generators.push_back(duration_generator_100);
    duration_generators.push_back(duration_generator_10);

    variate_generator< mt19937, lognormal_distribution<> > err_duration_generator_1000(mt19937(seed+4), lognormal_distribution<>(err_means[0], err_stdevs[0]));
    variate_generator< mt19937, lognormal_distribution<> > err_duration_generator_100(mt19937(seed+5), lognormal_distribution<>(err_means[1], err_stdevs[1]));
    variate_generator< mt19937, lognormal_distribution<> > err_duration_generator_10(mt19937(seed+6), lognormal_distribution<>(err_means[2], err_stdevs[2]));
    std::vector<variate_generator< mt19937, lognormal_distribution<> >> err_duration_generators;
    err_duration_generators.push_back(err_duration_generator_1000);
    err_duration_generators.push_back(err_duration_generator_100);
    err_duration_generators.push_back(err_duration_generator_10);

    for (size_t i = 0; i < 3; i++) {
        LOG(INFO) << i << ":Mean:" << means[i];
        LOG(INFO) << i << ":Stdevs:" << stdevs[i];
        LOG(INFO) << i << ":ErrMean:" << err_means[i];
        LOG(INFO) << i << ":ErrStdevs:" << err_stdevs[i];
        LOG(INFO) << i << ":Errors:" << errors[i];
    }

    // Write distribution sample to files
    /*
    std::ofstream log_outfile("log_duration.dat");
    std::ofstream emp_outfile("emp_duration.dat");
    for (size_t i = 0; i < 9317; i++) {
    //for (size_t i = 0; i < 13945; i++) {
    //for (size_t i = 0; i < 10670; i++) {
        log_outfile << duration_generator() << std::endl;
        emp_outfile << generator(error_generator, data[0]) << std::endl;
    }
    log_outfile.close();
    emp_outfile.close();
    */

    // Intiate
    size_t iterations = 100;
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
    filename << "dna_stats" << ".dat";
    outfile.open(filename.str());

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < iterations; j++) {
            double simulation_time = run_simulation(num_walks, samples_per_walk, num_workers, num_jobs_per_sample, quorum, errors[i], duration_generators[i], err_duration_generators[i], error_generator);
            simulation_times.push_back(simulation_time);
        }

        // Collect and print statistics.

        double sum_time = std::accumulate(simulation_times.begin(), simulation_times.end(), 0.0);
        double min_time = *std::min_element(simulation_times.begin(), simulation_times.end());
        double max_time = *std::max_element(simulation_times.begin(), simulation_times.end());
        double avg_time = sum_time / simulation_times.size();
        double sqr_time = std::inner_product(simulation_times.begin(), simulation_times.end(), simulation_times.begin(), 0.0);
        double std_time = std::sqrt(sqr_time/ simulation_times.size() - avg_time* avg_time);

        outfile << "------------------------------------------" << std::endl;
        outfile << "-------------------- " << i << " -------------------" << std::endl;
        outfile << "------------------------------------------" << std::endl;
        outfile << std::endl;

        outfile << "------------------------------------------" << std::endl;
        outfile << "---------------- RUN TIMES ---------------" << std::endl;
        outfile << "------------------------------------------" << std::endl;
        outfile << "Min run time: " << std::fixed <<  min_time << std::endl;
        outfile << "Max run time: " << std::fixed << max_time << std::endl;
        outfile << "Mean run time: " << std::fixed << avg_time << std::endl;
        outfile << "Standard deviation of run time: " << std::fixed << std_time << std::endl;
        outfile << std::endl;

    }
    outfile.close();
}
