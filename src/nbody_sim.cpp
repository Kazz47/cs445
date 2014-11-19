#include <iostream>
#include <fstream>
#include <queue>
#include <string>
#include <limits>
#include <cmath>

#include <glog/logging.h>

#define ONLY_EVENT_TYPE     0
#define NUMBER_EVENT_TYPES  2
#define NUMBER_TOPOLOGY_TYPES  2

const static unsigned int COMPUTE = 0;
const static unsigned int SEND_PACKET = 1;
const std::string event_names[NUMBER_EVENT_TYPES] = {
    "COMPUTE",
    "SEND_PACKET"
};

const static unsigned int RING = 0;
const static unsigned int GRID = 1;
const static unsigned int CUBE = 2;
const std::string topology_names[NUMBER_TOPOLOGY_TYPES] = {
    "RING",
    "GRID",
    "HYPERCUBE"
};

// Stats vectors
std::vector<double> simulation_times;

// Output File
std::ofstream outfile;

class Packet {
    public:
        size_t x;
        size_t y;
        size_t z;
        size_t size;
        size_t iteration;

        Packet(size_t x, size_t y, size_t z, size_t iteration) : x(x), y(y), z(z), iteration(iteration) {}

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

double run_simulation(size_t nodes, size_t data, double compute_time, double latency, char topology) {
    double simulation_time_s = 0;
    double previous_time_s = 0;

    std::priority_queue<Event*, std::vector<Event*>, CompareEvent> heap;

    // Put initial events in the heap and set worker error delays
    for (int node_id = 0; node_id < nodes; node_id++) {
        heap.push(new Event(simulation_time_s + node_id, COMPUTE));
    }

    // Put initial validator events in heap
    heap.push(new Event(simulation_time_s + 10, CHECK_VALID));
    heap.push(new Event(simulation_time_s + 86400, CHECK_WORKERS));

    // Reset Globals
    jobs_complete = 0;

    VLOG(2) << "Start Simulation...";

    while (!heap.empty() && !jobsDone()) {
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
                        //heap.push(new Event(simulation_time_s + generator(error_generator, data), RETURN_SUCCESS, current_client, next_job));
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
                current_job->successClient(current_event->client, job_queue, quorum);
                success_jobs.push(current_job);
                if (!worker_bans.at(current_client) && nextAvailableJob(current_client, available_jobs, next_job)) {
                    // Add another request event
                    if (error_generator() <= error_percent) {
                        heap.push(new Event(simulation_time_s + err_duration_generator(), RETURN_ERROR, current_client, next_job));
                    } else {
                        heap.push(new Event(simulation_time_s + duration_generator(), RETURN_SUCCESS, current_client, next_job));
                        //heap.push(new Event(simulation_time_s + generator(error_generator, data), RETURN_SUCCESS, current_client, next_job));
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
                        //heap.push(new Event(simulation_time_s + generator(error_generator, data), RETURN_SUCCESS, current_client, next_job));
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
    total_jobs = num_walks * samples_per_walk;

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
            double simulation_time = run_simulation(num_walks, samples_per_walk, num_workers, num_jobs_per_sample, quorum, errors[i], data[i], duration_generators[i], err_duration_generators[i], error_generator);
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
        outfile << "Min: " << std::fixed <<  min_time << std::endl;
        outfile << "Max: " << std::fixed << max_time << std::endl;
        outfile << "Mean: " << std::fixed << avg_time << std::endl;
        outfile << "Stdev: " << std::fixed << std_time << std::endl;
        outfile << std::endl;

    }
    outfile.close();
}
