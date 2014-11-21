#include <fstream>
#include <string>
#include <limits>
#include <random>
#include <algorithm>
#include <queue>
#include <vector>

#include <glog/logging.h>

// Stats vectors
std::vector<double> estimated_queue_lengths;

// Output File
std::ofstream outfile;

struct Packet {
    unsigned int type;
};

double run_simulation(double lambda, std::vector<double> mu_vec, size_t n) {
    VLOG(2) << "Start Simulation...";

    // Init RNG
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(0,1);

    double total_queue_lengths = 0;
    size_t arrivals = 0;
    std::queue<size_t> packet_queue;
    // Simulation steps
    while (arrivals < n) {
        double arrival = distribution(generator);
        double service = distribution(generator);

        // Arrival
        if (arrival <= lambda) {
            VLOG(2) << "Packet arrived";
            arrivals++;
            packet_queue.push(arrivals % mu_vec.size());
        }

        // Service
        if (!packet_queue.empty()) {
            if (service <= mu_vec.at(packet_queue.front())) {
                packet_queue.pop();
            }
        }
        total_queue_lengths += packet_queue.size();
    }

    double average_queue_length = total_queue_lengths / n;
    VLOG(1) << "The simulation ended with average queue length: " << average_queue_length;
    return average_queue_length;
}

int main(int argc, char **argv) {
    //Initialize Google Logging
    google::InitGoogleLogging(argv[0]);
    //Log to Stderr
    FLAGS_logtostderr = 1;

    // Intiate
    double lambda = 1.0;
    std::vector<double> mu = {0.5, 0.25};
    size_t n = 1000;
    size_t iterations = 1000;

    if (argc >= 2) {
        n = atoi(argv[1]);
    }
    if (argc >= 3) {
        lambda = atof(argv[2]);
    }
    if (argc >= 4) {
        mu.clear();
        for (int i = 3; i < argc; i++) {
            mu.push_back(atof(argv[i]));
        }
    }

    LOG(INFO) << "Lambda: " << lambda;
    for (int i = 0; i < mu.size(); i ++) {
        LOG(INFO) << "Mu: " << mu.at(i);
    }
    LOG(INFO) << "N: " << n;
    LOG(INFO) << "Iterations: " << iterations;

    // Open files
    std::stringstream filename;
    filename << "fifo_stats" << ".dat";
    outfile.open(filename.str());

    for (int i = 0; i < iterations; i++) {
        double estimated_queue_length = run_simulation(lambda, mu, n);
        estimated_queue_lengths.push_back(estimated_queue_length);
    }

    // Collect and print statistics.

    double sum_time = std::accumulate(estimated_queue_lengths.begin(), estimated_queue_lengths.end(), 0.0);
    double min_time = *std::min_element(estimated_queue_lengths.begin(), estimated_queue_lengths.end());
    double max_time = *std::max_element(estimated_queue_lengths.begin(), estimated_queue_lengths.end());
    double avg_time = sum_time / estimated_queue_lengths.size();
    double sqr_time = std::inner_product(estimated_queue_lengths.begin(), estimated_queue_lengths.end(), estimated_queue_lengths.begin(), 0.0);
    double std_time = std::sqrt(sqr_time/ estimated_queue_lengths.size() - avg_time* avg_time);

    outfile << "------------------------------------------" << std::endl;
    outfile << "-------------- QUEUE LENGTHS -------------" << std::endl;
    outfile << "------------------------------------------" << std::endl;
    outfile << "Min: " << std::fixed <<  min_time << std::endl;
    outfile << "Max: " << std::fixed << max_time << std::endl;
    outfile << "Mean: " << std::fixed << avg_time << std::endl;
    outfile << "Stdev: " << std::fixed << std_time << std::endl;
    outfile << std::endl;

    outfile.close();
}

