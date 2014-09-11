#include <iostream>
#include <vector>
#include <string>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>

using boost::variate_generator;
using boost::mt19937;
using boost::exponential_distribution;

int main(int argc, char **argv) {
    int seed = time(0);
    double distribution_mean = 3.0;
    variate_generator< mt19937, exponential_distribution<> > rand_generator(mt19937(seed), exponential_distribution<>(1/distribution_mean));
    std::cout << rand_generator() << std::endl;
}
