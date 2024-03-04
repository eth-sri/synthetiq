#include "randomhelper.h"

/**
 * @brief Constructs a RandomHelper object.
 * 
 * This constructor initializes the random number generator and the distribution for generating random numbers between 0 and 1.
 * It also initializes an empty map for integer distributions.
 */
RandomHelper::RandomHelper() 
{
    mt = std::mt19937(rd());
    distribution01 = std::uniform_real_distribution<double>(0.0, 1.0);
    int_ditributions = std::map<int, std::uniform_int_distribution<int>>();
}

/**
 * Generates a random number between 0 and 1.
 *
 * @return The generated random number.
 */
double RandomHelper::random01(){
    return distribution01(mt);
}

/**
 * Generates a random number within a specified range.
 *
 * @param min The minimum value of the range.
 * @param max The maximum value of the range.
 * @return A random number within the specified range.
 */
double RandomHelper::random_uniform(double min, double max) 
{ 
    return random01() * (max - min) + min;
}

/**
 * Generates a random integer between 0 and the specified maximum value (exclusive).
 * If a distribution for the specified maximum value already exists, it uses that distribution.
 * Otherwise, it creates a new uniform_int_distribution for the maximum value and stores it for future use.
 *
 * @param max The maximum value (exclusive) for the random integer.
 * @return A random integer between 0 and max-1.
 */
int RandomHelper::randomInt(int max) 
{
    auto distrib = int_ditributions.find(max);
    if (distrib != int_ditributions.end()){
        return distrib->second(mt);
    }
    std::uniform_int_distribution<int> distribution(0, max - 1);
    int_ditributions[max] = distribution;
    return distribution(mt);
}

/**
 * Sets the seed for the random number generator.
 * 
 * @param seed The seed value to initialize the random number generator.
 */
void RandomHelper::seed(int seed) {
    mt.seed(seed);
}