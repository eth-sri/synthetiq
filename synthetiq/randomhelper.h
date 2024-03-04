#ifndef DEF_RDMH
#define DEF_RDMH

#include <random>
#include <map>

class RandomHelper{
    public:
        RandomHelper();
        double random01();
        double random_uniform(double min, double max);
        int randomInt(int max); // [0, max[
        void seed(int seed);
        
    private:
        std::random_device rd;
        std::mt19937 mt;
        std::uniform_real_distribution<double> distribution01;
        std::map<int, std::uniform_int_distribution<int>> int_ditributions;
};

#endif