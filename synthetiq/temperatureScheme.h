#ifndef DEF_TEMP_SCHEME
#define DEF_TEMP_SCHEME
#include <vector>
#include <Eigen/Dense>
#include <memory>

class TemperatureScheme {
    public:
        TemperatureScheme() {};
        TemperatureScheme(TemperatureScheme const& other);
        TemperatureScheme(double temperature) : temperature(temperature), start_temperature(temperature) {};
        virtual void updateTemperature(int cur_step, int n_accepted, int nbElementsCircuit) = 0;
        virtual std::shared_ptr<TemperatureScheme> clone() = 0;
        void reset();
        double getTemperature();
    protected:
        double temperature;
        double start_temperature;

};

class ExponentialTemperatureScheme : public TemperatureScheme {
    public:
        ExponentialTemperatureScheme(double temperature=0.04, double n_steps_normalizer=100.0);
        ExponentialTemperatureScheme(ExponentialTemperatureScheme const& other);
        std::shared_ptr<TemperatureScheme> clone();
        void reset();
        void updateTemperature(int cur_step, int n_accepted, int nbElementsCircuit);
    private:
        double n_steps_normalizer;
};

#endif