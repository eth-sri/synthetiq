#include "temperatureScheme.h"
#include <iostream>

/**
 * @brief Copy constructor for the TemperatureScheme class.
 * 
 * @param other The TemperatureScheme object to be copied.
 */
TemperatureScheme::TemperatureScheme(TemperatureScheme const& other) {
    temperature = other.temperature;
    start_temperature = other.start_temperature;
}

/**
 * @brief Get the temperature value.
 * 
 * @return The temperature value.
 */
double TemperatureScheme::getTemperature() {
    return temperature;
}

/**
 * @brief Resets the temperature to its initial value.
 */
void TemperatureScheme::reset() {
    temperature = start_temperature;
}

/**
 * @brief Constructs an ExponentialTemperatureScheme object.
 * 
 * @param temperature The initial temperature value.
 * @param n_steps_normalizer The normalizing factor for the number of steps.
 */
ExponentialTemperatureScheme::ExponentialTemperatureScheme(double temperature, double n_steps_normalizer) : 
        TemperatureScheme(temperature), n_steps_normalizer(n_steps_normalizer) {

}

/**
 * @brief Copy constructor for ExponentialTemperatureScheme.
 * 
 * @param other The ExponentialTemperatureScheme object to be copied.
 */
ExponentialTemperatureScheme::ExponentialTemperatureScheme(ExponentialTemperatureScheme const& other) : TemperatureScheme(other) {
    n_steps_normalizer = other.n_steps_normalizer;
}

/**
 * @brief Returns a shared pointer to a clone of the TemperatureScheme object.
 * 
 * @return std::shared_ptr<TemperatureScheme> A shared pointer to a clone of the TemperatureScheme object.
 */
std::shared_ptr<TemperatureScheme> ExponentialTemperatureScheme::clone() {
    return std::make_shared<ExponentialTemperatureScheme>(*this);
}

/**
 * Updates the temperature based on the current step, number of accepted moves, number of rejected moves,
 * number of elements in the circuit, and number of qubits.
 *
 * @param cur_step The current step in the optimization process.
 * @param n_accepted The number of accepted moves so far.
 * @param cur_rejected The number of rejected moves in the current step.
 * @param nbElementsCircuit The number of elements in the circuit.
 * @param nb_qbs The number of qubits.
 */
void ExponentialTemperatureScheme::updateTemperature(int cur_step, int n_accepted, int nbElementsCircuit) {
    temperature = start_temperature * std::exp(-n_accepted / (nbElementsCircuit * n_steps_normalizer));
}