#include "gateScheme.h"
#include <iostream>

/**
 * @brief Constructs a GateScheme object with the specified parameters. This updates the circuit size bounds in the SA run
 * 
 * @param min_start_gates The minimum number of starting gates.
 * @param max_start_gates The maximum number of starting gates.
 * @param beta The beta value.
 * @param min_factor The minimum factor. (f_min)
 * @param max_factor The maximum factor. (f_max)
 */
GateScheme::GateScheme(int min_start_gates, int max_start_gates, double beta, double min_factor, double max_factor) : min_start_gates(min_start_gates), 
    max_start_gates(max_start_gates), beta(beta), min_factor(min_factor), max_factor(max_factor) {
        min_gates = min_start_gates;
        max_gates = max_start_gates;
        current_best_gates = max_start_gates;
        start_best_gates = max_start_gates;
}

/**
 * @brief Resets the gate scheme to its initial state.
 * 
 * This function sets the minimum and maximum number of gates to their initial values,
 * and resets the current best gates to the starting best gates.
 */
void GateScheme::reset() {
    min_gates = min_start_gates;
    max_gates = max_start_gates;
    current_best_gates = start_best_gates;
}

/**
 * @brief Sets the number of start best gates.
 * 
 * This function sets the number of start best gates for the GateScheme object.
 * 
 * @param gates The number of start best gates to set.
 */
void GateScheme::setStartBestGates(int gates) {
    start_best_gates = gates;
}

/**
 * Updates the gate scheme based on the number of gates found.
 * The gate scheme is adjusted to ensure that the number of gates falls within the specified range.
 *
 * @param gate_count_found The number of gates found.
 */
void GateScheme::update(int gate_count_found) {
    current_best_gates = std::min(current_best_gates, gate_count_found);

    if (max_factor * current_best_gates > max_gates) {
        int size_difference = (int) (max_factor * current_best_gates - max_gates);
        max_gates += std::max(1, (int) (beta * size_difference));
    } else {
        max_gates = (int) (max_factor * current_best_gates);
    }

    if (min_factor * current_best_gates > min_gates) {
        int size_difference = (int) (min_factor * current_best_gates - min_gates);
        min_gates += std::max(1, (int) (beta * size_difference));
    } else {
        min_gates = (int) (min_factor * current_best_gates);
    }
}

/**
 * Returns a random number of start gates within the specified range.
 *
 * @param rh The RandomHelper object used to generate random numbers.
 * @return The number of start gates.
 */
int GateScheme::getStartGates(RandomHelper& rh) {
    if (max_gates == 0) {
        return 1;
    }
    return rh.randomInt(max_gates - min_gates) + min_gates + 1;
}

/**
 * @brief Sets the value of beta for the GateScheme.
 * 
 * @param new_beta The new value of beta.
 */
void GateScheme::setBeta(double new_beta) {
    beta = new_beta;
}

/**
 * @brief Get the value of beta.
 * 
 * @return The value of beta.
 */
double GateScheme::getBeta() {
    return beta;
}

/**
 * Sets the maximum factor for the GateScheme.
 *
 * @param new_max_factor The new maximum factor to set.
 */
void GateScheme::setMaxFactor(double new_max_factor) {
    max_factor = new_max_factor;
}

/**
 * @brief Sets the minimum factor for the GateScheme.
 * 
 * This function sets the minimum factor for the GateScheme object.
 * The minimum factor determines the lower bound for the factor used in the GateScheme.
 * 
 * @param new_min_factor The new minimum factor value to be set.
 */
void GateScheme::setMinFactor(double new_min_factor) {
    min_factor = new_min_factor;
}

/**
 * @brief Sets the maximum number of start gates for the GateScheme.
 * 
 * This function updates the value of max_start_gates and resets the GateScheme.
 * 
 * @param new_max_start_gates The new maximum number of start gates.
 */
void GateScheme::setMaxStartGates(int new_max_start_gates) {
    max_start_gates = new_max_start_gates;
    reset();
}

/**
 * @brief Sets the minimum number of start gates for the GateScheme.
 * 
 * This function updates the value of min_start_gates and resets the GateScheme.
 * 
 * @param new_min_start_gates The new minimum number of start gates.
 */
void GateScheme::setMinStartGates(int new_min_start_gates) {
    min_start_gates = new_min_start_gates;
    reset();
}

/**
 * Returns the minimum factor of the gate scheme.
 *
 * @return The minimum factor of the gate scheme.
 */
double GateScheme::getMinFactor() {
    return min_factor;
}

/**
 * @brief Get the maximum factor of the gate scheme.
 * 
 * @return The maximum factor.
 */
double GateScheme::getMaxFactor() {
    return max_factor;
}

/**
 * Returns the minimum number of start gates for the GateScheme.
 *
 * @return The minimum number of start gates.
 */
int GateScheme::getMinStartGates() {
    return min_start_gates;
}

/**
 * Returns the maximum number of start gates allowed in the gate scheme.
 *
 * @return The maximum number of start gates.
 */
int GateScheme::getMaxStartGates() {
    return max_start_gates;
}

/**
 * @brief Get the minimum number of gates in the gate scheme.
 * 
 * @return The minimum number of gates.
 */
int GateScheme::getMinGates() {
    return min_gates;
}

/**
 * @brief Get the maximum number of gates allowed in the gate scheme.
 * 
 * @return The maximum number of gates.
 */
int GateScheme::getMaxGates() {
    return max_gates;
}

