#include "mcmc.h"
#include <math.h>
#include <sstream>
#include <omp.h>
#include <filesystem>

/**
 * @brief Default constructor for MCMCResult class.
 * Initializes the member variables with default values.
 */
MCMCResult::MCMCResult() {
    best_energy = 100000;
    best_eq = 0; 
}


/**
 * @brief Constructs a new MCMCResult object.
 * 
 * @param best A shared pointer to the best GateCircuit.
 */
MCMCResult::MCMCResult(std::shared_ptr<GateCircuit> best) {
    best_energy = 100000;
    best_eq = 0;
    circuit_best = best;
}
 
/**
 * Returns a string representation of the MCMCResult object.
 * The string includes information about the step reached, total steps, best circuit,
 * energy, equality cost, performance, number of no changes in matrix, number of no changes in circuit,
 * number of rejections, and number of no changes in probability.
 * 
 * @return A string representation of the MCMCResult object.
 */
std::string MCMCResult::print(){
    std::ostringstream oss;
    oss << circuit_best->print_qasm();
    return oss.str();
}

/**
 * @brief Constructs an instance of the MCMC class.
 * 
 * @param random_helper The RandomHelper object used for random number generation.
 * @param proba_id The probability of performing an identity operation.
 * @param enable_permutations Flag indicating whether permutations are enabled.
 */
MCMC::MCMC(RandomHelper& random_helper, double proba_id, bool enable_permutations):  random_helper(random_helper), enable_permutations(enable_permutations) {   
    equalitycomp = std::make_shared<FroebeniusCostComputer>(FroebeniusCostComputer());
    mutator = std::make_shared<Mutator>(Mutator(proba_id));
    temp_scheme = std::make_shared<ExponentialTemperatureScheme>(ExponentialTemperatureScheme());
}

/**
 * @brief Sets the equality computer for the MCMC algorithm.
 * 
 * @param e The shared pointer to the EqualityComputer object.
 */
void MCMC::set_eq_comp(std::shared_ptr<EqualityComputer> e){
    equalitycomp = e;
}

/**
 * @brief Sets the exact equality computer for the MCMC algorithm.
 * 
 * @param e The shared pointer to the EqualityComputer object.
 */
void MCMC::set_exact_eq_comp(std::shared_ptr<EqualityComputer> e){
    exact_eq_comp = e;
}

/**
 * @brief Sets the mutator for the MCMC algorithm.
 * 
 * @param m A shared pointer to the Mutator object.
 */
void MCMC::set_mutator(std::shared_ptr<Mutator> m){
    mutator = m;
}

/**
 * @brief Sets the temperature scheme for the MCMC algorithm.
 * 
 * @param t A shared pointer to a TemperatureScheme object.
 */
void MCMC::set_temp_scheme(std::shared_ptr<TemperatureScheme> t){
    temp_scheme = t;
}

/**
 * @brief Returns the equality computer used by the MCMC algorithm.
 * 
 * @return A shared pointer to the EqualityComputer object.
 */
std::shared_ptr<EqualityComputer> MCMC::getEqualityComputer() {
    return equalitycomp;
};

/**
 * @brief Returns the mutator object used by the MCMC algorithm.
 * 
 * @return A shared pointer to the Mutator object.
 */
std::shared_ptr<Mutator> MCMC::getMutator() {
    return mutator;
};

/**
 * @brief Returns the temperature scheme used in the MCMC algorithm.
 * 
 * @return A shared pointer to the TemperatureScheme object.
 */
std::shared_ptr<TemperatureScheme> MCMC::getTemperatureScheme() {
    return temp_scheme;
};

/**
 * Calculates the equality cost of a gate circuit with respect to a given QubitIndependentPartialMatrix.
 * 
 * @param circ The gate circuit to evaluate.
 * @param matrix_obj The QubitIndependentPartialMatrix to compare against.
 * @param circ_helper The CircuitHelper object containing circuit information.
 * @param log_val A boolean flag indicating whether to log the calculated value.
 * @return The normalized equality cost of the gate circuit.
 */
double MCMC::equalityCost(GateCircuit &circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& circ_helper, bool log_val){
    double equality_cost = 0.0;
    if (enable_permutations) {
        equality_cost = equalitycomp->normalizedEqualityCost(circ, matrix_obj, circ_helper); //~= 1
    } else {
        equality_cost = equalitycomp->normalizedEqualityCost(circ, matrix_obj.original, circ_helper);
    }
    return equality_cost;
}

/**
 * Calculates the equality cost of a gate circuit.
 *
 * This function calculates the equality cost of a gate circuit by comparing it with a given matrix object.
 * The equality cost is a measure of how closely the circuit matches the desired matrix.
 *
 * @param circ The gate circuit to calculate the equality cost for.
 * @param matrix_obj The matrix object to compare the circuit with.
 * @param circ_helper The circuit helper object used for circuit operations.
 * @param log_val A boolean flag indicating whether to log the equality cost value.
 * @return The equality cost of the gate circuit.
 */
double MCMC::equalityCost(GateCircuit &circ, PartialMatrix& matrix_obj, CircuitHelper& circ_helper, bool log_val){
    double equality_cost = equalitycomp->normalizedEqualityCost(circ, matrix_obj, circ_helper); //~= 1
    return equality_cost;
}

/**
 * Calculates the total energy cost for a given equilibrium cost and performance cost.
 *
 * @param eq_cost The equilibrium cost.
 * @param perf_cost The performance cost.
 * @return The total energy cost.
 */
double MCMC::getEnergy(double eq_cost){
    return eq_cost;
}

/**
 * Calculates the energy of a given gate circuit.
 * The energy is determined by the equality cost of the circuit.
 *
 * @param circ The gate circuit to calculate the energy for.
 * @param matrix_obj The QubitIndependentPartialMatrix object representing the circuit's matrix.
 * @param circ_helper The CircuitHelper object providing helper functions for the circuit.
 * @return The energy of the gate circuit.
 */
double MCMC::getEnergy(GateCircuit &circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& circ_helper){
    double equality_cost = equalityCost(circ, matrix_obj, circ_helper);
    return equality_cost;
}

/**
 * Corrects the result by ensuring qubit independence.
 * 
 * This function calculates the equality costs for each matrix in the QubitIndependentPartialMatrix object
 * and selects the matrix with the minimum cost. If the inverse_info flag for the selected matrix is true,
 * the circuit_best object is inverted. Then, the qubits of the circuit_best object are changed according
 * to the qbs_info of the selected matrix.
 * 
 * @param res The MCMCResult object to be corrected.
 * @param matrix_obj The QubitIndependentPartialMatrix object containing the matrices and associated information.
 * @param circ_helper The CircuitHelper object used for equality cost calculation.
 */
void MCMC::correctResultQubitIndependence(MCMCResult& res, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& circ_helper) {
    std::vector<double> equality_costs = {};
    for (int i = 0; i < matrix_obj.matrices.size(); i++) {
        equality_costs.push_back(exact_eq_comp->normalizedEqualityCost(*res.circuit_best, *matrix_obj.matrices[i], circ_helper));
    } 
    int min_element = std::min_element(equality_costs.begin(), equality_costs.end()) - equality_costs.begin();
    if (matrix_obj.inverse_info[min_element]) {
        res.circuit_best->invert();
    }
    res.circuit_best->changeQubits(matrix_obj.qbs_info[min_element]);
}

/**
 * Logs debug information about the MCMC process.
 * 
 * @param circ The GateCircuit object representing the current circuit.
 * @param matrix_obj The QubitIndependentPartialMatrix object representing the current matrix.
 * @param current_result The MCMCResult object containing the current MCMC result.
 * @param file The ofstream object used to write the debug information to a file.
 * @param ch The CircuitHelper object used to assist with circuit operations.
 */
void MCMC::log_debug_information(GateCircuit& circ, QubitIndependentPartialMatrix& matrix_obj, 
                                MCMCResult& current_result, std::ofstream& file, CircuitHelper& ch) {
    double eq_cost = equalityCost(circ, matrix_obj, ch);
    double temperature = temp_scheme->getTemperature();

    file << "Temperature: " << std::to_string(temperature) << std::endl;
    file << "Current Equality: " << std::to_string(eq_cost) << std::endl;
    file << "Best Energy: " << std::to_string(current_result.best_energy) << std::endl;
    file << "Best Equality: " << std::to_string(current_result.best_eq) << std::endl;
    file << "Current Circuit" << std::endl;
    file << circ.print_qasm() << std::endl;
    file << "Current Matrix" << std::endl;
    file << circ.toMatrix() << std::endl;
    file << " " << std::endl;
    file << " " << std::endl;
}

std::string MCMC::find_next_available_file(std::string& debug_folder) {
    int current_file = 0;
    while (std::filesystem::exists(debug_folder + std::to_string(current_file) + ".txt")) {
        current_file++;
    }
    return debug_folder + std::to_string(current_file) + ".txt";
}
