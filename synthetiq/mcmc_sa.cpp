#include "mcmc_sa.h"
#include <math.h>
#include <sstream>
#include <chrono>
#include <omp.h>


/**
 * @brief Constructor for the MCMC_Sa class.
 * 
 * @param random_helper The RandomHelper object used for generating random numbers.
 * @param proba_id The probability of performing an identity move.
 * @param factor_nb_steps The factor by which the number of steps is multiplied.
 * @param enable_permutations Flag indicating whether permutations are enabled.
 */
MCMC_Sa::MCMC_Sa(RandomHelper& random_helper, double proba_id, double factor_nb_steps, bool enable_permutations): 
                        MCMC(random_helper, proba_id, enable_permutations), factor_nb_steps(factor_nb_steps) {
}

/**
 * @brief Clones the MCMC object.
 * 
 * @return A shared pointer to the cloned MCMC object.
 */
std::shared_ptr<MCMC> MCMC_Sa::clone() {
    std::shared_ptr<MCMC> cloned = std::make_shared<MCMC_Sa>(*this);
    cloned->set_mutator(mutator->clone());
    cloned->set_temp_scheme(temp_scheme->clone());
    cloned->set_eq_comp(equalitycomp->clone());
    return cloned;
}

/**
 * Calculates the probability of a given gate circuit using the Metropolis-Hastings algorithm.
 * The probability is calculated based on the energy of the circuit, the given matrix object,
 * the circuit helper, and the temperature.
 *
 * @param circ The gate circuit to calculate the probability for.
 * @param matrix_obj The QubitIndependentPartialMatrix object representing the matrix.
 * @param circ_helper The CircuitHelper object providing circuit-related helper functions.
 * @param temperature The temperature parameter for the Metropolis-Hastings algorithm.
 * @return The probability of the gate circuit.
 */
double MCMC_Sa::getProbability(GateCircuit &circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& circ_helper, double temperature){
    double tot_cost = getEnergy(circ, matrix_obj, circ_helper);
    return std::exp(-tot_cost / temperature);
}

/**
 * Determines whether to accept a mutation based on the Metropolis criterion.
 * 
 * @param rd_val The random value used for acceptance probability calculation.
 * @param candidate_energy The energy of the candidate solution.
 * @param cur_energy The energy of the current solution.
 * @param temperature The current temperature of the system.
 * @return True if the mutation should be accepted, false otherwise.
 */
bool MCMC_Sa::acceptMutation(double rd_val, double candidate_energy, double cur_energy, double temperature){
    if(candidate_energy <= cur_energy)
        return true;
    double proba = std::exp(-(candidate_energy - cur_energy) / temperature);
    return rd_val <= proba;
}

/**
 * Calculates the index of the closest matrix in the given QubitIndependentPartialMatrix to the given GateCircuit.
 * The closeness is determined based on the equality cost between the matrices.
 *
 * @param matrix_obj The QubitIndependentPartialMatrix object containing the matrices to compare.
 * @param circ The GateCircuit object to compare against the matrices.
 * @param circ_helper The CircuitHelper object used for calculating the equality cost.
 * @return The index of the closest matrix in the QubitIndependentPartialMatrix.
 */
int MCMC_Sa::calculateClosestMatrix(QubitIndependentPartialMatrix& matrix_obj, GateCircuit& circ, CircuitHelper& circ_helper) {
    double minimum_cost = INFINITY;
    int best_index = 0;
    for (int i = 0; i < matrix_obj.matrices.size(); i++) {
        double cost = equalityCost(circ, *matrix_obj.matrices[i], circ_helper);
        if (cost < minimum_cost) {
            minimum_cost = cost;
            best_index = i;
        }
    }
    return best_index;
}

/**
 * Runs the MCMC algorithm using the Simulated Annealing variant.
 * 
 * @param matrix_obj The QubitIndependentPartialMatrix object representing the matrix.
 * @param init The initial GateCircuit object.
 * @param circ_helper The CircuitHelper object providing circuit-related helper functions.
 * @param debug Flag indicating whether to log debug information.
 * @param debug_folder The folder to log debug information to.
 * @return The MCMCResult object containing the results of the MCMC algorithm.
 */
MCMCResult MCMC_Sa::run(QubitIndependentPartialMatrix& matrix_obj, std::shared_ptr<GateCircuit> init, CircuitHelper& circ_helper, 
                         bool debug, std::string debug_folder) {
    MCMCResult res;
    std::string debug_file = "";
    std::ofstream file;
    if (debug) {
        debug_file = find_next_available_file(debug_folder);
        file.open(debug_file);
    }
    std::shared_ptr<GateCircuit> candidate_circuit = init->clone();
    res.best_eq = equalityCost(*init, matrix_obj, circ_helper); //~= 1    
    double cur_energy = getEnergy(res.best_eq);
    res.best_energy = cur_energy;
    res.circuit_best = candidate_circuit->clone();
    int n_accepted_mutations = 0;
    temp_scheme->reset();
    bool found = false;

    for (int cur_step = 0; cur_step < factor_nb_steps * init->nbElements(); cur_step++) {
        temp_scheme->updateTemperature(cur_step, n_accepted_mutations, init->nbElements());
        bool unchanged = mutator->mutate(*candidate_circuit, random_helper);
        if (unchanged) {
            continue;
        }

        double candidate_eq_cost = equalityCost(*candidate_circuit, matrix_obj, circ_helper);
        double candidate_energy = getEnergy(candidate_eq_cost);
        double u = random_helper.random01();
        if (acceptMutation(u, candidate_energy, cur_energy, temp_scheme->getTemperature())) { //candidate accepted
            
            n_accepted_mutations++;
            if (candidate_energy < res.best_energy)
            {   
                res.best_energy = candidate_energy;
                res.circuit_best = candidate_circuit->clone();
                res.best_eq = candidate_eq_cost; //~= 1
                if (exact_eq_comp->normalizedEqualityCost(*res.circuit_best, matrix_obj, circ_helper) < 1e-3) { // we assume cost is 0 for found and otherwise higher
                    found = true;
                    break;
                }
            }

            cur_energy = candidate_energy;
        } else {
            mutator->undo_mutation(*candidate_circuit);
        }
        if (debug) {
            log_debug_information(*candidate_circuit, matrix_obj, res, file, circ_helper);
        }
    }
    if (found) {
        correctResultQubitIndependence(res, matrix_obj, circ_helper);
    }
    
    if (debug) {
        file.close();
    }

    return res;
}