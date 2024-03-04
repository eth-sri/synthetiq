#include "randomCircuit.h"
#include "mutation.h"

/**
 * @brief Constructs a RandomCircuitGen object.
 * 
 * @param random_helper The RandomHelper object used for generating random numbers.
 * @param id_prob The probability of generating an identity gate.
 * @param ensure_non_id Flag indicating whether to ensure the generated circuit is non-identity.
 */
RandomCircuitGen::RandomCircuitGen(RandomHelper& random_helper, double id_prob, bool ensure_non_id): random_helper(random_helper), id_prob(id_prob), ensure_non_id(ensure_non_id) {}

/**
 * @brief Generates a random circuit with the specified number of gates and qubits.
 * 
 * @param nb_gates The number of gates in the circuit.
 * @param nb_qbs The number of qubits in the circuit.
 * @param ch The CircuitHelper object used for resynthesis.
 * @param mutator The Mutator object used for mutating the circuit.
 * @return A random GateCircuit object.
 */
GateCircuit RandomCircuitGen::randomGateCircuit(int nb_gates, int nb_qbs, CircuitHelper& ch, Mutator& mutator){
    int size = nb_gates;
    GateCircuit res(size, nb_qbs, ch);
    for(int g = 0; g < size; g++){
        res.mutate(random_helper, g, 1.0, id_prob, 0.5);
        if (ensure_non_id) {
            while (res.getGates()[g]->name == ch.id_gate->name) {
                res.mutate(random_helper, g, 1.0, id_prob, 0.5);
            }
        }
    }
    return res;
}

/**
 * Generates a random gate circuit.
 *
 * @param min_nb_gates The minimum number of gates in the circuit.
 * @param max_nb_gates The maximum number of gates in the circuit.
 * @param nb_qbs The number of qubits in the circuit.
 * @param ch The CircuitHelper object used for circuit generation.
 * @param mutator The Mutator object used for circuit mutation.
 * @return The generated gate circuit.
 */
GateCircuit RandomCircuitGen::randomGateCircuit(int min_nb_gates, int max_nb_gates, int nb_qbs, CircuitHelper& ch, Mutator& mutator){
    int size = min_nb_gates + random_helper.randomInt(max_nb_gates - min_nb_gates);
    return randomGateCircuit(size, nb_qbs, ch, mutator);
}