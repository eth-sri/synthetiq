#include "generator.h"
#include "resynthesis.h"
#include "temperatureScheme.h"
#include <filesystem>

/**
 * @brief Constructs a Generator object. The object generates random circuits and matrices.
 * 
 * @param append_mode Specifies whether the Generator should be in append mode or not.
 */
Generator::Generator (bool append_mode) : append_mode(append_mode) {

}


/**
 * Generates a random gate circuit.
 *
 * @param n_gate The number of gates in the circuit.
 * @param n_qubits The number of qubits in the circuit.
 * @param rg The RandomCircuitGen object used for generating random circuits.
 * @param ch The CircuitHelper object used for circuit manipulation.
 * @param rh The RandomHelper object used for random number generation.
 * @param mutator The Mutator object used for circuit mutation.
 * @param extra_gates The number of extra gates to add to the circuit.
 * @return A shared pointer to the generated GateCircuit object.
 */
std::shared_ptr<GateCircuit> Generator::generateCircuit(int n_gate, int n_qubits, 
                            RandomCircuitGen& rg, CircuitHelper& ch, RandomHelper& rh, Mutator& mutator, int extra_gates) {
    bool new_circuit_found = false;
    GateCircuit rand_circ = GateCircuit(ch);
    Resynthesize resynth = Resynthesize();
    while (!new_circuit_found) {
        rand_circ = rg.randomGateCircuit(n_gate + extra_gates, n_qubits, ch, mutator);
        std::shared_ptr<GateCircuit> rand_circ_ptr = std::make_shared<GateCircuit>(rand_circ);
        resynth.run(rand_circ_ptr, *rand_circ.getCircuitHelper());
        int nb_ref_gates = rand_circ_ptr->getNbNonIdGates(); 
        if (nb_ref_gates >= n_gate) {
            new_circuit_found = true;
        }
    }  
    return std::make_shared<GateCircuit>(rand_circ); 
}

/**
 * @brief Generates the matrix corresponding to a random circuit. It decides the indices at which the matrix is specified
 * 
 * @param circuit The circuit for which the matrix is to be generated.
 * @param type The type of matrix to be generated.
 * @param rg The RandomCircuitGen object used for generating random circuits.
 * @param ch The CircuitHelper object used for circuit manipulation.
 * @param rh The RandomHelper object used for random number generation.
 * @return The generated matrix.
 */
PartialMatrix Generator::generateMatrixFromCircuit(std::shared_ptr<GateCircuit> circuit, int type, 
                                                    RandomCircuitGen& rg, CircuitHelper& ch, RandomHelper& rh) {
    PartialMatrix matrix = PartialMatrix(*circuit, ch);
    int size = pow(2, circuit->nb_qbs);
    if (type == Isometry) {
        std::vector<int> cols_specified = {};
        while (cols_specified.size() == 0 || cols_specified.size() == size) {
            double random_prob = rh.random01();
            cols_specified = {};
            for (int i = 0; i < size; i++) {
                if (rh.random01() > random_prob) {
                    cols_specified.push_back(i);
                }
            }
        }
        BoolMatrix cover(size, size);
        cover.setConstant(size, size, false);
        for (int col : cols_specified) {
            cover.col(col) = Eigen::Matrix<bool, Eigen::Dynamic, 1>::Constant(size, true);
        }
        matrix = PartialMatrix(circuit->toMatrix(), cover);
    } else if (type == Random) {
        BoolMatrix cover(size, size);
        cover.setConstant(size, size, false);
        double random_prob = rh.random01();
        while (cover.count() == 0 || cover.count() == pow(size, 2)) {
            cover.setConstant(size, size, false);
            random_prob = rh.random01();
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    if (random_prob < rh.random01()) {
                        cover(i, j) = true;
                    }
                }
            }
        }
        matrix = PartialMatrix(circuit->toMatrix(), cover);
    }
    return matrix;
}
