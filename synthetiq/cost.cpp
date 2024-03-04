#include "cost.h"
#include "utils.h"
#include <Eigen/Dense>
#include <complex>


/**
 * Calculates the normalized equality cost for a given gate circuit, partial matrix, and circuit helper.
 * 
 * @param circ The gate circuit to calculate the cost for.
 * @param matrix_obj The partial matrix to calculate the cost for.
 * @param ch The circuit helper to assist in the cost calculation.
 * @return The normalized equality cost.
 */
double EqualityComputer::normalizedEqualityCost(GateCircuit& circ, PartialMatrix& matrix_obj, CircuitHelper& ch){
    return 0.0;
}

/**
 * Calculates the normalized equality cost for a given gate circuit, qubit independent partial matrix, and circuit helper.
 * The normalized equality cost is the minimum cost among all the partial matrices in the matrix object.
 *
 * @param circ The gate circuit.
 * @param matrix_obj The qubit independent partial matrix object.
 * @param ch The circuit helper.
 * @return The normalized equality cost.
 */
double EqualityComputer::normalizedEqualityCost(GateCircuit& circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& ch) {
    double cur_norm = INFINITY;
    for (std::shared_ptr<PartialMatrix> matrix : matrix_obj.matrices) {
        cur_norm = std::min(cur_norm, normalizedEqualityCost(circ, *matrix, ch));
    }
    return cur_norm;
}

/**
 * @brief Clones the ExactEqualityComputer object.
 * 
 * @return A shared pointer to a new ExactEqualityComputer object.
 */
std::shared_ptr<EqualityComputer> ExactEqualityComputer::clone() {
    return std::make_shared<ExactEqualityComputer>(*this);
}

/**
 * Calculates the normalized equality cost between a gate circuit and a partial matrix constraint. Returns 1 if equal, 0 otherwise.
 * 
 * @param circ The gate circuit to compare.
 * @param constraint The partial matrix constraint.
 * @param ch The CircuitHelper object.
 * @return The normalized equality cost.
 */
double ExactEqualityComputer::normalizedEqualityCost(GateCircuit& circ, PartialMatrix& constraint, CircuitHelper& ch){
	// note we assume all non specified elements of the constraint are set to 0
    // This part of the code implements E_{1, part} from the paper, but rewritten to increase performance.
    double normalization_cst = 0.0;
	// note we assume all non specified elements of the constraint are set to 0
	Eigen::MatrixXcd matr_circ = circ.toMatrix();
    double circ_size = 0.0;
    std::complex<double> conj = 0.0;
	int size = pow(2, circ.nb_qbs);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (constraint.cover(i, j)) {
				normalization_cst += 1;
                circ_size += pow(std::abs(matr_circ(i, j)), 2);
                conj += std::conj(constraint.matrix(i, j)) * matr_circ(i, j);
            }
		}
	}
	double distance = 1 / std::sqrt(2) * std::sqrt(std::max(0.0, constraint.squared_norm + circ_size - 2 * std::abs(conj))) / std::sqrt(std::sqrt(normalization_cst));

    if (distance <= tolerance) {
        return 0;
    }
    return 1;
}

/**
 * Calculates the normalized equality cost for a given gate circuit, qubit independent partial matrix, and circuit helper.
 * The normalized equality cost is the minimum cost among all the partial matrices in the matrix object.
 *
 * @param circ The gate circuit.
 * @param matrix_obj The qubit independent partial matrix object.
 * @param ch The circuit helper.
 * @return The normalized equality cost.
 */
double ExactEqualityComputer::normalizedEqualityCost(GateCircuit& circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& ch) {
    double cur_norm = INFINITY;
    for (std::shared_ptr<PartialMatrix> matrix : matrix_obj.matrices) {
        cur_norm = std::min(cur_norm, normalizedEqualityCost(circ, *matrix, ch));
    }
    return cur_norm;
}

/**
 * @brief Default constructor for the FroebeniusCostComputer class.
 * This class implements E_{1, part} from the paper, but rewritten for increased performance.
 */
FroebeniusCostComputer::FroebeniusCostComputer() {}

/**
 * @brief Clones the FroebeniusCostComputer object.
 * 
 * @return A shared pointer to a new FroebeniusCostComputer object that is a clone of the current object.
 */
std::shared_ptr<EqualityComputer> FroebeniusCostComputer::clone() {
    return std::make_shared<FroebeniusCostComputer>(*this);
}

/**
 * Calculates the normalized equality cost for a given gate circuit, partial matrix constraint, and circuit helper.
 * This function implements E_{1, part} from the paper, but rewritten for increased performance.
 *
 * @param circ The gate circuit.
 * @param constraint The partial matrix constraint.
 * @param ch The circuit helper.
 * @return The normalized equality cost.
 */
double FroebeniusCostComputer::normalizedEqualityCost(GateCircuit& circ, PartialMatrix& constraint, CircuitHelper& ch){
    double normalization_cst = 0.0;
	Eigen::MatrixXcd matr_circ = circ.toMatrix();
    double circ_size = 0.0;
    std::complex<double> conj = 0.0;
	int size = pow(2, circ.nb_qbs);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (constraint.cover(i, j)) {
				normalization_cst += 1;
                circ_size += pow(std::abs(matr_circ(i, j)), 2);
                conj += std::conj(constraint.matrix(i, j)) * matr_circ(i, j);
            }
		}
	}
    // max is for rounding errors
    double cost = std::sqrt(std::max(0.0, constraint.squared_norm + circ_size - 2 * std::abs(conj))) / std::sqrt(std::sqrt(normalization_cst));
	return cost;
}

/**
 * @brief Default constructor for the SimpleFroebeniusCostComputer class.
 * The class implements the equality cost if not rewritten as done in the paper
 */
SimpleFroebeniusCostComputer::SimpleFroebeniusCostComputer() {}

/**
 * @brief Clones the SimpleFroebeniusCostComputer object.
 * 
 * @return A shared pointer to a cloned SimpleFroebeniusCostComputer object.
 */
std::shared_ptr<EqualityComputer> SimpleFroebeniusCostComputer::clone() {
    return std::make_shared<SimpleFroebeniusCostComputer>(*this);
}

/**
 * Calculates the normalized equality cost for a given gate circuit, partial matrix constraint, and circuit helper.
 * The function implements the equality cost if not rewritten as done in the paper
 * @param circ The gate circuit.
 * @param constraint The partial matrix constraint.
 * @param ch The circuit helper.
 * @return The normalized equality cost.
 */
double SimpleFroebeniusCostComputer::normalizedEqualityCost(GateCircuit& circ, PartialMatrix& constraint, CircuitHelper& ch){
    double normalization_cst = 0.0;
	Eigen::MatrixXcd matr_circ = circ.toMatrix();
    double circ_size = 0.0;
    std::complex<double> conj = 0.0;
	int size = pow(2, circ.nb_qbs);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (constraint.cover(i, j)) {
				normalization_cst += 1;
                circ_size += pow(std::abs(matr_circ(i, j)), 2);
                conj += std::conj(constraint.matrix(i, j)) * matr_circ(i, j);
            }
		}
	}
	return std::sqrt(std::max(0.0, 1 - std::abs(conj) / std::sqrt(normalization_cst))) * std::sqrt(2);
}

/**
 * Constructor for the GatesSumComputer class. This calculates the cost based on the sum of the costs of the gates in the circuit.
 
 */
GatesSumComputer::GatesSumComputer() {

}

/**
 * @brief Constructs a GatesSumComputer object by copying another GatesSumComputer object.
 * 
 * @param other The GatesSumComputer object to be copied.
 */
GatesSumComputer::GatesSumComputer(GatesSumComputer const& other) {
}

/**
 * @brief Clones the GatesSumComputer object.
 * 
 * @return A shared pointer to a cloned GatesSumComputer object.
 */
std::shared_ptr<PerformanceComputer> GatesSumComputer::clone() {
    return std::make_shared<GatesSumComputer>(*this);
}

/**
 * Calculates the normalized performance cost of a gate circuit.
 * 
 * @param circ The gate circuit to calculate the cost for.
 * @param ch The circuit helper object.
 * @return The normalized performance cost of the gate circuit.
 */
double GatesSumComputer::normalizedPerformanceCost(GateCircuit& circ, CircuitHelper& ch){
    double performance_cost = circ.getCost();
    performance_cost /= (circ.performanceNormalizationCst() * ch.max_cost_basic);
    return performance_cost;
}