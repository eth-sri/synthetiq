#include "matrixGenerator.h"
#include "partialMatrix.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <sstream>
#include <string>

/**
 * @brief Adds ancilla qubits to a matrix.
 * 
 * @param matrix The matrix to which the ancilla qubits are added.
 * @param n_ancilla The number of ancilla qubits to add.
 * @param name The name of the matrix.
 * @return PartialMatrix The matrix with the ancilla qubits added.
 */
PartialMatrix MatrixGenerator::addAncilla(PartialMatrix& matrix, int n_ancilla, std::string name){
    int total_qubits = matrix.getNQubits() + n_ancilla;
    int size = pow(2, total_qubits);
    int size_matrix = pow(2, matrix.getNQubits());
    Eigen::MatrixXcd result_matrix(size, size);
    result_matrix.setZero();
    result_matrix.topLeftCorner(size_matrix, size_matrix) = matrix.getMatrix();
    BoolMatrix cover(size, size);
    cover.setConstant(size, size, false);
    cover.topLeftCorner(size_matrix, size_matrix) = matrix.getCover();
    cover.bottomLeftCorner(size_matrix, size - size_matrix).setConstant(true);
    return PartialMatrix(result_matrix, cover, name);
};

/**
 * @brief Adds dirty qubits to a matrix.
 * 
 * @param matrix The matrix to which the dirty qubits are added.
 * @param n_dirty The number of dirty qubits to add.
 * @param name The name of the matrix.
 * @return PartialMatrix The matrix with the dirty qubits added.
 */
PartialMatrix MatrixGenerator::addDirtyQubits(PartialMatrix& matrix, int n_dirty, std::string name){
    int total_qubits = matrix.getNQubits() + n_dirty;
    int size = pow(2, total_qubits);
    Eigen::MatrixXcd result_matrix(size, size);
    result_matrix.setZero();
    BoolMatrix cover(size, size);
    cover.setConstant(size, size, true);
    int size_identity = pow(2, n_dirty);
    Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(size_identity, size_identity);
    BoolMatrix identity_bool = BoolMatrix::Identity(size_identity, size_identity);
    result_matrix = Eigen::kroneckerProduct(identity, matrix.getMatrix());
    cover = Eigen::kroneckerProduct(identity_bool, matrix.getCover());
    return PartialMatrix(result_matrix, cover, name);
};
