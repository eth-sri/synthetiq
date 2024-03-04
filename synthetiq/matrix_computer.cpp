#include "matrix_computer.h"
#include "utils.h"
#include <sstream>
#include <string>
#include <iterator>
#include <filesystem>

/**
 * @brief Default constructor for the LinearMatrixComputer class. The linear matrix computer just computes the matrix of the circuit by multiplying the matrices of the gates.
 */
LinearMatrixComputer::LinearMatrixComputer() {
}

/**
 * Updates the matrix.
 * 
 * @param position The position at which the circuit was updated.
 * @param list_gates The list of gates to calculate the matrix from.
 */
void LinearMatrixComputer::updateMatrix(int position, std::vector<std::shared_ptr<Gate>> list_gates) {
    calculateMatrix(list_gates);
}

/**
 * Calculates the linear matrix representation of a quantum circuit.
 * 
 * @param list_gates The list of gates in the circuit.
 */
void LinearMatrixComputer::calculateMatrix(std::vector<std::shared_ptr<Gate>> list_gates) {
    int nb_qbs = list_gates[0]->nb_qbs;
    matrix = Eigen::MatrixXcd::Identity(pow(2, nb_qbs), pow(2, nb_qbs));
    for (int i = 0; i < list_gates.size(); i++) {
        matrix = list_gates[i]->matrix * matrix;
    }
}

/**
 * @brief Returns the matrix representation of the LinearMatrixComputer.
 * 
 * @return The matrix representation of the LinearMatrixComputer.
 */
Eigen::MatrixXcd LinearMatrixComputer::getMatrix() {
    return matrix;
}

/**
 * Constructor for the ChunkMatrixComputer class. This computer divides the circuit into chunks and computes the matrix of each chunk. 
 * The final matrix is the product of the matrices of the chunks. Only when a position in a specific chunk is updated, the matrix of that chunk is recomputed.
 */
ChunkMatrixComputer::ChunkMatrixComputer() {

}

/**
 * @brief Constructs a ChunkMatrixComputer object.
 * 
 * @param nb_gates The number of gates.
 * @param nb_qbs The number of qubits.
 */
ChunkMatrixComputer::ChunkMatrixComputer(int nb_gates, int nb_qbs): nb_gates(nb_gates), nb_qbs(nb_qbs) {
    initializeChunks(nb_gates, nb_qbs);
}

/**
 * Initializes the chunks for the ChunkMatrixComputer.
 * 
 * @param nb_gates The number of gates.
 * @param n_qubits The number of qubits.
 */
void ChunkMatrixComputer::initializeChunks(int nb_gates, int n_qubits) {
    nb_qbs = n_qubits;
    chunks = {};
    int nb_chunks = (int) std::sqrt(nb_gates);
    for (int i = 0; i < nb_chunks; i++) {
        chunks.push_back(Eigen::MatrixXcd::Identity(pow(2, nb_qbs), pow(2, nb_qbs)));
    }
}

/**
 * @brief Updates the matrix of the ChunkMatrixComputer using the given list of gates.
 * 
 * This function divides the list of gates into chunks and computes the matrix for each chunk.
 * The matrix for each chunk is computed by multiplying the gates' matrices in the chunk.
 * The computed matrices are stored in the chunks vector. Only when a position in a specific chunk is updated, the matrix of that chunk is recomputed.
 * 
 * @param position The position at which to update the matrix.
 * @param list_gates The list of gates to use for matrix computation.
 */
void ChunkMatrixComputer::updateMatrix(int position, std::vector<std::shared_ptr<Gate>> list_gates) {
    int chunk_factor = std::ceil(nb_gates / (double) chunks.size());
    int chunk = std::floor(position / chunk_factor);
    chunks[chunk] = Eigen::MatrixXcd::Identity(1 << nb_qbs, 1 << nb_qbs);
    for (int i = chunk * chunk_factor; i < std::min((chunk + 1) * chunk_factor, (int) list_gates.size()); i++) {
        chunks[chunk] = list_gates[i]->matrix * chunks[chunk];
    }
    matrix_computed = false;
}

/**
 * Calculates the matrix representation of a quantum circuit given a list of gates.
 * If the size of the list of gates is different from the number of gates in the circuit,
 * the chunks are initialized accordingly.
 * 
 * @param list_gates The list of gates in the circuit.
 */
void ChunkMatrixComputer::calculateMatrix(std::vector<std::shared_ptr<Gate>> list_gates) {
    if (list_gates.size() != nb_gates){ 
        initializeChunks(list_gates.size(), nb_qbs);
    }

    int chunk_factor = std::ceil(nb_gates / (double) chunks.size());
    for (int i = 0; i < chunks.size(); i++) {
        updateMatrix(i * chunk_factor, list_gates);
    }
}

/**
 * Gets the matrix representation of the ChunkMatrixComputer.
 *
 * @return The matrix representation of the ChunkMatrixComputer.
 */
Eigen::MatrixXcd ChunkMatrixComputer::getMatrix() {
    if (!matrix_computed) {
        matrix = Eigen::MatrixXcd::Identity(1 << nb_qbs, 1 << nb_qbs);
        for (auto chunk = chunks.begin(); chunk < chunks.end(); chunk++) {
            matrix = *chunk * matrix;
        }
        matrix_computed = true;
    }
    return matrix;
}

/**
 * Constructor for the BinaryMatrixComputer class. This class computes the matrix of the circuit by using a binary tree.
 */
BinaryMatrixComputer::BinaryMatrixComputer() {
    
}

/**
 * @brief Constructs a BinaryMatrixComputer object.
 * 
 * @param n_gates The number of gates.
 * @param nb_qbs The number of qubits.
 */
BinaryMatrixComputer::BinaryMatrixComputer(int n_gates, int nb_qbs): nb_gates(n_gates), nb_qbs(nb_qbs) {
    initializeTree(n_gates, nb_qbs);
}

/**
 * Initializes the tree data structure for the BinaryMatrixComputer.
 * 
 * @param n_gates The number of gates in the circuit.
 * @param n_qubits The number of qubits in the circuit.
 */
void BinaryMatrixComputer::initializeTree(int n_gates, int n_qubits) {
    nb_qbs = n_qubits;
    int depth = std::max((int) std::ceil(std::log2(n_gates)), 1);
    tree = {};
    for (int i = 0; i < depth; i++) {
        tree.push_back({});
        int size = std::ceil(n_gates / pow(2, i + 1));
        for (int j = 0; j < size; j++) {
            tree[i].push_back(Eigen::MatrixXcd::Identity(pow(2, nb_qbs), pow(2, nb_qbs)));
        }
    }
}

/**
 * Updates the binary matrix at the specified index with the given list of gates.
 * The binary matrix is updated based on the position of the gates in the list.
 * The updated matrix is stored in the tree data structure.
 *
 * @param i The index at which the circuit was changed
 * @param list_gates The list of gates.
 */
void BinaryMatrixComputer::updateMatrix(int i, std::vector<std::shared_ptr<Gate>> list_gates) {
    int starting_pos = tree[0].size() - 1 - i / 2;
    if (i % 2 == 0) {
        if (i == list_gates.size() - 1) {
            //note: the gates need to be turned around because the last gate is applied
            tree[0][starting_pos] = list_gates[i]->matrix;
        } else {
            tree[0][starting_pos] = list_gates[i + 1]->matrix * list_gates[i]->matrix;
        }
    } else {
        tree[0][starting_pos] = list_gates[i]->matrix * list_gates[i - 1]->matrix;
    }

    for (int current_depth = 1; current_depth < tree.size(); current_depth++) {
        int position = starting_pos/ pow(2, current_depth);
        int earlier_position = starting_pos / pow(2, current_depth - 1);
        if (earlier_position % 2 == 0) {
            if (earlier_position == tree[current_depth - 1].size() - 1) {
                tree[current_depth][position] = tree[current_depth - 1][earlier_position];
            } else {
                tree[current_depth][position] = tree[current_depth - 1][earlier_position] * tree[current_depth - 1][earlier_position + 1];
            }
        } else {
            tree[current_depth][position] = tree[current_depth - 1][earlier_position - 1] * tree[current_depth - 1][earlier_position];
        }
    }
}

/**
 * Calculates the matrix representation of a quantum circuit given a list of gates.
 * If the size of the list of gates is different from the expected number of gates, the tree is initialized.
 * The matrix calculation is performed by traversing the tree and applying the gate matrices.
 * 
 * @param list_gates The list of gates to be applied in the circuit.
 */
void BinaryMatrixComputer::calculateMatrix(std::vector<std::shared_ptr<Gate>> list_gates) {
    if (list_gates.size() != nb_gates) {
        initializeTree(list_gates.size(), nb_qbs);
    }

    for (int i = 0; i < list_gates.size(); i += 2) {
        if (i == list_gates.size() - 1) {
            tree[0][tree[0].size() - 1 - i / 2] = list_gates[i]->matrix;
        } else {
            tree[0][tree[0].size() - 1 - i / 2] = list_gates[i + 1]->matrix * list_gates[i]->matrix;
        }
    }
    
    for (int current_depth = 1; current_depth < tree.size(); current_depth++) {
        for (int earlier_position = 0; earlier_position < tree[current_depth - 1].size(); earlier_position += 2) {
            int position = earlier_position / 2;
            if (earlier_position == tree[current_depth - 1].size() - 1) {
                tree[current_depth][position] = tree[current_depth - 1][earlier_position];
            } else {
                tree[current_depth][position] = tree[current_depth - 1][earlier_position] * tree[current_depth - 1][earlier_position + 1];
            }
        }
    }
}

/**
 * @brief Returns the matrix representation of the BinaryMatrixComputer.
 * 
 * This function returns the matrix representation of the BinaryMatrixComputer
 * as an Eigen::MatrixXcd object.
 * 
 * @return The matrix representation of the BinaryMatrixComputer.
 */
Eigen::MatrixXcd BinaryMatrixComputer::getMatrix() {
    return tree[tree.size() - 1][0];
}
