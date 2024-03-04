#include <Eigen/Dense>
#include "utils.h"
#include <iostream>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/Core>
#include <sstream>
#include <string>
#include <iterator>
#include <regex>


/**
 * @brief Converts a decimal number to a binary representation.
 * 
 * @param number The decimal number to convert.
 * @param size The size of the binary representation.
 * @return std::vector<int> The binary representation of the decimal number.
 */
std::vector<int> Utils::dec2bin(int number, int size) {
    std::vector<int> output(size);
    int d = size - 1;
    while (number > 0 && d >= 0) {
        if (number % 2 == 1) {
            output[d] = 1;
            number -= 1;
        } else {
            output[d] = 0;
        }
        d -= 1;
        number /= 2;
    }
    return output;
}

/**
 * Converts a binary number represented as a vector of integers to its decimal equivalent.
 * 
 * @param bin The binary number to be converted.
 * @return The decimal equivalent of the binary number.
 */
int Utils::bin2dec(std::vector<int> bin) {
    int output = 0;
    for (int i = bin.size() - 1; i >= 0; i--) {
        output += pow(2, bin.size() - i - 1) * bin[i];
    }
    return output;
};

/**
 * Computes the row change for a matrix with qubits in a specific order.
 * 
 * @param qbs The qubits in the order they appear in the matrix.
 * @param nb_qbs The number of qubits in the matrix.
 */
std::vector<int> Utils::getRowChangeQubitIndependent(std::vector<int> qbs, int nb_qbs) {
    std::vector<int> row_change = {};
    for (int i = 0; i < pow(2, nb_qbs); i++) {
            std::vector<int> bin = dec2bin(i, nb_qbs);
            std::vector<int> to_bin(nb_qbs);
            for (int j = 0; j < nb_qbs; j++) {
                to_bin[nb_qbs - 1 - qbs[j]] = bin[nb_qbs - j - 1];
            }
            int to_dec = Utils::bin2dec(to_bin);
            row_change.push_back(to_dec);
    }
    return row_change;
}

/**
 * Increases the number of qubits of a matrix by adding identity matrices to the left and right.
 * 
 * @param matrix The matrix to be increased.
 * @param nb_qbs The number of qubits to increase to.
 * @return The increased matrix.
 */
Eigen::MatrixXcd Utils::increaseQubits(Eigen::MatrixXcd& matrix, int nb_qbs) {
    int nb_extra_qbs = nb_qbs - std::log2(matrix.rows());
    Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(pow(2, nb_extra_qbs), pow(2, nb_extra_qbs));
    return Eigen::kroneckerProduct(identity, matrix).eval();
}

/**
 * Changes the order of the qubits in a matrix.
 * 
 * @param matrix The matrix to be changed.
 * @param qbs_change The order of the qubits in the matrix.
 * @return The changed matrix.
 */
Eigen::MatrixXcd Utils::changeQubits(Eigen::MatrixXcd& matrix, std::vector<int> qbs_change) {
    Eigen::MatrixXcd new_matrix = matrix;
    int nb_qbs = std::log2(matrix.rows());
    std::vector<int> change = Utils::getRowChangeQubitIndependent(qbs_change, nb_qbs);
    for (int i = 0; i < pow(2, nb_qbs); i++) {
        //Note: there is a mismatch between this one and the one used by PartialMatrix
        new_matrix.row(change[i]) = matrix.row(i);
    }
    Eigen::MatrixXcd new_new_matrix = new_matrix;
    for (int i = 0; i < pow(2, nb_qbs); i++) {
        new_new_matrix.col(change[i]) = new_matrix.col(i);
    }
    return new_new_matrix;
}

/**
 * @brief Compares two complex-valued matrices for equality.
 * 
 * This function checks if two matrices have the same dimensions and if the absolute difference
 * between corresponding elements is within a certain epsilon value.
 * 
 * @param matrix1 The first matrix to compare.
 * @param matrix2 The second matrix to compare.
 * @return True if the matrices are equal, false otherwise.
 */
bool Utils::matricesEqual(Eigen::MatrixXcd& matrix1, Eigen::MatrixXcd& matrix2) {
    if (matrix1.rows() != matrix2.rows() || matrix1.cols() != matrix2.cols()) {
        return false;
    }
    for (int i = 0; i < matrix1.rows(); i++) {
        for (int j = 0; j < matrix1.cols(); j++) {
            if (std::abs(matrix1(i, j) - matrix2(i, j)) > epsilon) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Reads a line from a qasm file and splits it into a vector of strings.
 * 
 * @param line The line to be split.
 * @return The vector of strings.
 */
std::vector<std::string> Utils::readQasmLineBase(std::string line) {
    std::stringstream string_stream(line);
    std::istream_iterator<std::string> begin(string_stream);
    std::istream_iterator<std::string> end;
    std::vector<std::string> splitted_line(begin, end);
    std::vector<std::string> totally_splitted_line = {};
    // just some magic to parse the string properly, with and without spaces etc.
    for (std::string word : splitted_line) {
        if (word.find(',') > 0 && word.find(",") != std::string::npos) {
            std::stringstream ss(word);
            std::string str;
            while (getline(ss, str, ',')) {
                totally_splitted_line.push_back(str);
            }
        } else {
            totally_splitted_line.push_back(word);
        }
    }
    return totally_splitted_line;
}

/**
 * @brief Reads the name from a QASM line.
 * 
 * This function extracts the name from a QASM line and returns it as a string.
 * 
 * @param line The QASM line to extract the name from.
 * @return The extracted name as a string.
 */
std::string Utils::readNameQasmLine(std::string line) {
    return Utils::readQasmLineBase(line)[0];
}

/**
 * @brief Reads the acting qubits from a QASM line and returns them as a vector of integers.
 * 
 * @param line The QASM line to read.
 * @return std::vector<int> The acting qubits as a vector of integers.
 */
std::vector<int> Utils::readActingQubitsQasmLine(std::string line) {
    std::vector<std::string> splitted_line = Utils::readQasmLineBase(line);
    std::vector<int> acting_qubits_gate = {};
    for (int i = 1; i < splitted_line.size(); i++) {
        acting_qubits_gate.push_back(std::stoi(std::regex_replace(splitted_line[i], std::regex(R"([\D])"), "")));
    }
    return acting_qubits_gate;
}

/**
 * Reads the number of qubits from a QASM line.
 * 
 * @param line The QASM line to parse.
 * @return The number of qubits extracted from the line.
 */
int Utils::readNbQbsQasmLine(std::string line) {
    return std::stoi(std::regex_replace(line, std::regex(R"([\D])"), ""));
}

/**
 * @brief Computes the trace of the conjugate product of two matrices.
 * 
 * @param matrix1 The first matrix.
 * @param matrix2 The second matrix.
 */
std::complex<double> Utils::traceConjugateProduct(Eigen::MatrixXcd& matrix1, Eigen::MatrixXcd& matrix2) {
    std::complex<double> output = 0.0;
    for (int i = 0; i < matrix1.rows(); i++) {
        for (int j = 0; j < matrix1.cols(); j++) {
            output += std::conj(matrix1(i, j)) * matrix2(i, j);
        }
    }
    return output;
};