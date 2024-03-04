#include "partialMatrix.h"
#include "circuit.h"
#include "utils.h"
#include <sstream>
#include <string>

/**
 * @brief Default constructor for the BaseConstraint class.
 */
BaseConstraint::BaseConstraint() {};

/**
 * Returns the number of qubits in the base constraint.
 *
 * @return The number of qubits.
 */
int BaseConstraint::getNQubits() {
    return n_qubits;
};

/**
 * Calculates the initialization of the PartialMatrix.
 * This function sets the matrix to 0 if the cover is false at that position,
 * and calculates the squared norm of the matrix.
 */
void PartialMatrix::calculateInitialization () {
    matrix = cover.select(matrix, 0);
    array = matrix.array();
    squared_norm = std::abs(Utils::traceConjugateProduct(matrix, matrix));
}

/**
 * @brief Default constructor for the PartialMatrix class.
 */
PartialMatrix::PartialMatrix() {};

/**
 * @brief Constructs a PartialMatrix object.
 * 
 * This constructor initializes a PartialMatrix object with the given parameters.
 * 
 * @param matrix The Eigen::MatrixXcd representing the matrix.
 * @param cover The BoolMatrix representing the cover.
 * @param name The name of the PartialMatrix.
 */
PartialMatrix::PartialMatrix(Eigen::MatrixXcd matrix, 
                                     BoolMatrix cover, std::string name) 
                                     : matrix(matrix), cover(cover), name(name), BaseConstraint(std::log2(matrix.rows())) {
    calculateInitialization();
};

/**
 * @brief Constructs a PartialMatrix object.
 * 
 * @param circuit The GateCircuit object.
 * @param ch The CircuitHelper object.
 * @param name The name of the PartialMatrix.
 */
PartialMatrix::PartialMatrix(GateCircuit& circuit, CircuitHelper& ch, std::string name) : name(name), BaseConstraint(circuit.nb_qbs) {
    matrix = circuit.toMatrix();
    cover = BoolMatrix::Constant(pow(2, circuit.nb_qbs), pow(2, circuit.nb_qbs), true);
    calculateInitialization();
};

/**
 * @brief Copy constructor for PartialMatrix class.
 * 
 * This constructor creates a new PartialMatrix object by copying the contents of another PartialMatrix object.
 * It initializes the matrix, cover, and name members by copying the corresponding members from the other object.
 * It also calculates the initialization of the new object.
 * 
 * @param other The PartialMatrix object to be copied.
 */
PartialMatrix::PartialMatrix(PartialMatrix const& other) : BaseConstraint(other.n_qubits) {
    matrix = other.matrix;
    cover = other.cover;
    name = other.name;
    calculateInitialization();
}

/**
 * @brief Constructs a PartialMatrix object from a given file.
 * 
 * @param filename The name of the file containing the matrix data.
 */
PartialMatrix::PartialMatrix(std::string filename) {
    std::ifstream file;
    file.open(filename);
    file >> name; 
    file >> n_qubits;
    matrix = Eigen::MatrixXcd::Constant(std::pow(2, n_qubits), std::pow(2, n_qubits), 0);
    for (int i = 0; i < std::pow(2, n_qubits); i++) {
        for (int j = 0; j < std::pow(2, n_qubits); j++)  {
            file >> matrix(i, j);
        }
    }

    cover = BoolMatrix::Constant(std::pow(2, n_qubits), std::pow(2, n_qubits), false);
    for (int i = 0; i < std::pow(2, n_qubits); i++) {
        for (int j = 0; j < std::pow(2, n_qubits); j++)  {
            file >> cover(i, j);
        }
    }
    file.close();
    calculateInitialization();
};

/**
 * @brief Get the matrix stored in the PartialMatrix object.
 * 
 * @return Eigen::MatrixXcd The matrix stored in the PartialMatrix object.
 */
Eigen::MatrixXcd PartialMatrix::getMatrix() {
    return matrix;
}

/**
 * @brief Get the Eigen::ArrayXXcd object.
 * 
 * @return Eigen::ArrayXXcd The Eigen array of complex doubles.
 */
Eigen::ArrayXXcd PartialMatrix::getArray() {
    return array;
}

/**
 * @brief Returns the cover matrix.
 * 
 * @return The cover matrix.
 */
BoolMatrix PartialMatrix::getCover() {
    return cover;
}

/**
 * Returns the number of constraints in the PartialMatrix.
 *
 * @return The number of constraints.
 */
int PartialMatrix::nConstraints() {
    return cover.count();
}

/**
 * Writes the PartialMatrix object to a file.
 * 
 * @param filename The name of the file to write to.
 */
void PartialMatrix::write(std::string filename) {
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n");
    std::ofstream file;
    file.open(filename);
    file << name << std::endl;
    file << getNQubits() << std::endl;
    file << getMatrix().format(CSVFormat) << std::endl;
    file << getCover().format(CSVFormat) << std::endl;
    file.close();
}

/**
 * @brief Constructs a QubitIndependentPartialMatrix object.
 * 
 * @param matrix The original PartialMatrix object.
 * @param use_independent_qbs Flag indicating whether to use independent qubits.
 * @param use_inverse Flag indicating whether to use the inverse of the matrix.
 */
QubitIndependentPartialMatrix::QubitIndependentPartialMatrix(PartialMatrix& matrix, bool use_independent_qbs, bool use_inverse): 
                use_independent_qbs(use_independent_qbs), use_inverse(use_inverse) {
    original = PartialMatrix(matrix);
    calculateQubitIndependent();
}

/**
 * @brief Constructs a QubitIndependentPartialMatrix object.
 * 
 * @param matrix The input matrix.
 * @param cover The boolean matrix representing the qubit coverage.
 * @param name The name of the partial matrix.
 * @param use_independent_qbs Flag indicating whether to use independent qubits.
 * @param use_inverse Flag indicating whether to use the inverse of the matrix.
 */
QubitIndependentPartialMatrix::QubitIndependentPartialMatrix(Eigen::MatrixXcd matrix, BoolMatrix cover, std::string name, 
                bool use_independent_qbs, bool use_inverse): use_inverse(use_inverse), use_independent_qbs(use_independent_qbs) {
    original = PartialMatrix(matrix, cover, name);
    calculateQubitIndependent();
}
/**
 * @brief Constructs a QubitIndependentPartialMatrix object.
 * 
 * @param filename The filename of the partial matrix.
 * @param use_independent_qbs Flag indicating whether to use independent qubits.
 * @param use_inverse Flag indicating whether to use the inverse of the partial matrix.
 */
QubitIndependentPartialMatrix::QubitIndependentPartialMatrix(std::string filename, bool use_independent_qbs, bool use_inverse) : 
                use_inverse(use_inverse), use_independent_qbs(use_independent_qbs) {
    original = PartialMatrix(filename);
    calculateQubitIndependent();
}

/**
 * @brief Copy constructor for QubitIndependentPartialMatrix.
 * 
 * This constructor creates a new QubitIndependentPartialMatrix object by copying the contents of another QubitIndependentPartialMatrix object.
 * 
 * @param other The QubitIndependentPartialMatrix object to be copied.
 */
QubitIndependentPartialMatrix::QubitIndependentPartialMatrix(QubitIndependentPartialMatrix const& other) : 
        use_inverse(other.use_inverse), use_independent_qbs(other.use_independent_qbs) {
    original = PartialMatrix(other.original);
    calculateQubitIndependent();
}

/**
 * @brief Constructs a QubitIndependentPartialMatrix object.
 * 
 * @param circuit The GateCircuit object.
 * @param ch The CircuitHelper object.
 * @param name The name of the partial matrix.
 * @param use_independent_qbs Flag indicating whether to use independent qubits.
 * @param use_inverse Flag indicating whether to use the inverse of the partial matrix.
 */
QubitIndependentPartialMatrix::QubitIndependentPartialMatrix(GateCircuit& circuit, CircuitHelper& ch, std::string name, 
            bool use_independent_qbs, bool use_inverse) : use_inverse(use_inverse), use_independent_qbs(use_independent_qbs) {
    original = PartialMatrix(circuit, ch, name);
    calculateQubitIndependent();
}

/**
 * Adds a PartialMatrix to the QubitIndependentPartialMatrix.
 * 
 * @param matrix The PartialMatrix to be added.
 * @return True if the matrix was successfully added, false otherwise.
 */
bool QubitIndependentPartialMatrix::addMatrix(PartialMatrix& matrix) {
    for (std::shared_ptr<PartialMatrix> matr : matrices) {
        if ((matr->matrix - matrix.matrix).norm() < 1e-6 && (matr->cover.array() != matrix.cover.array()).count() == 0) {
            return false;
        }
    }
    return true;
}

/**
 * Calculates the qubit-independent partial matrices.
 * This function iterates through all possible permutations of qubits and generates partial matrices based on the permutation.
 * The generated partial matrices are added to the matrices vector if they meet certain conditions.
 * If use_independent_qbs is true, the partial matrices are added regardless of the size of the matrices vector.
 * If use_inverse is true, the inverse of each generated partial matrix is also added to the matrices vector.
 */
void QubitIndependentPartialMatrix::calculateQubitIndependent() {
    std::vector<int> qbs = {};
    qbs_info = {};
    matrices = {};
    for (int i = 0; i < original.getNQubits(); i++) {
        qbs.push_back(i);
    }

    do {
        Eigen::MatrixXcd new_matrix = original.matrix;
        BoolMatrix new_cover = original.cover;
        std::vector<int> change = Utils::getRowChangeQubitIndependent(qbs, getNQubits());
        for (int i = 0; i < pow(2, getNQubits()); i++) {
            new_matrix.row(i) = original.matrix.row(change[i]);
            new_cover.row(i) = original.cover.row(change[i]);
        }
        Eigen::MatrixXcd new_new_matrix = new_matrix;
        BoolMatrix new_new_cover = new_cover;
        for (int i = 0; i < pow(2, getNQubits()); i++) {
            new_new_matrix.col(i) = new_matrix.col(change[i]);
            new_new_cover.col(i) = new_cover.col(change[i]);
        }
        PartialMatrix matr = PartialMatrix(new_new_matrix, new_new_cover);
        
        if (addMatrix(matr) && (use_independent_qbs || matrices.size() == 0)) {
            std::shared_ptr<PartialMatrix> ptr_matr = std::make_shared<PartialMatrix>(matr);
            matrices.push_back(ptr_matr);
            qbs_info.push_back(qbs);
            inverse_info.push_back(false);
        }
        Eigen::MatrixXcd inverse_matr = new_new_matrix.conjugate();
        BoolMatrix inverse_cover = new_new_cover.transpose();
        PartialMatrix matr2 = PartialMatrix(inverse_matr, inverse_cover);
        if (addMatrix(matr2) && use_inverse) {
            std::shared_ptr<PartialMatrix> ptr_matr2 = std::make_shared<PartialMatrix>(matr2);
            matrices.push_back(ptr_matr2);
            qbs_info.push_back(qbs);
            inverse_info.push_back(true);
        }
        
    } while (std::next_permutation(qbs.begin(), qbs.end()));
}

/**
 * Returns the number of qubits in the original matrix.
 *
 * @return The number of qubits in the original matrix.
 */
int QubitIndependentPartialMatrix::getNQubits() {
    return original.n_qubits;
};

/**
 * @brief Clones the QubitIndependentPartialMatrix object.
 * 
 * @return A shared pointer to the cloned QubitIndependentPartialMatrix object.
 */
std::shared_ptr<QubitIndependentPartialMatrix> QubitIndependentPartialMatrix::clone() {
    return std::make_shared<QubitIndependentPartialMatrix>(QubitIndependentPartialMatrix(original, use_independent_qbs, use_inverse));
}