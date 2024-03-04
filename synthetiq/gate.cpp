#include "gate.h"
#include "utils.h"
#include <sstream>
#include <string>
#include <iterator>
#include <filesystem>

/**
 * Default constructor for the Gate class.
 */
Gate::Gate() {

}

/**
 * @brief Constructs a Gate object.
 * 
 * @param name The name of the gate.
 * @param matrix The matrix representation of the gate.
 * @param acting_qubits The qubits on which the gate acts.
 * @param cost The cost associated with the gate.
 */
Gate::Gate(std::string name, Eigen::MatrixXcd matrix, std::vector<int> acting_qubits, double cost) : 
                name(name), matrix(matrix), acting_qubits(acting_qubits), nb_qbs((int) std::log2(matrix.rows())), cost(cost) {

}

/**
 * @brief Converts the Gate object to a string representation.
 * 
 * @return The string representation of the Gate object.
 */
std::string Gate::print() {
    std::string out = name;
    for (int i = 0; i < acting_qubits.size(); i++) {
        out += " " + std::to_string(acting_qubits[i]);
    }
    return out;
}

/**
 * Compares the current Gate object with another Gate object for equality.
 * Two gates are considered equal if they have the same name and the same acting qubits.
 * 
 * @param other_gate The Gate object to compare with.
 * @return True if the gates are equal, false otherwise.
 */
bool Gate::equals(Gate& other_gate) {
    if (other_gate.name != name || other_gate.acting_qubits.size() != acting_qubits.size()) {
        return false;
    }
    for (int i = 0; i < acting_qubits.size(); i++) {
        if (acting_qubits[i] != other_gate.acting_qubits[i]) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Constructs a BasicGate object from a file. A basic gate is a gate that is not decomposed into other gates.
 * 
 * @param filename The name of the file containing the gate information.
 */
BasicGate::BasicGate(std::string filename) {
    std::ifstream file;
    file.open(filename);
    file >> name; 
    file >> nb_qbs;
    file >> cost;
    acting_qubits = {};
    std::string line;
    file.ignore();
    std::getline(file, line);
    std::istringstream string_stream(line);
    int current_int;
    while (string_stream >> current_int) {
        acting_qubits.push_back(current_int);
    }

    matrix = Eigen::MatrixXcd::Constant(std::pow(2, nb_qbs), std::pow(2, nb_qbs), 0);
    for (int i = 0; i < std::pow(2, nb_qbs); i++) {
        for (int j = 0; j < std::pow(2, nb_qbs); j++)  {
            file >> matrix(i, j);
        }
    }
    file.close();
}

/**
 * @brief Constructor for the BasicGate class.
 * 
 * @param name The name of the gate.
 * @param matrix The matrix representation of the gate.
 * @param acting_qubits The qubits on which the gate acts.
 * @param cost The cost associated with the gate.
 */
BasicGate::BasicGate(std::string name, Eigen::MatrixXcd matrix, std::vector<int> acting_qubits, double cost) :
    Gate(name, matrix, acting_qubits, cost) {
}

/**
 * @brief Checks if the gate is a basic gate.
 * 
 * @return true if the gate is a basic gate, false otherwise.
 */
bool BasicGate::isBasicGate() {
    return true;
}

/**
 * @brief A vector of shared pointers to Gate objects.
 */
std::vector<std::shared_ptr<Gate>> BasicGate::decomposeInBasicGates() {
    return {nullptr};
}

/**
 * Finds the correct gate from a list of allowed gates based on the gate name and acting qubits.
 *
 * @param allowed_gates A vector of shared pointers to Gate objects representing the allowed gates.
 * @param gate_name The name of the gate to find.
 * @param acting_qbs_gate A vector of integers representing the acting qubits of the gate to find.
 * @return A shared pointer to the correct Gate object if found, nullptr otherwise.
 */
std::shared_ptr<Gate> Gate::findCorrectGate(std::vector<std::shared_ptr<Gate>> allowed_gates, std::string gate_name, std::vector<int> acting_qbs_gate) {
    for (std::shared_ptr<Gate> allowed_gate : allowed_gates) {
        if (gate_name == allowed_gate->name && acting_qbs_gate.size() == allowed_gate->acting_qubits.size()) {
            bool good_qbs = true;
            for (int i = 0; i < acting_qbs_gate.size(); i++) {
                if (acting_qbs_gate[i] != allowed_gate->acting_qubits[i]) {
                    good_qbs = false;
                    break;
                }
            }
            if (good_qbs) {
                return allowed_gate;
            }
        }
    }
    return nullptr;
}

/**
 * @brief Constructs a CompositeGate object. A composite gate is a gate that is decomposed into other gates.
 * 
 * @param filename The name of the file to read the gate information from.
 * @param allowed_gates A vector of shared pointers to allowed gate objects.
 */
CompositeGate::CompositeGate(std::string filename, std::vector<std::shared_ptr<Gate>> allowed_gates) {
    if (filename.substr(filename.find_last_of(".") + 1) == "qasm") {
        readFromFileQasm(filename, allowed_gates);
    } else {
        readFromFileTxT(filename, allowed_gates);
    }
}

/**
 * @brief Constructs a CompositeGate object.
 * 
 * @param name The name of the gate.
 * @param matrix The matrix representation of the gate.
 * @param acting_qubits The qubits on which the gate acts.
 * @param cost The cost associated with the gate.
 * @param decomposition The decomposition of the gate into a sequence of simpler gates.
 */
CompositeGate::CompositeGate(std::string name, Eigen::MatrixXcd matrix, std::vector<int> acting_qubits, double cost, 
                             std::vector<std::shared_ptr<Gate>> decomposition) : Gate(name, matrix, acting_qubits, cost), decomposition(decomposition) {
}

/**
 * Reads a composite gate from a text file.
 * 
 * @param filename The name of the file to read from.
 * @param allowed_gates The vector of allowed gates.
 */
void CompositeGate::readFromFileTxT(std::string filename, std::vector<std::shared_ptr<Gate>> allowed_gates) {
    std::ifstream file;
    file.open(filename);
    file >> name;
    file >> nb_qbs;
    if (allowed_gates.size() == 0 || allowed_gates[0]->nb_qbs < nb_qbs) {
        return;
    }
    cost = 0;
    decomposition = {};
    nb_qbs = allowed_gates[0]->nb_qbs;
    matrix = Eigen::MatrixXcd::Identity(pow(2, nb_qbs), pow(2, nb_qbs));
    acting_qubits = {};
    std::string line;
    file.ignore();
    std::getline(file, line);
    std::istringstream string_stream(line);
    int current_int;
    while (string_stream >> current_int) {
        acting_qubits.push_back(current_int);
    }

    while (std::getline(file, line)) {
        bool found_line = false;
        if (line.size() == 0) {
            continue;
        }
        std::string gate_name = Utils::readNameQasmLine(line);
        std::vector<int> acting_qubits_gate = Utils::readActingQubitsQasmLine(line);
        std::shared_ptr<Gate> gate_to_add = findCorrectGate(allowed_gates, gate_name, acting_qubits_gate);
        if (gate_to_add) {
                decomposition.push_back(gate_to_add);
                matrix = gate_to_add->matrix * matrix;
                cost += gate_to_add->cost;
        } else {
            throw std::invalid_argument("Gate in file cannot be constructed. File: " + filename + ". Problem with line: " + line);
        }
    }
}

/**
 * Reads a quantum circuit from a QASM file and constructs the composite gate.
 * 
 * @param filename The path to the QASM file.
 * @param allowed_gates The vector of allowed gates for the composite gate.
 */
void CompositeGate::readFromFileQasm(std::string filename, std::vector<std::shared_ptr<Gate>> allowed_gates) {
    std::ifstream file;
    file.open(filename);
    name = filename.substr(filename.find_last_of("/"), filename.find_last_of("."));
    std::string line;
    file.ignore();
    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);

    nb_qbs = Utils::readNbQbsQasmLine(line);

    if (allowed_gates.size() == 0 || allowed_gates[0]->nb_qbs < nb_qbs) {
        return;
    }

    acting_qubits = {};
    for (int i = 0; i < nb_qbs; i++) {
        acting_qubits.push_back(i);
    }
    
    cost = 0;
    decomposition = {};
    nb_qbs = allowed_gates[0]->nb_qbs;
    matrix = Eigen::MatrixXcd::Identity(pow(2, nb_qbs), pow(2, nb_qbs));

    while (std::getline(file, line)) {
        bool found_line = false;
        if (line.size() == 0) {
            continue;
        }
        std::string gate_name = Utils::readNameQasmLine(line);
        std::vector<int> acting_qubits_gate = Utils::readActingQubitsQasmLine(line);
        std::shared_ptr<Gate> gate_to_add = findCorrectGate(allowed_gates, gate_name, acting_qubits_gate);
        if (gate_to_add) {
                decomposition.push_back(gate_to_add);
                matrix = gate_to_add->matrix * matrix;
                cost += gate_to_add->cost;
        } else {
            throw std::invalid_argument("Gate in file cannot be constructed. File: " + filename + ". Problem with line: " + line);
        }
    }
}

/**
 * @brief Checks if the composite gate is a basic gate.
 * 
 * @return true if the composite gate is a basic gate, false otherwise.
 */
bool CompositeGate::isBasicGate() {
    return false;
}

/**
 * @brief A vector of shared pointers to Gate objects.
 */
std::vector<std::shared_ptr<Gate>> CompositeGate::decomposeInBasicGates() {
    return decomposition;
}

/**
 * @brief Changes the qubits in the decomposition of the composite gate.
 * 
 * This function modifies the acting qubits and decomposition of the composite gate
 * based on the provided qubits change and allowed gates. It updates the matrix
 * of the composite gate accordingly.
 * 
 * @param qbs_change The vector of qubits change, specifying the new qubit indices
 * for the acting qubits.
 * @param allowed_gates The vector of allowed gates that can be used for decomposition.
 * @throws std::invalid_argument If a gate in the decomposition cannot be qubit changed.
 */
void CompositeGate::changeQubitsDecomposition(std::vector<int> qbs_change, std::vector<std::shared_ptr<Gate>> allowed_gates) {
    matrix = Eigen::MatrixXcd::Identity(pow(2, nb_qbs), pow(2, nb_qbs));
    for (int i = 0; i < acting_qubits.size(); i++) {
        acting_qubits[i] = qbs_change[acting_qubits[i]];
    }
    for (int i = 0; i < decomposition.size(); i++) {
        std::vector<int> new_acting_qbs = {};
        for (int j = 0; j < decomposition[i]->acting_qubits.size(); j++) {
            new_acting_qbs.push_back(qbs_change[decomposition[i]->acting_qubits[j]]);
        }
        std::shared_ptr<Gate> new_decomp_gate = findCorrectGate(allowed_gates, decomposition[i]->name, new_acting_qbs);
        if (new_decomp_gate) {
            decomposition[i] = new_decomp_gate;
        } else {
            throw std::invalid_argument("Gate in file cannot be qubit changed. Name: " + name);
        }
        matrix = new_decomp_gate->matrix * matrix;
    }
}