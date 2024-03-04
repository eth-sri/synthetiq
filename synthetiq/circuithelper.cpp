#include "circuithelper.h"
#include "utils.h"
#include <sstream>
#include <string>
#include <iterator>
#include <filesystem>

/**
 * Checks if a gate with the given name and acting qubits is already present in the readable gates. This is used to avoid the same gate being read and stored multiple times.
 * 
 * @param name The name of the gate.
 * @param acting_qubits The vector of acting qubits.
 * @return True if the gate is already present, false otherwise.
 */
bool CircuitHelper::isAlreadyPresent(std::string name, std::vector<int> acting_qubits) {
    std::shared_ptr<Gate> gate = Gate::findCorrectGate(readable_gates, name, acting_qubits);
    if (gate) {
        return true;
    }
    return false;
}

/**
 * @brief Constructs a CircuitHelper object.
 * 
 * @param nb_qbs The number of qubits in the circuit.
 * @param basic_gate_folder The folder path containing the basic gate files.
 * @param composite_gate_folder The folder path containing the composite gate files.
 * @param read_gate_folder The folder path containing the read gate files.
 */
CircuitHelper::CircuitHelper(int nb_qbs, std::string basic_gate_folder, std::string composite_gate_folder, std::string read_gate_folder) : nb_qbs(nb_qbs), 
    basic_gate_folder(basic_gate_folder), composite_gate_folder(composite_gate_folder), read_gate_folder(read_gate_folder) {
    id_gate = std::make_shared<BasicGate>(BasicGate("id", Eigen::MatrixXcd::Identity(pow(2, nb_qbs), pow(2, nb_qbs)), {0}, 0.0));
    all_gates.push_back(id_gate);
    readable_gates.push_back(id_gate);
    max_cost_basic = 0;
    if (std::filesystem::is_directory(basic_gate_folder)) {
        readBasicGateFolder(basic_gate_folder);
    }
    
    max_cost_all = max_cost_basic;
    if (std::filesystem::is_directory(composite_gate_folder)) {
        readCompositeGateFolder(composite_gate_folder);
    }

    if (std::filesystem::is_directory(read_gate_folder)) {
        readCompositeGateFolder(read_gate_folder, true);
    }
}

/**
 * Reads the basic gate folder and populates the basic_gates vector with the gates found.
 * 
 * @param folder The path to the folder containing the basic gates.
 */
void CircuitHelper::readBasicGateFolder(std::string folder) {
    for (const auto & file : std::filesystem::directory_iterator(folder)) {
        BasicGate gate = BasicGate(file.path().string());
        if (gate.nb_qbs > nb_qbs) {
            continue;
        }

        Eigen::MatrixXcd resized_matrix = Utils::increaseQubits(gate.matrix, nb_qbs);
        BasicGate gate_correct_qbs = BasicGate(gate.name, resized_matrix, gate.acting_qubits, gate.cost);
        if (gate.cost > max_cost_basic) {
            max_cost_basic = gate.cost;
        }
        std::vector<int> qbs = {};
        for (int i = 0; i < nb_qbs; i++) {
            qbs.push_back(i);
        }
        basic_gates_by_name.push_back(std::vector<std::shared_ptr<Gate>>());
        do {
            std::vector<int> new_acting_qbs = {};
            for (int i = 0; i < gate_correct_qbs.acting_qubits.size(); i++) {
                new_acting_qbs.push_back(qbs[gate_correct_qbs.acting_qubits[i]]);
            }
            Eigen::MatrixXcd new_matrix = Utils::changeQubits(gate_correct_qbs.matrix, qbs);
            if (!isAlreadyPresent(gate_correct_qbs.name, new_acting_qbs)) {
                BasicGate new_gate = BasicGate(gate_correct_qbs.name, new_matrix, new_acting_qbs, gate_correct_qbs.cost);
                std::shared_ptr<Gate> gate_ptr = std::make_shared<BasicGate>(new_gate);
                basic_gates_by_name[basic_gates_by_name.size() - 1].push_back(gate_ptr);
                basic_gates.push_back(gate_ptr);
                all_gates.push_back(gate_ptr);
                readable_gates.push_back(gate_ptr);
            }
        } while (std::next_permutation(qbs.begin(), qbs.end()));
    }
    
}

/**
 * Reads the composite gate folder and processes each file.
 * 
 * @param folder The path to the folder containing the composite gate files.
 * @param read_folder Flag indicating whether to read the entire folder or not.
 */
void CircuitHelper::readCompositeGateFolder(std::string folder, bool read_folder) {
    for (const auto & file : std::filesystem::directory_iterator(folder)) {
        CompositeGate gate = CompositeGate(file.path().string(), basic_gates);
        if (gate.nb_qbs > nb_qbs) {
            continue;
        }
        if (gate.cost > max_cost_all && !read_folder) {
            max_cost_all = gate.cost;
        }
        std::vector<int> qbs = {};
        for (int i = 0; i < nb_qbs; i++) {
            qbs.push_back(i);
        }
        if (!read_folder) {
            composite_gates_by_name.push_back(std::vector<std::shared_ptr<Gate>>());
        }
        do {
            std::vector<int> new_acting_qbs = {};
            for (int i = 0; i < gate.acting_qubits.size(); i++) {
                new_acting_qbs.push_back(qbs[gate.acting_qubits[0]]);
            }
            if (!isAlreadyPresent(gate.name, new_acting_qbs)) {
                CompositeGate new_gate = CompositeGate(gate);
                new_gate.changeQubitsDecomposition(qbs, basic_gates);
                std::shared_ptr<Gate> gate_ptr = std::make_shared<CompositeGate>(new_gate);
                
                if (!read_folder) {
                    composite_gates_by_name[composite_gates_by_name.size() - 1].push_back(gate_ptr);
                    composite_gates.push_back(gate_ptr);
                    all_gates.push_back(gate_ptr);
                }
                readable_gates.push_back(gate_ptr);
            }
        } while (std::next_permutation(qbs.begin(), qbs.end()));
    }
}

/**
 * @brief Inverts a gate.
 * 
 * @param gate The gate to invert.
 */
std::shared_ptr<Gate> CircuitHelper::invertGate(std::shared_ptr<Gate> gate) {
    for (std::shared_ptr<Gate> other_gate : all_gates) {
        if (gate->matrix == other_gate->matrix.conjugate()) {
            return other_gate;
        }
    }
    return nullptr;
}