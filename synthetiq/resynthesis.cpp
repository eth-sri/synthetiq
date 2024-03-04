#include "resynthesis.h"
#include "utils.h"
#include <sstream>

/**
 * Runs the resynthesis algorithm on the given gate circuit.
 * 
 * @param circuit The gate circuit to be resynthesized.
 * @param ch The CircuitHelper object used for resynthesis.
 * @param second_time Flag indicating if it is the second time running the algorithm.
 */
void Resynthesize::run(std::shared_ptr<GateCircuit> circuit, CircuitHelper& ch, bool second_time) {

    circuit->stopMatrixComputer();
    int gate_index = 0;
    while(gate_index < circuit->nbElements() - 1) {
        int changed = change(circuit, gate_index, ch, second_time);
        if (changed < 0) {
            gate_index += std::max(changed, -gate_index);
        } else if (changed >= 0) {
            gate_index += changed;
        }
    }
    if (!second_time) {
        circuit->rotate();
        run(circuit, ch, true);
        circuit->rotate();
    }
    circuit->startMatrixComputer();
};

/**
 * Calculates the number of gates before the given index that commute with the gate at the given index.
 * 
 * @param index The index of the gate to check.
 * @param circuit A shared pointer to the GateCircuit object.
 * @return The number of gates before the given index that commute with the gate at the given index.
 */
int Resynthesize::NbBeforeCommutes(int index, std::shared_ptr<GateCircuit> circuit) {
    int nb_commutes = 0;
    for (int i = index - 1; i >= 0; i--) {
        if (commutes(circuit->getGates()[index], circuit->getGates()[i])) {
            nb_commutes += 1;
        }
        else {
            break;
        }
    }
    return nb_commutes;
}

/**
 * Determines the priority change between two gates in a circuit.
 * The priority change is based on several criteria, including the number of commutes before each gate,
 * the number of acting qubits, and the gate name.
 * 
 * @param index The index of the first gate.
 * @param circuit The gate circuit.
 * @param second_time Indicates whether it is the second time the priority change is being evaluated.
 * @return True if the priority change is higher, False otherwise.
 */
bool Resynthesize::priorityChange(int index, std::shared_ptr<GateCircuit> circuit, bool second_time) {
    if (watch_depth && second_time) {
        std::shared_ptr<Gate> gate1 = circuit->getGates()[index];
        std::shared_ptr<Gate> gate2 = circuit->getGates()[index + 1];
        int original_depth = circuit->getDepth(depth_gates);
        circuit->placeGateAt(index, gate2);
        circuit->placeGateAt(index + 1, gate1);
        int new_depth = circuit->getDepth(depth_gates);
        circuit->placeGateAt(index, gate1);
        circuit->placeGateAt(index + 1, gate2);
        if (new_depth > original_depth) {
            return false;
        } else if (new_depth < original_depth) {
            return true;
        }
    }
    
    int commutes1 = NbBeforeCommutes(index, circuit);
    int commutes2 = NbBeforeCommutes(index + 1, circuit);
    std::shared_ptr<Gate> gate1 = circuit->getGates()[index];
    std::shared_ptr<Gate> gate2 = circuit->getGates()[index + 1];
    if (commutes1 < commutes2 - 1) {
        return true;
    } else if (commutes1 > commutes2 - 1) {
        return false;
    } else if (gate1->acting_qubits.size() > gate2->acting_qubits.size()) {
        return true;
    } else if (gate1->acting_qubits.size() < gate2->acting_qubits.size()) {
        return false;
    } else if (gate1->name > gate2->name) {
        return true;
    } else if (gate1->name < gate2->name) {
        return false;
    }

    for (int i = 0; i < gate1->acting_qubits.size(); i++) {
        if (gate1->acting_qubits[i] > gate2->acting_qubits[i]) {
            return true;
        } else if (gate1->acting_qubits[i] < gate2->acting_qubits[i]) {
            return false;
        }
    }
    return false;
}

/**
 * Checks if two gates commute with each other.
 * 
 * @param gate1 The first gate.
 * @param gate2 The second gate.
 * @return True if the gates commute, false otherwise.
 */
bool Resynthesize::commutes(std::shared_ptr<Gate> gate1, std::shared_ptr<Gate> gate2) {
    return gate1->matrix * gate2->matrix == gate2->matrix * gate1->matrix;
};

/**
 * @brief Changes the order of gates in a circuit based on certain conditions.
 * 
 * This function is responsible for reordering gates in a circuit based on specific conditions.
 * It takes a shared pointer to a GateCircuit object, the index of the gate to be changed,
 * a CircuitHelper object, and a boolean flag indicating whether it is the second time the function is called.
 * 
 * @param circuit A shared pointer to a GateCircuit object representing the circuit.
 * @param gate_index The index of the gate to be changed.
 * @param ch A reference to a CircuitHelper object.
 * @param second_time A boolean flag indicating whether it is the second time the function is called.
 * @return An integer value indicating how many steps to retrace in the circuit before applying the next change.
 */
int Resynthesize::change(std::shared_ptr<GateCircuit> circuit, int gate_index, CircuitHelper& ch, bool second_time) {
    std::shared_ptr<Gate> gate1 = circuit->getGates()[gate_index];
    std::shared_ptr<Gate> gate2 = circuit->getGates()[gate_index + 1];
    if (gate2->name == ch.id_gate->name) {
        return 1;
    }
    Eigen::MatrixXcd mult = gate2->matrix;
    double cost = gate2->cost;
    if (!second_time) {
        for (int n_extra_gates = 0; n_extra_gates < max_gate_mult - 1; n_extra_gates++) {
            if (gate_index - n_extra_gates < 0) {
                break;
            } else if (circuit->getGates()[gate_index - n_extra_gates]->name == ch.id_gate->name) {
                continue;
            }
            mult = circuit->getGates()[gate_index - n_extra_gates]->matrix * mult;
            cost += circuit->getGates()[gate_index - n_extra_gates]->cost;
            for (int i = 0; i < ch.all_gates.size(); i++) {
                if (cost > ch.all_gates[i]->cost && Utils::matricesEqual(mult, ch.all_gates[i]->matrix)) {
                    circuit->placeGateAt(gate_index - n_extra_gates, ch.all_gates[i]);
                    for (int j = -n_extra_gates + 1; j < 2; j++) {
                        circuit->placeGateAt(gate_index + j, ch.id_gate);
                    }
                    return -1 - n_extra_gates;
                }
            }
        }
    }
    
    if (gate1->name == ch.id_gate->name && gate2->name != ch.id_gate->name) {
        circuit->placeGateAt(gate_index, gate2);
        circuit->placeGateAt(gate_index + 1, gate1);
        return -1;
    }

    if (gate1->name != ch.id_gate->name && gate2->name != ch.id_gate->name && commutes(gate1, gate2)) {
        if (priorityChange(gate_index, circuit, second_time)) {
            circuit->placeGateAt(gate_index, gate2);
            circuit->placeGateAt(gate_index + 1, gate1);
            return -1;
        } 
    }
    return 1;
};
