#include "circuit.h"
#include "utils.h"
#include <sstream>
#include <string>
#include <iterator>
#include <filesystem>


/**
 * @brief Constructs a GateCircuit object.
 * 
 * @param circ_helper The CircuitHelper object.
 * @param matrix_computer_type The type of MatrixComputer to use.
 */
GateCircuit::GateCircuit(CircuitHelper& circ_helper, MatrixComputerType matrix_computer_type) : matrix_computer_type(matrix_computer_type) {
    ch = std::make_shared<CircuitHelper>(circ_helper);
}

/**
 * @brief Constructs a GateCircuit object.
 * 
 * @param nb_g The number of gates in the circuit.
 * @param nb_q The number of qubits in the circuit.
 * @param circ_helper The CircuitHelper object.
 * @param matrix_computer_type The type of MatrixComputer to use.
 */
GateCircuit::GateCircuit(int nb_g, int nb_q, CircuitHelper& circ_helper, MatrixComputerType matrix_computer_type) : matrix_computer_type(matrix_computer_type) {
    nb_qbs = nb_q;
    ch = std::make_shared<CircuitHelper>(circ_helper);
    list_gates = {};
    for(int cur_g = 0; cur_g < nb_g; cur_g++){
        list_gates.push_back(ch->id_gate);
    }
    initializeMatrixComputer();
    matrixComputer->calculateMatrix(list_gates);
    calculateCost();
}

/**
 * @brief Copy constructor for GateCircuit class.
 * 
 * @param other The GateCircuit object to be copied.
 */
GateCircuit::GateCircuit(GateCircuit const& other): nb_qbs(other.nb_qbs), ch(other.ch) {
    list_gates = other.list_gates;
    old_gate = other.old_gate;
    new_gate = other.new_gate;
    nb_non_id_gates = other.nb_non_id_gates;
    cost = other.cost;
    pos_mutation = other.pos_mutation;
    matrix_computer_type = other.matrix_computer_type;
    initializeMatrixComputer();
    matrixComputer->calculateMatrix(list_gates);
}

/**
 * @brief Clones the GateCircuit object.
 * 
 * @return A shared pointer to the cloned GateCircuit object.
 */
std::shared_ptr<GateCircuit> GateCircuit::clone() {
    return std::make_shared<GateCircuit>(*this);
}

/**
 * @brief Constructs a GateCircuit object.
 * 
 * @param gates The vector of shared pointers to Gate objects representing the gates in the circuit.
 * @param nb_q The number of qubits in the circuit.
 * @param circ_helper The CircuitHelper object used for circuit operations.
 * @param matrix_computer_type The type of MatrixComputer used for matrix calculations.
 */
GateCircuit::GateCircuit(std::vector<std::shared_ptr<Gate>> gates, int nb_q, CircuitHelper& circ_helper, MatrixComputerType matrix_computer_type): 
        nb_qbs(nb_q), list_gates(gates), matrix_computer_type(matrix_computer_type) {
    ch = std::make_shared<CircuitHelper>(circ_helper);
    initializeMatrixComputer();
    matrixComputer->calculateMatrix(gates);
    calculateCost();
}

/**
 * Returns the number of elements in the GateCircuit.
 *
 * @return The number of elements in the GateCircuit.
 */
int GateCircuit::nbElements(){
    return list_gates.size();
}

/**
 * Reads a gate from a QASM input line and adds it to the gate circuit.
 * 
 * @param line The QASM input line containing the gate information.
 * @throws std::invalid_argument if the gate in the line cannot be constructed.
 */
void GateCircuit::readGateFromQasmInputLine(std::string line) {
    if(line.size() == 0) //empty line
        return;
    std::string gate_name = Utils::readNameQasmLine(line);
    std::vector<int> acting_qubits_gate = Utils::readActingQubitsQasmLine(line);
    std::shared_ptr<Gate> gate_to_add = Gate::findCorrectGate(ch->readable_gates, gate_name, acting_qubits_gate);
    if (gate_to_add) {
            list_gates.push_back(gate_to_add);
    } else {
        throw std::invalid_argument("Gate in line cannot be constructed. Line: " + line);
    }
}

/**
 * Reads a gate circuit from an input file.
 * 
 * @param filename The name of the input file.
 */
void GateCircuit::readFromInput(std::string filename){ //cin stream
    std::ifstream file;
    file.open(filename);
    std::string line;
    std::regex words_regex("[^\\s\\[\\],]+");

    // ignore first two lines
    std::getline(file, line);
    std::getline(file, line);
    //3rd line is nb qb: qreg qubits[nb]
    std::getline(file, line);
    auto words_begin = std::sregex_iterator(line.begin(), line.end(), words_regex);
    words_begin++;
    words_begin++;
    nb_qbs = atoi((*words_begin).str().c_str());
    if (nb_qbs != ch->nb_qbs) {
        ch = std::make_shared<CircuitHelper>(CircuitHelper(nb_qbs, ch->basic_gate_folder, ch->composite_gate_folder, ch->read_gate_folder));
    }
    while (std::getline(file, line)) {
        readGateFromQasmInputLine(line);
    }
    initializeMatrixComputer();
    matrixComputer->calculateMatrix(list_gates);
    calculateCost();
}

/**
 * Reads a gate circuit from a non-Qasm input file.
 * 
 * @param filename The name of the input file.
 * @param nb_l_to_ignore The number of lines to ignore at the beginning of the file.
 * @param nb_qbs_in_circ The number of qubits in the circuit.
 */
void GateCircuit::readFromInputNonQuasm(std::string filename, int nb_l_to_ignore, int nb_qbs_in_circ){
    std::ifstream file;
    file.open(filename);
    std::string line;
    nb_qbs = nb_qbs_in_circ;
    if (nb_qbs != ch->nb_qbs) {
        ch = std::make_shared<CircuitHelper>(CircuitHelper(nb_qbs, ch->basic_gate_folder, ch->composite_gate_folder, ch->read_gate_folder));
    }

    for(int i = 0; i < nb_l_to_ignore; i++)
        std::getline(file, line);

    while (std::getline(file, line)){
        readGateFromQasmInputLine(line);
    }
    initializeMatrixComputer();
    matrixComputer->calculateMatrix(list_gates);
    calculateCost();
}

/**
 * Reads the gate circuit from an input file stream.
 * 
 * @param file The input file stream to read from.
 */
void GateCircuit::readFromIfStream(std::ifstream& file){
    std::string line;
    std::regex words_regex("[^\\s\\[\\],]+");

    std::getline(file, line);
    std::istringstream buffer(line);
    std::vector<std::string> splitOnSpace((std::istream_iterator<std::string>(buffer)), 
                                 std::istream_iterator<std::string>());
    nb_qbs = atoi(splitOnSpace[2].c_str());
    int nb_gates = atoi(splitOnSpace[6].c_str());
    if (nb_qbs != ch->nb_qbs) {
        ch = std::make_shared<CircuitHelper>(CircuitHelper(nb_qbs, ch->basic_gate_folder, ch->composite_gate_folder, ch->read_gate_folder));
    }

    for(int i = 0; i < nb_gates; i++){
        std::getline(file, line);
        readGateFromQasmInputLine(line);
    }
    initializeMatrixComputer();
    matrixComputer->calculateMatrix(list_gates);
    calculateCost();
}

/**
 * @brief Returns a string representation of the GateCircuit in QASM format.
 * 
 * This function generates a QASM representation of the GateCircuit object,
 * including the necessary QASM header and gate definitions.
 * 
 * @return A string containing the QASM representation of the GateCircuit.
 */
std::string GateCircuit::print_qasm(){
    std::ostringstream oss;
    oss << "OPENQASM 2.0;\n" << "include \"qelib1.inc\";\n" << "qreg qubits[" << std::to_string(nb_qbs) << "];\n";
    std::vector<std::string> gate_names = {"h", "p", "pdg", "t", "tdg", "id", "cx", "x", "y", "z"}; 
    for(size_t i = 0; i < list_gates.size(); i++)
        if(list_gates[i]->name != ch->id_gate->name) {
            oss << list_gates[i]->name;
            for (int qb = 0; qb < list_gates[i]->acting_qubits.size(); qb++) {
                oss << " qubits[" << std::to_string(list_gates[i]->acting_qubits[qb]) << "]";
                if (qb != list_gates[i]->acting_qubits.size() - 1) {
                    oss << ",";
                }
            }                
            oss << ";\n";
        }
    return oss.str();
}

/**
 * Initializes the matrix computer based on the type specified in `matrix_computer_type`.
 * If `matrix_computer_type` is Linear, a LinearMatrixComputer is created and assigned to `matrixComputer`.
 * If `matrix_computer_type` is Chunk, a ChunkMatrixComputer is created with the specified number of gates and qubits, and assigned to `matrixComputer`.
 * If `matrix_computer_type` is Binary, a BinaryMatrixComputer is created with the specified number of gates and qubits, and assigned to `matrixComputer`.
 */
void GateCircuit::initializeMatrixComputer() {
    if (matrix_computer_type == Linear) {
        matrixComputer = std::make_shared<LinearMatrixComputer>(LinearMatrixComputer());
    } else if (matrix_computer_type == Chunk) {
        matrixComputer = std::make_shared<ChunkMatrixComputer>(ChunkMatrixComputer(list_gates.size(), nb_qbs));
    } else if (matrix_computer_type == Binary) {
        matrixComputer = std::make_shared<BinaryMatrixComputer>(BinaryMatrixComputer(list_gates.size(), nb_qbs));
    }
}

/**
 * Returns the count of gates in the circuit that match the given gate names.
 *
 * @param gate_names A vector of gate names to search for.
 * @return The count of gates that match the given gate names.
 */
int GateCircuit::getCount(std::vector<std::string> gate_names) {
    int count = 0;
    for (std::shared_ptr<Gate> gate : list_gates) {
        if (std::find(gate_names.begin(), gate_names.end(), gate->name) != gate_names.end()) {
            count += 1;
        }
    }
    return count;
}

/**
 * Calculates the depth of the gate circuit based on the given gate names.
 * The depth of a circuit is defined as the maximum number of layers of gates that need to be applied
 * in order to execute all the gates in the circuit.
 *
 * @param gate_names A vector of gate names for which the depth needs to be calculated.
 * @return The depth of the circuit.
 */
int GateCircuit::getDepth(std::vector<std::string> gate_names) {
    std::vector<int> depths(nb_qbs);
    for (std::shared_ptr<Gate> gate : list_gates) {
        bool find = std::find(gate_names.begin(), gate_names.end(), gate->name) != gate_names.end();
        int max_depth_acting = 0;
        for (int index : gate->acting_qubits) {
            max_depth_acting = std::max(max_depth_acting, depths[index]);
        }

        for (int index : gate->acting_qubits) {
            if (find) {
                depths[index] = max_depth_acting + 1;
            } else {
                depths[index] = max_depth_acting;
            }
        }
    }
    return *std::max_element(depths.begin(), depths.end());
}

/**
 * @brief Replaces the gate at the specified position with a new gate.
 * 
 * It updates the cost of the circuit based on the new gate and the gate being replaced.
 * The new gate is assigned to the specified position in the list of gates.
 * If the matrix computer is currently calculating the matrix (this is turned of in the simplification pass for more efficiency), it updates the matrix with the new gate.
 * 
 * @param position The position of the gate to be replaced.
 * @param gate The new gate to be placed at the specified position.
 */
void GateCircuit::placeGateAt(int position, std::shared_ptr<Gate> gate){ //replaces gate at pos with new gate, assumes qbs of size 2
    updateCost(gate, list_gates[position]);
    new_gate = gate;
    old_gate = list_gates[position];
    list_gates[position] = new_gate;
    if (calculating_matrix_computer) {
        matrixComputer->updateMatrix(position, list_gates);
    }
    
    pos_mutation = position;
}

/**
 * Expands composite gates in the circuit.
 * This function replaces composite gates with their equivalent basic gates.
 * It modifies the list of gates in the circuit and recalculates the cost.
 */
void GateCircuit::expandCompositeGates() {
    stopMatrixComputer();
    std::vector<std::shared_ptr<Gate>> new_list_gates;
    for (int i = 0; i < list_gates.size(); i++) {
        if (list_gates[i]->isBasicGate()) {
            new_list_gates.push_back(list_gates[i]);
        } else {
            for (std::shared_ptr<Gate> gate : list_gates[i]->decomposeInBasicGates()) {
                new_list_gates.push_back(gate);
            }
        }
    }
    list_gates = new_list_gates;
    calculateCost();
    startMatrixComputer();
}

/**
 * Calculates the cost of the gate circuit.
 * The cost is calculated by summing up the costs of all the gates in the circuit.
 * Additionally, it counts the number of non-identity gates in the circuit.
 */
void GateCircuit::calculateCost() {
    cost = 0.0;
    nb_non_id_gates = 0;
    for (std::shared_ptr<Gate> gate : list_gates) {
        cost += gate->cost;
        if (gate->name != ch->id_gate->name) {
            nb_non_id_gates += 1;
        }
    }
}

/**
 * Updates the cost of the gate circuit by subtracting the cost of the old gate and adding the cost of the new gate.
 * Additionally, it adjusts the count of non-identity gates based on whether the old or new gate is an identity gate.
 * 
 * @param new_gate The new gate to be added to the circuit.
 * @param old_gate The old gate to be removed from the circuit.
 */
void GateCircuit::updateCost(std::shared_ptr<Gate> new_gate, std::shared_ptr<Gate> old_gate) {
    cost += new_gate->cost - old_gate->cost;
    if (new_gate->name == ch->id_gate->name) {
        nb_non_id_gates -= 1;
    }
    if (old_gate->name == ch->id_gate->name) {
        nb_non_id_gates += 1;
    }
}

/**
 * Returns the number of non-identity gates in the circuit.
 *
 * @return The number of non-identity gates.
 */
int GateCircuit::getNbNonIdGates() {
    return nb_non_id_gates;
}

/**
 * Returns the cost of the gate circuit.
 *
 * @return The cost of the gate circuit.
 */
double GateCircuit::getCost() {
    return cost;
}

/**
 * @brief Converts the GateCircuit to a matrix representation.
 * 
 * @return The matrix representation of the GateCircuit.
 */
Eigen::MatrixXcd GateCircuit::toMatrix(){
    return matrixComputer->getMatrix();
}

/**
 * Stops the matrix computer for the gate circuit.
 */
void GateCircuit::stopMatrixComputer() {
    calculating_matrix_computer = false;
}

/**
 * Starts the matrix computer to calculate the matrix for the gate circuit.
 */
void GateCircuit::startMatrixComputer() {
    calculating_matrix_computer = true;
    matrixComputer->calculateMatrix(list_gates);
}

/**
 * @brief Returns the list of gates in the circuit.
 * 
 * @return The list of gates in the circuit.
 */
const std::vector<std::shared_ptr<Gate>> GateCircuit::getGates(){
    return list_gates;
}

/**
 * @brief Get the CircuitHelper object.
 * 
 * @return std::shared_ptr<CircuitHelper> The CircuitHelper object.
 */
std::shared_ptr<CircuitHelper> GateCircuit::getCircuitHelper() {
    return ch;
}

/**
 * Calculates the performance normalization constant for the GateCircuit.
 * The performance normalization constant is equal to the number of elements in the circuit.
 *
 * @return The performance normalization constant.
 */
int GateCircuit::performanceNormalizationCst(){
    return nbElements();
}

/**
 * Mutates a gate in the circuit at the specified position.
 * 
 * @param random_helper The random number generator helper.
 * @param position The position of the gate to mutate.
 * @param proportional_prob The proportional probability for composite gates in the circuit, P_comp in the paper.
 * @param proba_id The probability of selecting an identity gate.
 * @return True if the gate remains unchanged after mutation, false otherwise.
 */
bool GateCircuit::mutate(RandomHelper& random_helper, int position, double proportional_prob, double proba_id, double proba_name){
    bool unchanged = false;
    //create potential mutation
    double mutation = random_helper.random01();
    old_gate = list_gates[position];
    pos_mutation = position;
    // create candidate
    if (mutation < proba_id){
        new_gate = ch->id_gate;
    }
    else if (random_helper.random01() < proba_name) {
        double prob_basic_gate = (ch->basic_gates_by_name.size()) / (ch->basic_gates_by_name.size() + proportional_prob * ch->composite_gates_by_name.size());
        if (random_helper.random01() < prob_basic_gate) {
            int random1 = random_helper.randomInt(ch->basic_gates_by_name.size());
            int random2 = random_helper.randomInt(ch->basic_gates_by_name[random1].size());
            new_gate = ch->basic_gates_by_name[random1][random2];
        } else {
            int random1 = random_helper.randomInt(ch->composite_gates_by_name.size());
            int random2 = random_helper.randomInt(ch->composite_gates_by_name[random1].size());
            new_gate = ch->composite_gates_by_name[random1][random2];
        }
    }
    else {
        double prob_basic_gate = (ch->basic_gates.size()) / (ch->basic_gates.size() + proportional_prob * ch->composite_gates.size());
        if (random_helper.random01() < prob_basic_gate) {
            new_gate = ch->basic_gates[random_helper.randomInt(ch->basic_gates.size())];
        } else {
            new_gate = ch->composite_gates[random_helper.randomInt(ch->composite_gates.size())];
        }
    }
    if (old_gate->equals(*new_gate)) {
        unchanged = true;
    } else {
        placeGateAt(position, new_gate);
    }
    return unchanged;
}

/**
 * Mutates the GateCircuit by randomly selecting a position and applying a mutation at that position.
 * 
 * @param random_helper The RandomHelper object used for generating random numbers.
 * @param proportional_prob The probability of applying a proportional mutation.
 * @param proba_id The probability of applying an identity mutation.
 * @return True if the mutation was successful, false otherwise.
 */
bool GateCircuit::mutate(RandomHelper& random_helper, double proportional_prob, double proba_id, double proba_name){
    int position = random_helper.randomInt(nbElements());
    return mutate(random_helper, position, proportional_prob, proba_id, proba_name);
}


/**
 * Undoes the previous mutation by placing the original gate back at the specified position.
 */
void GateCircuit::undoMutation() {
    placeGateAt(pos_mutation, old_gate);
}
/**
 * Applies the mutation again to another GateCircuit.
 * 
 * @param other The GateCircuit to apply the mutation to.
 */
void GateCircuit::applyMutationAgain(GateCircuit& other){
    other.placeGateAt(pos_mutation, new_gate);
}

/**
 * @brief Changes the order of qubits in the gate circuit.
 * 
 * This function rearranges the qubits in the gate circuit according to the given order.
 * It stops the matrix computer, iterates through each gate in the circuit, and updates the acting qubits
 * based on the new order. It then finds the correct gate with the updated acting qubits and places it
 * at the corresponding position in the circuit. Finally, it restarts the matrix computer.
 * 
 * @param qbs_order The new order of qubits.
 */
void GateCircuit::changeQubits(std::vector<int> qbs_order) {
    stopMatrixComputer();
    for (int i = 0; i < nbElements(); i++) {
        std::vector<int> acting_qubits = {};
        for (int j = 0; j < list_gates[i]->acting_qubits.size(); j++) {
            if (list_gates[i]->name != ch->id_gate->name) {
                acting_qubits.push_back(qbs_order[list_gates[i]->acting_qubits[j]]);
            } else {
                acting_qubits.push_back(list_gates[i]->acting_qubits[j]);
            }
        }
        std::shared_ptr<Gate> gate = Gate::findCorrectGate(ch->all_gates, list_gates[i]->name, acting_qubits);
        placeGateAt(i, gate);

    }
    startMatrixComputer();
}

/**
 * Inverts the gate circuit by inverting each gate in reverse order.
 * This method stops the matrix computer, creates a new list of inverted gates,
 * and then places the inverted gates back into the circuit in the original order.
 * Finally, it restarts the matrix computer.
 */
void GateCircuit::invert() {
    stopMatrixComputer();
    std::vector<std::shared_ptr<Gate>> new_list_gates;
    for (int i = nbElements() - 1; i >= 0; i--){
        std::shared_ptr<Gate> inverted_gate = ch->invertGate(list_gates[i]);
        new_list_gates.push_back(inverted_gate);
    }
    for (int i = 0; i < nbElements(); i++) {
        placeGateAt(i, new_list_gates[i]);
    }
    startMatrixComputer();
};

/**
 * Rotates the gates in the circuit.
 * This function reverses the order of the gates in the circuit.
 * It stops the matrix computer, creates a copy of the list of gates,
 * and then places the gates in reverse order in the circuit.
 * Finally, it starts the matrix computer again.
 */
void GateCircuit::rotate() {
    stopMatrixComputer();
    std::vector<std::shared_ptr<Gate>> list_gates_copy = list_gates;
    for (int i = nbElements() - 1; i >= 0; i--){
        placeGateAt(nbElements() - 1 - i, list_gates_copy[i]);
    }
    startMatrixComputer();
};
