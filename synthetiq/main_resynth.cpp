#include "algo.h"
#include "circuit.h"
#include "resynthesis.h"

/**
 * @brief Applies the simplification pass of the algorithm
 * 
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line arguments.
 * @return int The exit status of the program.
 */
int main(int argc, char* argv[]) {
    Parser parser = Parser();
    parser.parse(argc, argv);
    CircuitHelper ch = CircuitHelper(1, parser.base_gate_folder + parser.gate_set, parser.base_gate_folder + parser.composite_gate_folder);
    GateCircuit circuit_obj = GateCircuit(ch);
    circuit_obj.readFromInput(parser.base_input_folder + parser.input_name);
    ch = *circuit_obj.getCircuitHelper();
    std::shared_ptr<GateCircuit> circuit = std::make_shared<GateCircuit>(circuit_obj);
    Resynthesize resynth = Resynthesize(12, true);
    int gates_before = circuit->getNbNonIdGates();
    int tcount_before = circuit->getCount({"t", "tdg"});
    int tdepth_before = circuit->getDepth({"t", "tdg"});
    resynth.run(circuit, ch);
    int gates_after = circuit->getNbNonIdGates();
    int tcount_after = circuit->getCount({"t", "tdg"});
    int tdepth_after = circuit->getDepth({"t", "tdg"});
    std::cout << "Gates: " << gates_before << " -> " << gates_after << std::endl;
    std::cout << "T-count: " << tcount_before << " -> " << tcount_after << std::endl;
    std::cout << "T-depth: " << tdepth_before << " -> " << tdepth_after << std::endl;
}
