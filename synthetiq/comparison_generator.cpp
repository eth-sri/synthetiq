#include "circuit.h"
#include "generator.h"
#include "resynthesis.h"
#include "algo.h"

#include <omp.h>
#include <iostream>
#include <vector>
#include <thread>
#include <Eigen/Dense>
#include <filesystem>
#include <map>

int count_files(std::string directory, int nb_qbs, int n_gates, Type type) {
    std::string sub_filename = std::to_string(n_gates) + "_" + std::to_string(type) + "_";
    int count = 0;
    for (const auto & entry : std::filesystem::directory_iterator(directory)) {
        if (entry.path().string().find(sub_filename) != std::string::npos) {
            count += 1;
        }
    }
    return count;
}


/**
 * @brief The main function that generates and saves quantum circuits based on different parameters.
 * 
 * This function iterates over different combinations of gateset, number of qubits, number of gates, and type.
 * For each combination, it generates a random circuit and calculates its optimal T-count and depth using an algorithm.
 * If the generated circuit meets certain criteria, it saves the circuit to a file.
 */
int main(){
    std::vector<int> gates = {10};
    int extra_gates_allowed = 0;
    for (std::string gateset : {"CliffordT"}) {
        for (int n_qbs : {2,3,4}) {
            for (int n_gate : gates) {
                for (Type type : {Full, Isometry, Random}) {
                    std::string directory = "data/input/61" + std::to_string(n_qbs) +"qbs/";

                    Mutator mutator = Mutator(0.0, 1.0);
                    CircuitHelper ch = CircuitHelper(n_qbs, "data/gates/" + gateset);
                    RandomHelper rh = RandomHelper();
                    rh.seed(42);
                    RandomCircuitGen rg = RandomCircuitGen(rh, 0.0, true);

                    Generator comp = Generator();
                    int number_per_gate = 10;

                    std::cout << "Working on qbs " << n_qbs << std::endl;
                    std::cout << "Working on set " << gateset << std::endl;
                    std::cout << "Working on gates " << n_gate << std::endl;
                    std::cout << "Working on type " << type << std::endl;
                    int i = count_files(directory, n_qbs, n_gate, type);
                    std::cout << i << std::endl;
                    while (i < number_per_gate) {
                        bool good_one = true;

                        Algorithm algo = Algorithm();
                        algo.parser.n_threads = 128;
                        algo.parser.n_found_stop = 1;
                        algo.parser.optimal_gatecount = n_gate - 1;
                        algo.parser.time_allowed = 100;
                        algo.parser.save_any_circuit = false;
                        
                        std::shared_ptr<GateCircuit> original_circ = comp.generateCircuit(n_gate, n_qbs, rg, ch, rh, mutator, extra_gates_allowed);
                        PartialMatrix original_matrix = comp.generateMatrixFromCircuit(original_circ, type, rg, ch, rh);

                        while (original_matrix.getCover().count() == 0) {
                            original_circ = comp.generateCircuit(n_gate, n_qbs, rg, ch, rh, mutator, extra_gates_allowed);
                            original_matrix = comp.generateMatrixFromCircuit(original_circ, type, rg, ch, rh);
                        }

                        int optimal_t_count = original_circ->getCount({"t", "tdg"});
                        int optimal_depth = original_circ->getDepth({"t", "tdg"});

                        QubitIndependentPartialMatrix matrix = QubitIndependentPartialMatrix(original_matrix, true, false);

                        std::cout << "starting run" << std::endl;
                        std::map<std::string, double> result = algo.run_inner_loop(std::make_shared<QubitIndependentPartialMatrix>(matrix));
                        good_one = (result["tcount"] != -1) && (result["gatecount"] >= n_gate);
                        optimal_t_count = std::min(optimal_t_count, (int) result["tcount"]);
                        optimal_depth = std::min(optimal_depth, (int) result["depth"]);
                        
                        if (good_one) {
                            std::string filename = directory + std::to_string(n_gate) + "_" + std::to_string(type) + "_" + std::to_string(optimal_t_count) + "_" + std::to_string(optimal_depth) + "_" + std::to_string(i) + ".txt";
                            std::cout << filename << std::endl;
                            matrix.original.write(filename);
                            i++;
                        }
                    }
                }
            }
        }
    }
}
