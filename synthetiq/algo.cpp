#include "algo.h"

#include "circuit.h"
#include "mcmc_sa.h"
#include "cost.h"
#include "randomhelper.h"
#include "randomCircuit.h"
#include "resynthesis.h"
#include "partialMatrix.h"
#include "gateScheme.h"
#include "matrixGenerator.h"

#include <omp.h>
#include <iostream>
#include <vector>
#include <map>
#include <thread>
#include <Eigen/Dense>
#include <filesystem>
#include <sstream>

Parser::Parser() {

}

/**
 * Parses the command line arguments and sets the corresponding member variables.
 * 
 * @param argc The number of command line arguments.
 * @param argv The array of command line arguments.
 */
void Parser::parse(int argc, char* argv[]) {
    input_name = argv[1];
    format = input_name.substr(input_name.rfind(".") + 1, input_name.size());
    // this is default parameter for output_folder
    output_folder = input_name.substr(0, input_name.rfind(".")) + "/";

    for (int arg = 2; arg < argc; arg++) {
        if (std::string(argv[arg]) == "--threads"  || std::string(argv[arg]) == "-h") {
            n_threads = std::stoi(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--save") {
            save_times = true;
        } else if (std::string(argv[arg]) == "--time" || std::string(argv[arg]) == "-t") {
                time_allowed = std::stod(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--circuits" || std::string(argv[arg]) == "-c") {
                n_found_stop = std::stoi(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--gates"  || std::string(argv[arg]) == "-g") {
            int gates = std::stoi(argv[arg + 1]);
            gateScheme.setMinStartGates(gateScheme.getMinFactor() * gates);
            gateScheme.setMaxStartGates(gateScheme.getMaxFactor() * gates);
            gateScheme.setStartBestGates(gates);
            gateScheme.reset();
        } else if (std::string(argv[arg]) == "--output" || std::string(argv[arg]) == "-o") {
            output_folder = argv[arg + 1];
            output_folder += "/";
        } else if (std::string(argv[arg]) == "--format" || std::string(argv[arg]) == "-f") {
            format = argv[arg + 1];
        } else if (std::string(argv[arg]) == "--ancilla" || std::string(argv[arg]) == "-a") {
            n_ancillas = std::stoi(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--dirty" || std::string(argv[arg]) == "-d") {
            n_dirty = std::stoi(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--expand" || std::string(argv[arg]) == "-e") {
            expand_composite = true;
        } else if (std::string(argv[arg]) == "--qubit_independence" || std::string(argv[arg]) == "-q") {
            qubit_independent = false;
        } else if (std::string(argv[arg]) == "--tcount" || std::string(argv[arg]) == "-tc") {
            optimal_tcount = std::stoi(argv[arg + 1]);
        }  else if (std::string(argv[arg]) == "--tdepth" || std::string(argv[arg]) == "-td") {
            optimal_tdepth = std::stoi(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--gate-count" || std::string(argv[arg]) == "-gc") {
            optimal_gatecount = std::stoi(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--cost-required" || std::string(argv[arg]) == "-cr") {
            cost_required = std::stod(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--gate-set" || std::string(argv[arg]) == "-gs") {
            gate_set = argv[arg + 1];
        } else if (std::string(argv[arg]) == "--composite-gates" || std::string(argv[arg]) == "-cg") {
            composite_gate_folder = argv[arg + 1];
        } else if (std::string(argv[arg]) == "--save-all" || std::string(argv[arg]) == "-sa") {
            save_all_circuits = true;
        } else if (std::string(argv[arg]) == "--start-temp" || std::string(argv[arg]) == "-st") {
            start_temp_base = std::stod(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--epsilon" || std::string(argv[arg]) == "-eps") {
            epsilon = std::stod(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--beta") {
            double beta = std::stod(argv[arg + 1]);
            gateScheme.setBeta(beta);
        } else if (std::string(argv[arg]) == "--fmin") {
            double fmin = std::stod(argv[arg + 1]);
            gateScheme.setMinFactor(fmin);
        } else if (std::string(argv[arg]) == "--fmax") {
            double fmax = std::stod(argv[arg + 1]);
            gateScheme.setMaxFactor(fmax);
        } else if (std::string(argv[arg]) == "--pcomp") {
            pcomp = std::stod(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--pid") {
            pid = std::stod(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--no-perms") {
            enable_permutations = false;
        } else if (std::string(argv[arg]) == "--no-resynth") {
            do_resynth = false;
        } else if (std::string(argv[arg]) == "--simple") {
            simple_cost = true;
        } else if (std::string(argv[arg]) == "--n-norm") {
            n_norm = std::stod(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--iterations-factor") {
            iterations_factor = std::stoi(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--gs-no-update") {
            update_gate_scheme = false;
        } else if (std::string(argv[arg]) == "--n-start-gates") {
            int n_start_gates = std::stoi(argv[arg + 1]);
            gateScheme.setMinStartGates(n_start_gates);
            gateScheme.setMaxStartGates(n_start_gates + 1);
            update_gate_scheme = false;
        } else if (std::string(argv[arg]) == "--times-file") {
            times_file = argv[arg + 1];
        } else if (std::string(argv[arg]) == "--optimization-number") {
            optimization_numb = std::stoi(argv[arg + 1]);
        } else if (std::string(argv[arg]) == "--no-save-circuits") {
            save_any_circuit = false;
        } else if (std::string(argv[arg]) == "--no-optimize-depth") {
            optimize_depth = false;
        } else if (std::string(argv[arg]) == "--inverse-independent") {
            inverse_independent = true;
        } else if (std::string(argv[arg]) == "--depth-gates") {
            std::string depth_gates_string = argv[arg + 1];
            depth_gates = {};
            std::stringstream ss(depth_gates_string);
            std::string token;

            while (std::getline(ss, token, ',')) {
                depth_gates.push_back(token);
            }
        }
    }

    // print all the parameters
    std::cout << "Input file: " << input_name << std::endl;
    std::cout << "Output folder: " << output_folder << std::endl;
    std::cout << "Format: " << format << std::endl;
    std::cout << "Gate set: " << gate_set << std::endl;
    std::cout << "Composite gate folder: " << composite_gate_folder << std::endl;
    std::cout << "Number of threads: " << n_threads << std::endl;
    std::cout << "Number of circuits to find: " << n_found_stop << std::endl;
    std::cout << "Number of ancillas: " << n_ancillas << std::endl;
    std::cout << "Number of dirty qubits: " << n_dirty << std::endl;
    std::cout << "Optimal t-count: " << optimal_tcount << std::endl;
    std::cout << "Optimal t-depth: " << optimal_tdepth << std::endl;
    std::cout << "Optimal gate count: " << optimal_gatecount << std::endl;
    std::cout << "Cost required: " << cost_required << std::endl;
    std::cout << "Time allowed: " << time_allowed << std::endl;
    std::cout << "Start temp base: " << start_temp_base << std::endl;
    std::cout << "Epsilon: " << epsilon << std::endl;
    std::cout << "Beta: " << gateScheme.getBeta() << std::endl;
    std::cout << "Min factor: " << gateScheme.getMinFactor() << std::endl;
    std::cout << "Max factor: " << gateScheme.getMaxFactor() << std::endl;
    std::cout << "Pcomp: " << pcomp << std::endl;
    std::cout << "Pid: " << pid << std::endl;
    std::cout << "Enable permutations: " << enable_permutations << std::endl;
    std::cout << "Do resynth: " << do_resynth << std::endl;
    std::cout << "Simple cost: " << simple_cost << std::endl;
    std::cout << "N norm: " << n_norm << std::endl;
    std::cout << "Iterations factor: " << iterations_factor << std::endl;
    std::cout << "Update gate scheme: " << update_gate_scheme << std::endl;
    std::cout << "Save times: " << save_times << std::endl;
    std::cout << "Expand composite: " << expand_composite << std::endl;
    std::cout << "Qubit independent: " << qubit_independent << std::endl;
    std::cout << "Save all circuits: " << save_all_circuits << std::endl;
    std::cout << "Save any circuit: " << save_any_circuit << std::endl;
    std::cout << "Times file: " << times_file << std::endl;
    std::cout << "Optimization number: " << optimization_numb << std::endl;
    std::cout << "Optimize depth: " << optimize_depth << std::endl;
    std::cout << "Inverse independent: " << inverse_independent << std::endl;
    std::cout << "Depth gates: ";
    for (int i = 0; i < depth_gates.size(); i++) {
        std::cout << depth_gates[i] << " ";
    }
    std::cout << std::endl;

    total_output_folder = base_output_folder + output_folder;
}


/**
 * Creates the output folder for the parser.
 * This function iterates through the total_output_folder path and creates each folder in the path if it does not exist.
 */
void Parser::createOutputFolder() {
    int found_slash = total_output_folder.find("/");
    while (found_slash != std::string::npos) {
        std::string current_folder = total_output_folder.substr(0, found_slash);
        if (!std::filesystem::is_directory(current_folder)) {
            std::filesystem::create_directory(current_folder);
        }
        found_slash = total_output_folder.find("/", found_slash + 1);
    }
}

/**
 * @brief Constructor for the Algorithm class.
 * 
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line arguments.
 */
Algorithm::Algorithm(int argc, char* argv[]) {
    parser.parse(argc, argv);
    parser.createOutputFolder();
    t1_total = std::chrono::high_resolution_clock::now();
    t2_total = std::chrono::high_resolution_clock::now();
    time_taken_total = 0.0;
}

/**
 * @brief Default constructor for the Algorithm class.
 * Initializes the t1_total and t2_total variables with the current time,
 * and sets the time_taken_total variable to 0.0.
 */
Algorithm::Algorithm() {
    t1_total = std::chrono::high_resolution_clock::now();
    t2_total = std::chrono::high_resolution_clock::now();
    time_taken_total = 0.0;
}

/**
 * @brief Initializes a shared pointer to a QubitIndependentPartialMatrix object.
 * 
 * This function creates a QubitIndependentPartialMatrix object and returns a shared pointer to it.
 * The QubitIndependentPartialMatrix is initialized with the provided original matrix, qubit independence flag, and inverse independence flag.
 * 
 * @param original The original matrix used to initialize the QubitIndependentPartialMatrix.
 * @param qubit_independent A flag indicating whether the QubitIndependentPartialMatrix should use qubit independence for faster discovery.
 * @param inverse_independent A flag indicating whether the QubitIndependentPartialMatrix should use inverse independence for faster discovery.
 * @return std::shared_ptr<QubitIndependentPartialMatrix> A shared pointer to the initialized QubitIndependentPartialMatrix object.
 */
std::shared_ptr<QubitIndependentPartialMatrix> Algorithm::initMatrix() {
    PartialMatrix original;
    MatrixGenerator matrix_gen = MatrixGenerator();
    CircuitHelper ch = CircuitHelper(1, parser.base_gate_folder + parser.gate_set, parser.base_gate_folder + parser.composite_gate_folder);
    if (parser.format == "qasm") {
        GateCircuit circuit_obj = GateCircuit(ch);
        circuit_obj.readFromInput(parser.base_input_folder + parser.input_name);
        ch = *circuit_obj.getCircuitHelper();
        original = PartialMatrix(std::ref(circuit_obj), ch);
    } else {
        original = PartialMatrix(parser.base_input_folder + parser.input_name);
    }     
    if (parser.n_ancillas > 0) {
        original = matrix_gen.addAncilla(original, parser.n_ancillas);
    }   
    if (parser.n_dirty > 0) {
        original = matrix_gen.addDirtyQubits(original, parser.n_dirty);
    }
    QubitIndependentPartialMatrix matrix = QubitIndependentPartialMatrix(original, parser.qubit_independent, parser.inverse_independent);
    return std::make_shared<QubitIndependentPartialMatrix>(matrix);
}

/**
 * @brief Runs the inner loop of the algorithm.
 * 
 * This function runs the inner loop of the algorithm, which is the loop that runs the MCMC algorithm until a circuit is found, or until the time is out.
 * The function returns a map of the best t-count, t-depth, gate count, and cost found in the inner loop.
 * 
 * @param input_matrix A shared pointer to the QubitIndependentPartialMatrix object to run the inner loop on.
 * @param run The number of times the inner loop has been run so far (used for seeding).
 */
std::map<std::string, double> Algorithm::run_inner_loop(std::shared_ptr<QubitIndependentPartialMatrix> input_matrix, int run) {
    bool stop_inner = false;
    std::map<std::string, double> output_map;
    output_map["tcount"] = -1;
    output_map["tdepth"] = -1;
    output_map["gatecount"] = -1;
    output_map["cost"] = -1;
    auto t_start_inner = std::chrono::high_resolution_clock::now();
    auto t_end_inner = std::chrono::high_resolution_clock::now();
    double time_taken_inner = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t_end_inner-t_start_inner).count();
    omp_set_num_threads(parser.n_threads);
    #pragma omp parallel num_threads(parser.n_threads) 
    {   
        int id = omp_get_thread_num();
        RandomHelper random_helper = RandomHelper();
        random_helper.seed(id + run * parser.n_threads);
        std::shared_ptr<ExactEqualityComputer> exact_comp = std::make_shared<ExactEqualityComputer>(ExactEqualityComputer(parser.epsilon));
        RandomCircuitGen random_gen = RandomCircuitGen(random_helper, parser.pid);
        Resynthesize resynth = Resynthesize(parser.optimization_numb, parser.optimize_depth);
        std::shared_ptr<QubitIndependentPartialMatrix> matrix = input_matrix->clone();
        CircuitHelper ch = CircuitHelper(matrix->getNQubits(), parser.base_gate_folder + parser.gate_set, parser.base_gate_folder + parser.composite_gate_folder);
        GatesSumComputer perf_comp = GatesSumComputer();
        MCMC_Sa algo2 = MCMC_Sa(random_helper, parser.pid, parser.iterations_factor * matrix->getNQubits(), parser.enable_permutations);
        ExponentialTemperatureScheme temp_scheme = ExponentialTemperatureScheme(parser.start_temp_base / std::sqrt(pow(2.0, matrix->getNQubits())), parser.n_norm);
        algo2.set_temp_scheme(std::make_shared<ExponentialTemperatureScheme>(temp_scheme));
        
        if (!parser.simple_cost) {
            std::shared_ptr<FroebeniusCostComputer> froeb = std::make_shared<FroebeniusCostComputer>(FroebeniusCostComputer());
            algo2.set_eq_comp(froeb);
        } else {
            std::shared_ptr<SimpleFroebeniusCostComputer> froeb = std::make_shared<SimpleFroebeniusCostComputer>(SimpleFroebeniusCostComputer());
            algo2.set_eq_comp(froeb);
        }
        algo2.set_exact_eq_comp(exact_comp);
        Mutator mutator = Mutator(parser.pid, parser.pcomp);
        algo2.set_mutator(std::make_shared<Mutator>(mutator));
        while (time_taken_total < parser.time_allowed && n_found_so_far < parser.n_found_stop && !stop_inner) {
            n_runs += 1;
            std::shared_ptr<GateCircuit> circ_init;
            int startGates = parser.gateScheme.getStartGates(random_helper);
            circ_init = std::make_shared<GateCircuit>(random_gen.randomGateCircuit(startGates, matrix->original.getNQubits(), ch, *algo2.getMutator()));
            MCMCResult res = algo2.run(*matrix, circ_init, ch, false);
            bool found = exact_comp->normalizedEqualityCost(*res.circuit_best, *matrix, ch) < 1e-3;
            
            if (found) {
                successful_runs += 1;
                std::shared_ptr<GateCircuit> best = res.circuit_best;
                int gates_before = best->getNbNonIdGates();
                int tcount_before = best->getCount({"t", "tdg"});
                int tdepth_before = best->getDepth({"t", "tdg"});

                if (parser.do_resynth) {
                    resynth.run(best, ch);
                }

                if (parser.update_gate_scheme) {
                    parser.gateScheme.update(best->getNbNonIdGates());
                }

                if (parser.expand_composite) {
                    best->expandCompositeGates();
                    if (parser.do_resynth) {
                        resynth.run(best, ch);
                    }
                }

                int gates_after = best->getNbNonIdGates();
                int tcount_after = best->getCount({"t", "tdg"});
                int tdepth_after = best->getDepth({"t", "tdg"});

                if (output_map["tcount"] == -1) {
                    output_map["tcount"] = tcount_after;
                    output_map["tdepth"] = tdepth_after;
                    output_map["gatecount"] = gates_after;
                    output_map["cost"] = best->getCost();
                } else {
                    output_map["tcount"] = std::min(output_map["tcount"], (double) tcount_after);
                    output_map["tdepth"] = std::min(output_map["tdepth"], (double) tdepth_after);
                    output_map["gatecount"] = std::min(output_map["gatecount"], (double) gates_after);
                    output_map["cost"] = std::min(output_map["cost"], best->getCost());
                }
                
                bool save_found_tcount = true;
                bool save_found_tdepth = true;
                bool save_found_gatecount = true;
                bool save_found_cost = true;

                if (parser.optimal_tcount > -1) {
                    int tcount = best->getCount({"t", "tdg"});
                    save_found_tcount = tcount <= parser.optimal_tcount && !stop_inner;
                } 
                if (parser.optimal_tdepth > -1) {
                    int tdepth = best->getDepth({"t", "tdg"});
                    save_found_tdepth = tdepth <= parser.optimal_tdepth  && !stop_inner;
                }

                if (parser.optimal_gatecount > -1) {
                    int gatecount = best->getNbNonIdGates();
                    save_found_gatecount = gatecount <= parser.optimal_gatecount && !stop_inner;
                }

                if (parser.cost_required > -1) {
                    save_found_cost = best->getCost() <= parser.cost_required && !stop_inner;
                }

                bool save_found = save_found_tcount && save_found_tdepth && save_found_gatecount && save_found_cost;

                if ((parser.optimal_tdepth > -1 || parser.optimal_tcount > -1 || parser.cost_required > -1 || parser.optimal_gatecount > -1) && save_found) {
                    stop_inner = true;
                    optimal_runs += 1;
                }

                if (save_found) {
                    t_end_inner = std::chrono::high_resolution_clock::now();
                    time_taken_inner = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t_end_inner-t_start_inner).count();
                    times.push_back(time_taken_inner);
                    std::cout << parser.total_output_folder << " " << time_taken_inner << " " << n_found_so_far + 1 << std::endl;
                    std::cout << "Gates: " << gates_before << " -> " << gates_after << std::endl;
                    std::cout << "T-count: " << tcount_before << " -> " << tcount_after << std::endl;
                    std::cout << "T-depth: " << tdepth_before << " -> " << tdepth_after << std::endl;
                    t_start_inner = std::chrono::high_resolution_clock::now();
                    n_found_so_far += 1;
                }

                if ((save_found || parser.save_all_circuits) && parser.save_any_circuit) {
                    double performance = best->getCost();
                    int count = best->getCount(parser.depth_gates);
                    int depth = best->getDepth(parser.depth_gates);
                    std::string filename = std::to_string(performance) + "-" + std::to_string(count) + "-" + std::to_string(depth) + "-" + std::to_string(id) + "-" + std::to_string(n_found_so_far) + ".qasm";
                    std::ofstream myfile;
                    myfile.open(parser.total_output_folder + filename);
                    myfile <<  best->print_qasm() << std::endl;
                    myfile.close();
                }
            }
            if (id == 0) {
                t2_total = std::chrono::high_resolution_clock::now();
                time_taken_total = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t2_total-t1_total).count();
            }
        }   
    }

    return output_map;
}

/**
 * Runs the algorithm until the time limit or the number of desired solutions is reached.
 */
void Algorithm::run() {
    int run = 0;
    while (time_taken_total < parser.time_allowed && n_found_so_far < parser.n_found_stop) {
        std::shared_ptr<QubitIndependentPartialMatrix> matrix = initMatrix();
        run_inner_loop(matrix, run);
        t2_total = std::chrono::high_resolution_clock::now();
        time_taken_total = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t2_total-t1_total).count();
        parser.gateScheme.reset();
        run += 1;

    }
    
    if (parser.save_times) {
        std::ofstream file;
        file.open(parser.times_file, std::ios_base::app);
        double mean_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
        double sq_sum = std::inner_product(times.begin(), times.end(), times.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / times.size() - mean_time * mean_time);
        file << parser.total_output_folder << "\t" << mean_time << "\t" << stdev << "\t" << n_runs << "\t" << successful_runs << "\t" << optimal_runs << std::endl;
        file.close();
    }
}