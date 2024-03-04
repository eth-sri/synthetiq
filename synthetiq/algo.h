#ifndef DEF_ALGO
#define DEF_ALGO

#include <chrono>
#include <memory>
#include <map>

#include "gateScheme.h"
#include "partialMatrix.h"

class Parser {
    public:
        Parser();
        void parse(int argc, char* argv[]);
        void createOutputFolder();
       
        std::string base_input_folder = "data/input/";
        std::string base_output_folder = "data/output/";
        std::string base_gate_folder = "data/gates/";

        std::string total_output_folder = "";
        std::string input_name;
        std::string output_folder;
        std::string format;
        std::string gate_set = "CliffordT";
        std::string composite_gate_folder = "composite_gates";
        std::string times_file = "data/times.csv";

        int n_threads = 1;
        int n_found_stop = 10;
        int n_ancillas = 0;
        int n_dirty = 0;
        int optimal_tcount = -1;
        int optimal_tdepth = -1;
        int optimal_gatecount = -1;
        double cost_required = -1;
        
        double time_allowed = 100.0;
        double start_temp_base = 0.1;
        double epsilon = 1e-6;

        bool save_times = false;
        bool expand_composite = false;
        bool qubit_independent = true;
        bool inverse_independent = false;
        bool save_all_circuits = false;
        bool save_any_circuit = true;
        
        double pid = 0.3;
        double pcomp = 0.2;
        bool enable_permutations = true;
        bool do_resynth = true;
        bool simple_cost = false;

        double n_norm = 80.0;
        int iterations_factor = 40;

        bool update_gate_scheme = true;

        int optimization_numb = 12;
        bool optimize_depth = true;

        std::vector<std::string> depth_gates = {"t", "tdg"};

        GateScheme gateScheme = GateScheme();

};

class Algorithm {
    public:
        Algorithm(int argc, char* argv[]);
        Algorithm();
        void run();
        std::map<std::string, double> run_inner_loop(std::shared_ptr<QubitIndependentPartialMatrix> matrix, int run=0);
        std::shared_ptr<QubitIndependentPartialMatrix> initMatrix();

        Parser parser = Parser();
        int n_found_so_far = 0;
        int n_runs = 0;
        int successful_runs = 0;
        int optimal_runs = 0;
        std::chrono::time_point<std::chrono::high_resolution_clock> t1_total = std::chrono::high_resolution_clock::now();
        std::chrono::time_point<std::chrono::high_resolution_clock> t2_total = std::chrono::high_resolution_clock::now();
        double time_taken_total;
        std::vector<double> times;
};


#endif