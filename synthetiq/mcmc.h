#ifndef DEF_MCMC
#define DEF_MCMC
#include "circuit.h"
#include "cost.h"
#include "randomhelper.h"
#include "mutation.h"
#include "partialMatrix.h"
#include "temperatureScheme.h"
#include <map>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <memory>
#include <vector>

enum MonitoredValType {Equality, Performance, Probability, AcceptedEquality, Energy, Equality_froebenius, Equality_entanglement};
const int nbTypesMonitorable = 5;

class MCMCResult{
    public:
        std::shared_ptr<GateCircuit> circuit_best;
        double best_energy;
        double best_eq; //from best proba hein

        MCMCResult();
        MCMCResult(std::shared_ptr<GateCircuit> circuitBest);
        std::string print();
};

class MCMC{
    public:
        //consts for cost computation
        
        MCMC(RandomHelper& random_helper, double proba_id=0.2, bool enable_permutations=true);
        virtual std::shared_ptr<MCMC> clone() = 0;
    //    ~MCMC();

        void set_eq_comp(std::shared_ptr<EqualityComputer>e); //to change perf from default
        void set_exact_eq_comp(std::shared_ptr<EqualityComputer>e); //to change perf from default
        void set_mutator(std::shared_ptr<Mutator> m); //change mutations from default
        void set_temp_scheme(std::shared_ptr<TemperatureScheme> t);
        std::shared_ptr<EqualityComputer> getEqualityComputer();
        std::shared_ptr<Mutator> getMutator();
        std::shared_ptr<TemperatureScheme> getTemperatureScheme();

        virtual MCMCResult run(QubitIndependentPartialMatrix& matrix_obj, std::shared_ptr<GateCircuit> init, CircuitHelper& circ_helper, bool debug=false, std::string debug_folder="debug/") = 0;
        double getEnergy(GateCircuit &circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& circ_helper);
        double getEnergy(double eq_cost);
        double equalityCost(GateCircuit &circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& circ_helper, bool log_val = false);
        double equalityCost(GateCircuit &circ, PartialMatrix& matrix_obj, CircuitHelper& circ_helper, bool log_val = false);
        bool acceptMutation(double rd_val, double candidate_energy, double cur_energy, double temperature);
        void correctResultQubitIndependence(MCMCResult& res, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& circ_helper);
        void log_debug_information(GateCircuit& circ, QubitIndependentPartialMatrix& matrix_obj, MCMCResult& current_result, std::ofstream& file, CircuitHelper& ch);
        std::string find_next_available_file(std::string& debug_folder);

    protected:
        RandomHelper& random_helper;
        std::shared_ptr<EqualityComputer> equalitycomp;
        std::shared_ptr<Mutator> mutator;
        std::shared_ptr<TemperatureScheme> temp_scheme;
        std::shared_ptr<EqualityComputer> exact_eq_comp;
        bool enable_permutations = true;
};

#endif