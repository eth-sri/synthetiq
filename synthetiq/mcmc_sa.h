#ifndef DEF_MCMC_SA
#define DEF_MCMC_SA
#include "circuit.h"
#include "cost.h"
#include "randomhelper.h"
#include "mutation.h"
#include "mcmc.h"
#include "partialMatrix.h"
#include "temperatureScheme.h"
#include <map>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <memory>
#include <vector>


class MCMC_Sa : public MCMC{
    public:
        //usceq = relevance of equality vs performance: score = usceq * equalityScore + perfScore
        MCMC_Sa(RandomHelper& random_helper, double proba_id=0.1, double factor_nb_steps=375.0, bool enable_permutations=true);

        MCMCResult run(QubitIndependentPartialMatrix& matrix_obj, std::shared_ptr<GateCircuit> init, CircuitHelper& circ_helper, bool debug=false, std::string debug_folder="debug/");
        std::shared_ptr<MCMC> clone();
        double getProbability(GateCircuit& circ, QubitIndependentPartialMatrix& matrix_obj, CircuitHelper& circ_helper, double temperature);
        bool acceptMutation(double rd_val, double candidate_energy, double cur_energy, double temperature);
        int calculateClosestMatrix(QubitIndependentPartialMatrix& matrix_obj, GateCircuit& circ, CircuitHelper& circ_helper);
    private:
        const double factor_nb_steps;
};

#endif