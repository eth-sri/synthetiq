#ifndef DEF_GENERATOR
#define DEF_GENERATOR
#include <Eigen/Dense>
#include <iostream>
#include "mcmc.h"
#include "mcmc_sa.h"
#include "partialMatrix.h"
#include "randomCircuit.h"
#include "randomhelper.h"

enum Type {Full = 0, Isometry = 1, Random = 2};
class Generator{
    public:
        Generator(bool append_mode=false);
        std::shared_ptr<GateCircuit> generateCircuit(int n_gate, int n_qubits, 
                                      RandomCircuitGen& rg, CircuitHelper& ch, RandomHelper& rh, Mutator& mutator, int extra_gates=0);
        PartialMatrix generateMatrixFromCircuit(std::shared_ptr<GateCircuit> circuit, int type, 
                                      RandomCircuitGen& rg, CircuitHelper& ch, RandomHelper& rh);
    private:
     bool append_mode;
};

#endif