#ifndef DEF_RANDOMCIRCUIT
#define DEF_RANDOMCIRCUIT

#include "randomhelper.h"
#include "circuit.h"
#include "mutation.h"

class RandomCircuitGen{
    public:
        RandomCircuitGen(RandomHelper& random_helper, double id_prob = 0.7, bool ensure_non_id = false);
        GateCircuit randomGateCircuit(int nb_gates, int nb_qbs, CircuitHelper& ch, Mutator& mutator);
        GateCircuit randomGateCircuit(int min_nb_gates, int max_nb_gates, int nb_qbs, CircuitHelper& ch, Mutator& mutator);

    private:
        RandomHelper& random_helper;
        double id_prob = 0.7;
        bool ensure_non_id = false;
};

#endif