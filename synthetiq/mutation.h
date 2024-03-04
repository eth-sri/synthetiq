#ifndef DEF_MUT
#define DEF_MUT

#include "circuit.h"
#include "randomhelper.h"
#include <Eigen/Dense>

class Mutator{
    //usage: call mutate which picks a random mutaiton and applies it to candidate. Then either call undo mutation to get candidate back to what it was 
    //or call apply_mutation again to apply the same mutation to another circuit
    public:
        Mutator(double proba_id=0.2, double proportional_prob=0.2, double proba_name=0.5): proba_id(proba_id), proportional_prob(proportional_prob), proba_name(proba_name) {};
        Mutator(Mutator const& other) : proba_id(other.proba_id), proportional_prob(other.proportional_prob) {};
        std::shared_ptr<Mutator> clone();
        bool mutate(GateCircuit& candidate, RandomHelper& rh); 
        bool mutate_at_pos(GateCircuit& candidate, int position, RandomHelper& rh); 
        void undo_mutation(GateCircuit& candidate);
        void apply_mutation_again(GateCircuit& copyMutationFrom, GateCircuit& applyMutationTo);
    protected:
        double proba_id;
        double proportional_prob;
        double proba_name;
};

#endif