#include "mutation.h"
#include <vector>

/**
 * Undoes the last mutation performed on the given GateCircuit.
 * 
 * @param candidate The GateCircuit to undo the mutation on.
 */
void Mutator::undo_mutation(GateCircuit& candidate){
    candidate.undoMutation();
}

/**
 * Mutates the given GateCircuit using the provided RandomHelper.
 * 
 * @param candidate The GateCircuit to be mutated.
 * @param random_helper The RandomHelper used for generating random mutations.
 * @return True if the mutation was successful, false otherwise.
 */
bool Mutator::mutate(GateCircuit &candidate, RandomHelper& random_helper){
    return candidate.mutate(random_helper, proportional_prob, proba_id, proba_name);
}

/**
 * Mutates the GateCircuit at the specified position.
 * 
 * @param candidate The GateCircuit to be mutated.
 * @param position The position at which the mutation should occur.
 * @param rh The RandomHelper object used for generating random numbers.
 * @return True if the mutation was successful, false otherwise.
 */
bool Mutator::mutate_at_pos(GateCircuit& candidate, int position, RandomHelper& rh) {
    return candidate.mutate(rh, position, proportional_prob, proba_id, proba_name);
}

/**
 * Applies the mutation again from one GateCircuit to another.
 * 
 * @param copyMutationFrom The GateCircuit to copy the mutation from.
 * @param applyMutationTo The GateCircuit to apply the mutation to.
 */
void Mutator::apply_mutation_again(GateCircuit& copyMutationFrom, GateCircuit& applyMutationTo){
    copyMutationFrom.applyMutationAgain(applyMutationTo);
}

/**
 * @brief Clones the Mutator object.
 * 
 * @return A shared pointer to a new Mutator object that is a copy of the current object.
 */
std::shared_ptr<Mutator> Mutator::clone() {
    return std::make_shared<Mutator>(*this);
}
