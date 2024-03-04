
import itertools
import cirq
from typing import List

if __name__ == "qsyn.util.rule_utils": 
    from qsyn.util.state_search_utils import *
    from qsyn.util.utils import NC, BOOL, ENTANGLE, PHASING, possible_moments_of_qregsize
else :
    from util.state_search_utils import *
    from util.utils import NC, BOOL, ENTANGLE, PHASING, possible_moments_of_qregsize

def syntatic_separated_component_gates(g_li : List[cirq.Gate]) :
    separated_cg = list()
    for g in g_li:
        if g.num_qubits()==1 :
            separated_cg.append(g)
        elif isinstance(g, cirq.ControlledGate):
            separated_cg.append(g.sub_gate)
        elif g == cirq.CZ:
            separated_cg.append(cirq.Z)
        elif g == cirq.CNOT or g == cirq.CX:
            separated_cg.append(cirq.X)
        else :
            separated_cg.append(g)
    return separated_cg



# possibily removing redudant rules

#check C_RULE_ID - OPEN-AND-CLOSE_B is valid
def check_C_RULE_OPEN_AND_CLOSE_B(sep_cg : List[cirq.Gate]):
    nc_cgs = [ g for g in sep_cg  if NC in get_action_prop_scheme_of_gate(g)]
    bool_cgs  = [ g for g in sep_cg if g not in nc_cgs and BOOL in get_action_prop_scheme_of_gate(g)]

    if not bool_cgs or not nc_cgs : 
        return False

    test_qc_size = max([ g.num_qubits() for g in nc_cgs + bool_cgs])
    
    moments_from_nc_cgs = posible_moments_by_qregsize_component_gates(test_qc_size, nc_cgs)
    moments_from_bool_cgs = posible_moments_by_qregsize_component_gates(test_qc_size, bool_cgs)
    for moments_to_feed_qc in itertools.product(moments_from_nc_cgs,moments_from_bool_cgs,moments_from_nc_cgs):
        qc = cirq.Circuit(moments_to_feed_qc)
        if is_classical_operation(cirq.unitary(qc)):
            return True

#check C_RULE_ID - OPEN-AND-CLOSE_P is valid
def check_C_RULE_OPEN_AND_CLOSE_P(sep_cg : List[cirq.Gate]):
    nc_cgs = [ g for g in sep_cg  if NC in get_action_prop_scheme_of_gate(g)]
    p_cgs  = [ g for g in sep_cg if g not in nc_cgs and PHASING in get_action_prop_scheme_of_gate(g)]
    if nc_cgs == [] or p_cgs == [] :
        return False
    test_qc_size = max([ g.num_qubits() for g in nc_cgs + p_cgs])
    if not p_cgs or not nc_cgs : 
        return False
    test_qc_size = max([ g.num_qubits() for g in nc_cgs + p_cgs])
    moments_from_nc_cgs = posible_moments_by_qregsize_component_gates(test_qc_size, nc_cgs)
    moments_from_p_cgs = posible_moments_by_qregsize_component_gates(test_qc_size, p_cgs)
    for moments_to_feed_qc in itertools.product(moments_from_nc_cgs,moments_from_p_cgs,moments_from_nc_cgs):
        qc = cirq.Circuit(moments_to_feed_qc)
        if is_classical_operation(cirq.unitary(qc)):
            return True
        else : 
            print(qc)

#check P_RULE_ID - OPEN-AND-CLOSE_B is valid
def check_P_RULE_OPEN_AND_CLOSE_B(sep_cg : List[cirq.Gate]):
    nc_cgs = [ g for g in sep_cg  if NC in get_action_prop_scheme_of_gate(g)]
    bool_cgs  = [ g for g in sep_cg if g not in nc_cgs and BOOL in get_action_prop_scheme_of_gate(g)]
    if not bool_cgs or not nc_cgs : 
        return False
    test_qc_size = max([ g.num_qubits() for g in nc_cgs + bool_cgs])
    moments_from_nc_cgs = posible_moments_by_qregsize_component_gates(test_qc_size, nc_cgs)
    moments_from_bool_cgs = posible_moments_by_qregsize_component_gates(test_qc_size, bool_cgs)
    for moments_to_feed_qc in itertools.product(moments_from_nc_cgs,moments_from_bool_cgs,moments_from_nc_cgs):
        qc = cirq.Circuit(moments_to_feed_qc)
        if is_phasing_operation(cirq.unitary(qc)):
            return True
        # else: 
        #     print(qc)
        #     print(is_phasing_operation(cirq.unitary(qc)))
        #     print(is_phasing_operation(cirq.Z))

        #     print(cirq.unitary(qc))
        #     print(cirq.unitary(cirq.Z))

def posible_moments_by_qregsize_component_gates(qreg_size : int, cg : List[cirq.Gate]) -> List[cirq.Moment]: 
    qreg = cirq.LineQubit.range(qreg_size)
    valid_num_qubits = list(set([g.num_qubits() for g in cg]))
    possible_moments = possible_moments_of_qregsize(qreg_size, valid_num_qubits)
    to_return = list()
    for moment in possible_moments:
        blocks_for_moments = list()
        for block in moment:
            gate_ops_for_block = list()
            q_tobe_applied = [qreg[q_i] for q_i in block[1]]
            for perm_q in itertools.permutations(q_tobe_applied, r = len(q_tobe_applied)):
                for g in cg:
                    if g.num_qubits() == len(perm_q):
                        gate_ops_for_block.append(g(*perm_q))
            blocks_for_moments.append(gate_ops_for_block)
        
        for concrete_moment_li in itertools.product(*blocks_for_moments):
            res_concret_moment = cirq.Moment(concrete_moment_li)
            to_return.append(res_concret_moment)

    return to_return