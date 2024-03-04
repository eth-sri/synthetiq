if __name__ == "qsyn.util.state_search_utils" :
    from qsyn.state import MyProperty, MyMoment, MomentBasedState, ACTION_PROP
    from qsyn.util.utils import *

else :
    from state import MyProperty, MyMoment, MomentBasedState, ACTION_PROP
    from util.utils import *

from os import dup
from typing import Tuple, Dict, List, Union, Set
import cirq
from itertools import permutations, combinations, product
import numpy as np

def qubits_giving_difference(in_arr : np.array, out_arr : np.array) -> set:
    in_arr_binary = classic_sv_to_binarys(in_arr)
    out_arr_binary = classic_sv_to_binarys(out_arr)
    if (not in_arr_binary) or (not out_arr_binary):
        return None
    diff_qbits = set()
    for i in range(len(in_arr_binary)):
        if in_arr_binary[i] != out_arr_binary[i]:
            diff_qbits.add(i)
    return diff_qbits


class InvalidState(Exception):
    def __init__(self):
        super().__init__('Invalid State.')


## Moment Insertion ##
def is_appendable(moment: MyMoment, block: MyProperty):
    block_qubits = block.qubits_on()
    if len(moment.blocks) == 0:
        return True
    for block_in_moment in moment.blocks:
        if set(block_in_moment.qubits_on()).intersection(set(block_qubits)):
            return False
    return True

def blocks_to_opt_moments(blocks) -> List[MyMoment]:
    if len(blocks) == 0:
        return None
    init_list = list()
    init_moment = MyMoment(blocks=init_list)
    gen_moments = list()
    curr_moment_tobe_filled = init_moment
    len_blocks = len(blocks)
    for idx, block in enumerate(blocks):
        if is_appendable(curr_moment_tobe_filled, block):
            curr_moment_tobe_filled.blocks.append(block)
        else:
            curr_moment_tobe_filled.blocks.sort(key= lambda s : str(s))
            gen_moments.append(curr_moment_tobe_filled)
            empty_li = list()
            curr_moment_tobe_filled = MyMoment(blocks=empty_li)
            curr_moment_tobe_filled.blocks.append(block)
        if idx == len_blocks - 1:
            curr_moment_tobe_filled.blocks.sort(key= lambda s : str(s))
            gen_moments.append(curr_moment_tobe_filled)
    if len(gen_moments) == 0:
        curr_moment_tobe_filled.blocks.sort(key= lambda s : str(s))
        gen_moments.append(curr_moment_tobe_filled)
    return gen_moments


def connected_lists_of_pair(qubits_tobe_entgled : Union[List[int], Tuple[int]], superposed_qubits : Union[List[int], Tuple[int]], inseparable_gate_qubit_num : List[int]) -> List[Tuple[int,int]]:
    assert set(superposed_qubits).issubset(set(qubits_tobe_entgled))
    to_return = list()
    # for r in inseparable_gate_qubit_num:
    #     list_of_pairs = [x for x in combinations(qubits_tobe_entgled, r=r)] # assume component_gate.num_qubits() <= 2
    #     for li_pair in product(list_of_pairs ,repeat=len(qubits_tobe_entgled)- r + 1):
    #         if len(li_pair) == len(set(li_pair)):
    #             if check_will_entangling(li_pair, qubits_tobe_entgled, superposed_qubits) :
    #                     to_return.append(li_pair)

    duplicated = list()
    num_ins_gate_app = len(qubits_tobe_entgled)- 1 #e.g if len(qubits_tobe_entgled)=4 then 3 application of CNOT is enough
    if inseparable_gate_qubit_num.count(max(inseparable_gate_qubit_num) )>=2  and max(inseparable_gate_qubit_num)>=3  : # if contain ins gate of more than 2 qubit, we need less ins gate
        num_ins_gate_app -=1
    # print(num_ins_gate_app)
    inseparable_gate_qubit_num_set = list(set(inseparable_gate_qubit_num))
    for i in range(1, num_ins_gate_app + 1):
        list_of_pairs_mixed = [x for r in inseparable_gate_qubit_num_set for x in combinations(qubits_tobe_entgled, r=r)  ] 
        prod_res = [ x for x in product(list_of_pairs_mixed ,repeat=i)]
        for li_pair in prod_res:
            if len(li_pair) == len(set(li_pair))  and li_pair not in duplicated:
                duplicated.append(li_pair)
                if check_will_entangling(li_pair, qubits_tobe_entgled, superposed_qubits) :
                    to_return.append(li_pair)
    return to_return

def connect_components(li_of_set): 
    # https://stackoverflow.com/questions/54673308/how-to-merge-sets-which-have-intersections-connected-components-algorithm
    pool = set(map(frozenset, li_of_set))
    groups = []
    while pool:
        groups.append(set(pool.pop()))
        while True:
            for candidate in pool:
                if groups[-1] & candidate:
                    groups[-1] |= candidate
                    pool.remove(candidate)
                    break
            else:
                break
    return groups

def check_will_entangling(li_pair : List[Tuple[int]], qubits_tobe_entgled : List[int], superposed_qubits : List[int]) -> bool:
    set_of_entagled_subregister = list()
    if not set(li_pair[0]).intersection(set(superposed_qubits)): 
        return False
    # superposed_qubits = q_sp
    # li_pair = q_1, .., q_n
    for x in li_pair:
        if set(x).intersection(set(superposed_qubits)):
            temp  = set(x)
            set_of_entagled_subregister.append(temp)
        else :
            for entangled_subregister in set_of_entagled_subregister:
                if set(x).intersection(entangled_subregister):
                    entangled_subregister.update(x)
                    break
    set_of_entagled_subregister = connect_components(set_of_entagled_subregister)
    return len(set_of_entagled_subregister) == 1 and set(set_of_entagled_subregister[0]) == set(qubits_tobe_entgled)



def get_action_prop_scheme_of_gate(gate : cirq.Gate) :
    if gate.num_qubits() == 1 :
        # entanglement inducement checking is essential in order 
        if is_nonclassical_operation(gate):
            return (NC,)
        elif is_phasing_operation(gate):
            return (PHASING,)
        elif is_classical_operation(gate):
            return (BOOL,)
    else :
        if isinstance(gate, cirq.ControlledGate):
            tuple_builder = [ BOOL for _ in range(gate.num_controls())]
            sub_gate = gate.sub_gate 
            tuple_builder += list(get_action_prop_scheme_of_gate(sub_gate))
            return tuple(tuple_builder)
        elif gate == cirq.CX or gate == cirq.CNOT :
            # raise NotImplementedError("gate for cirq.CX, cirq.CNOT are not implemented yet") 
            return (BOOL, BOOL)
        elif gate == cirq.CZ  : 
            return (BOOL, PHASING)
        elif is_classical_operation(gate):
            tuple_builder = [ BOOL for i in range(gate.num_qubits())]
            return tuple(tuple_builder)
        else :
            tuple_builder = list()
            # some_result_in_computational_basis = cirq.unitary(gate)[:,-1]
            some_result_in_computational_basis = derive_representative_res_in_comp_basis(gate)
            try : 
                for i in range(gate.num_qubits()):
                    curr_qubit_state = sub_state_vector(some_result_in_computational_basis, keep_indices=(i,), atol=1e-04)
                    tuple_builder.append(get_state_property(curr_qubit_state))
                return tuple(tuple_builder)
            except EntangledStateError:
                return (ENTANGLE, )
                

def get_action_property_of_gate_op(gate_op : cirq.Operation) -> MyProperty: 
    # gate is either single qubit gate or mutliple qubit gate
    prop_builder = dict()
    gate = gate_op.gate
    if gate.num_qubits() == 1 :
        if is_nonclassical_operation(gate):
            prop_builder[(gate_op.qubits[0].x,)] = NC
        elif is_phasing_operation(gate):
            prop_builder[(gate_op.qubits[0].x,)] = PHASING
        elif is_classical_operation(gate):
            prop_builder[(gate_op.qubits[0].x,)] = BOOL
        else:
            prop_builder[(gate_op.qubits[0].x,)] = ENTANGLE
        return MyProperty(action_property=prop_builder)
    else :  # multiple qubit gate
        if isinstance(gate, cirq.ControlledGate):
            sub_gate = gate.sub_gate
            targ_qubits = [ q for q in gate_op.qubits if q not in gate_op.controls ]
            for cq in gate_op.controls :
                prop_builder[(cq.x,)] = BOOL
            sub_gate_myprop = get_action_property_of_gate_op(sub_gate(*targ_qubits))
            new_dict = {**prop_builder, **sub_gate_myprop.action_property}
            return MyProperty(action_property=new_dict)
        elif gate == cirq.CX or gate == cirq.CNOT or is_classical_operation(gate):
            for q in gate_op.qubits:
                prop_builder[(q.x,)] = BOOL
            return MyProperty(action_property=prop_builder) 
        else : 
            some_result_in_computational_basis = derive_representative_res_in_comp_basis(gate)
            # some_result_in_computational_basis = cirq.unitary(gate)[:,-1]
            for i, q in enumerate(gate_op.qubits):
                curr_qubit_state = sub_state_vector(some_result_in_computational_basis, keep_indices=(i,), atol=1e-04)
                prop_builder[(q.x,)] = get_state_property(curr_qubit_state)
            return MyProperty(action_property=prop_builder)
                
