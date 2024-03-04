from __future__ import annotations
from cmath import phase
import itertools
from itertools import product, permutations, combinations
from warnings import catch_warnings
import cirq
from cirq import Simulator
from cirq import work
from cirq.ops.moment import Moment
import numpy as np
from typing import Tuple, Dict, List, Union

if __name__ == "state":
    from util.utils import *
    from synthesis import Synthesis
    from synthesis_spec.specification import Spec
else :
    from qsyn.util.utils import *
    from qsyn.synthesis import Synthesis
    from qsyn.synthesis_spec.specification import Spec
from cirq.linalg.transformations import EntangledStateError

NC = "NoClassical"
BOOL = "Bool"
PHASING = "PHASING"
ENTANGLE = "Entangle"

CONTROL = "Controls"
TARGET = "Targs"

MUX_PROP = [CONTROL, TARGET]
ACTION_LOCAL = [NC, BOOL, PHASING]
ACTION_PROP = ACTION_LOCAL + [ENTANGLE] + [LOSS_INFO]


def OVERALL_ACTION_RELATION(a: str, b: str):
    if (a == ENTANGLE and b in ACTION_PROP) or (b == ENTANGLE and a in ACTION_PROP):
        return ENTANGLE
    if (a == NC and b in ACTION_LOCAL) or (b == NC and a in ACTION_LOCAL):
        return NC
    if (a == PHASING and b in [BOOL, PHASING]) or (b == PHASING and a in [BOOL, PHASING]):
        return PHASING
    if (a == LOSS_INFO):
        if b == NC :
            return NC
        else : return LOSS_INFO
    if (b == LOSS_INFO):
        if a == NC :
            return NC
        else : return LOSS_INFO
    else:
        return BOOL


class MyProperty:

    def __init__(self, mux_property=None, action_property=None) -> None:
        # assert isinstance(mux_property, dict) or isinstance(action_property, dict)
        self.mux_property = mux_property
        self.action_property = action_property

    def __repr__(self):
        return str(self)

    def __str__(self):
        str_builder = "Block : "
        if self.mux_property:
            str_builder += "<mux>" + str(self.mux_property)
        if self.action_property:
            str_builder += "<act>" + str(self.action_property)

        return str_builder

    def copy(self):
        if self.mux_property:
            copied_mux_prop = self.mux_property.copy()
        else:
            copied_mux_prop = None
        if self.action_property:
            copied_act_prop = self.action_property.copy()
        else:
            copied_act_prop = None
        return MyProperty(copied_mux_prop, copied_act_prop)

    def qubits_on(self) -> Set[int]:
        qubits = list()
        if self.mux_property:
            for value in self.mux_property.values():
                qubits += value
        if self.action_property:
            for key in self.action_property.keys():
                qubits += key
        return set(qubits)

    def custom_str(self):
        string_builder = "Block: "
        # for key in self.mux_property.keys():
        #     string_builder += f"\t{key} ∣-> {self.mux_property[key]}\n"
        if self.action_property:
            sorted_one= sorted(self.action_property.items())
            for key, val in sorted_one:
                string_builder += f"{key}:{val},"

        return string_builder

    
    def evaluate(self, ss, component_gates, working_qubits=None) -> List[cirq.Operation]:
        # TODO : `evaluate` rotuine currently only works for action_property = EMPTY(i.e, deals only with mux_prop)
        if self.action_property != None:
            return self.evaluate_act_prop(ss)
    
    def evaluate_act_prop(self, ss) -> List[cirq.Operation]:
        working_qubits = ss.working_qubit
        num_qubit_block = len(self.qubits_on())
        evaluated_gate_operations = list()
        if len(self.action_property.keys()) == 1: #if it is single-qubit block
            q = list(self.action_property.keys())[0] 
            prop = self.action_property[q] 
            if len(q) >= 2: #signle entangling
                assert prop == ENTANGLE
                # keys = list(self.action_property.keys())
                for gate in ss.spec.component_gates:
                    if gate.num_qubits() == len(q) and is_entangling_operation(gate):
                        qubits_on_gate  = list(range(gate.num_qubits()))
                        for x in itertools.permutations(qubits_on_gate, r = len(qubits_on_gate)):
                            gate_qubit_apps = [ working_qubits[q[idx]] for idx in x]
                            evaluated_gate_operations.append(gate(*gate_qubit_apps))
            else : # single local
                for gate in ss.spec.component_gates:
                    try:
                        if gate.num_qubits() == len(q) and check_action_property_of_gate(gate, prop):
                            qubits_to_be_applied = [working_qubits[i] for i in q]
                            evaluated_gate_operations.append(gate(*qubits_to_be_applied))
                    except AssertionError :
                        print("Something happened in state eval")
                        print("Gate",gate)
                        print("q", q)
                        print("state", self) 
                        exit()
        else:
            for gate in ss.spec.component_gates:
                if num_qubit_block == gate.num_qubits():
                    if is_multiplexor(gate):
                        raise NotImplementedError("How did we get here?")
                    elif gate == cirq.CNOT or gate == cirq.CX: #hardcoded
                        can_be_applied_as_cqubit = [working_qubits[key[0]] for key, val in self.action_property.items() if val == BOOL or val == None]
                        can_be_applied_as_tqubit = [working_qubits[key[0]] for key, val in self.action_property.items() if val == BOOL or val == None]
                        for c, t in itertools.product(can_be_applied_as_cqubit, can_be_applied_as_tqubit):
                            if c != t:
                                evaluated_gate_operations.append(gate(*[c, t]))
                    elif gate == cirq.CZ :  #hardcoded
                        already_checked = set() # CZ(a,b) = CZ(b,a), CZ commute in qubit order
                        can_be_applied_as_cqubit = [working_qubits[key[0]] for key, val in self.action_property.items() if val == BOOL or val == None]
                        can_be_applied_as_tqubit = [working_qubits[key[0]] for key, val in self.action_property.items() if val == PHASING or val == None]
                        for c, t in itertools.product(can_be_applied_as_cqubit, can_be_applied_as_tqubit):
                            if c != t and  frozenset((c,t)) not in already_checked :
                                evaluated_gate_operations.append(gate(*[c, t]))
                                already_checked.add(frozenset((c,t)))
                    elif gate == cirq.CSWAP:
                        already_checked = set()
                        can_be_applied_as_cqubit = [working_qubits[key[0]] for key, val in self.action_property.items() if val == BOOL or val == None]
                        can_be_applied_as_tqubit_one = [working_qubits[key[0]] for key, val in self.action_property.items() if val == PHASING or val == None]
                        can_be_applied_as_tqubit_another = [working_qubits[key[0]] for key, val in self.action_property.items() if val == PHASING or val == None]
                        
                        for c, t, t_pri in itertools.product(can_be_applied_as_cqubit, can_be_applied_as_tqubit_one, can_be_applied_as_tqubit_another):
                            if c!=t and t!=t_pri and c!= t_pri  and frozenset((c,t,t_pri)) not in already_checked:
                                evaluated_gate_operations.append(gate(*[c,t,t_pri]))
                                already_checked.add(frozenset((c,t,t_pri)))

                    elif isinstance(gate, cirq.ControlledGate):
                        can_be_applied_as_cqubit = [working_qubits[key[0]] for key, val in self.action_property.items() if val == BOOL or val == None]
                        possible_cqubits = [x for x in itertools.product(can_be_applied_as_cqubit, repeat=gate.num_controls())]
                        if gate.sub_gate.num_qubits() == 1:
                            can_be_applied_as_tqubit = [working_qubits[key[0]] for key, val in self.action_property.items() if check_action_property_of_gate(gate.sub_gate, val) or val == None]
                            for c, t in itertools.product(possible_cqubits, can_be_applied_as_tqubit):
                                flatten = [x for x in c]
                                flatten += [t]
                                if  len(set(c)) == len(c) and len(set(flatten)) == len(flatten):
                                    evaluated_gate_operations.append(gate(*flatten))
                        else:
                            can_be_applied_as_tqubit = [working_qubits[key[0]] for key, val in self.action_property.items() if check_action_property_of_gate(gate.sub_gate, val) or val == None]
                            raise NotImplementedError("To Do")
                    else:
                        is_none =  [ v == None for _, v in self.action_property.items()]
                        if all(is_none):
                            qubits_on_gate  = list(range(gate.num_qubits()))
                            keys = list(self.action_property.keys())
                            for x in itertools.permutations(qubits_on_gate, r = len(qubits_on_gate)):
                                gate_qubit_apps = [ working_qubits[keys[idx][0]] for idx in x]
                                evaluated_gate_operations.append(gate(*gate_qubit_apps))
                        else:
                            #if the gate is non-single-qubit gate and it is not even a controlled gate, we need to check its property manually
                            # like QFT? iSWAP?
                            if gate in ss.component_prior[BOOL]: # TODO : Replace by prior computation
                                match_bool_li  = [ val == BOOL for val in self.action_property.values()]
                                if all(match_bool_li) : 
                                    qubits_on_gate  = list(range(gate.num_qubits()))
                                    keys = list(self.action_property.keys())
                                    for x in itertools.permutations(qubits_on_gate, r = len(qubits_on_gate)):
                                        for idx in x :
                                            gate_qubit_apps = [ working_qubits[keys[idx][0]] for idx in x]
                                        evaluated_gate_operations.append(gate(*gate_qubit_apps))
                            elif gate in ss.component_prior[PHASING]: #for like iSWAP # TODO : Replace by prior computation
                                keys = list(self.action_property.keys())
                                qubits_on_gate  = list(range(gate.num_qubits()))
                                for x in itertools.permutations(qubits_on_gate, r = len(qubits_on_gate)):
                                    is_bool_or_phase = [ self.action_property[keys[idx]] in [BOOL, PHASING ]for idx in x ]
                                    is_exist_phase = [ self.action_property[keys[idx]] == PHASING for idx in x ]
                                    if all(is_bool_or_phase) and any(is_exist_phase):
                                        for idx in x :
                                            gate_qubit_apps = [ working_qubits[keys[idx][0]] for idx in x]
                                        evaluated_gate_operations.append(gate(*gate_qubit_apps))
                            else :
                                # rep_res = derive_representative_res_in_comp_basis(gate)
                                # some_result_in_computational_basis = cirq.unitary(gate)[:,-1]
                                some_result_in_computational_basis = ss.representative_basis_sv_res[gate]
                                gate_sv_result_prop_seq = list()
                                try : 
                                    for i in range(gate.num_qubits()):
                                        curr_qubit_state = sub_state_vector(some_result_in_computational_basis, keep_indices=(i,), atol=1e-04)
                                        gate_sv_result_prop_seq.append(get_state_property(curr_qubit_state)) #TODO : Hardcoded
                                    qubits_on_gate  = list(range(gate.num_qubits()))
                                    keys = list(self.action_property.keys())
                                    for x in itertools.permutations(qubits_on_gate, r = len(qubits_on_gate)):
                                        match_bool_li = list()
                                        for idx in x :
                                            if gate_sv_result_prop_seq[idx] == NC:
                                                match_bool_li.append(self.action_property[keys[idx]] == gate_sv_result_prop_seq[idx])
                                            else : 
                                                match_bool_li.append((self.action_property[keys[idx]] in [PHASING, BOOL] ))
                                        if all(match_bool_li) :
                                            gate_qubit_apps = [ working_qubits[keys[idx][0]] for idx in x]
                                            evaluated_gate_operations.append(gate(*gate_qubit_apps))
                                except EntangledStateError:  # if so, assume given quantum operator makes full entanglement on gate.num_qbuits() of qubits
                                    continue
        return evaluated_gate_operations

    def is_multiplexor(self):
        if not self.mux_property:
            return False
        if len(self.mux_property[CONTROL]) != 0:
            return True
        return False

    def is_local_action(self):
        return len(self.action_property.keys()) == 1 and len(list(self.action_property.keys())[0]) == 1

    def is_empty_block(self):
        if self.action_property == None and self.mux_property== None : 
            return True
        elif len(self.action_property.keys()) == 0 :
            return True
        else: return False

    def no_entanglement(self):
        for key, val in self.action_property.items():
            if val == ENTANGLE :
                return False
        return True

    
    def contains_none(self):
        for key, val in self.action_property.items():
            if val == None :
                return True
        return False

    def overall_action(self):
        if not self.action_property:
            raise Exception(f"Overall action not defined for {str(self)}")
        curr_oa = BOOL
        for act_val in self.action_property.values():
            curr_oa = OVERALL_ACTION_RELATION(act_val, curr_oa)
        return curr_oa


class MyMoment:

    def __init__(self, blocks: List = None):
        if not blocks:
            self.blocks = []
        else:
            # check if input blocks compose valid moment
            collected_qubits = list()
            for block in blocks:
                if isinstance(block, MyProperty):
                    collected_qubits.append(block.qubits_on())
                elif isinstance(block, cirq.Operation):
                    collected_qubits.append(set([q.x for q in block.qubits]))
                else:
                    raise TypeError(f"Invalid block for moment : {blocks}")
            if not pairwise_disjoint(collected_qubits):
                raise ValueError(f"Invalid moment : {blocks}")
            # blocks.sort(key = lambda s : str(s))
            self.blocks = blocks

    def __str__(self): 
        temp_blocks = self.blocks.copy()
        temp_blocks.sort(key = lambda s : str(s))
        string_builder = "[" 
        for block in temp_blocks :
            if isinstance(block, MyProperty):
                string_builder +=  block.custom_str() 
            else :
                string_builder += f" {str(block)}"
        string_builder +="]"
        return string_builder
        # return str(self.blocks)

    def __repr__(self) -> str:
        return str(self)

    def copy(self):
        return MyMoment(blocks=self.blocks.copy())

    def is_block_applicable(self, block: MyProperty):
        return self.qubits_on().intersection(block.qubits_on())

    def is_abstraction_only(self) -> bool:
        for block in self.blocks:
            if not isinstance(block, MyProperty):
                return False
        return True

    def evaluate(self, ss, component_gates, working_qubits):
        to_apply_product = list()
        for block in self.blocks:
            try:
                if isinstance(block, MyProperty):
                    to_apply_product.append(block.evaluate(ss,component_gates, working_qubits))
                if isinstance(block, cirq.Operation):
                    to_apply_product.append([block])
            except AttributeError as e:
                print("Error in evaluating a moment")
                print(block)
                print(e)
                exit()
        return [x for x in product(*to_apply_product)]

    def qubits_on(self) -> Set[int]:
        qubits_on = set()
        for block in self.blocks:
            if isinstance(block, MyProperty):
                qubits_on = qubits_on.union(block.qubits_on())
            elif isinstance(block, cirq.Operation):
                qubits_on = qubits_on.union(set([q.x for q in block.qubits]))
        return qubits_on

    def is_block_appendable(self, prop_block: MyProperty):
        return (not self.qubits_on().intersection(prop_block.qubits_on()))

    def overall_action(self) -> str:
        # TODO : Redefine here
        curr_oa = BOOL
        for block in self.blocks:
            if isinstance(block, cirq.Operation):
                raise Exception("Not defined")
            else:
                curr_oa = OVERALL_ACTION_RELATION(
                    block.overall_action(), curr_oa)
        return curr_oa

    def is_all_entangling_block(self) -> bool :
        for block in self.blocks :
            if isinstance(block, cirq.Operation):
                raise Exception("Not defined")
            elif not (len(block.action_property.keys()) == 1 and list(block.action_property.values())[0] == ENTANGLE):
                return False
        return True
                

    def is_empty(self) -> bool :
        if self.blocks == None or len(self.blocks) == 0 :
            return True
        else : return False

class MomentBasedState:

    def __init__(self, moments: List[MyMoment] = None):
        if not moments:
            self.moments = []
        else:
            self.moments = moments

    def append(self, moment):
        self.moments.append(moment)

    def __repr__(self):
        return str(self)

    def __str__(self):
        # copied_moments = copy.deepcopy(self.moments)
        if len(self.moments) != 0:
            str_builder = ""
            for idx, m in enumerate(self.moments):
                str_builder += f"Moment #{idx} \n"
                str_builder += "\t" + "-" + str(m) + "\n"
            return str_builder
        else:
            return "EMPTY STATE"

    def copy(self):
        return MomentBasedState(moments=self.moments.copy())

    def evaluate(self,ss, component_gates, working_qubits):
        to_product = list()
        to_return = list()
        for moment in self.moments:
            assert isinstance(moment, MyMoment)
            to_product.append(moment.evaluate(ss, component_gates, working_qubits))
        for x in product(*to_product):
            gen_qc = cirq.Circuit()
            flatten = [z for y in x for z in y]
            gen_qc.append(flatten)
            to_return.append(gen_qc)
        return to_return

    def is_concrete_state(self):
        for moment in self.moments:
            for block in moment.blocks:
                if isinstance(block, MyProperty):
                    return False
        return True

    def is_empty_state(self):
        if not self.moments:
            return True
        return len(self.moments) == 0

    def contains_none(self):
        for moment in self.moments:
            for block in moment.blocks:
                if isinstance(block, MyProperty) and block.contains_none() :
                    return True
        return False

def from_qc_to_moment_based_state(qc: cirq.Circuit) -> MomentBasedState:
    moment_bsaed_state = MomentBasedState()

    for idx, qc_m in enumerate(qc):
        blocks_for_moment = [gate_app for gate_app in qc_m]
        moment_bsaed_state.append(MyMoment(blocks=blocks_for_moment))

    return moment_bsaed_state

## Iinit Abstraction of given (concrete V) ##


# ==============================================================================================================================
# Given V circuit from phase 1, this gives intiial trials of abstractionof $V$ to feed on further state search at Phase 2
# - qc : the fed V circuit from phase1
# - synthesis_problem : to get necessary information such as component priors
# ==============================================================================================================================

def initial_abstracted_state(qc: cirq.Circuit, synthesis_problem: Synthesis) -> List[MomentBasedState]:
    component_pior = synthesis_problem.component_prior
    # spec = synthesis_problem.spec
    returned_abstract_states = []
    moment_state_rep = from_qc_to_moment_based_state(qc)  # temporary
    for idx, m in enumerate(moment_state_rep.moments):
        abstracted_mb_state = moment_state_rep.copy()
        blocks_for_abs_moment = []
        for block in m.blocks:
            if isinstance(block, cirq.Operation):
                prop = gate_op_property(block, component_pior)
            else:
                raise Exception(
                    f"Here block {block} of type {type(block)} should be some concrete Gate Operation")
            blocks_for_abs_moment.append(prop)
        abstracted_mb_state.moments[idx] = MyMoment(
            blocks=blocks_for_abs_moment)
        returned_abstract_states.append(abstracted_mb_state)
    return returned_abstract_states


def possible_moments_of_properties(qreg_size: int, minimum_valid_num_qubits: int) -> List[MyMoment]:
    # result_moments = list()
    # res = possible_moments_of_qregsize(
    #     qreg_size=qreg_size, valid_num_qubits=set(list(range(minimum_valid_num_qubits, qreg_size+1))))

    return possible_moments_of_mux_properties(qreg_size=qreg_size, minimum_valid_num_qubits=minimum_valid_num_qubits)


def possible_moments_of_act_properties(qreg_size: int, minimum_valid_num_qubits: int) -> List[MyMoment]:
    qreg = list(range(qreg_size))
    result_moments = list()
    res = possible_moments_of_qregsize(qreg_size=qreg_size, valid_num_qubits=set(list(range(minimum_valid_num_qubits, qreg_size + 1))))
    for m in res:
        to_producut_for_possible_blocks = list()
        for block in m:
            possible_prop_association_for_block = list()
            partitions_of_qubit_in_block = [x for x in partitions_of_list(block[1])]
            for partition in partitions_of_qubit_in_block:
                possible_assignments = list()
                for partition_elt in partition:
                    if len(partition_elt) > 1:
                        possible_assignments.append([ENTANGLE])
                    else:
                        possible_assignments.append(ACTION_LOCAL)
                for x in itertools.product(*possible_assignments, repeat=1):
                    action_prop_builder = dict()
                    for partition_elt, assignment in zip(partition, x):
                        action_prop_builder[partition_elt] = assignment
                    temp = MyProperty(action_property=action_prop_builder)
                    possible_prop_association_for_block.append(temp)
            to_producut_for_possible_blocks.append(
                possible_prop_association_for_block)
        for x in itertools.product(*to_producut_for_possible_blocks, repeat=1):
            blocks_of_moment = list(x)
            generated_moment = MyMoment(blocks=blocks_of_moment)
            result_moments.append(generated_moment)
    return result_moments

def possible_act_property_blocks(qreg_size: int, minimum_valid_num_qubits: int ):
    qubits_in_int = list(range(qreg_size))
    to_return = list()
    for i in range(minimum_valid_num_qubits, qreg_size+1): 
        for x in itertools.combinations(qubits_in_int, r= i):
            for partition in partitions_of_list(x):
                possible_assignments = list()
                for partition_elt in partition:
                    if len(partition_elt) > 1 : 
                        possible_assignments.append([ENTANGLE])
                    else :
                        possible_assignments.append(ACTION_LOCAL)
                for x in itertools.product(*possible_assignments, repeat=1):
                    action_prop_builder = dict()
                    for partition_elt, assignment in zip(partition, x):
                        action_prop_builder[partition_elt] = assignment
                    temp = MyProperty(action_property=action_prop_builder)
                    to_return.append(temp)
    return to_return
    

def possible_moments_of_mux_properties(qreg_size: int, minimum_valid_num_qubits: int) -> List[MyMoment]:
    #  minimum_valid_num_qubits = min(component_gate.num_qubits)
    result_moments = list()
    res = possible_moments_of_qregsize(
        qreg_size=qreg_size, valid_num_qubits=set(list(range(minimum_valid_num_qubits, qreg_size + 1))))

    for m in res:
        to_producut_for_possible_blocks = []
        for block in m:
            possible_prop_association_for_block = list()
            # cases for mux_prop
            cases_for_mux_prop = list()
            already_tracked = list()
            for part_one, part_two in partition_to_two_sub(block[1]):
                if not (set([part_one, part_two]) in already_tracked):
                    if len(part_two) != 0:
                        cases_for_mux_prop.append(dict(
                            [(CONTROL, part_one), (TARGET, part_two)]))
                    if len(part_one) != 0:
                        cases_for_mux_prop.append(dict(
                            [(CONTROL, part_two), (TARGET, part_one)]))
                    already_tracked.append(set([part_one, part_two]))
            for mux_prop in cases_for_mux_prop:
                temp = MyProperty(mux_property=mux_prop)
                possible_prop_association_for_block.append(
                    temp)
            to_producut_for_possible_blocks.append(
                possible_prop_association_for_block)
        for x in itertools.product(*to_producut_for_possible_blocks, repeat=1):
            blocks_of_moment = list(x)
            generated_moment = MyMoment(blocks=blocks_of_moment)
            result_moments.append(generated_moment)
    return result_moments


# Rules from state -> state
def gate_op_property(gate_app: cirq.GateOperation, component_prior) -> MyProperty:
    mux_prop_builder = dict()
    action_prop_builder = dict()
    if gate_app.gate in component_prior["C_PHASING_GATES"]:
        mux_prop_builder[CONTROL] = tuple([q_i.x for q_i in gate_app.controls])
        mux_prop_builder[TARGET] = tuple([
            q_i.x for q_i in gate_app.qubits if q_i not in gate_app.controls]
        )
    res = MyProperty(mux_property=mux_prop_builder,
                     action_property=action_prop_builder)
    return res


# ================================
# Given specific block within a state(=abstracted QC),
# It merely spans by width
# ..., MUX(q_c, q_t), ... ====>  ...,MUX(q_c, q_t),MUX(q_c, q_t),...
# ================================
def multiplexor_span_in_width(state: MomentBasedState, moment_idx: int, block_idx: None) -> MyMoment:
    # If moment_loc is given while block_loc is None, we assume that len((state.moments[moment_idx].blocks))==1
    new_generated_state = state.copy()
    if not ((not block_idx) and len(state.moments[moment_idx].blocks) == 1):
        raise Exception(
            "if block idx is given as none, we automatically assume that # of block on the moment is 1")
    if not block_idx:
        block_to_deal = state.moments[moment_idx].blocks[0]
    else:
        block_to_deal = state.moments[moment_idx].blocks[block_idx]
    new_block = block_to_deal.copy()
    moment_to_insert = MyMoment(blocks=[new_block])
    new_generated_state.moments.insert(moment_idx, moment_to_insert)
    return new_generated_state

# ================================================================================================================================
# Target or Control Qubit remove
# ================================================================================================================================


def mutiplexor_remove_control(state: MomentBasedState, moment_idx: int, block_idx=None, control_qubit_to_remove=None) -> MyMoment:
    new_generated_state = state.copy()
    if not ((not block_idx) and len(state.moments[moment_idx].blocks) == 1):
        raise Exception(
            "if block idx is given as none, we automatically assume that # of block on the moment is 1")
    if not block_idx:
        block_to_deal = state.moments[moment_idx].blocks[0]
        block_idx = 0
    else:
        block_to_deal = state.moments[moment_idx].blocks[block_idx]

    if control_qubit_to_remove not in block_to_deal.mux_property[CONTROL]:
        raise Exception(
            "The control qubit to remove is not valid for the current property block")

    # generate new prop block
    new_block = block_to_deal.copy()
    new_tup = tuple([x for x in new_block.mux_property[CONTROL]
                     if not x == control_qubit_to_remove])
    new_block.mux_property[CONTROL] = new_tup

    new_generated_state.moments[moment_idx].blocks[block_idx] = new_block
    return new_generated_state


def mutiplexor_remove_targ(state: MomentBasedState, moment_idx: int, block_idx=None, targ_qubit_to_remove=None) -> MyMoment:
    new_generated_state = state.copy()
    if not ((not block_idx) and len(state.moments[moment_idx].blocks) == 1):
        raise Exception(
            "if block idx is given as none, we automatically assume that # of block on the moment is 1")
    if not block_idx:
        block_to_deal = state.moments[moment_idx].blocks[0]
        block_idx = 0

    else:
        block_to_deal = state.moments[moment_idx].blocks[block_idx]

    if targ_qubit_to_remove not in block_to_deal.mux_property[CONTROL]:
        raise Exception(
            "The control qubit to remove is not valid for the current property block")
    new_block = block_to_deal.copy()

    new_block = block_to_deal.copy()
    new_tup = tuple([x for x in new_block.mux_property[TARGET]
                     if not x == targ_qubit_to_remove])
    new_block.mux_property[TARGET] = new_tup
    return new_generated_state

# ================================================================================================================================
# Decompose into block operation that when merged in sequentially, satisfies(although not soundly) the given block's mux properties
# ================================================================================================================================


def distance_between_prop(targ_prop: MyProperty, another_prop: MyProperty) -> Union[float, int]:
    if (not targ_prop.mux_property) and (not targ_prop.action_property):
        return float('-inf')
    # For now, asssume targ_prop and another_prop only has mux_prop
    dist = 0
    if set(targ_prop.mux_property[CONTROL]).issubset(another_prop.mux_property[CONTROL]):
        dist += 1
    if set(targ_prop.mux_property[TARGET]).issubset(another_prop.mux_property[TARGET]):
        dist += 1
    return dist


def distance_between_abs_moment(moment: MyMoment, another_moment: MyMoment) -> Union[float, int]:
    # Assume blocks of both `moment` and `another_moment` has property for mux_prop only
    dist = 0
    for block_moment in moment.blocks:
        for block_another_moment in another_moment.blocks:
            dist += distance_between_prop(block_moment, block_another_moment)
    return dist
