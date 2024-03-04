import math
from itertools import product, permutations
from cirq import Simulator
from cirq.ops.pauli_gates import X
import numpy as np
from typing import Tuple, Dict, List, Union
from numpy.core.fromnumeric import resize
from timeit import default_timer as timer, repeat
if __name__ == "guide":
    from util.utils import *
    from state import *
    from rules import *
    from score import *
    from state_search import *
    from synthesis import Synthesis
    from synthesis_spec.specification import Spec

else :
    from qsyn.util.utils import *
    from qsyn.state import *
    from qsyn.rules import *
    from qsyn.score import *
    from qsyn.state_search import *
    from qsyn.synthesis import Synthesis
    from qsyn.synthesis_spec.specification import Spec
from cirq.value import state_vector_to_probabilities

## features for differ analysis
IDENTICAL = "IDENTICAL"


class IdentityCheckException(Exception):   
    def __init__(self):
        super().__init__('Given two sv to query attribute difference was actually identical!')


def sv_differ_by(phi_sv: np.array, psi_sv: np.array) -> Tuple[Union[NC, BOOL, PHASING, ENTANGLE], str]:
    #str may be some additional info for differ analysis, may be differ by qubit analysis?
    # psi_sv is the output!!
    
    assert len(phi_sv) == len(psi_sv)
    entgle_condi_phi = entangling_condition(phi_sv) #set of frozensets
    entgle_condi_psi = entangling_condition(psi_sv) #set of frozensets
    if entgle_condi_phi != entgle_condi_psi:
        return ENTANGLE, (entgle_condi_phi, entgle_condi_psi)
    prop_dist_of_phi_sv = state_vector_to_probabilities(phi_sv)
    prop_dist_of_psi_sv = state_vector_to_probabilities(psi_sv)
    prop_dist_of_phi_sv_nonzero = list(prop_dist_of_phi_sv[np.nonzero(prop_dist_of_phi_sv)])
    prop_dist_of_psi_sv_nonzero = list(prop_dist_of_psi_sv[np.nonzero(prop_dist_of_psi_sv)])
    prop_dist_of_phi_sv_nonzero.sort()
    prop_dist_of_psi_sv_nonzero.sort()
    prop_dist_of_phi_sv_nonzero= np.array(prop_dist_of_phi_sv_nonzero)
    prop_dist_of_psi_sv_nonzero= np.array(prop_dist_of_psi_sv_nonzero)
    if len(np.nonzero(phi_sv)[0]) !=  len(np.nonzero(psi_sv)[0]) or (not np.allclose(prop_dist_of_phi_sv_nonzero, prop_dist_of_psi_sv_nonzero, atol = 1e-04)):
        info_for_NC = set()
        for x in entgle_condi_phi :
            if len(x) == 1:
                qubit_local = list(x)[0]
                qubit_local_sv_of_phi = sub_state_vector( phi_sv, keep_indices=[qubit_local], atol=1e-04)
                qubit_local_sv_of_psi = sub_state_vector( psi_sv, keep_indices=[qubit_local], atol=1e-04)
                qubit_local_sv_of_phi = list(np.nonzero(qubit_local_sv_of_phi))
                qubit_local_sv_of_psi = list(np.nonzero(qubit_local_sv_of_psi))
                if len(qubit_local_sv_of_phi[0]) != len(qubit_local_sv_of_psi[0]):
                    info_for_NC.add(qubit_local)
        if len(info_for_NC) != 0 :
            return NC, info_for_NC
        else :
            return NC, None
    elif len(np.nonzero(phi_sv)[0]) ==  len(np.nonzero(psi_sv)[0]) and amplitudes_same_upto_order(phi_sv, psi_sv,atol=1e-04):
        if len(prop_dist_of_phi_sv_nonzero) == 1 and len(prop_dist_of_psi_sv_nonzero)==1:
            qbits_giving_diff = qubits_giving_difference(phi_sv,psi_sv)
            if qbits_giving_diff :
                return BOOL, qbits_giving_diff
            else : 
                return BOOL, None
        else : ## hardcodeed for later useage
            return BOOL, None
    else : 
        return PHASING, None
    
def will_induce_entgle(moment, source_entanglements, target_entanglements) -> bool: 
    # assert isinstance(moment, MyMoment)
    curr_entanglement_status = [set(x) for x in source_entanglements]
    for block in moment.blocks:
        for key, val in block.action_property.items():
            if val == ENTANGLE :
                curr_entanglement_status.append(set(key))
    curr_entanglement_status = [ frozenset(x) for x in connect_components(curr_entanglement_status)]
    curr_entanglement_status = set(curr_entanglement_status)
    assert type(curr_entanglement_status) == type(target_entanglements)
    return curr_entanglement_status == target_entanglements
        

def biased_moment_set_of_property(differ_by, additional_info, ss: Synthesis) -> List[List[Moment]]:
    to_return = list()
    if differ_by == PHASING:
        for x in ss.rules.possible_moments_of_act_properties:
            if x.overall_action() == PHASING:
                to_return.append(x)
    elif differ_by == BOOL:
        assert additional_info == None or isinstance(additional_info, set)
        if additional_info:
            for x in ss.rules.possible_moments_of_act_properties:
                if x.overall_action() == BOOL and len(additional_info.intersection(x.qubits_on())) != 0  : 
                    to_return.append(x)
            to_return.sort(key=lambda s : (not additional_info.issubset(s.qubits_on()) , - len(s.qubits_on()) )  )
        else :
            for x in ss.rules.possible_moments_of_act_properties:
                to_return.append(x)
    elif differ_by == ENTANGLE: #should be somehow ordered
        soruce_entgle, targ_entgle = additional_info
        flag_source_entangled = False
        if len(soruce_entgle) < len(targ_entgle):
            flag_source_entangled = True
            temp = soruce_entgle
            soruce_entgle = targ_entgle
            targ_entgle = temp
        for x in ss.rules.possible_moments_of_act_properties:
            if x.is_all_entangling_block()  and will_induce_entgle(x, soruce_entgle, targ_entgle ):
                to_return.append(x)
    elif differ_by == NC:
        for x in ss.rules.possible_moments_of_act_properties:
            if x.overall_action() == NC:
                to_return.append(x)
    else:
        raise Exception("Difference_by value must be on of ACTION_PROP")
    return to_return


def guided_append_rule_generation(V_state: MomentBasedState, ss: Synthesis, targ_attribute : Tuple) -> List[Tuple[str, dict]]: #return list of rules
    generated_rule = list()
    curr_qc = V_state.evaluate(ss, component_gates=ss.spec.component_gates, working_qubits=ss.working_qubit)
    # differs = list()
    if V_state.is_empty_state():
        curr_qc = cirq.Circuit()
        # for io_pair in ss.spec.spec_object.get_io_pairs():
            # in_nparr, out_nparr = io_pair
            # if (in_nparr == out_nparr).all() : 
                # continue
            # sv_differ_by_res = sv_differ_by( in_nparr, out_nparr)
            # if sv_differ_by_res not in differs:
                # differs.append(sv_differ_by_res)    
    elif len(curr_qc) != 1:
        print("Why this happend? in `guided_append_rule_generation`")
        exit()
    else:
        curr_qc = curr_qc[0]
        # for io_pair in ss.spec.spec_object.get_io_pairs():
        #     in_nparr, out_nparr = io_pair
        #     res = ss.simulator.simulate(curr_qc, qubit_order=ss.working_qubit, initial_state=in_nparr)
        #     differs.append(sv_differ_by(res.final_state_vector, out_nparr))

    for differy_by, additional_info in [targ_attribute]: #guide into appen rules
        for x in biased_moment_set_of_property(differy_by, additional_info, ss=ss):
            if len(curr_qc.moments) == 0:
                param_builder = dict([("state", from_qc_to_moment_based_state(curr_qc)), ("moment_to_append", x), ("append_position", (-1, 0))])
                generated_rule.append(("append", param_builder))
            else:
                for idx in [len(curr_qc.moments) -1]:
                # for idx in [-1, len(curr_qc.moments) - 1]:
                    param_builder = dict([("state", from_qc_to_moment_based_state(curr_qc)), ("moment_to_append", x), ("append_position", (idx, idx + 1))])
                    generated_rule.append(("append", param_builder))
    return generated_rule


def prop_block_distance(block_A : MyProperty, block_B : MyProperty) -> int:
    dist = 0
    if not set(block_A.action_property.keys()).issubset(set((block_B.action_property.keys()))):
        return +1000
    else:
        for key, val in block_A.action_property.items():
            if val != block_B.action_property[key]:
                dist += 100
        return dist

# def has_room_for_component_gate(state, ss):
#     possible_num_gates = [ g.num_qubits() for g in ss.spec.component_gates ]
#     for moment_idx, moment in enumerate(state.moments):
#         for block_idx, block in enumerate(moment.blocks):
#             if isinstance(block, MyProperty) and len(block.qubits_on()) in possible_num_gates:
#                 return True
#     return False

def num_of_unfit_block(state, ss) -> List :
    # counting number of blocks that is not fittable on component gates
    cnt = 0 
    for moment_idx, moment in enumerate(state.moments):
        for block_idx, block in enumerate(moment.blocks):
            if isinstance(block, MyProperty) and (len(block.qubits_on()) not in ss.valid_num_qubits):
                cnt += 1 
            elif isinstance(block, MyProperty) :
                act_prop_schem_match_booli = [sorted(list(block.action_property.values()), key= lambda s : str(s)) == sorted(act_prop_schem, key= lambda s : str(s)) 
                                              for act_prop_schem in ss.act_prop_schems_of_component_gate ]
                if not any(act_prop_schem_match_booli) :
                    cnt+=1
    return cnt

def list_of_unfit_block(state : MomentBasedState, ss ) -> List :
    # reutrn list of unfit block in [ (moment_idx, block_idx) ] 
    to_return = list()
    for moment_idx, moment in enumerate(state.moments):
        for block_idx, block in enumerate(moment.blocks):
            if isinstance(block, MyProperty) and not block.contains_none():
                if (len(block.qubits_on()) not in ss.valid_num_qubits):
                    to_return.append((moment_idx, block_idx))
                else : 
                    act_prop_schem_match_booli = [sorted(list(block.action_property.values()), key= lambda s : str(s)) == sorted(act_prop_schem, key= lambda s : str(s)) 
                                                for act_prop_schem in ss.act_prop_schems_of_component_gate ]
                    if not any(act_prop_schem_match_booli) :
                        to_return.append((moment_idx, block_idx))
    return to_return

def plausible_decompose_rule_for_further_decompose(state : MomentBasedState, ss ) -> List[Tuple[str, dict]]:
    generated_rule = list()
    li_of_unfit_block = list_of_unfit_block(state,ss)
    if len(li_of_unfit_block) >= 2 :
        return list()
    else :
        for (moment_idx, block_idx) in li_of_unfit_block : 
            targ_block = state.moments[moment_idx].blocks[block_idx]
            if isinstance(targ_block, MyProperty) and targ_block.action_property:
                valid_rule_ids = ss.rules.valid_rule_id_for_block(block=targ_block)
                if not valid_rule_ids:
                    params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", None), ("block_idx", block_idx)])
                    generated_rule.append(("prop_decompose", params))
                else :
                    for rule_id in valid_rule_ids:
                        params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", rule_id), ("block_idx", block_idx)])
                        generated_rule.append(("prop_decompose", params))
        return generated_rule


def plausible_decompose_rule_for_intensive_search(state : MomentBasedState, num_unfit_block : int,  ss) -> List[Tuple[str,dict]]: 
    generated_rule = list()
    if num_unfit_block >= 2 :
        return list()
    elif num_unfit_block == 0 :
        for moment_idx, moment in enumerate(state.moments):
            for block_idx, block in enumerate(moment.blocks):
                if isinstance(block, MyProperty) and block.action_property:
                    valid_rule_ids = ss.rules.valid_rule_id_for_block(block=block)
                    if not valid_rule_ids: # it is multi qubit decomposition
                        params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", None), ("block_idx", block_idx)])
                        generated_rule.append(("prop_decompose", params))
                    else:
                        for rule_id in valid_rule_ids:
                            params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", rule_id), ("block_idx", block_idx)])
                            generated_rule.append(("prop_decompose", params))
        return generated_rule
    else : # unfit block is only one, num_unfit_block == 1
        for moment_idx, moment in enumerate(state.moments):
            for block_idx, block in enumerate(moment.blocks):
                if isinstance(block, MyProperty) and (len(block.qubits_on()) not in ss.valid_num_qubits):
                    valid_rule_ids = ss.rules.valid_rule_id_for_block(block=block)
                    if not valid_rule_ids: # it is multi qubit decomposition
                        params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", None), ("block_idx", block_idx)])
                        generated_rule.append(("prop_decompose", params))
                    else:
                        for rule_id in valid_rule_ids:
                            params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", rule_id), ("block_idx", block_idx)])
                            generated_rule.append(("prop_decompose", params))                    
                    return generated_rule
    return generated_rule

def guided_replace_rule_generation(V_state: MomentBasedState, ss: Synthesis) -> List[Tuple[str,dict]]: 
    param_rule_with_score = list()
    generated_rule = list()
    curr_qc = V_state.evaluate(ss, component_gates=ss.spec.component_gates, working_qubits=ss.working_qubit) 
    assert len(curr_qc) == 1
    for moment_idx, moment in enumerate(V_state.moments):
        for block_idx, block in enumerate(moment.blocks): # here the block is an concrete quantum operation
            act_prop_of_block = get_action_property_of_gate_op(block)  
            for to_replace in ss.rules.possible_act_property_blocks :
                # if to_replace.overall_action() == PHASING:
                #     print(set(act_prop_of_block.action_property.items()))
                #     print(set(to_replace.action_property.items()))
                #     print(to_replace.overall_action())
                #     print(set(act_prop_of_block.action_property.items()).issubset(set(to_replace.action_property.items())))
                targ_overall_action = act_prop_of_block.overall_action()
                if targ_overall_action == LOSS_INFO:
                    if to_replace.overall_action() in [PHASING, BOOL] and set(act_prop_of_block.action_property.items()) <= set(to_replace.action_property.items()) :
                        param_builder = {
                            "state"        : V_state,
                            "replace_to"   : to_replace,
                            "moment_idx"   : moment_idx,
                            "block_idx"    : block_idx
                        }
                        param_rule_with_score.append( (param_builder,   ( prop_block_distance(act_prop_of_block, to_replace)   , len(act_prop_of_block.action_property.keys())    ))  ) 
                else :
                    if act_prop_of_block.overall_action() == to_replace.overall_action() and set(act_prop_of_block.action_property.items()) <= set(to_replace.action_property.items()) :
                        param_builder = {
                            "state"        : V_state,
                            "replace_to"   : to_replace,
                            "moment_idx"   : moment_idx,
                            "block_idx"    : block_idx
                        }
                        param_rule_with_score.append( (param_builder,   ( prop_block_distance(act_prop_of_block, to_replace)   , len(act_prop_of_block.action_property.keys())    ))  ) 
        param_rule_with_score.sort( key= lambda s : s[-1])
    generated_rule = [ ("replace", r)  for r, i in param_rule_with_score ]
    return generated_rule
    
def guided_rule_generation(V_state: MomentBasedState, ss: Synthesis) -> List[Tuple[str, dict]]: #return list of rules
    generated_rule = list()
    curr_qc = V_state.evaluate(component_gates=ss.spec.component_gates, working_qubits=ss.working_qubit)
    differs = list()
    if V_state.is_empty_state():
        curr_qc = cirq.Circuit()
        for io_pair in ss.spec.spec_object.get_io_pairs():
            in_nparr, out_nparr = io_pair
            differs.append(sv_differ_by( in_nparr, out_nparr) ) 
    elif curr_qc == None or len(curr_qc) == 0 :
        # state that with gates can be fiiled
        return ss.rules.prop_decompose_rule_space(state=V_state)
    elif len(curr_qc) > 1:
        raise ValueError("Value Error Here")
        # print("For guidance V_state must be single concrete q circuit.")
        # print("Else, it just returns the state's rule space")
        # print("Tinkle. No Guide Here")
        # return ss.rules.rule_space(state=V_state)
    else:
        curr_qc = curr_qc[0]
        for io_pair in ss.spec.spec_object.get_io_pairs():
            in_nparr, out_nparr = io_pair
            res = ss.simulator.simulate(curr_qc, qubit_order=ss.working_qubit, initial_state=in_nparr)
            differs.append(sv_differ_by(res.final_state_vector, out_nparr))


    for differy_by, additional_info in differs: #guide into appen rules
        for x in biased_moment_set_of_property(differy_by, additional_info, ss=ss):
            if len(curr_qc.moments) == 0: ## if initial state
                param_builder = dict([("state", from_qc_to_moment_based_state(curr_qc)), ("moment_to_append", x), ("append_position", (-1, 0))])
                generated_rule.append(("append", param_builder))
            else:
                for idx in [-1, len(curr_qc.moments) - 1]:
                    param_builder = dict([("state", from_qc_to_moment_based_state(curr_qc)), ("moment_to_append", x), ("append_position", (idx, idx + 1))])
                    generated_rule.append(("append", param_builder))

    return generated_rule

def init_targ_attribute(ss : Synthesis, initial_V_qc : Union[cirq.Circuit, MomentBasedState]) -> Tuple[Union[ENTANGLE, NC, BOOL, PHASING], str]: 
    if isinstance(initial_V_qc, MomentBasedState):
        if initial_V_qc.is_empty_state():
            initial_V_qc = None
        else : 
            eval_res = initial_V_qc.evaluate(ss, component_gates=ss.spec.component_gates, working_qubits=ss.working_qubit)
            assert len(eval_res) == 1
            initial_V_qc = eval_res[0]
    differs = list() # list of tuple[attribute-id, additional_info]
    for in_nparr, out_nparr in ss.spec.spec_object.get_io_pairs():
        if initial_V_qc:
            res_by_init_Vqc = ss.simulator.simulate(initial_V_qc, qubit_order=ss.working_qubit, initial_state=in_nparr).final_state_vector
        else:
            res_by_init_Vqc = in_nparr
        if not np.allclose(res_by_init_Vqc, out_nparr):
            differs.append(sv_differ_by(res_by_init_Vqc,  out_nparr))
    differs = handle_multiple_targ_attribute(differs)
    print("Initdiffers")
    print(differs)
    return representative_targ_attribute(differs)


def handle_multiple_targ_attribute(differs : List)-> List:
    # to deal with frozen set
    if len(differs) == 1:
        return differs
    res = list()
    for x in differs :
        if x not in res :
            res.append(x)
    return res

def set_rep_of_differs(differs : List)-> List:
    if len(differs)==1 : 
        return differs
    res = list()
    for x in differs:
        if x not in res and x!= IDENTICAL:
            res.append(x) 
    return res

def order_of_differ_for_rep(differ_by : Tuple):
    if differ_by[0] == ENTANGLE : 
        return 3
    if differ_by[0] == BOOL :
        return 0
    if differ_by[0] == PHASING :
        return 1
    if differ_by[0] == NC :
        return 2
    
    print(f"order computation for {differ_by}")
    exit()

def index_of_non_identical(differs):
    for idx, item in enumerate(differs):
        if item != IDENTICAL : return idx     

def representative_targ_attribute(differs : List) -> List :
    if len(differs) == 1 : 
        return differs[0]
    if differs == [IDENTICAL] :
        return list()
    idx_of_first_non_identical = index_of_non_identical(differs)
    curr_representative_att = differs[idx_of_first_non_identical]

    for differ in differs[idx_of_first_non_identical:]:
        if (differ != IDENTICAL 
            and curr_representative_att[0] == BOOL 
            and differ[0]==BOOL):
            if curr_representative_att[1] and differ[1]:
                if curr_representative_att[1]  <= differ[1]:
                    curr_representative_att  = differ
        elif (differ != IDENTICAL 
            and order_of_differ_for_rep(differ) > order_of_differ_for_rep(curr_representative_att)):
            curr_representative_att = differ
    assert curr_representative_att != None
    return curr_representative_att

# def att_filled_in_score(differs, targ_att) :
#     score = 0
#     for differ in differs :
#         if differ == IDENTICAL:
#             score +=1
#         if targ_att[0] == ENTANGLE :
#             if differ[0] != ENTANGLE:
#                 score+=1
#         elif targ_att[0] == NC :
#             if differ[0] not in [ENTANGLE, NC]:
#                 score+=1
#         elif targ_att[0] == PHASING:
#             if differ[0] not in [ENTANGLE, NC]: 
#                 score+=1
#         elif targ_att[0] == BOOL:
#             if differ[0] not in [ENTANGLE, NC, PHASING]: 
#                 score+=1
#     return score


def is_strict_decrease(targ_attribute, new_targ_att):
    if targ_attribute[0] == ENTANGLE :
        return new_targ_att[0] != ENTANGLE
    elif targ_attribute[0] == NC :
        return new_targ_att[0] not in [ENTANGLE, NC]
    elif targ_attribute[0] == PHASING :
        return new_targ_att[0] not in [ENTANGLE, NC, PHASING]
    elif targ_attribute[0] == BOOL :
        return new_targ_att[0] not in [ENTANGLE, NC, PHASING]
    elif targ_attribute == IDENTICAL:
        return new_targ_att == IDENTICAL
    else :
        print(f"is_strict_decrease check routine {targ_attribute}, {new_targ_att}")
        exit()

def attribute_fill_in_critertion_multiple_IO(ss : Synthesis, targ_attribute : str, 
                                            currently_found_V_qc :cirq.Circuit, 
                                            stacked_module : List[Tuple[str,...]], 
                                            input_wise_stacked_modules,
                                            input_wise_targ_att,
                                            action_of_qc):
    # new_input_wise_stacked_module  = input_wise_stacked_modules.copy()
    differs = list()
    att_fill_in_score =  0
    curr_qc_input_wise_targ_att = dict() 
    for in_nparr, out_nparr in ss.spec.spec_object.get_io_pairs():
        res_by_curr_Vqc   = action_of_qc[tuple(in_nparr)]
        if not np.allclose(res_by_curr_Vqc, out_nparr) :
            sv_differ = sv_differ_by(res_by_curr_Vqc, out_nparr)
            # new_input_wise_stacked_module[basis_rep(in_nparr)].append(sv_differ)
            differs.append(sv_differ)
            curr_qc_input_wise_targ_att[basis_rep(in_nparr)] = sv_differ
            if is_strict_decrease(input_wise_targ_att[basis_rep(in_nparr)], sv_differ):
                att_fill_in_score+=1
        else :
            differs.append(IDENTICAL)
            curr_qc_input_wise_targ_att[basis_rep(in_nparr)] = IDENTICAL
            if is_strict_decrease(input_wise_targ_att[basis_rep(in_nparr)], IDENTICAL):
                att_fill_in_score+=1

    reprsentative_differ = representative_targ_attribute(differs) 
    # aditional fill-in criterion


    if att_fill_in_score >= math.floor(ss.spec.num_of_ios/2):
        if True:
            updated_stacked_module = stacked_module.copy()
            updated_stacked_module.append(targ_attribute)
            return updated_stacked_module, True, reprsentative_differ, None, curr_qc_input_wise_targ_att
        else :
            return None, False, None, None, None
    else :
        return stacked_module.copy(), False, None, None, None
    

def derive_next_targ_attribute_for_no_prunungs(ss : Synthesis, 
                                        action_of_qc : dict): 
    if ss.is_state_preparation: 
        in_nparr, out_nparr  =  ss.spec.spec_object.get_io_pairs()[0]
        res_by_curr_Vqc = action_of_qc[tuple(in_nparr)]
        if not np.allclose(res_by_curr_Vqc, out_nparr) :
            sv_differ = sv_differ_by(res_by_curr_Vqc, out_nparr)
        else :
            raise Exception("Should not be happendend")
        return sv_differ, None
    else :
        input_wise_att_diff_builder = dict()
        differs = list()
        for in_nparr, out_nparr in ss.spec.spec_object.get_io_pairs():
            res_by_curr_Vqc = action_of_qc[tuple(in_nparr)]
            if not np.allclose(res_by_curr_Vqc, out_nparr) :
                sv_differ = sv_differ_by(res_by_curr_Vqc, out_nparr)
                # new_input_wise_stacked_module[basis_rep(in_nparr)].append(sv_differ)
                differs.append(sv_differ)
                input_wise_att_diff_builder[str(basis_rep(in_nparr))]  = sv_differ
            else :
                differs.append(IDENTICAL)
                input_wise_att_diff_builder[str(basis_rep(in_nparr))]  = IDENTICAL
        reprsentative_differ = representative_targ_attribute(differs) 
        return reprsentative_differ, input_wise_att_diff_builder

def attribute_fill_in_criterion_state_prep(ss : Synthesis, targ_attribute : str, 
                                        currently_found_V_qc :cirq.Circuit, 
                                        stacked_module : List[Tuple[str,...]], 
                                        action_of_qc : dict,
                                        module,
                                        action_of_prev_C): 
    in_nparr, out_nparr  =  ss.spec.spec_object.get_io_pairs()[0]
    res_by_curr_Vqc = action_of_qc[tuple(in_nparr)]
    if not np.allclose(res_by_curr_Vqc, out_nparr) :
        sv_differ = sv_differ_by(res_by_curr_Vqc, out_nparr)
    else :
        raise Exception("Should not be happendend")
    # TEST : THIS IS INDEED AN BOTTLENECK
    # QC_in = action_of_prev_C[tuple(in_nparr)]
    # if not np.allclose(QC_in,res_by_curr_Vqc):
    #     attribute_of_module = sv_differ_by(QC_in, res_by_curr_Vqc)
    # end_time = time.time()
    # print("TIME", end_time-start_time)
    reprsentative_differ = sv_differ
    # if ss.is_state_preparation : # if the problem is state_prep 
    if targ_attribute[0] == ENTANGLE :
        if  reprsentative_differ[0] != ENTANGLE:
            updated_stacked_module = stacked_module.copy()
            updated_stacked_module.append(targ_attribute)
            return updated_stacked_module, True, reprsentative_differ, None 
    elif targ_attribute[0] == NC :
        if reprsentative_differ[0] not in [ENTANGLE, NC]:
            updated_stacked_module = stacked_module.copy()
            updated_stacked_module.append(targ_attribute)
            return updated_stacked_module, True, reprsentative_differ, None
    elif targ_attribute[0] == PHASING :
        if reprsentative_differ[0] == BOOL:
        # if reprsentative_differ[0] not in [ENTANGLE, NC, PHASING]: #i.e, if it is BOOL
            updated_stacked_module = stacked_module.copy()
            updated_stacked_module.append(targ_attribute)
            return updated_stacked_module, True, reprsentative_differ, None
    elif targ_attribute[0] == BOOL :
        if reprsentative_differ[0] == BOOL:
        # if reprsentative_differ[0] not in [ENTANGLE, NC, PHASING] : #i.e, if it is BOOL
            updated_stacked_module = stacked_module.copy()
            updated_stacked_module.append(targ_attribute)
            return updated_stacked_module, True, reprsentative_differ, None
    else :
        exit()
    return stacked_module.copy(), False, None, None # base case




def concrete_attribute_fill_in_critertion_multiple_IO(ss : Synthesis,
                                             targ_attribute : str, 
                                            currently_found_V_qc :cirq.Circuit, 
                                            stacked_module : List[Tuple[str,...]], 
                                            input_wise_stacked_modules,
                                            input_wise_targ_att,
                                            action_of_qc,
                                            module : cirq.Circuit,
                                            action_of_prev_C):
    # input : QC
    differs = list()
    att_fill_in_score =  0
    curr_qc_input_wise_targ_att = dict()  #  QC ;M 
    for in_nparr, out_nparr in ss.spec.spec_object.get_io_pairs():
        # module_sv = ss.simulator.simulate(module, qubit_order=ss.working_qubit, initial_state=in_nparr).final_state_vector
        #THIS becomes MC|in>
        res_by_curr_Vqc = action_of_qc[tuple(in_nparr)] 
        # THIS becomes C|in> 
        action_of_C  = action_of_prev_C[tuple(in_nparr)]
        #GOAL : To calculate Att_{C|in>}(M) = MC|in> \ominus C|in>

        if not np.allclose(res_by_curr_Vqc, out_nparr) :
            sv_differ = sv_differ_by(res_by_curr_Vqc, out_nparr) #required later module stack ups
            differs.append(sv_differ)
            curr_qc_input_wise_targ_att[basis_rep(in_nparr)] = sv_differ
            if not np.allclose(action_of_C,res_by_curr_Vqc):
                attribute_of_module = sv_differ_by(action_of_C, res_by_curr_Vqc) #MC|in> \ominus C|in>
            else : 
                continue
            if attribute_of_module[0] == input_wise_targ_att[basis_rep(in_nparr)][0] :
                att_fill_in_score+=1
        else :
            differs.append(IDENTICAL)
            curr_qc_input_wise_targ_att[basis_rep(in_nparr)] = IDENTICAL
            if np.allclose(res_by_curr_Vqc, out_nparr):
                att_fill_in_score+=1
            
    reprsentative_differ = representative_targ_attribute(differs) 

    if att_fill_in_score >= math.floor(ss.spec.num_of_ios/2):
        if True:
        # if flag_for_inputwise_targ_crit == True:
        # if  targ_attribute[0] != reprsentative_differ[0]: #TODO fitted for sqrt_X toffoli
            updated_stacked_module = stacked_module.copy()
            updated_stacked_module.append(targ_attribute)
            return updated_stacked_module, True, reprsentative_differ, None, curr_qc_input_wise_targ_att
        else :
            return None, False, None, None, None
    else :
        return stacked_module.copy(), False, None, None, None



def concrete_attribute_fill_in_criterion_state_prep_with(ss : Synthesis,
                                                        targ_attribute ,
                                                        currently_found_V_qc :cirq.Circuit, 
                                                        stacked_module : List[Tuple[str,...]], 
                                                        action_of_qc : dict,
                                                        module : cirq.Circuit,
                                                        action_of_prev_C ):

    in_nparr, out_nparr  =  ss.spec.spec_object.get_io_pairs()[0]
    #MQC |in> is as res_by_curr_Vqc
    res_by_curr_Vqc = action_of_qc[tuple(in_nparr)]
    if not np.allclose(res_by_curr_Vqc, out_nparr) :
        sv_differ = sv_differ_by(out_nparr, res_by_curr_Vqc)
    else :
        return None, True, None, None
        # raise Exception("It is identity! ")
    reprsentative_differ = sv_differ
    # module_sv = ss.simulator.simulate(module, qubit_order=ss.working_qubit, initial_state=res_by_curr_Vqc).final_state_vector
    # QC|in> is as follwong
    QC_in = action_of_prev_C[tuple(in_nparr)]
        
    if not np.allclose(QC_in,res_by_curr_Vqc):
        attribute_of_module = sv_differ_by(QC_in, res_by_curr_Vqc)
    else :
        return stacked_module.copy(), False, None, None
    if attribute_of_module[0] == targ_attribute[0]: 
        if targ_attribute[0] == ENTANGLE :
            if attribute_of_module[1] == targ_attribute[1]:
                updated_stacked_module = stacked_module.copy()
                updated_stacked_module.append(targ_attribute)
                return updated_stacked_module, True, reprsentative_differ, None
        else :
            updated_stacked_module = stacked_module.copy()
            updated_stacked_module.append(targ_attribute)
            return updated_stacked_module, True, reprsentative_differ, None
    else :
        return stacked_module.copy(), False, None, None
    return stacked_module.copy(), False, None, None