import cirq
import numpy as np
from typing import Tuple, Dict, List, Union

if __name__ == "state_search": 
    from util.utils import *
    from state import *
    from rules import *
    from guide import *
    from score import *
    from synthesis import Synthesis
    from synthesis_spec.specification import Spec
    from module_restrain import *
else : 
    from qsyn.util.utils import *
    from qsyn.state import *
    from qsyn.rules import *
    from qsyn.guide import *
    from qsyn.score import *
    from qsyn.synthesis import Synthesis
    from qsyn.synthesis_spec.specification import Spec
    from qsyn.module_restrain import *
from timeit import default_timer as timer, repeat




class SearchStateQc():
    # upper_state, rule, curr_state, stacked_module, curr_targ_differ
    def __init__(self, ss , upper_state, applied_rule, curr_qc, action_of_qc, concrete_curr_qc, stacked_modules, targ_diff, 
                input_wise_stacked_module = None,
                input_wise_targ_att = None ,
                seq_of_att_of_stacked_moduels = None,
                no_pruning = False ) -> None:
        self.ss = ss
        self.is_state_preparation = ss.is_state_preparation
        self.upper_state = upper_state
        self.applied_rule = applied_rule
        self.curr_qc = curr_qc
        self.action_of_qc = action_of_qc

        bool_li_is_classic = list()
        for _, curr_out_arr in action_of_qc.items() :
            bool_li_is_classic.append(is_classic_state(curr_out_arr, atol =1e-05))
        self.is_curr_out_arr_all_classic = all(bool_li_is_classic)
        self.concrete_curr_qc=concrete_curr_qc
        self.stacked_modules = stacked_modules
        self.no_pruning = no_pruning
        if not self.is_state_preparation :
            assert seq_of_att_of_stacked_moduels != None and isinstance(seq_of_att_of_stacked_moduels, list)
        self.seq_of_att_of_stacked_moduels = seq_of_att_of_stacked_moduels

        self.targ_diff  = targ_diff
        
        if self.is_state_preparation : 
            self.input_wise_stacked_module = None
            self.input_wise_targ_att = None
        else : 
            # assert input_wise_stacked_module
            self.input_wise_stacked_module = input_wise_stacked_module
            self.input_wise_targ_att = input_wise_targ_att
            already_added_att = list()
            att_counter = { BOOL : 0, PHASING : 0, NC :0, ENTANGLE : 0}
            for att_val in (self.input_wise_targ_att.values()):
                if att_val!= IDENTICAL :
                    already_added_att.append(att_val[0])
                    att_counter[att_val[0]]  = att_counter[att_val[0]]+1
            max_num_of_cnt = max(att_counter.values())
            atts_of_max_count =list(filter( lambda x :att_counter[x]==max_num_of_cnt ,   list(att_counter.keys())   ))
            self.atts_of_max_count = atts_of_max_count

        self.cost = self.cost_val()
        assert self.cost or self.cost ==0 

    def cost_val(self):
        # smaller cost higher priority
        cost = None
        if self.is_state_preparation:
            if self.stacked_modules : 
                prefered_seq_val  = len(self.stacked_modules)>=2 and self.targ_diff[0] in [BOOL]
                if prefered_seq_val : return (-10000,)
            cost = (len(self.stacked_modules)
                    ,count_gate(self.concrete_curr_qc)
                    ,score_of_differ(self.targ_diff)
                    )
            return cost
        elif (not self.is_state_preparation) :
            # if self.targ_diff[0] by implementation it is the only non-idential attribute gap among inputs
            # memo : not in [s[0] for s in self.stacked_modules]
            # if False:
            if self.targ_diff[0] == BOOL and len(self.stacked_modules)>=4 and self.targ_diff[0]  :
                cost = (-10000,)
            else : 
                # computing cost of attributes
                # also conpute cost of score = check number of IDENTITY in self.input_wise_targ_att
                cost_omega = 0
                score_E = 0
                score_for_non_triv =0 
                for input_spec, att in self.input_wise_targ_att.items():
                    if att == IDENTICAL:
                        if input_spec in self.ss.spec.spec_object.get_identity_input_specs():
                            score_E += 1 # place holder can be adjusted later
                        else :
                            score_E += 1
                            score_for_non_triv  +=1

                # cal cost_omega
                for att in self.atts_of_max_count:
                    if att in [ENTANGLE, NC, PHASING, BOOL]:
                        cost_omega+= score_of_differ_for_multi_IO((att, None))
                    else : raise Exception(f"Invalt Att val {att} ")

                num_of_non_iden_spec= len(self.ss.spec.spec_object.get_identity_input_specs())
                cost = (len(self.stacked_modules)
                        ,   count_gate(self.concrete_curr_qc)
                        , - score_for_non_triv
                        ,  cost_omega
                        # , not (score_for_non_triv >= 1)
                        # , -score_E
                        # ,not (score_E >= self.ss.spec.num_of_ios/4)
                        # , self.targ_diff[0] == ENTANGLE
                        # ,set(self.atts_of_max_count).intersection(set(self.seq_of_att_of_stacked_moduels)) != set()
                        # ,self.targ_diff[0]  in [s[0] for s in self.stacked_modules]
                        # ,score_E
                        )
            return cost

    def __eq__(self, __o: object) -> bool:
        return self.concrete_curr_qc == __o.concrete_curr_qc

    def __hash__(self) -> int:
        return hash(str(self.concrete_curr_qc))

class StateSearchSynthesis(Synthesis):

    def is_involution(self, circuit_to_check: cirq.Circuit, gate_op_to_check: cirq.ops.GateOperation):
        if len(circuit_to_check) == 0:
            return False

        gate_op = None
        gate_op_to_check_qubits = sorted( list([q.x for q in gate_op_to_check.qubits]))
        rep_qubit_index = gate_op_to_check_qubits[0]
        assert isinstance(rep_qubit_index, int)
        for moment_loc in reversed(range(len(circuit_to_check))):
            gate_opt_on_same_qubit = circuit_to_check.operation_at(self.working_qubit[rep_qubit_index], moment_loc)
            if gate_opt_on_same_qubit:
                gate_op = gate_opt_on_same_qubit
                moment_idx_of_gate_op = moment_loc
                break
            else:
                continue
        if not gate_op:
            return False
        elif (gate_op == gate_op_to_check
            and gate_op.gate in self.component_prior["INVOLUTIONARY"]):
            return True
        elif (gate_op.qubits == gate_op_to_check.qubits):
            if gate_op.gate in self.component_prior["INVERSES"].keys():
                if (gate_op_to_check.gate == self.component_prior["INVERSES"][gate_op.gate]):
                    return True
            if gate_op_to_check.gate in self.component_prior["IDENTITY_N"].keys() and str(gate_op_to_check) == str(gate_op):
                iden_N = self.component_prior["IDENTITY_N"][gate_op_to_check.gate]
                collected_moment_idx = list()
                for moment_loc in reversed(range(len(circuit_to_check))):
                    temp_gate_app = circuit_to_check.operation_at(gate_op_to_check.qubits[0], moment_loc)
                    if temp_gate_app == gate_op_to_check : collected_moment_idx  = [moment_loc] + collected_moment_idx
                if len(collected_moment_idx) == iden_N-1 and check_continuity(collected_moment_idx) and moment_idx_of_gate_op in collected_moment_idx:
                    return True
        else:
            return False

    def __init__(self, spec: Spec, 
                    prepare_naive_module_gen : bool, 
                    entangle_template : bool,
                    module_gate_num : int):
        super().__init__(spec)
        self.MAX_MODULE_NUM = 5
        print(f"MAX MODULE NUM is : {self.MAX_MODULE_NUM}")
        self.cache_for_state_eval = dict()
        io_specs = self.spec.spec_object.get_io_pairs()
        self.num_io_spec = len(io_specs)
        self.num_critical_ios = 0
        for io_pair in io_specs:
            in_arr, out_arr = io_pair
            if is_identity_mapping(in_arr, out_arr):
                self.num_critical_ios += 1

        self.component_prior = component_priors(components=self._spec.component_gates, is_att_guide=True)
            
        self.act_prop_schems_of_component_gate = set([get_action_prop_scheme_of_gate(g) for g in self.spec.component_gates])
        # self.act_prop_schems_of_component_gate = scheme
        self.is_loss_info_exist = [LOSS_INFO in scheme for scheme in self.act_prop_schems_of_component_gate]
        
        self.min_sp_gate_qubit_num = min([  g.num_qubits() for g in self.component_prior[NC]]  + [  g.num_qubits() for g in self.component_prior[ENTANGLE]])
        self.rules = Rules(spec=self.spec, 
                            do_rule_prune= not any(self.is_loss_info_exist),
                            min_sp_gate_qubit_num = self.min_sp_gate_qubit_num )

        
        self.valid_num_qubits = [g.num_qubits() for g in self.spec.component_gates]
        self.min_valid_num_qubits = min(self.valid_num_qubits)
        self.is_state_preparation = len(io_specs) == 1
        self.representative_basis_sv_res = dict()  
        for gate in self.spec.component_gates :
            self.representative_basis_sv_res[gate] = derive_representative_res_in_comp_basis(gate)

        if prepare_naive_module_gen :
            self.module_qc_collection_by_level = {
                0 : [cirq.Circuit()],
                1 : list() ,
                2 : list()   ,
                3 : list()   ,
                4 : list()   ,
                5 : list()   ,
                6 : list()    ,
                7 : list()                    
                }
            self.module_qc_collection_by_att_level = {
                ENTANGLE : { 
                    0 : [cirq.Circuit()],
                    1 : [],
                    2 : [],
                    3 : [],
                    4 : []
                },
                NC :  { 
                    0 : [cirq.Circuit()],
                    1 : [],
                    2 : [],
                    3 : [],
                    4 : []
                },
                PHASING :  { 
                    0 : [cirq.Circuit()],
                    1 : [],
                    2 : [],
                    3 : [],
                    4 : []
                },
                BOOL :  { 
                    0 : [cirq.Circuit()],
                    1 : [],
                    2 : [],
                    3 : [],
                    4 : []
                }
            }
            print("prior module generating")
            print(f"module_gate_num is {module_gate_num}")
            self.module_gate_num= module_gate_num
            i = 3
            for i in range(1, i + 1):
                assert self.module_qc_collection_by_level[i] == []
                for work_circuit in self.module_qc_collection_by_level[i-1]:
                    for gate_op in self.gate_operations:
                        new_circuit = work_circuit.copy()
                        if not self.is_involution(new_circuit, gate_op):
                            new_circuit.append(gate_op)
                            self.module_qc_collection_by_level[i].append(new_circuit)
                print(f"prior module gen at level {i} done")
            print("prior module gen done")

    def do_qc_eval_criterion(self, mbs): 
        for m in mbs.moments : 
            for b in m.blocks :
                if isinstance(b, MyProperty):
                    if ENTANGLE in  b.action_property.values():
                        return True
                    if len(b.action_property.keys()) not in self.valid_num_qubits or len(b.action_property.keys()) < self.min_valid_num_qubits  :
                        return False
                    elif not b.contains_none() and not (any(self.is_loss_info_exist)) :
                        act_prop_schem_match_booli  = [sorted(list(b.action_property.values()), key= lambda s : str(s)) == sorted(act_prop_schem, key= lambda s : str(s)) for act_prop_schem in self.act_prop_schems_of_component_gate]
                        if not any(act_prop_schem_match_booli):
                            return False
        return True
        

    def modular_search(self, 
                                 init_qc: Union[MomentBasedState, cirq.Circuit], 
                                 rule_depth=2,
                                 stacked_module = None,
                                 concrete_criterion = False,
                                 naive_module_gen = False, 
                                 entangle_template = False,
                                 no_pruning = False,
                                 no_mutate = False,
                                 start_time = None,
                                 timeout = 300) ->Union[None, cirq.Circuit]:
        init_qc_as_concrete = None
        print("Modular Search state search!")
        print(f"concrete_criterion : {concrete_criterion}")
        print(f"naive_module_gen : {naive_module_gen}")
        if not stacked_module :
            stacked_module = list()
        if isinstance(init_qc, cirq.Circuit):
            init_V_qc_as_state = from_qc_to_moment_based_state(init_qc)
            init_qc_as_concrete = init_qc
        elif init_qc == None:  # if it is initial state
            init_V_qc_as_state = MomentBasedState(moments=None)
            init_qc_as_concrete = cirq.Circuit()
        elif isinstance(init_qc, MomentBasedState) and (not init_qc.is_empty_state()) :
            init_V_qc_as_state = init_qc
            init_qc_eval_res  = init_qc.evaluate(self, self.spec.component_gates, self.working_qubit)
            assert len(init_qc_eval_res) == 1
            init_qc_as_concrete = init_qc_eval_res[0]
        else:
            print("Initial input of V_qc for state search is not concrete state. Did you implemented?")
            exit()
        
        print("**Initial State is**")
        print(init_V_qc_as_state)
        print("********************")
        work_list_of_qc_state = list() # this wokring list should be list of concrete qcs
        already_checked_qc_state = set()
        cnt = 0  
        init_targ_attribute_val = init_targ_attribute(self, init_V_qc_as_state ) # Attribute Difference val

        ### Initialization 
        if self.is_state_preparation:
            in_arr, out_arr =  self.spec.spec_object.get_io_pairs()[0]
            action_builder = { tuple(in_arr) : in_arr }
            init_search_state_qc = SearchStateQc(
                                ss = self,
                                upper_state = None,
                                applied_rule= None,
                                curr_qc = init_V_qc_as_state,
                                action_of_qc= action_builder,
                                concrete_curr_qc= init_qc_as_concrete,
                                stacked_modules= stacked_module,
                                targ_diff = init_targ_attribute_val
            )
        else : 
            init_input_wise_stacked_module_builder = dict()
            init_input_wise_targ_att = dict()
            action_builder = dict()
            for in_arr, out_arr in self.spec.spec_object.get_io_pairs():
                init_input_wise_stacked_module_builder[basis_rep(in_arr)] = list()
                action_builder[tuple(in_arr)] = in_arr
                if not np.allclose(in_arr, out_arr, atol=1e-04):
                    init_input_wise_targ_att[basis_rep(in_arr)] = sv_differ_by(in_arr, out_arr)
                else :
                    init_input_wise_targ_att[basis_rep(in_arr)] = IDENTICAL
            init_search_state_qc = SearchStateQc(
                                ss = self,
                                upper_state = None,
                                applied_rule= None,
                                curr_qc = init_V_qc_as_state,
                                action_of_qc=action_builder,
                                concrete_curr_qc= init_qc_as_concrete,
                                stacked_modules= stacked_module,
                                targ_diff = init_targ_attribute_val,
                                input_wise_stacked_module=init_input_wise_stacked_module_builder,
                                input_wise_targ_att = init_input_wise_targ_att,
                                seq_of_att_of_stacked_moduels = list()
            )
        ##### initialization done
        dict_num_of_passed_module = dict()
        work_list_of_qc_state.append(init_search_state_qc)
        while True:
            cnt += 1
            if cnt % 100 == 0:
                print(f"Cnt,,{cnt}")
            if len(work_list_of_qc_state) == 0:
                print(f"Not Found. State Search is Over.")
                return None, None
            if len(work_list_of_qc_state)<=10:
                # if work_list_of_qc_state is small enough, duplicate removal is fast and do so
                # else could be costly,,,but heuristcly set the `small` to be of length less than 10
                work_list_of_qc_state = list(set(work_list_of_qc_state))
                work_list_of_qc_state.sort(key=lambda w : str(w.concrete_curr_qc))
            
            work_list_of_qc_state.sort(key = lambda s : s.cost)     # sorting custom       
            
            curr_state = work_list_of_qc_state.pop(0)
            dict_num_of_passed_module[str(curr_state)] = 0
            assert isinstance(curr_state, SearchStateQc)
            already_checked_qc_state.add(str(curr_state.curr_qc))

            if len(curr_state.stacked_modules) <= self.MAX_MODULE_NUM-1:
                print("===================")
                print("curr state\n")
                print(curr_state.curr_qc)
                print("curr_targ_differ", curr_state.targ_diff)
                print("check_order : stacked, targ_diff", curr_state.stacked_modules, curr_state.targ_diff)
                print("cost", curr_state.cost)
                # print("ord value", curr_state.cost)

                # flag_at_least_one_module_filled = False
                if naive_module_gen and ( curr_state.targ_diff[0] not in [ENTANGLE] or  (curr_state.targ_diff[0] == ENTANGLE and not entangle_template) ) :
                    generator_based = False
                    curr_module_pool = self.init_naive_module_gen(curr_qc=curr_state.curr_qc, targ_diff = curr_state.targ_diff)
                    if curr_state.is_state_preparation:
                        curr_module_pool = self.optimize_module_candidate_for_naive_gen(curr_state, curr_module_pool,  curr_state.targ_diff)
                    else :
                        module_pool_collector = list()
                        print("Choose Best Count Attribute!")
                        print(curr_state.input_wise_targ_att)
                        print(curr_state.atts_of_max_count)
                        for att_of_max_count in curr_state.atts_of_max_count:
                            if att_of_max_count != PHASING :
                                best_att_of_max_count = (att_of_max_count, set())
                                for att_val in curr_state.input_wise_targ_att.values():
                                    if att_val[0] == att_of_max_count :
                                        if att_of_max_count == ENTANGLE:
                                            if best_att_of_max_count[1] == set():
                                                to_comapre_att_info = ()
                                            else:
                                                to_comapre_att_info = best_att_of_max_count[1]
                                            if best_att_of_max_count[0] == att_val[0] and to_comapre_att_info <= att_val[1]:
                                                    best_att_of_max_count = att_val
                                        else : 
                                            if best_att_of_max_count[1] and att_val[1] :
                                                if best_att_of_max_count[0] == att_val[0] and best_att_of_max_count[1] <= att_val[1]:
                                                    best_att_of_max_count = att_val
                                            else :
                                                if best_att_of_max_count[0] == att_val[0] :
                                                    best_att_of_max_count = att_val
                            else :
                                best_att_of_max_count = (PHASING, None)
                            # temp_collect = [(m, best_att_of_max_count[0]) for  m in self.optimize_module_candidate_for_naive_gen(curr_module_pool, best_att_of_max_count)]
                            module_pool_collector +=  [(m, best_att_of_max_count[0]) for  m in self.optimize_module_candidate_for_naive_gen(curr_state, curr_module_pool,  best_att_of_max_count)]
                        curr_module_pool = module_pool_collector

                elif not naive_module_gen :
                    raise NotImplementedError("Given option '--naive_module_gen' for current usage. This branch was only for holder of later development")

                elif naive_module_gen and curr_state.targ_diff[0]  == ENTANGLE and entangle_template :
                    # generating ENTANGLING MODULES
                    generator_based = True
                    curr_module_pool = self.Gamma(curr_qc=curr_state.curr_qc, targ_diff = curr_state.targ_diff, max_level=rule_depth)
                    print(curr_state.input_wise_targ_att)
                while True:
                    curr_time = timer()
                    if curr_time - start_time > timeout:
                        print("TIMEOUT")
                        print(f"TIMEOUT was {timeout}")
                        return None, None
                    if generator_based :
                        module = next(curr_module_pool, None)
                        if not self.is_state_preparation :
                            modules_input_indept_att= ENTANGLE
                            # currently, if generator based than module is entanlging
                    else :
                        try :
                            item_of_module = curr_module_pool.pop(0)
                            if isinstance(item_of_module, tuple) and (not self.is_state_preparation)  : # handle multio problem
                                module = item_of_module[0]
                                modules_input_indept_att= item_of_module[1]
                            elif isinstance(item_of_module, cirq.Circuit) and self.is_state_preparation:
                                module = item_of_module
                            else :
                                raise Exception(f"Invalid Case of handling item_of_module {item_of_module}")
                        except IndexError :
                            module = None

                    if not naive_module_gen and module == None :
                        print("** cnt_num_of_passed_module **", dict_num_of_passed_module[str(curr_state)])
                        break
                    if naive_module_gen and module == None :
                        break
                    
                    qc = curr_state.concrete_curr_qc.copy()
                    qc.append([ y for x  in module for y in x ])
                    is_passed, action_of_qc = self.test_synthesis_att_guided(circuit_to_check=qc, phase_equiv=self.spec.equiv_phase)
                    # action of QC;M
                    if is_passed and not concrete_criterion: 
                        print("MODULES were ", curr_state.stacked_modules + [curr_state.targ_diff])

                        return qc, curr_state.stacked_modules + [curr_state.targ_diff]
                    if not concrete_criterion and (not no_pruning):
                        if self.is_state_preparation : 
                            updated_stacked_module, is_attribute_filled, next_target_differ, _ =   attribute_fill_in_criterion_state_prep(self, 
                                                                                                                                        curr_state.targ_diff,   
                                                                                                                                        qc,
                                                                                                                                        curr_state.stacked_modules,
                                                                                                                                        action_of_qc,
                                                                                                                                        module,
                                                                                                                                        curr_state.action_of_qc)
                        else : 
                            updated_stacked_module, is_attribute_filled, next_target_differ, new_input_wise_stacked_modules, new_input_wise_targ_att = attribute_fill_in_critertion_multiple_IO(self, 
                                                                                                                                                    curr_state.targ_diff, 
                                                                                                                                                    qc, 
                                                                                                                                                    curr_state.stacked_modules,
                                                                                                                                                    curr_state.input_wise_stacked_module,
                                                                                                                                                    curr_state.input_wise_targ_att,
                                                                                                                                                    action_of_qc)
                    elif concrete_criterion and (not no_pruning) :
                        if self.is_state_preparation :
                            updated_stacked_module, is_attribute_filled, next_target_differ, _ =   concrete_attribute_fill_in_criterion_state_prep_with(self, 
                                                                                                                                        curr_state.targ_diff,   
                                                                                                                                        qc,
                                                                                                                                        curr_state.stacked_modules,
                                                                                                                                        action_of_qc,
                                                                                                                                        module,
                                                                                                                                        curr_state.action_of_qc)
                        else : 
                            new_seq_of_att = curr_state.seq_of_att_of_stacked_moduels + [modules_input_indept_att]
                            updated_stacked_module, is_attribute_filled, next_target_differ, new_input_wise_stacked_modules, new_input_wise_targ_att = concrete_attribute_fill_in_critertion_multiple_IO(self, 
                                                                                                                        curr_state.targ_diff, 
                                                                                                                        qc, 
                                                                                                                        curr_state.stacked_modules,
                                                                                                                        curr_state.input_wise_stacked_module,
                                                                                                                        curr_state.input_wise_targ_att,
                                                                                                                        action_of_qc,
                                                                                                                        module,
                                                                                                                        curr_state.action_of_qc)
                    if no_pruning :
                        is_attribute_filled = None
                        if self.is_state_preparation :
                                updated_stacked_module, is_attribute_filled, new_input_wise_stacked_modules = curr_state.stacked_modules + [curr_state.targ_diff], None, None, 
                        else: # then multi-io problem
                                updated_stacked_module, is_attribute_filled, new_input_wise_stacked_modules, new_input_wise_targ_att = curr_state.stacked_modules + [curr_state.targ_diff], None, None, None
                        next_target_differ, new_input_wise_targ_att = derive_next_targ_attribute_for_no_prunungs(ss = self,action_of_qc =action_of_qc )
                    if no_pruning or is_attribute_filled : 
                        if is_passed and concrete_criterion:
                            print("MODULES were ", curr_state.stacked_modules + [curr_state.targ_diff])
                            return qc
                        dict_num_of_passed_module[str(curr_state)] += 1
                        if self.is_state_preparation : 
                            to_append_search_state = SearchStateQc(
                                ss = self,
                                upper_state=None,
                                applied_rule = None,
                                curr_qc = from_qc_to_moment_based_state(qc),
                                action_of_qc = action_of_qc,
                                concrete_curr_qc= qc,
                                stacked_modules=updated_stacked_module,
                                targ_diff=next_target_differ,
                                no_pruning= no_pruning
                            )
                        else : 
                            new_seq_of_att = curr_state.seq_of_att_of_stacked_moduels + [modules_input_indept_att]
                            to_append_search_state = SearchStateQc(
                                ss = self,
                                upper_state=None,
                                applied_rule = None,
                                curr_qc = from_qc_to_moment_based_state(qc),
                                action_of_qc = action_of_qc,
                                concrete_curr_qc= qc,
                                stacked_modules=updated_stacked_module,
                                targ_diff=next_target_differ,
                                input_wise_stacked_module=new_input_wise_stacked_modules,
                                input_wise_targ_att=new_input_wise_targ_att,
                                seq_of_att_of_stacked_moduels= new_seq_of_att,
                                no_pruning= no_pruning
                            )
                        if len(to_append_search_state.stacked_modules) <= self.MAX_MODULE_NUM-1:
                            work_list_of_qc_state.append(to_append_search_state)
                print("** cnt_num_of_passed_module after seeing all modules in pool **",  dict_num_of_passed_module[str(curr_state)])


    def module_of_inseparables(self, num_gate:int):
        ins_gates = [ g for g in self.spec.component_gates if g.num_qubits()>=2]
        to_return = list()
        module_collect_by_level = {
            0 : [cirq.Circuit()],
            1: list(),
            2: list(),
            3: list()
        }
        i=3
        for i in range(1, i+1):
            for work_circuit in module_collect_by_level[i-1]:        
                for gate_op in self.gate_operations:
                    if gate_op.gate in ins_gates:
                        new_circuit = work_circuit.copy()
                        if not self.is_involution(new_circuit, gate_op):
                            new_circuit.append(gate_op)
                            module_collect_by_level[i].append(new_circuit)
                            to_return.append(new_circuit)
        return to_return 
        
        # input()

    def Gamma(self, curr_qc : MomentBasedState, targ_diff, max_level : int): # the module generation
        gen_mbs = list()  #state, level
        one_level = [(MomentBasedState(moments = [x]), 1)for x in biased_moment_set_of_property(targ_diff[0], targ_diff[1], self) ]
        gen_mbs += one_level
        while gen_mbs:
            mbs, level_mbs = gen_mbs.pop(0)
            if  self.do_qc_eval_criterion(mbs) :
                eval_res_qcs = mbs.evaluate(self, self.spec.component_gates, self.working_qubit)
                for res_mod in eval_res_qcs:
                    yield res_mod
            if level_mbs+1 <= max_level:
                for rule_mode, rule_param in self.rules.prop_decompose_rule_space(mbs):
                    decomposed = self.rules.run_rule(mode=rule_mode, params=rule_param)
                    if isinstance(decomposed, list):
                        for z in decomposed :
                            gen_mbs.append((z, level_mbs+1))
                    elif isinstance(decomposed, MomentBasedState):
                        gen_mbs.append((decomposed, level_mbs+1))
                    elif decomposed == None:
                        continue
                    else:
                        raise Exception(f"invalid value for rule_param {decomposed}")

    def init_naive_module_gen(self, curr_qc : MomentBasedState, targ_diff): # the module generation
        to_return = list()
        max_gate_num_for_module = None
        if targ_diff[0] == NC:
            max_gate_num_for_module  = self.module_gate_num -1
        if targ_diff[0] == BOOL and self.module_gate_num >=4:
            max_gate_num_for_module  = self.module_gate_num -1
            print(max_gate_num_for_module)
        else :
            max_gate_num_for_module  = self.module_gate_num
        for key in list(sorted(list(self.module_qc_collection_by_level.keys()))):
            if key != 0  and key <= max_gate_num_for_module:
                 to_return += self.module_qc_collection_by_level[key]
        return to_return


    def optimize_module_candidate_for_naive_gen(self, curr_state, modules, targ_diff):
        print(f"GEN for {targ_diff}")
        assert targ_diff[0] in [ NC, PHASING, BOOL]
        to_return = list()
        if targ_diff[0] == NC :
            collect = list()
            for m in modules:
                if restrain_to_sp_modules(self, m, targ_diff):
                    collect.append(m)
            print("GEN DONE", len(collect))
            return collect
        elif targ_diff[0] == PHASING :
            collect = list()
            for m in modules:
                if restrain_to_phase_modules(self,  m, targ_diff):
                    collect.append(m)
            print("GEN DONE", len(collect))
            return collect

        elif targ_diff[0] == BOOL:
            collect =   list()
            for m in modules :
                if restrain_to_bool_modules(self,curr_state, m, targ_diff):
                    collect.append(m) 
            # generate +1 level by appending soem gate operations on m
            if self.module_gate_num >=4:
                plus_one_level_collects = list()
                plus_one_level_collects += collect
                print("Further Evolve Start")
                for m in reversed(collect) :
                    if count_gate(m)>=3:
                        plus_one_level_collects += constrained_module_evolve(self, m)
                    elif count_gate(m)==2:
                        break
                print("GEN DONE", len(plus_one_level_collects))

                if len(plus_one_level_collects) >= 30000:
                    print("the generated module is too many,,,maybe module of smaller length will be enough")
                    # collect.reverse() # try non trivial one first
                    return collect
                else:
                    # plus_one_level_collects.reverse()
                    return plus_one_level_collects
            else:
                print("GEN DONE", len(collect))
                # collect.reverse()
                return collect


def score_of_differ(differ_by : str):
    # heuristic for cost-guidance
    # set integer to be Entangle < PHASING= BOOL < NC for single IO case
    # according to our decreasing modularization, Ent>NC,Phase,BOOL order modules
    # fill Entangle first, others later
    # NC is too big, prioritize BOOL and PHASE
    if differ_by[0] == ENTANGLE :  
        return -1
    if differ_by[0] == PHASING :
        return 0
    if differ_by[0] == BOOL :
        return 0
    if differ_by[0] == NC :
        return 2
    print(f"order computation for {differ_by} is invalid. exit program.")
    exit()

def score_of_differ_for_multi_IO(differ_by : str):
    # heuristic for cost-guidance
    # no specific order assumption in modularizaiton in multi IO case
    # order by smallest $M_\omega$
    # for multi IO case
    if differ_by[0] == ENTANGLE :  
        return 2
    if differ_by[0] == PHASING :
        return 0
    if differ_by[0] == BOOL :
        return 0
    if differ_by[0] == NC :
        return 1
    print(f"order computation for {differ_by} is invalid. exit program.")
    exit()
