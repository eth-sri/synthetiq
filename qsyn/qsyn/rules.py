from xml.dom import NotSupportedErr
from cirq.ops import moment
if __name__ == "rules": 
    from state_search import *
    from state import *
    from set_synthesis import SPEC_OBJ, load_benchmark
    from util.utils import *
    from util.state_search_utils import *
    from util.rule_utils import *
else :
    from qsyn.state_search import *
    from qsyn.state import *
    from qsyn.set_synthesis import SPEC_OBJ, load_benchmark
    from qsyn.util.utils import *
    from qsyn.util.state_search_utils import *
    from qsyn.util.rule_utils import *
from itertools import product
import cirq
import pprint

SPAN = "span"
OPEN_AND_CLOSE_B = "open_and_close_B"
OPEN_AND_CLOSE_P = "open_and_close_P"
OPEN_AND_CLOSE_NC = "open_and_close_NC"
ADD = "add"
TRIPLE_SPAN = "triple_span"

class InvalidDecomposeException(Exception):
    def __init__(self, msg):
        super().__init__(msg)

class Rules:    

    '''
    Rules are either 'State -> List of States' or 'State -> State'
    '''

    def __init__(self, spec=Spec, do_rule_prune = True, min_sp_gate_qubit_num : int = None) -> None:
        self.queried_states = list()
        self.spec = spec
        gate_num_qubits = [gate.num_qubits() for gate in self.spec.component_gates]
        minimum_valid_num_qubits = min(gate_num_qubits)
        self.possible_moments_of_properties     = possible_moments_of_properties(qreg_size=self.spec.qreg_size, minimum_valid_num_qubits=minimum_valid_num_qubits)
        self.possible_moments_of_mux_properties = possible_moments_of_mux_properties(qreg_size=self.spec.qreg_size, minimum_valid_num_qubits=minimum_valid_num_qubits)
        self.possible_moments_of_act_properties = possible_moments_of_act_properties(qreg_size=self.spec.qreg_size, minimum_valid_num_qubits=minimum_valid_num_qubits)
        self.possible_act_property_blocks       = possible_act_property_blocks(qreg_size=self.spec.qreg_size, minimum_valid_num_qubits=minimum_valid_num_qubits)
        # self.inseparable_gate_qubit_num = min([gate.num_qubits() for gate in self.spec.component_gates if gate.num_qubits()>=2])
        self.inseparable_gate_qubit_num = [ i for i in gate_num_qubits if i>=2 ]
        self.min_sp_gate_qubit_num = min_sp_gate_qubit_num
        # Rule Selection Routine
        print("Rule Selection Rotuine")
        self.C_RULE_ID =  [TRIPLE_SPAN, SPAN, OPEN_AND_CLOSE_B, OPEN_AND_CLOSE_P, OPEN_AND_CLOSE_NC] 
        # should check OPEN_AND_CLOSE_B, OPEN_AND_CLOSE_P, OPEN_AND_CLOSE_NC is not redundant
        self.P_RULE_ID =  [SPAN, OPEN_AND_CLOSE_B]   
        # should check OPEN_AND_CLOSE_B is not redundant
        self.NC_RULE_ID = [ADD]
        # if redudant, remove it

        print("C_RULE_ID", self.C_RULE_ID)
        print("P_RULE_ID", self.P_RULE_ID)
        print("NC_RULE_ID", self.NC_RULE_ID)
        
    def valid_rule_id_for_block(self, block: MyProperty) -> List[str]:
        if len(block.action_property.keys()) == 1:
            prop_to_decompose = list(block.action_property.values())[0]
            if prop_to_decompose == BOOL:
                return self.C_RULE_ID
            if prop_to_decompose == PHASING:
                return self.P_RULE_ID
            if prop_to_decompose == NC:
                return self.NC_RULE_ID
        else:
            # TODO : HardCoded
            return None

    def rule_to_string(self, rule_mode: str, rule_params: dict):
        pp = pprint.PrettyPrinter(indent=4)
        return f"Rule : {rule_mode}" + "\n" + "Params : " + pp.pformat(rule_params)

    def run_rule(self, mode: str, params: dict) -> List[MomentBasedState]:
        if mode == "replace":
            return self.replace_into_prop_block_rule(**params)
        elif mode == "append":
            return self.add_props_in_moment_rule(**params)
        elif mode == "mux_decompose":
            return self.muitplexor_decompose(**params)
        elif mode == "prop_decompose":
            return self.property_decompose(**params)
        else:
            raise Exception("Rule not run")

    def print_rule_spce(self, state: MomentBasedState) -> List[Tuple[str, dict]]: #printing whole rule space
        rules = self.rule_space(state)
        for rule_mode, rule_param in rules:
            print(self.rule_to_string(rule_mode=rule_mode, rule_params=rule_param))


    def append_rule_space(self, state: MomentBasedState) -> List[Tuple[str, dict]]:
        applicable_add_block_rules = list()
        for idx in [-1, len(state.moments) - 1]:
            for moment in self.possible_moments_of_act_properties:
                params = dict([("state", state),
                               ("moment_to_append", moment),
                               ("append_position", (idx, idx + 1))])
                applicable_add_block_rules.append(("append", params))
        return applicable_add_block_rules

    def prop_decompose_rule_space(self, state: MomentBasedState) -> List[Tuple[str, dict]]:
        assert(isinstance(state, MomentBasedState))
        applicable_prop_decomposition_rule = list()
        for moment_idx, moment in enumerate(state.moments):
            for block_idx, block in enumerate(moment.blocks):
                if isinstance(block, MyProperty) and block.action_property:
                    valid_rule_ids = self.valid_rule_id_for_block(block=block)
                    if not valid_rule_ids:
                        params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", None), ("block_idx", block_idx)])
                        applicable_prop_decomposition_rule.append(("prop_decompose", params))
                    else:
                        for rule_id in valid_rule_ids:
                            params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", rule_id), ("block_idx", block_idx)])
                            applicable_prop_decomposition_rule.append(("prop_decompose", params))
        return applicable_prop_decomposition_rule

    def replace_rule_space(self, state : MomentBasedState) -> List[Tuple[str, dict]]:
        applicable_replace_rule_space = list()
        for moment_idx, moment in enumerate(state.moments):
            for block_idx, block in enumerate(moment.blocks):
                 if isinstance(block, cirq.Operation):
                     for replace_to in self.possible_act_property_blocks :
                        params = {
                            "state"        : state,
                            "replace_to"   : replace_to,
                            "moment_idx"   : moment_idx,
                            "block_idx"    : block_idx
                        }
                        applicable_replace_rule_space.append(("replace", params))
        return applicable_replace_rule_space

    def rule_space(self, state: MomentBasedState) -> List[Tuple[str, dict]]:
        raise NotSupportedErr("rule_space feature not supported anymore?")
        assert isinstance(state, MomentBasedState)
        applicable_abstract_rules = list()
        # set of all possible paramters of abstrac_rule will defined set of all possible abstract_rules on given state
        # hence we define the space as set of args(or parameters) instead of set(class) of functional values

        # # 1. block abstraction
        # for idx in range(0, len(state.moments)):
        #     for moment in self.possible_moments_of_properties:
        #         params = dict([("state", state), ("abstract_by",
        #                       moment), ("abstract_idx", idx)])
        #         applicable_abstract_rules.append(("abstraction", params))

        # 2. block append - append by mux prop & act prop
        applicable_add_block_rules = list()
        for idx in [-1, len(state.moments) - 1]:
            # for idx in range(-1, len(state.moments)):
            # for moment in self.possible_moments_of_mux_properties:
            #     params = dict([("state", state),
            #                    ("moment_to_append", moment),
            #                    ("append_position", (idx, idx + 1))])
            #     applicable_add_block_rules.append(("append", params))
            for moment in self.possible_moments_of_act_properties:
                params = dict([("state", state),
                               ("moment_to_append", moment),
                               ("append_position", (idx, idx + 1))])
                applicable_add_block_rules.append(("append", params))
        applicable_multiplexor_decomposition_rule = list()

        # # 3 -1. decompose in multiplexor
        # for moment_idx, moment in enumerate(state.moments):
        #     for block_idx, block in enumerate(moment.blocks):
        #         if isinstance(block, MyProperty) and block.is_multiplexor():
        #             for num_decompose in range(1, 3 + 1):
        #                 params = dict([("state", state), ("moment_idx", moment_idx),
        #                               ("block_idx", block_idx), ("num_decompose", num_decompose)])
        #                 applicable_multiplexor_decomposition_rule.append(
        #                     ("mux_decompose", params)
        #                 )
        # if given state includs some block with property of moment with property ....

        # 3-2 decompose in property
        applicable_prop_decomposition_rule = list()
        for moment_idx, moment in enumerate(state.moments):
            for block_idx, block in enumerate(moment.blocks):
                if isinstance(block, MyProperty) and block.action_property:
                    valid_rule_ids = self.valid_rule_id_for_block(block=block)
                    if not valid_rule_ids:
                        params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", None), ("block_idx", block_idx)])
                        applicable_prop_decomposition_rule.append(("prop_decompose", params))
                    else:
                        for rule_id in valid_rule_ids:
                            params = dict([("state", state), ("moment_idx", moment_idx), ("rule_id", rule_id), ("block_idx", block_idx)])
                            applicable_prop_decomposition_rule.append(("prop_decompose", params))

        return applicable_prop_decomposition_rule + applicable_abstract_rules + applicable_add_block_rules + applicable_multiplexor_decomposition_rule

    ## Rule Runners ##

    def replace_into_prop_block_rule(self, state: MomentBasedState, replace_to : MyProperty, moment_idx: int, block_idx=None) -> MomentBasedState: # res is not list 
        to_return = list()
        if ((not isinstance(block_idx, int)) and len(state.moments[moment_idx].blocks) != 1):
            raise InvalidDecomposeException(
                "if block idx is given as none, we automatically assume that # of block on the moment is 1"
                + f"The current state to decompose was  \n {str(state)} \n and moment_idx was {moment_idx}, block idx was {block_idx}")
        if block_idx == None:
            block_to_deal = state.moments[moment_idx].blocks[0]
            block_idx = 0
        else:
            block_to_deal = state.moments[moment_idx].blocks[block_idx]
        moments_for_new_state = list()
        for enum_idx_mm, enum_moment in enumerate(state.moments): 
            if moment_idx != enum_idx_mm :
                moments_for_new_state.append(enum_moment)
            else :
                moment_possibly_for_lo_blocs_builder = list()
                new_moment_builder = list()
                new_moment_builder.append(replace_to)
                lo_blocks = [  enum_block   for enum_idx_block, enum_block in enumerate(enum_moment.blocks) if enum_idx_block != block_idx]
                moment_possibly_for_lo_blocs_builder = lo_blocks
                moment_possibly_for_lo_blocks = MyMoment(blocks=moment_possibly_for_lo_blocs_builder)
                if moment_possibly_for_lo_blocks.qubits_on().intersection(replace_to.qubits_on()):
                    if len(moment_possibly_for_lo_blocs_builder) > 0 :
                        new_moment = MyMoment(blocks=new_moment_builder)
                        moments_for_new_state.append(new_moment)
                        moments_for_new_state.append(moment_possibly_for_lo_blocks)
                    else :
                        new_moment = MyMoment(blocks=new_moment_builder)
                        moments_for_new_state.append(new_moment)

                else : 
                    new_moment_builder+= lo_blocks       
                    new_moment = MyMoment(blocks=new_moment_builder)
                    moments_for_new_state.append(new_moment)

        new_state = MomentBasedState(moments=moments_for_new_state)
        return new_state

    def add_prop_block_rule(self, state: MomentBasedState, prop_block_to_append=MyProperty, append_position=Tuple[int, int]):
        # ====================
        # 2. Append Rules
        # ====================
        '''
        append_position (i,j) : decides where to append the given prop block with in the circuit and i+1==j must hold
        '''
        raise NotImplementedError("add_prop_block_rule")
        new_state = state.copy()
        l_moment, r_moment = append_position
        assert l_moment < r_moment
        if state.moments[l_moment].is_block_appendable(prop_block_to_append):
            new_state.moments[l_moment].append(prop_block_to_append)
        elif state.moments[r_moment].is_block_appendable(prop_block_to_append):
            new_state.moments[r_moment].append(prop_block_to_append)
        else:
            new_moment = MyMoment(blocks=[prop_block_to_append])
            new_state.moments.insert(r_moment, new_moment)
        return new_state

    def add_props_in_moment_rule(self, state: MomentBasedState, moment_to_append=MyMoment, append_position=Tuple[int, int]):
        '''
        append_position (i,j) : decides where to append the given prop block with in the circuit and i+1==j must hold
        '''
        new_state = state.copy()
        l_moment, r_moment = append_position
        assert l_moment < r_moment and r_moment <= len(state.moments)
        new_state.moments.insert(r_moment, moment_to_append)
        return new_state

    # Running Multiplexor decompose rules
    def muitplexor_decompose(self, state: MomentBasedState, moment_idx: int, block_idx=None, num_decompose=2) -> List[MomentBasedState]:
        '''
        (moment idx, block idx) decideds the real position of block withinthe state(=abstracted circuit) to give decomposition on
        '''
        to_return = list()  # will be list of block
        if ((not isinstance(block_idx, int)) and len(state.moments[moment_idx].blocks) != 1):
            raise InvalidDecomposeException(
                "if block idx is given as none, we automatically assume that # of block on the moment is 1"
                + f"The current state to decompose was  \n {str(state)} \n and moment_idx was {moment_idx}, block idx was {block_idx}")
        if block_idx == None:
            block_to_deal = state.moments[moment_idx].blocks[0]
            block_idx = 0
        else:
            block_to_deal = state.moments[moment_idx].blocks[block_idx]

        # gather combinations of ctrl, targ qubits.
        control_combis = [x for i in range(1, num_decompose + 1) for x in combinations(block_to_deal.mux_property[CONTROL], r=i)]
        targ_combis = [x for i in range(1, num_decompose + 1) for x in combinations(block_to_deal.mux_property[TARGET], r=i)]

        possible_control_decomposes = [x for x in product(control_combis, repeat=num_decompose)]

        possible_target_decomposes = [x for x in product(targ_combis, repeat=num_decompose)]
        res_decomposes = list()
        for con, targ in product(possible_control_decomposes, possible_target_decomposes):
            list_of_decomposes = list()
            for i in range(num_decompose):
                mux_prop_builder = dict()
                mux_prop_builder[CONTROL] = con[i]
                mux_prop_builder[TARGET] = targ[i]
                curr_res_block = MyProperty(mux_property=mux_prop_builder)
                curr_res_moment = MyMoment(blocks=[curr_res_block])
                list_of_decomposes.append(curr_res_moment)
            res_decomposes.append(list_of_decomposes)

        for decompose_res in res_decomposes:
            if len(state.moments[moment_idx].blocks) == 1:
                newly_generated_li_of_moment = list()
                for idx, moment in enumerate(state.moments):
                    assert isinstance(moment, MyMoment)
                    if idx == moment_idx:
                        newly_generated_li_of_moment += decompose_res
                    else:
                        newly_generated_li_of_moment.append(moment)
                new_state = MomentBasedState(
                    moments=newly_generated_li_of_moment)
                to_return.append(new_state)
            else:
                # insertion routine
                newly_generated_li_of_moment = list()
                for idx, moment in enumerate(state.moments):
                    assert isinstance(moment, MyMoment)
                    if idx == moment_idx:
                        preserved_blocks = [block for l, block in enumerate(moment.blocks) if l != block_idx]
                        intial_moment_of_newly_spanend = MyMoment(blocks=(preserved_blocks + decompose_res[0].blocks))
                        newly_generated_li_of_moment.append(intial_moment_of_newly_spanend)
                        for i in range(1, len(decompose_res)):
                            newly_generated_li_of_moment.append(decompose_res[i])
                    else:
                        newly_generated_li_of_moment.append(moment)
                new_state = MomentBasedState(moments=newly_generated_li_of_moment)
                to_return.append(new_state)
        if len(to_return) == 0:
            raise Exception("The returned list of states after multiplexor decompose rule should not be zero" +
                            f"\n The env was as follows. \n state: {state}, \n moment_idx : {moment_idx}, \n block_idx :{block_idx}, \n num_decompose : {num_decompose}")
        return to_return

    # Running Property decompose rules
    def property_decompose(self, state: MomentBasedState, moment_idx: int, rule_id: str = None, block_idx=None) -> List[MomentBasedState]:
        to_return = list()
        if ((not isinstance(block_idx, int)) and len(state.moments[moment_idx].blocks) != 1):
            raise InvalidDecomposeException(
                "if block idx is given as none, we automatically assume that # of block on the moment is 1"
                + f"The current state to decompose was  \n {str(state)} \n and moment_idx was {moment_idx}, block idx was {block_idx}")
        if block_idx == None:
            block_to_deal = state.moments[moment_idx].blocks[0]
            block_idx = 0
        else:
            block_to_deal = state.moments[moment_idx].blocks[block_idx]

        if len(block_to_deal.action_property.keys()) == 1:  # if we decompose only for single length block
            try:
                to_gen_in_moment, res_seq_of_blocks = self.prop_decompose_of_single_qubit_block(block=block_to_deal, rule_id=rule_id)
            except TypeError as e:
                print(f"Type error occuredThe block_to_deal was {block_to_deal}")
                exit()
            if not res_seq_of_blocks:
                return None
            if to_gen_in_moment:
                for blocks in res_seq_of_blocks:
                    decmp_into_li_moments = blocks_to_opt_moments(blocks=blocks)  # a list of moments
                    moments_for_new_mbs = list()
                    # copy.deepcopy(state.moments)
                    for idx, moment in enumerate(state.moments):
                        if idx == moment_idx:
                            moments_for_new_mbs += decmp_into_li_moments
                            # may do deepcopy
                            leftover_blocks = [block for blk_idx_in_mm, block in enumerate(moment.blocks) if blk_idx_in_mm != block_idx] # blk = block
                            if not len(leftover_blocks) == 0:
                                mm_for_leftover_blocks = MyMoment(blocks=leftover_blocks)
                                moments_for_new_mbs.append(mm_for_leftover_blocks)
                        else:
                            moments_for_new_mbs.append(moment)
                    new_mbs = MomentBasedState(moments=moments_for_new_mbs)
                    to_return.append(new_mbs)
            else:  # BOOL, PHASE, NC Decompose?
                decmp_into_li_moments = blocks_to_opt_moments(blocks=res_seq_of_blocks)
                moments_for_new_mbs = list()
                for idx, moment in enumerate(state.moments):
                    if idx == moment_idx:
                        moments_for_new_mbs += decmp_into_li_moments
                        leftover_blocks = [block for blk_idx_in_mm, block in enumerate(moment.blocks) if blk_idx_in_mm != block_idx]
                        if len(leftover_blocks) != 0 :
                            mm_for_leftover_blocks = MyMoment(blocks=leftover_blocks)
                            moments_for_new_mbs.append(mm_for_leftover_blocks)
                    else:
                        moments_for_new_mbs.append(moment)
                new_mbs = MomentBasedState(moments=moments_for_new_mbs)
                to_return.append(new_mbs)
            return to_return

        elif len(block_to_deal.action_property.keys()) > 1:
            list_of_spanned_moments = self.prop_decompose_by_block(block=block_to_deal) #return list of list of moments
        
            if not list_of_spanned_moments:
                return None
            for spanned_moments in list_of_spanned_moments: 
                newly_gen_li_of_moments = list()
                for idx, moment in enumerate(state.moments):
                    if idx == moment_idx:
                        newly_gen_li_of_moments += spanned_moments
                        leftover_blocks = [block for blk_idx_in_mm, block in enumerate(moment.blocks) if blk_idx_in_mm != block_idx]
                        if len(leftover_blocks) != 0 :
                            mm_for_leftover_blocks = MyMoment(blocks=leftover_blocks)
                            newly_gen_li_of_moments.append(mm_for_leftover_blocks)
                    else:
                        newly_gen_li_of_moments.append(state.moments[idx])

                res_state = MomentBasedState(moments=newly_gen_li_of_moments)
                to_return.append(res_state)
            return to_return
        elif len(block_to_deal.action_property.keys()) in [3,4]:
            return None
        else:
            print(f"Something wrong when decomposition. The block was {block_to_deal}")
            exit()

    def prop_decompose_of_single_qubit_block(self, block: MyProperty, rule_id: str) -> Union[List[MyProperty], List[List[MyProperty]]]:
        prop_to_decompose = list(block.action_property.values())[0]
        if prop_to_decompose == BOOL:
            return False, prop_decompose_single_bool(block=block, rule_id=rule_id)  # HardCoded
        elif prop_to_decompose == PHASING:
            return False, prop_decompose_single_phase(block=block, rule_id=rule_id)  # TODO
        elif prop_to_decompose == NC:
            return False, property_decompose_single_nc(block=block, rule_id=rule_id)  # TODO
        elif prop_to_decompose == ENTANGLE:
            return True, prop_decompose_single_entangle(block=block,
                                                        inseparable_gate_qubit_num = self.inseparable_gate_qubit_num,
                                                        min_sp_gate_qubit_num=self.min_sp_gate_qubit_num)
        elif prop_to_decompose == None :
            return False, list()
        else:
            raise Exception(f"Property to decompose \"{prop_to_decompose}\" is some what nonvalidate one. ")
    
    def prop_decompose_by_block_aux(self, block : MyProperty) -> Union[List[MyProperty], List[List[MyProperty]]]:
        # decompose a block by atomic rules of single-qubit block decompositions
        # input `block` is a temporary (partial) block that is choosen to be decompose from the main block
        to_return = list()
        to_product = list()
        for key, val in block.action_property.items() :
            temp_block_builder = {key : val}
            temp_block = MyProperty(action_property = temp_block_builder)
            valid_rules_for_temp_single_block = self.valid_rule_id_for_block(block = temp_block)
            if valid_rules_for_temp_single_block == None :
                valid_rules_for_temp_single_block = [None]
            curr_possible_single_block_decomposes = list()
            for rule_id in valid_rules_for_temp_single_block :
                to_gen_in_moment, res = self.prop_decompose_of_single_qubit_block(block=temp_block, rule_id = rule_id)
                if res == None or len(res)==0 :
                    # print("why res here is none")
                    # print(temp_block)
                    return list()
                if isinstance(res[0], list):
                    for decomposed_blocks  in res :
                        print("asdfdsaf") ###????
                elif isinstance(res[0], MyProperty):
                    res = [res]
                for decomposed_blocks  in res :
                    # decomposed blocks are either len two or three 
                    if len(decomposed_blocks) == 2:
                        left_non_append = decomposed_blocks.copy()
                        right_non_append = decomposed_blocks.copy()
                        center_non_append = decomposed_blocks.copy()
                        left_non_append = [None]  + left_non_append 
                        right_non_append = right_non_append + [None]
                        center_non_append = [ decomposed_blocks[0], None, decomposed_blocks[1]]
                        curr_possible_single_block_decomposes.append(left_non_append)
                        curr_possible_single_block_decomposes.append(right_non_append)
                        curr_possible_single_block_decomposes.append(center_non_append)
                    elif len(decomposed_blocks) == 3 : 
                        curr_possible_single_block_decomposes.append(decomposed_blocks)
            to_product.append(curr_possible_single_block_decomposes)

        for x in itertools.product(*to_product):
            for_first_moment  = list()
            for_second_moment = list()
            for_third_moment  = list()
            for y in x : 
                for_first_moment.append(y[0])
                for_second_moment.append(y[1])
                for_third_moment.append(y[2])
            block_for_first_moment_builder = dict()
            for x in for_first_moment :
                if x != None :
                    x_keys = list(x.action_property.keys())
                    block_for_first_moment_builder[ x_keys[0]] = x.action_property[ x_keys[0]]
            block_for_first_moment = MyProperty(action_property=block_for_first_moment_builder)

            block_for_second_moment_builder = dict()
            for x in for_second_moment :
                if x != None :
                    x_keys = list(x.action_property.keys())
                    block_for_second_moment_builder[ x_keys[0]] = x.action_property[ x_keys[0]]

            block_for_second_moment = MyProperty(action_property=block_for_second_moment_builder)

            block_for_third_moment_builder = dict()
            for x in for_third_moment :
                if x != None :
                    x_keys = list(x.action_property.keys())
                    block_for_third_moment_builder[ x_keys[0]] = x.action_property[ x_keys[0]]
            block_for_third_moment = MyProperty(action_property=block_for_third_moment_builder)

            to_return.append([block_for_first_moment,
                              block_for_second_moment,
                              block_for_third_moment ])
        return to_return

    def prop_decompose_by_block(self, block: MyProperty) -> List[List[Moment]]:
        to_return = list()
        if block.action_property == None:
            return None
        elif len(block.action_property.keys()) > 1 and block.no_entanglement():
            keys = list(block.action_property.keys())
            for will_decompose in get_all_choices_of_qreg(keys) :
                hold_qubits = [ k  for k in keys if not set(k).issubset(set(will_decompose))]
                hold_qubits_li_of_int = [ hq[0] for hq in hold_qubits ]
                hold_qubit_three_partitions = into_n_partitions(hold_qubits_li_of_int, n=3)
                temp_block_to_decompose_builder =  dict([((x,) , block.action_property[(x,)]) for x in will_decompose])
                temp_block_to_decompose = MyProperty(action_property=temp_block_to_decompose_builder)
                duplication_checker = set()
                for x in self.prop_decompose_by_block_aux(temp_block_to_decompose): 
                    for hold_qubit_partition in hold_qubit_three_partitions:
                        assert len(hold_qubit_partition)==3
                        first, second, third = x[0], x[1], x[2] # each of them is block
                        first_hold_q, scond_hold_q, third_hold_q = hold_qubit_partition
                        if len(hold_qubits) != 0 :
                            first_with_hold_qubits_builder  = dict()
                            second_with_hold_qubits_builder  = dict()
                            third_with_hold_qubits_builder  = dict()

                            for hq in first_hold_q : first_with_hold_qubits_builder[(hq,)]  = block.action_property[(hq,)]
                            for hq in scond_hold_q : second_with_hold_qubits_builder[(hq,)] = block.action_property[(hq,)] 
                            for hq in third_hold_q : third_with_hold_qubits_builder[(hq,)]  = block.action_property[(hq,)]

                            for key in first.action_property.keys()  : first_with_hold_qubits_builder[key]  = first.action_property[key]
                            for key in second.action_property.keys() : second_with_hold_qubits_builder[key] = second.action_property[key]
                            for key in third.action_property.keys()  : third_with_hold_qubits_builder[key]  = third.action_property[key]
                            first_with_hold_qubits  = MyProperty(action_property= first_with_hold_qubits_builder )
                            second_with_hold_qubits = MyProperty(action_property= second_with_hold_qubits_builder)
                            third_with_hold_qubits  = MyProperty(action_property= third_with_hold_qubits_builder )

                            gathered_moments = list()
                            for y in  [first_with_hold_qubits, second_with_hold_qubits, third_with_hold_qubits]:
                                if not y.is_empty_block():
                                    gathered_moments.append(y)
                            opt_gathered_moments = blocks_to_opt_moments(gathered_moments)
                            if str(opt_gathered_moments) not in duplication_checker : 
                                to_return.append(opt_gathered_moments)
                                duplication_checker.add(str(opt_gathered_moments))

                        else :
                            gathered_moments = list()
                            for y in [first, second, third]:
                                if not y.is_empty_block():
                                    gathered_moments.append(y)
                            opt_gathered_moments = blocks_to_opt_moments(gathered_moments)
                            if str(opt_gathered_moments) not in duplication_checker : 
                                to_return.append(opt_gathered_moments)
                                duplication_checker.add(str(opt_gathered_moments))
        return to_return



# Single Decomposition
def prop_decompose_single_bool(block: MyProperty, rule_id: str) -> List[MyProperty]: 
    # sq for single qubit
    qubit_position = list(block.action_property.keys())[0][0]
    assert len(block.action_property.keys()) == 1 and block.action_property[(qubit_position,)] == BOOL
    if rule_id == SPAN:
        one_c = block.copy()
        two_c = block.copy()
        return [one_c, two_c]
    if rule_id == TRIPLE_SPAN :
        one_c = block.copy()
        two_c = block.copy()
        three_c = block.copy()
        return [one_c, two_c,three_c]
    elif rule_id == OPEN_AND_CLOSE_B:
        one_nc_builder = dict([((qubit_position,), NC)])
        two_nc_builder = dict([((qubit_position,), NC)])
        one_nc = MyProperty(action_property=one_nc_builder)
        two_nc = MyProperty(action_property=two_nc_builder)
        center_c_builder = dict([((qubit_position,), BOOL)])
        center_c = MyProperty(action_property=center_c_builder)
        return [one_nc, center_c, two_nc]
    elif rule_id == OPEN_AND_CLOSE_P:
        one_nc_builder = dict([((qubit_position,), NC)])
        two_nc_builder = dict([((qubit_position,), NC)])
        one_nc = MyProperty(action_property=one_nc_builder)
        two_nc = MyProperty(action_property=two_nc_builder)
        center_c_builder = dict([((qubit_position,), PHASING)])
        center_c = MyProperty(action_property=center_c_builder)
        return [one_nc, center_c, two_nc]
    elif rule_id == OPEN_AND_CLOSE_NC :
        one_nc_builder = dict([((qubit_position,), NC)])
        two_nc_builder = dict([((qubit_position,), NC)])
        one_nc = MyProperty(action_property=one_nc_builder)
        two_nc = MyProperty(action_property=two_nc_builder)
        return [one_nc, two_nc]
    else:
        raise Exception(f"Rule of named [{rule_id}]  for classical decomposition is not valid")


def prop_decompose_single_entangle(block: MyProperty, inseparable_gate_qubit_num, min_sp_gate_qubit_num) -> List[List[MyProperty]]:
    to_return = list()
    qubits_to_be_entangled = set(list(block.action_property.keys())[0])
    inseparable_gate_qubit_num = [ i for i in inseparable_gate_qubit_num if i <= len(qubits_to_be_entangled)]
    # inseparable_gate_qubit_num = list((inseparable_gate_qubit_num))
    for i in range(1, len(qubits_to_be_entangled) + 1):
        for qubits_to_be_superposed in itertools.combinations(qubits_to_be_entangled, r=i):
            # superpose_block_builder = {(qbit_superposed,): NC}  # to be used as one of result
            to_superpose_blocks = list()
            # if min_sp_gate_qubit_num == 1:
            for q in qubits_to_be_superposed :
                superpose_block_builder = dict([ ((q,) , NC)  ])
                superpose_block = MyProperty(action_property=superpose_block_builder)
                to_superpose_blocks.append(superpose_block)

            for conn_li_pair in connected_lists_of_pair(qubits_to_be_entangled, qubits_to_be_superposed, inseparable_gate_qubit_num):
                # this charaterize one single state to be produced
                to_append = to_superpose_blocks.copy()
                for pair in conn_li_pair :
                    entnagling_block_builder = dict( [((q,) , None)  for q in pair])
                    entnagling_block = MyProperty(action_property=entnagling_block_builder)
                    to_append.append(entnagling_block)
                to_return.append(to_append)

            if  min_sp_gate_qubit_num >=2 : 
                to_superpose_blocks = list()
                if len(qubits_to_be_superposed)  ==  min_sp_gate_qubit_num :
                    superpose_block_builder = dict()
                    for q in qubits_to_be_superposed:
                        superpose_block_builder[(q, )] = NC
                    superpose_block = MyProperty(action_property=superpose_block_builder)
                    to_superpose_blocks.append(superpose_block)
                    for conn_li_pair in connected_lists_of_pair(qubits_to_be_entangled, qubits_to_be_superposed, inseparable_gate_qubit_num):
                        to_append = to_superpose_blocks.copy()
                        for pair in conn_li_pair :
                            entnagling_block_builder = dict( [((q,) , None)  for q in pair])
                            entnagling_block = MyProperty(action_property=entnagling_block_builder)
                            to_append.append(entnagling_block)
                        to_return.append(to_append)
    return to_return

def property_decompose_single_nc(block: MyProperty, rule_id : str):
    # Decompose of Single NC
    qubit_position = list(block.action_property.keys())[0][0]
    assert len(block.action_property.keys()) == 1 and block.action_property[(qubit_position,)] == NC
    if rule_id == ADD :
        one_nc = block.copy()
        to_add_builder = dict([((qubit_position,), None)])
        to_add = MyProperty(action_property=to_add_builder)
        return [to_add, one_nc]
    else :
        raise Exception(f"Rule of named [{rule_id}] for nc decomposition is not valid")

def prop_decompose_single_phase(block: MyProperty,rule_id: str):
    qubit_position = list(block.action_property.keys())[0][0]
    assert len(block.action_property.keys()) == 1 and block.action_property[(qubit_position,)] == PHASING
    if rule_id == SPAN:
        one_c = block.copy()
        two_c = block.copy()
        return [one_c, two_c]
    elif rule_id == OPEN_AND_CLOSE_B:
        one_nc_builder = dict([((qubit_position,), NC)])
        two_nc_builder = dict([((qubit_position,), NC)])
        one_nc = MyProperty(action_property=one_nc_builder)
        two_nc = MyProperty(action_property=two_nc_builder)
        center_c_builder = dict([((qubit_position,), BOOL)])
        center_c = MyProperty(action_property=center_c_builder)
        return [one_nc, center_c, two_nc]
    else:
        raise Exception(f"Rule of named [{rule_id}]  for phasing decomposition is not valid")

