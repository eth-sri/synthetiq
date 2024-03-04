from types import new_class
import cirq
import numpy as np
if __name__ == "qsyn.synthesis":
    from qsyn.synthesis_spec.specification import Spec
    from qsyn.util.utils import component_priors, check_continuity
    from qsyn.util.utils import normalize_gate_operations
else:
    from synthesis_spec.specification import Spec
    from util.utils import component_priors,check_continuity
    from util.utils import normalize_gate_operations

from typing import List, Tuple, Union
from itertools import permutations, product, groupby
from cirq import Simulator
import pprint
from timeit import default_timer as timer, repeat


TOLERANCE = 1e-04

class NotFoundCircuitException(Exception):
    def __init__(self, msg: str) -> None:
        self.msg = msg

    def __str__(self) -> str:
        return "The target circuit does not exist within the search space"


class Synthesis():
    def __init__(self, spec: Spec):
        self._spec = spec
        self._qreg = cirq.LineQubit.range(spec.qreg_size)
        self.gate_operations: List[cirq.GateOperation] = self._generate_spanned_gate_set()

        
        print("Do Gate Operations Normalization")
        temp = normalize_gate_operations(self.gate_operations, self._qreg)
        print(f"Reducde to {len(self.gate_operations)}->{len(temp)}")
        self.gate_operations = temp

        # print(self.gate_operations)
        # print(len(self.gate_operations))
        # normalize_gate_operations(self.gate_operations, self._qreg)
        # input()


        self.gate_operations.sort(key = lambda s : str(s), reverse=True)
        self.simulator = Simulator()

    def _generate_spanned_gate_set(self) -> List[cirq.GateOperation]:
        spanned_get_set = []
        if not self.spec.component_gates:
            print("The Component Set is Empty")
            return []
        self.spec.component_gates.sort(key = lambda s : str(s))
        for component_gate in self.spec.component_gates:
            for qubit_to_apply in list(permutations(self.working_qubit, component_gate.num_qubits())):
                spanned_get_set.append(component_gate(*qubit_to_apply))

        return list(set(spanned_get_set)) # for some operations (eg, iSWAP ) gate normalization could occured

    @property
    def spec(self):
        return self._spec

    @property
    def spanned_set(self) -> List[cirq.GateOperation]:
        return self.gate_operations

    @property
    def working_qubit(self):
        return self._qreg

    ####
    # Synthesis Validation Check Functions
    ####
    def test_circuit_equivalence(self, circuit_to_check: cirq.Circuit) -> bool:
        try:
            cirq.testing.assert_allclose_up_to_global_phase(
                circuit_to_check.unitary(),
                self.spec.spec_object.unitary(), 
                atol=TOLERANCE
            )
            return True
        except AssertionError:
            return False

    def test_iopair_phase_equiv_act_cache(self, circuit_to_check: cirq.Circuit) -> bool:
        built_dict = dict()
        flag = True
        for io_pair in self.spec.spec_object.get_io_pairs():
            in_nparr, out_nparr = io_pair
            simulator = Simulator()
            res = simulator.simulate(circuit_to_check, qubit_order=self.working_qubit, initial_state=in_nparr)
            built_dict[tuple(in_nparr)]  = res.final_state_vector
            if not cirq.linalg.allclose_up_to_global_phase(out_nparr, res.final_state_vector, atol=TOLERANCE):
                flag = False
        return flag, built_dict

    def test_iopair_strict_equiv_act_cache(self, circuit_to_check):
        built_dict = dict()
        flag = True
        for io_pair in self.spec.spec_object.get_io_pairs():
            in_nparr, out_nparr = io_pair
            simulator = Simulator()
            res = simulator.simulate(circuit_to_check, qubit_order=self.working_qubit, initial_state=in_nparr)
            built_dict[tuple(in_nparr)]  = res.final_state_vector
            if not np.allclose(res.final_state_vector, out_nparr, atol=TOLERANCE):
                flag = False
        return flag, built_dict

    def test_iopair_phase_equiv(self, circuit_to_check: cirq.Circuit) -> bool:
        for io_pair in self.spec.spec_object.get_io_pairs():
            in_nparr, out_nparr = io_pair
            simulator = Simulator()
            res = simulator.simulate(circuit_to_check, qubit_order=self.working_qubit, initial_state=in_nparr)
            if not cirq.linalg.allclose_up_to_global_phase(out_nparr, res.final_state_vector, atol=TOLERANCE):
                return False
        return True

    def test_iopair_strict_equiv(self, circuit_to_check):
        for io_pair in self.spec.spec_object.get_io_pairs():
            in_nparr, out_nparr = io_pair
            simulator = Simulator()
            res = simulator.simulate(circuit_to_check, qubit_order=self.working_qubit, initial_state=in_nparr)
            if not np.allclose(res.final_state_vector, out_nparr, atol=TOLERANCE):
                return False
        return True

    
    def test_unitary_phase_equiv(self, circuit_to_check):
        # try:
        #     cirq.testing.assert_allclose_up_to_global_phase(
        #         circuit_to_check.unitary(),
        #         self.spec.spec_object,
        #         atol=TOLERANCE
        #     )
        #     return True
        # except AssertionError:
        #     return False
        raise NotImplementedError
    def test_unitary_strict_equiv(self, circuit_to_check):
        raise NotImplementedError

    def test_synthesis(self, circuit_to_check: cirq.Circuit, phase_equiv: bool) -> bool:
        if self.spec.type_of_spec == self.spec.PART_SYN_SPEC:
            if phase_equiv:
                return self.test_iopair_phase_equiv(circuit_to_check)
            else:
                return self.test_iopair_strict_equiv(circuit_to_check)
        else:
            raise NotImplementedError()

    def test_synthesis_att_guided(self, circuit_to_check: cirq.Circuit, phase_equiv: bool) -> bool:
        if self.spec.type_of_spec == self.spec.PART_SYN_SPEC:
            if phase_equiv:
                return self.test_iopair_phase_equiv_act_cache(circuit_to_check)
            else:
                return self.test_iopair_strict_equiv_act_cache(circuit_to_check)
        # else:
        #     if phase_equiv:
        #         return self.test_unitary_phase_equiv(circuit_to_check)
        #     else:
        #         return self.test_unitary_strict_equiv(circuit_to_check)



class BFSSynthesis(Synthesis):

    def __init__(self,
                 spec: Spec,):
        super().__init__(spec=spec)
        self.component_prior = component_priors(components=self._spec.component_gates, is_att_guide=False)
    def get_search_space_size(self) -> int:
        sum = 0
        len_mu = len(self.spanned_set)
        for num_instruct in range(1, self.spec.get_max_inst_number()+1):
            sum += len_mu ** num_instruct
        return sum

    ####
    # Synthesis aux functions
    ####

    def is_involution(self, circuit_to_check: cirq.Circuit, gate_op_to_check: cirq.ops.GateOperation):
        if len(circuit_to_check) == 0:
            return False
        gate_op = None
        gate_op_to_check_qubits = sorted( list([q.x for q in gate_op_to_check.qubits]))
        rep_qubit_index = gate_op_to_check_qubits[0]
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
        # if not (gate_op == gate_op_to_check):
        #     return False
        # elif (gate_op == gate_op_to_check
        #       and np.allclose(np.matmul(cirq.unitary(gate_op.gate), cirq.unitary(gate_op_to_check.gate)),
        #                       np.identity(cirq.unitary(gate_op.gate).shape[0]))):
        elif (gate_op == gate_op_to_check
              and gate_op.gate in self.component_prior["INVOLUTIONARY"]):
            return True
        elif (gate_op.qubits == gate_op_to_check.qubits):
            if gate_op.gate in self.component_prior["INVERSES"].keys():
                if (gate_op_to_check.gate == self.component_prior["INVERSES"][gate_op.gate]):
                    # print("PING!")
                    # print(circuit_to_check)
                    # print(gate_op)
                    # print(gate_op_to_check)
                    # input()
                    return True
            if gate_op_to_check.gate in self.component_prior["IDENTITY_N"].keys() and str(gate_op_to_check) == str(gate_op):
                
                # print(circuit_to_check)
                # print(circuit_to_check[0])
                # print(gate_op)
                # print(gate_op_to_check)
                # print(self.component_prior["IDENTITY_N"])
                iden_N = self.component_prior["IDENTITY_N"][gate_op_to_check.gate]
                collected_moment_idx = list()
                for moment_loc in reversed(range(len(circuit_to_check))):
                    temp_gate_app = circuit_to_check.operation_at(gate_op_to_check.qubits[0], moment_loc)
                    if temp_gate_app == gate_op_to_check : collected_moment_idx  = [moment_loc] + collected_moment_idx
                if len(collected_moment_idx) == iden_N-1 and check_continuity(collected_moment_idx) and moment_idx_of_gate_op in collected_moment_idx:
                    return True
        else:
            return False

    def is_identitical_in_matrix(self, circuit_to_check, gate_op_to_check, workset_to_check, workset_to_check_str):
        constructed_ahead = circuit_to_check.copy()
        constructed_ahead.append(gate_op_to_check)
        if constructed_ahead is None:
            raise ValueError(f"constructed_ahead is {constructed_ahead}")

        for worked_circuit in workset_to_check:
            if worked_circuit.all_qubits() == constructed_ahead.all_qubits() and np.allclose(cirq.unitary(constructed_ahead), cirq.unitary(worked_circuit), atol=1e-04):
                return True
        return False

    def is_identitical(self, circuit_to_check, gate_op_to_check, workset_to_check, workset_to_check_str):
        constructed_ahead = circuit_to_check.copy()
        constructed_ahead.append(gate_op_to_check)
        if constructed_ahead is None:
            raise ValueError(f"constructed_ahead is {constructed_ahead}")
        # str_to_check = re.sub(r"[\n\t\s]*", "", repr(constructed_ahead))
        str_to_check = repr(constructed_ahead)
        # res = (str_to_check in workset_to_check_str)
        res = constructed_ahead in workset_to_check
        return res

    ####
    # Enumerative Synthesis Algorithms
    ####
    def prune_bfs_enumeration_synthesis(self, timeout: int = 3600, involution_prune: bool = False, identity_prune: bool = False) -> Tuple[Union[str, int], cirq.Circuit]:
        start_time = timer()
        # input()
        # logging.info(f"[Pruning Mode] involution_prune : {involution_prune}, identity_prune : {identity_prune}")
        print(f"[Pruning Mode] involution_prune : {involution_prune}, identity_prune : {identity_prune}")

        cnt = 0
        identity_prune_cnt = 0
        involution_prune_cnt = 0

        do_prune = involution_prune or identity_prune
        workset_circuits = [(cirq.Circuit())]

        for inst_num in range(0, self.spec.get_max_inst_number()):
            # logging.info(
                # f"Now in Depth {inst_num}. Which means, I am generating of which depth is {inst_num+1} ")
            # logging.info(
                # f"The Length of working set in current level is : {len(workset_circuits)}")
            print(
                f"Now in Depth {inst_num}. Which means, I am generating of which depth is {inst_num+1} ")
            print(
                f"The Length of working set in current level is : {len(workset_circuits)}")

            new_workset_circuits = []
            new_workset_circuits_str = []
            for work_circuit in workset_circuits:
                for gate_op in self.gate_operations:
                    # some debung infos
                    cnt += 1
                    # check timeout
                    if cnt % 1000 == 0:
                        current_time = timer()
                        if current_time - start_time > timeout:
                            # logging.info("TIME OUT!")
                            return "TIME OUT", None
                    # if cnt % 10000 == 0:
                        # logging.info(f"Progressing..{cnt}")
                    
                    synthesis_cadndiate = work_circuit.copy()
                    if do_prune:
                        if involution_prune and self.is_involution(synthesis_cadndiate, gate_op):
                            involution_prune_cnt += 1
                            continue  # Skip Doing Test
                        elif identity_prune and self.is_identitical(synthesis_cadndiate, gate_op, new_workset_circuits, new_workset_circuits_str):
                            # elif identity_prune and self.is_identitical_in_matrix(synthesis_cadndiate, gate_op, new_workset_circuits,new_workset_circuits_str) :
                            identity_prune_cnt += 1
                            continue  # Skip Doing Test
                        else:
                            synthesis_cadndiate.append(gate_op)
                            new_workset_circuits.append(synthesis_cadndiate)
                            str_to_append = repr(synthesis_cadndiate)
                            # str_to_append = re.sub(r"[\n\t\s]*", "", repr(synthesis_cadndiate))
                            new_workset_circuits_str.append(str_to_append)
                    else:
                        synthesis_cadndiate.append(gate_op)
                        new_workset_circuits.append(synthesis_cadndiate)
                        # str_to_append = re.sub(r"[\n\t\s]*", "", repr(synthesis_cadndiate))
                        str_to_append = repr(synthesis_cadndiate)
                        new_workset_circuits_str.append(str_to_append)
                    if self.test_synthesis(synthesis_cadndiate, self.spec.equiv_phase):
                        elapsed_time = timer() - start_time
                        # logging.info("Synthesis Done.")
                        # logging.info("Synthesized Circuit is" +
                                    #  "\n" + f"{str(synthesis_cadndiate)}")
                        # logging.info(f"Time : {elapsed_time}")
                        # logging.info(f"Searched Space : {cnt+1}")
                        print(
                            f"Involution Pruned Count : {involution_prune_cnt}")
                        print(
                            f"Identity Pruned Count : {identity_prune_cnt}")

                        return elapsed_time, synthesis_cadndiate
            workset_circuits = new_workset_circuits
        raise NotFoundCircuitException("Not Found!")

    # def naive_bfs_enumerate_synthesis(self):
    #     cnt = 0
    #     for num_instruct in range(1, self.spec.get_max_inst_number()+1):
    #         for instructions_to_apply in list(product(self.gate_operations, repeat=num_instruct)):
    #             synthesized = cirq.Circuit()
    #             synthesized.append(list(instructions_to_apply))
    #             if self.test_synthesis(synthesized, self.spec.equiv_phase):
    #                 logging.info("Synthesis Done.")
    #                 print(synthesized)
    #                 return synthesized
    #             cnt += 1
    #             if cnt % 10000 == 0:
    #                 logging.info(f"Progressing..{cnt}")
    #     logging.info("==NOT FOUND==")
