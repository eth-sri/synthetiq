import itertools
import cirq
import numpy as np
import math
from cirq.qis import density_matrix_from_state_vector
from cirq.linalg import sub_state_vector
from typing import Union, List, Tuple, Dict, Set
from itertools import chain, combinations, product
from cirq.linalg.transformations import EntangledStateError
# print(__name__)
if __name__ == "util.utils":
    from util.gate_property_checker import *
else : 
    from qsyn.util.gate_property_checker import *

NC = "NoClassical"
BOOL = "Bool"
PHASING = "PHASING"
ENTANGLE = "Entangle"

CONTROL = "Controls"
TARGET = "Targs"

# LOSS_INFO = [BOOL, PHASING]
LOSS_INFO = "LOSS_INFO"
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    from itertools import chain, combinations
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

def partitions_of_list(s):
    # [A,B,C,D] => [[A],[B,C,D]] , [[A,B], [C,D]], ...
    if len(s) > 0:
        for i in range(1, len(s) + 1):
            first, rest = s[:i], s[i:]
            for p in partitions_of_list(rest):
                yield [first] + p
    else:
        yield []

def partition_to_two_sub(s):
    res_partitions = list()
    for i in range(0, math.floor(len(s) / 2) + 1):
        if i == 0:
            res_partitions.append((tuple(), s))
        else:
            for x in combinations(s, r=i):
                elts_of_x = [k for k in x]
                lo_elts_of_x = [k for k in s if k not in elts_of_x]
                res_partitions.append((tuple(elts_of_x), tuple(lo_elts_of_x)))
    return res_partitions

def component_priors(components: cirq.Gate, is_att_guide = False):
    # TODO : Optimize Here
    prior = dict()
    bool_gates = []
    phasing_gates = []
    nc_gates = [] # nc = superposition (meaning non_classical)
    entangling_gate = []
    valid_num_qubits = set()

    for gate in components:
        # if isinstance(gate, cirq.ControlledGate):
        #     if is_nonclassical_operation(cirq.unitary(gate.sub_gate)):
        #         nc_gates.append(gate.sub_gate)
        #     if is_phasing_operation(cirq.unitary(gate.sub_gate)):
        #         phasing_gates.append(gate.sub_gate)
            #     nc_gates.append(gate.sub_gate)
            # if is_entangling_operation_for_unitary(cirq.unitary(gate.sub_gate)):
            #     entangling_gate.append(gate.sub_gate)
            # elif is_nonclassical_operation(cirq.unitary(gate.sub_gate)):
            #     nc_gates.append(gate.sub_gate)
            # elif is_phasing_operation(cirq.unitary(gate.sub_gate)):
            #     phasing_gates.append(gate.sub_gate)
            # elif is_classical_operation(cirq.unitary(gate.sub_gate)):
            #     bool_gates.append(gate.sub_gate)

        valid_num_qubits.add(gate.num_qubits())
        if is_entangling_operation_for_unitary(cirq.unitary(gate)):
            entangling_gate.append(gate)
        elif is_nonclassical_operation(cirq.unitary(gate)):
            nc_gates.append(gate)
        elif is_phasing_operation(cirq.unitary(gate)):
            phasing_gates.append(gate)
        elif is_classical_operation(cirq.unitary(gate)):
            bool_gates.append(gate)
    prior[NC]       = nc_gates
    prior[BOOL]     = bool_gates
    prior[PHASING]  = phasing_gates
    prior[ENTANGLE] = entangling_gate


    print("bool",prior[BOOL]    ) 
    print("phasing",prior[PHASING] ) 
    print("sp",prior[NC]      ) 
    print("entangle",prior[ENTANGLE]) 



    prior["VALID_NUM_QUBITS"] = valid_num_qubits
    prior["INVOLUTIONARY"] = []
    prior["INVERSES"] = dict()
    for gate in components:
        if cirq.linalg.allclose_up_to_global_phase(np.matmul(cirq.unitary(gate), cirq.unitary(gate)), np.identity(cirq.unitary(gate).shape[0])):
            prior["INVOLUTIONARY"].append(gate)

    for gate_A in components:
        for gate_B in components:
            if (gate_A != gate_B and gate_A.num_qubits() == gate_B.num_qubits() and np.allclose(np.matmul(cirq.unitary(gate_A), cirq.unitary(gate_B)), np.identity(cirq.unitary(gate_A).shape[0]))):
                prior["INVERSES"][gate_A] = gate_B

    prior["IDENTITY_N"] = dict()
    for gate_A in components:
        flag = False
        if gate_A not in prior["INVOLUTIONARY"]:
            for i in range(3,11):
                if cirq.linalg.allclose_up_to_global_phase(np.linalg.matrix_power(cirq.unitary(gate), i), np.identity(cirq.unitary(gate).shape[0])):
                    prior["IDENTITY_N"][gate_A] = i
                    flag=True
                if flag :
                    break
    print(prior["INVOLUTIONARY"])
    print(prior["IDENTITY_N"])
    print(prior["IDENTITY_N"].keys())
    print(prior["INVERSES"])
    return prior

def get_state_property(state_vec : np.array): 
    # import warnings
    # warnings.warn('get_state_property is hardocded yet')
    if is_nonclassical_state(state_vec):
        return NC
    else : return LOSS_INFO
  
######### ## state prop ## #########
def is_nonclassical_state(state_vec : np.array):
    return len(np.nonzero(state_vec)[0]) > 1

def is_classic_state(state_vec: np.array, atol=None):
    if atol is None:
        atol = 1e-07
    assert len(state_vec.shape) == 1
    for i in range(state_vec.shape[0]):
        if abs(state_vec[i] - 0) < atol:
            continue
        elif abs(state_vec[i] - 1) < atol:
            return True
    return False

def is_entangled_state(state_vec: np.array, subsys_A):
    dens_mat = density_matrix_from_state_vector(state_vec, indices=subsys_A)  # Check the concrete action of density_matrix_from_state_vector
    res_mat = np.linalg.matrix_power(dens_mat, 2)
    return not np.allclose(1, np.trace(res_mat), atol=1e-07)

##############################


def basis_rep(state_vec: np.array) -> str:
    return cirq.qis.dirac_notation(state_vec)


def is_same_amp_type(amp1, amp2) -> bool:
    return np.allclose(np.isreal(amp1), np.isreal(amp2))


def is_hermitian(mat: np.ndarray, atol=None):
    if atol is None:
        atol = 1e-04
    return np.allclose(mat, mat.H, atol=atol)



def split(array, nrows, ncols):
    r, h = array.shape
    return (array.reshape(h // nrows, nrows, -1, ncols)
                 .swapaxes(1, 2)
                 .reshape(-1, nrows, ncols))


def is_entangling_operation(to_check: Union[cirq.Circuit, np.array, cirq.Gate], M: int = None, N: int = None) -> bool:

    assert isinstance(to_check, cirq.Gate)
    for i in range(2**to_check.num_qubits()):
        try:
            for j in range(to_check.num_qubits()):
                sub_state_vector(cirq.unitary(to_check)[:,i], keep_indices=(j,), atol=1e-04)
        except EntangledStateError:
            return True

    return False


class InvalidMomentCheck(Exception):  
    def __init__(self):
        super().__init__('Invalid moment check.')


def is_to_exclude(qc: cirq.Circuit, component_prior):
    if len(qc) == 1:
        return False
    involutionary_gates = component_prior["INVOLUTIONARY"]
    for idx, moment in enumerate(qc):
        if idx == len(qc) - 1:
            return False

        for res_gate_op in moment.operations:
            try:
                if res_gate_op.gate not in involutionary_gates:
                    raise InvalidMomentCheck()
                next_moment_idx = qc.next_moment_operating_on(
                    (res_gate_op.qubits), idx + 1)
                if not next_moment_idx:
                    raise InvalidMomentCheck()
                next_moment_gate_op = qc[next_moment_idx].operation_at(
                    res_gate_op.qubits[0])
                if next_moment_gate_op == res_gate_op:
                    return True
            except InvalidMomentCheck as e:
                pass
    # TODO : [Exclude QC Checker] Enhancement, Exclude QC that contains sequene that is in inverse relation
    return False


def partitions(n, I=1):
    # From https://stackoverflow.com/questions/10035752/elegant-python-code-for-integer-partitioning
    yield (n,)
    for i in range(I, n // 2 + 1):
        for p in partitions(n - i, i):
            yield (i,) + p


def pairwise_disjoint(sets):
    union = set().union(*sets)
    n = sum(len(u) for u in sets)
    return n == len(union)


# def is_classical_sv(sv: np.array):
#     return is_nonclassical_state(sv)


def classic_sv_to_binarys(sv: np.array):
    in_dirac = cirq.qis.dirac_notation(sv)
    binary_str = in_dirac.replace("⟩", "").replace("|", "")
    if not "j" in binary_str:
        return list(binary_str)
    else : return None

def chunks(li, n):
    # https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    for i in range(0, len(li), n):
        yield li[i:i + n]


def possible_moments_of_qregsize(qreg_size: int, valid_num_qubits: Set[int]):
    full_qreg = tuple(range(qreg_size))
    res_moments = []
    collected_permparti = set()
    for i in range(1, qreg_size + 1):
        for partition in partitions(i):
            if not set(partition).issubset(valid_num_qubits):
                continue
            collected_permparti.add(partition)
    for qbitnum_placement in sorted(list(collected_permparti)):
        already_coverd = []
        qbit_placements = []
        for k in qbitnum_placement:
            qbit_placements.append(
                set([l for l in combinations(range(qreg_size), r=k)]))
        for state_placements in product(*qbit_placements):
            if set([l for l in state_placements]) in already_coverd:
                continue
            if pairwise_disjoint(state_placements):
                state = [('Generic', tuple(l)) for l in state_placements]
                res_moments.append(state)
                already_coverd.append(set([l for l in state_placements]))
    return res_moments


    
def is_multiplexor(gate: cirq.Gate):
    # TODO : Multiplexor Detector
    # import warnings
    # warnings.warn('is_mutiplexor ftn is not yet implemented')
    return False


def is_identity_mapping(in_arr: np.array, out__arr: np.array):
    return np.allclose(in_arr, out__arr, atol=1e-04)


def entangling_condition(sv: np.array):
    num_qubits = int(math.log2(len(sv)))
    res = rec_entangling_condition(sv=sv, sys=list(
        range(num_qubits)), collected_partis=list())
    res = [frozenset(r) for r in res]
    return set(res)


def to_fed_subsv_cal(sys, subsys):
    to_return = list()
    for idx, x in enumerate(sys):
        if x in subsys:
            to_return.append(idx)
    return to_return


def rec_entangling_condition(sv: np.array, sys: List[int], collected_partis: List[List[int]]) -> List[List[int]]:
    num_qubits = int(math.log2(len(sv)))
    max_parti = math.floor(num_qubits / 2)
    already_checked = list()
    # check one-and-others, two-and-other, dotdotdot..
    for i in range(1, max_parti + 1):
        # checking i-and-others
        subsys_As = [x for x in combinations(sys, i)]
        for subsys_A in subsys_As:
            try:
                subsys_B = list(set(sys) - set(subsys_A))
                if subsys_B not in already_checked and subsys_A not in already_checked:
                    res_sub_sv_A = sub_state_vector(
                        sv, keep_indices=list(to_fed_subsv_cal(sys, subsys_A)), atol=1e-04)
                    res_sub_sv_B = sub_state_vector(
                        sv, keep_indices=list(to_fed_subsv_cal(sys, subsys_B)), atol=1e-04)
                    return (rec_entangling_condition(sv=res_sub_sv_A, sys=list(subsys_A), collected_partis=collected_partis.copy())
                            + rec_entangling_condition(sv=res_sub_sv_B, sys=list(subsys_B), collected_partis=collected_partis.copy()))
                already_checked.append(subsys_A)
                already_checked.append(subsys_B)
            except EntangledStateError:
                pass
    collected_partis.append(sys)
    return collected_partis

def amplitudes_same_upto_order(arr1, arr2, atol) -> bool:
    #numpy different
    arr1_nonzeros = (arr1[np.nonzero(arr1)])
    arr2_nonzeros = (arr2[np.nonzero(arr2)])
    indicies = [i for  i in range(len(arr1_nonzeros))]
    
    temp_bool_li =   [not ((arr1_nonzeros[i] > 0) ^ (arr2_nonzeros[i] > 0))   for  i in range(len(arr1_nonzeros))    ]

    arr1_all_real =  [  np.conjugate(x) == x  for x in  arr1_nonzeros   ] 
    arr2_all_real =  [  np.conjugate(x) == x  for x in  arr2_nonzeros   ] 
    arr_1_plus_cnt = 0
    arr_2_plus_cnt = 0
    for i in range(len(arr1_nonzeros)):
        if np.conjugate(arr1_nonzeros[i]) == arr1_nonzeros[i]:
            if arr1_nonzeros[i] > 0 :
                arr_1_plus_cnt+=1
        if np.conjugate(arr2_nonzeros[i]) == arr2_nonzeros[i]:
            if arr2_nonzeros[i] > 0 :
                arr_2_plus_cnt+=1

    if all(arr1_all_real) and all(arr2_all_real) and arr_1_plus_cnt != arr_2_plus_cnt:
        # print(arr1_nonzeros)
        # print(arr2_nonzeros)
        # print(arr_1_plus_cnt)
        # print(arr_2_plus_cnt)
        # input()
        return False
    if all(arr1_all_real) and all(arr2_all_real) and arr_1_plus_cnt == arr_2_plus_cnt and len(arr1_nonzeros) == len(arr2_nonzeros):
        return True
    elif (all(arr1_all_real) and not all(arr2_all_real)) or  (not all(arr1_all_real) and  all(arr2_all_real)):
        return False 
    # input()
    for x in itertools.permutations(arr1_nonzeros ,len(indicies)):
        if np.allclose(np.array(x), arr2_nonzeros, atol=atol) :return True
    return False


def get_all_choices_of_qreg(qreg : List[int]):
    qreg = [  i for q in qreg for  i in q  ]
    res = [x for i in range(1,len(qreg)+1) for x in itertools.combinations(qreg, r=i)]
    return res
    
def into_n_partitions(li : List, n : int):
    # generate & check based :(
    to_return = list()
    choices_of_li = [ x for i in range(0, len(li) +1) for x in itertools.combinations(li, r = i) ] 
    

    for y in itertools.product(choices_of_li, repeat= n) :
        flattened = list()
        for x in y:
            flattened+=list(x)
        if set(flattened) == set(li) and len(flattened) == len(li):
            to_return.append(y)
    return to_return


# should be implemented in here due to cyclic import
def check_action_property_of_gate(gate, act_prop):
    # for single qubit gate only
    # assert gate.num_qubits() == 1 
    if act_prop == ENTANGLE:
        return is_entangling_operation(gate)
    elif act_prop == PHASING:
        return is_phasing_operation(gate)
    elif act_prop == BOOL:
        return is_classical_operation(gate)
    elif act_prop == None:
        return True
    else:
        return is_nonclassical_operation(gate)

def test_io(ss , qc : cirq.Circuit):
    for in_nparr, out_nparr in ss.spec.spec_object.get_io_pairs():
        res_by_curr_Vqc = ss.simulator.simulate(qc, qubit_order=ss.working_qubit, initial_state=in_nparr).final_state_vector
        print(f"Curr action : {basis_rep(in_nparr)}->{basis_rep(res_by_curr_Vqc)}")
        # print(entangling_condition(res_by_curr_Vqc))

def check_continuity(my_list):
    #https://stackoverflow.com/questions/48596542/how-to-check-all-the-integers-in-the-list-are-continuous
    return all(a+1==b for a, b in zip(my_list, my_list[1:]))

def derive_representative_res_in_comp_basis(gate : cirq.Gate):
    unitary = cirq.unitary(gate)
    bool_li = [is_classic_state( unitary[:,i]) for i in range(0, gate.num_qubits()) ]
    for idx, b in reversed(list(enumerate(bool_li))):
        if not b :
            return unitary[:,idx]
    return unitary[:,-1]


def operations_at_q(module, q):
    res = list()
    for i in range(len(module)):
        op = module.operation_at(q,i)
        if op:
            res.append(op)
    return res

def qubits_of_sp_in_entangling_module(module, component_prior):
    res = list()
    for q in module.all_qubits():
        op = module.operation_at(q,0)
        if op and op.gate in component_prior[NC]:
            res.append(q)
    return res


def count_gate(qc : cirq.Circuit)->int :
    count = 0
    for moment in qc:
        for gate_op in moment:
            count+=1
    return count


def gate_op_in_list_by_unitary(gate_op, li_of_gate_op,working_qubits) -> bool :
    
    for other_op in li_of_gate_op :
        if other_op.gate != gate_op.gate:
            continue
        elif set(other_op.qubits) != set(gate_op.qubits):
            continue
        else:
            other_op_qc = cirq.Circuit()
            gate_op_qc = cirq.Circuit()
            other_op_qc.append(cirq.IdentityGate(len(working_qubits))(*working_qubits))
            other_op_qc.append(other_op)
            gate_op_qc.append(cirq.IdentityGate(len(working_qubits))(*working_qubits))
            gate_op_qc.append(gate_op)
            if np.allclose(cirq.unitary(gate_op_qc), cirq.unitary(other_op_qc),atol=1e-04):
                return True
    return False
    
def normalize_gate_operations(gate_operations : List[cirq.Operation], working_qubits):
    normalized_ones = list()
    for gate_op in gate_operations :
        if not gate_op_in_list_by_unitary(gate_op,normalized_ones,working_qubits):
            normalized_ones.append(gate_op)
    return normalized_ones




#  def check_qubits_of_sp_position(module : cirq.Circuit, component_prior):
#     res = list()
#     for q in module.all_qubits():
#         print(operations_at_q(module,q))
#         for op in operations_at_q(module,q):
#             if op.gate in component_prior[NC]:
#                 res.append(q)
#                 break
#     print(module)
#     print(res)
#     input()