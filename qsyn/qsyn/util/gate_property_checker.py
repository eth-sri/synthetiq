import cirq
import numpy as np
from cirq.qis import density_matrix_from_state_vector
from cirq.linalg import sub_state_vector
from typing import Union, List, Tuple, Dict, Set
from itertools import chain, combinations, product
from cirq.linalg.transformations import EntangledStateError
import math

def is_classical_operation(to_check: Union[cirq.Circuit, cirq.Gate, np.array]):
    if isinstance(to_check, cirq.Circuit) or isinstance(to_check, cirq.Gate):
        to_check = cirq.unitary(to_check)

    for i in range(to_check.shape[0]):  # column
        for j in range(to_check.shape[1]):  # row
            curr_elt = to_check[j][i]
            if curr_elt == 1 or np.isclose(curr_elt, 1, atol=1e-04):
                break
            elif not (curr_elt == 0 or np.isclose(curr_elt, 0, atol=1e-04)):
                return False
    return True

def is_phasing_operation(to_check: Union[cirq.Gate, np.array]):
    if is_classical_operation(to_check):  # TODO : Optimize?
        return False
    
    if isinstance(to_check, cirq.Gate):
        to_check = cirq.unitary(to_check)
    for i in range(to_check.shape[0]):  # column
        for j in range(to_check.shape[1]):  # row
            curr_elt = to_check[j][i]
            if not(np.absolute(curr_elt) == 0 or np.isclose( np.absolute(curr_elt) , 1, atol=1e-04)):
                return False 
    return True

def is_nonclassical_operation(to_check: Union[cirq.Gate, np.array]):
    if isinstance(to_check, cirq.Gate):
        to_check = cirq.unitary(to_check)
    for i in range(to_check.shape[0]):  # column
        for j in range(to_check.shape[1]):  # row
            curr_elt = to_check[j][i]
            if (curr_elt != 0 and curr_elt != 1 and not(np.allclose(np.absolute(curr_elt),1,atol=1e-04 ))and  not( np.allclose(curr_elt, 1, atol=1e-04)      ) ): 
                return True
    return False



def is_entangling_operation_for_unitary(to_check: Union[cirq.Circuit, np.array, cirq.Gate]) -> bool:
    assert isinstance(to_check, np.ndarray)
    
    shape = to_check.shape
    num_qubit = math.log2(shape[0])
    for i in range(shape[0]):
        try:
            for j in range(int(num_qubit)):
                sub_state_vector((to_check)[:,i], keep_indices=(j,), atol=1e-04)
        except EntangledStateError:
            return True

    return False
