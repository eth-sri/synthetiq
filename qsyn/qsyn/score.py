import itertools
from itertools import product, permutations
import cirq
from cirq import Simulator
import numpy as np
from typing import Tuple, Dict, List, Union
from numpy.core.fromnumeric import resize

if __name__ == "score":
    from util.utils import *
    from state import *
    from rules import *
    from synthesis import Synthesis
    from synthesis_spec.specification import Spec
else :
    from qsyn.util.utils import *
    from qsyn.state import *
    from qsyn.rules import *
    from qsyn.synthesis import Synthesis
    from qsyn.synthesis_spec.specification import Spec
import warnings
# ===================================================
# Score for Generic States, for generin state ranking
# ===================================================


# ==========================================================
# Score for V qc, for feeding the appropriate one to Phase-2
# ==========================================================

def score_E(ss: Synthesis, qc_to_check: Union[cirq.Circuit, MomentBasedState], is_state_prep) -> Tuple[int, List]:
    # ss for state search

    if isinstance(qc_to_check , MomentBasedState):
        eval_res = qc_to_check.evaluate(ss, component_gates=ss.spec.component_gates, working_qubits=ss.working_qubit)
        assert len(eval_res) == 1
        qc_to_check = eval_res[0]

    if is_state_prep == False:
        return score_E_unitary_syn(ss, qc_to_check)
    elif is_state_prep == True:
        raise NotImplementedError
        # return score_E_state_prep(ss, qc_to_check)
    else:
        raise Exception("is_state_prep bool must be set")


def score_E_unitary_syn(ss: Synthesis, qc_to_check: cirq.Circuit):
    violated_io_s = list()
    # violated_io_s is list of tuple as (outputspec, currently_evaluated_output_by_input_spec which must be violated )
    score_E = 0

    score_identity_io = 0
    score_criticial_io = 0

    simulator = Simulator()
    io_specs = ss.spec.spec_object.get_io_pairs()
    for io_pair in io_specs:
        in_nparr, out_nparr = io_pair
        res = simulator.simulate(
            qc_to_check, qubit_order=ss.working_qubit, initial_state=in_nparr)
        if is_identity_mapping(in_nparr, out_nparr):
            weight = 5
        else:
            weight = 10
    
        if ss.spec.equiv_phase :
            if cirq.linalg.allclose_up_to_global_phase(res.final_state_vector, out_nparr, atol=1e-04):
                score_E += (1 * weight)
            else : violated_io_s.append((in_nparr, out_nparr, res.final_state_vector))
        else :
            if np.allclose(res.final_state_vector, out_nparr, atol=1e-04):
                score_E += (1 * weight)
            else : violated_io_s.append((in_nparr, out_nparr, res.final_state_vector))

    return score_E, violated_io_s


def score_E_state_prep(ss: Synthesis, qc_to_check: cirq.Circuit):
    violated_io_s = list()
    # violated_io_s is list of tuple as (input_spec, currently_evaluated_output_by_input_spec which must be violated )

    score_E = 0
    simulator = Simulator()
    io_specs = ss.spec.spec_object.get_io_pairs()
    in_nparr, out_nparr = io_specs[0]
    res = simulator.simulate(
        qc_to_check, qubit_order=ss.working_qubit, initial_state=in_nparr)
    res = res.final_state_vector.round(8)
    warnings.warn('score_E for ')
    if np.count_nonzero(res) == np.count_nonzero(out_nparr):
        score_E += 1
        if is_same_amp_type(res, out_nparr):
            score_E += 1
            if np.allclose(res[res != 0], out_nparr[out_nparr != 0]):
                score_E += 1
    return score_E, violated_io_s


def is_identity_mapping(in_arr: np.array, out__arr: np.array):
    return np.allclose(in_arr, out__arr, atol=1e-04)


def sim_measure(arr_one, arr_two):
    return 0
