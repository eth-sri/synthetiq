import cirq
import logging
import time
import json
import os
import numpy as np
import pprint
from typing import List, Tuple, Union
import argparse
import math
if __name__ == "set_synthesis":
    from synthesis_spec.specification import Spec
    from synthesis_spec.utils import IOpairs
    from synthesis import Synthesis, BFSSynthesis
else:
    from .synthesis_spec.specification import Spec
    from .synthesis_spec.utils import IOpairs
    from .synthesis import Synthesis, BFSSynthesis

from timeit import default_timer as timer


class NMR_INTERACT(cirq.Gate):
    def __init__(self):
        super(NMR_INTERACT, self)

    def _num_qubits_(self):
        return 2

    def _unitary_(self):
        return (np.identity(4) + np.kron(cirq.unitary(cirq.Z), cirq.unitary(cirq.Z)) * (0. + 1.j)) * (1 / np.sqrt(2))

    def _circuit_diagram_info_(self, args):
        return "NMR", "NMR"


class BELL(cirq.Gate):
    def __init__(self):
        super(BELL, self)

    def _num_qubits_(self):
        return 2

    def _decompose_(self, qubits):
        a, b = qubits
        yield cirq.H(a)
        yield cirq.CNOT(a, b)

    def _circuit_diagram_info_(self, args):
        return ["BELL"] * self.num_qubits()



class GHZThree(cirq.Gate):
    def __init__(self):
        super(GHZThree, self)

    def _num_qubits_(self):
        return 3

    def _decompose_(self, qubits):
        a, b,c = qubits
        yield cirq.H(a)
        yield cirq.CNOT(a, b)
        yield cirq.CNOT(b, c)


    def _circuit_diagram_info_(self, args):
        return ["GHZ"] * self.num_qubits()


##JSON_ID's##
ID = "ID"
SPEC_TYPE = "Spec Type"
SPEC = "Spec"
QREG_SIZE = "qreg_size"
MAX_INST_ALLOW = "max_inst_allowed"
SPEC_OBJ = "spec_object"
COMPONENTS = "components"
EQUIV_PHASE = "Equiv Phase"

z = np.exp(1j * np.pi / 4) # e^(iπ/4)

# Your two-qubit entangling gate as numpy array:
basis_unitary = 0.5 * np.array([
    [1j*z, z, z, -1j*z],
    [1/z, 1j/z, -1j/z, 1/z],
    [1/z, -1j/z, 1j/z, 1/z],
    [-1j*z, z, z, 1j*z]
])


##Component Gates##
TO_CIRQ_GATE = {
    "H": cirq.H,
    "CZ": cirq.CZ,
    "SCZ": cirq.CZ ** 0.5,
    "U": cirq.ops.MatrixGate(basis_unitary),
    "CX": cirq.CX,
    "CH": cirq.ops.ControlledGate(cirq.H),
    "CNOT": cirq.CNOT,
    "iSWAP": cirq.ISWAP,
    "X": cirq.X,
    "Y": cirq.Y,
    "X^1/2": cirq.X**0.5,
    "C_X^1/2": cirq.ops.ControlledGate(cirq.X**0.5),
    "ANTI_C_X^1/2": (cirq.X**0.5).controlled(num_controls=1, control_values=(0,)),
    "C_(X^1/2_REV)": cirq.ops.ControlledGate(cirq.X**(-0.5)),
    "Y^1/2": cirq.Y**0.5,
    "exp(2pi/3)" : cirq.ZPowGate(exponent= 2 /3 , global_shift = 0),
    "Z^1/2": cirq.Z**0.5,
    "SWAP": cirq.SWAP,
    "T": cirq.T,
    "S": cirq.S,
    "T_REV": cirq.T ** -1,
    "S_REV": cirq.S ** -1,
    "Y": cirq.Y,
    "Z": cirq.Z,
    "QFT_1_W/O_REV": cirq.ops.QuantumFourierTransformGate(1, without_reverse=True),
    "QFT_1_W/O_REV_INV": cirq.ops.QuantumFourierTransformGate(1, without_reverse=True) ** -1,
    "QFT_2_W/O_REV": cirq.ops.QuantumFourierTransformGate(2, without_reverse=True),
    "QFT_2_W_REV": cirq.ops.QuantumFourierTransformGate(2, without_reverse=False),
    "QFT_2_W/O_REV_INV": cirq.ops.QuantumFourierTransformGate(2, without_reverse=True) ** -1,
    "QFT_3": cirq.ops.QuantumFourierTransformGate(3),
    "QFT_4": cirq.ops.QuantumFourierTransformGate(4),
    "R_2": cirq.ZPowGate(exponent=2 / (2**2)),  # S
    "R_1": cirq.ZPowGate(exponent=2 / (2**1)),  # Z
    "C_R_2": cirq.ops.ControlledGate(
        sub_gate=cirq.ZPowGate(exponent=2 / (2**2)) #FOR DRAPER CS
    ),
    "CS": cirq.ops.ControlledGate(
        sub_gate=cirq.ZPowGate(exponent=2 / (2**2)) #FOR DRAPER CS
    ),
    "CT": cirq.ops.ControlledGate(sub_gate = cirq.T),
    "C_R_1": cirq.ops.ControlledGate(
        sub_gate=cirq.ZPowGate(exponent=2 / (2**1)) #FOR DRAPER CZ
    ),
    "RZ_pi/2": cirq.rz(math.pi/2),
    "sqrt_SWAP": cirq.SWAP ** 0.5,
    "NMR_INTERACT": NMR_INTERACT(),
    "RZ(-pi/8)": cirq.rz(-1/8 * math.pi),
    "RZ(-pi/4)": cirq.rz(-math.pi * 1/4),
    "RZ(pi/4)": cirq.rz((math.pi * 1/4)),
    "RY(-pi)": cirq.ry(-math.pi),
    "RY(pi)": cirq.ry(math.pi),
    "Ry_for_W": cirq.ry(2 * np.arccos(1 / math.sqrt(3))),  # concide with $G_{1/3}$ gate. (1/sqrt(3), -sqrt(2/3), sqrt(2/3), 1/sqrt(3))
    "G2": cirq.ry(2 * np.arccos(1 / math.sqrt(2))),  
    "CG2": (cirq.ry(2 * np.arccos(1 / math.sqrt(2)))).controlled(num_controls=1),  
    "G3": cirq.ry(2 * np.arccos(1 / math.sqrt(3))),  # concide with $G_{1/3}$ gate. (1/sqrt(3), -sqrt(2/3), sqrt(2/3), 1/sqrt(3))
    "CG3": (cirq.ry(2 * np.arccos(1 / math.sqrt(3)))).controlled(num_controls=1),  # concide with $G_{1/3}$ gate. (1/sqrt(3), -sqrt(2/3), sqrt(2/3), 1/sqrt(3))
    "CG_2/5": (cirq.ry(2 * np.arccos(math.sqrt(2) / math.sqrt(5)))).controlled(num_controls=1),  
    "G4": cirq.ry(2 * np.arccos(1 / math.sqrt(4))),  
    "G_2/5": cirq.ry(2 * np.arccos(math.sqrt(2) / math.sqrt(5))),  
    "G_2/12": cirq.ry(2 * np.arccos(math.sqrt(2) / math.sqrt(12))),  # 
    "G_9/10": cirq.ry(2 * np.arccos(math.sqrt(9) / math.sqrt(10))),  # 
    "CG2": (cirq.ry(2 * np.arccos(1 / math.sqrt(2)))).controlled(num_controls=1),  
    "CG_9/10": (cirq.ry(2 * np.arccos(math.sqrt(9) / math.sqrt(10)))).controlled(num_controls=1),  
    "G5": cirq.ry(2 * np.arccos(1 / math.sqrt(5))),  # 
    "CG5": (cirq.ry(2 * np.arccos(1 / math.sqrt(5)))).controlled(num_controls=1),  
    "G6": cirq.ry(2 * np.arccos(1 / math.sqrt(6))),  
    "G7": cirq.ry(2 * np.arccos(1 / math.sqrt(7))),  
    "G_{2/3}": cirq.ry(2 * np.arccos(math.sqrt(2) / math.sqrt(3))),  
    "G_1/3": cirq.ry(2 * np.arccos(1 / math.sqrt(3))),
    "G_{1/10}": cirq.ry(2 * np.arccos(math.sqrt(1) / math.sqrt(10))),  
    "G_{5/9}": cirq.ry(2 * np.arccos(math.sqrt(5) / math.sqrt(9))),  
    "CG_{1/10}": (cirq.ry(2 * np.arccos(math.sqrt(1) / math.sqrt(10)))).controlled(num_controls=1),  
    "CG_{5/9}": (cirq.ry(2 * np.arccos(math.sqrt(5) / math.sqrt(9)))).controlled(num_controls=1),  
    "Rx_pi/2" : cirq.rx( math.pi/2 ),
    "Rx_-pi/2" : cirq.rx( -math.pi/2 ),
    "Ry_pi/2"  : cirq.ry(  math.pi/2 ),
    "Ry_pi/4"  : cirq.ry(  math.pi/4 ),
    "Ry_-pi/4"  : cirq.ry( -math.pi/4 ),
    "Toffoli"  : cirq.CCNOT,
    "BELL" : BELL(),
    "GHZ3" : GHZThree(),
    "CSWAP" : cirq.CSWAP,
    "CCCNOT" : (cirq.X).controlled(num_controls=3, control_values=(1,1,1)),
    "NCNOT" : (cirq.X).controlled(num_controls=1, control_values=(0,)),
    "NCNCX" : (cirq.X).controlled(num_controls=2, control_values=(0,0)),
}


def load_gate_components(components: str) -> List[cirq.Gate]:
    to_return = []
    components = components.replace(" ", "").split(",")
    for gate in components:
        to_return.append(TO_CIRQ_GATE[gate])

    return to_return


def load_benchmark(bench_id: str) -> Spec:
    dirname = os.path.dirname(__file__)

    with open(dirname + "/../benchmarks/" + bench_id + ".json") as json_file:
        json_data = json.load(json_file)

    if json_data[SPEC_TYPE] == Spec.PART_SYN_SPEC:
        ##load io-puts##
        targ_obj = IOpairs(io_pairs=[])

        # https://stackoverflow.com/questions/4841782/python-constructor-and-default-value
        for inout in json_data[SPEC][SPEC_OBJ]:
            input_str, output_str = inout["input"], inout["output"]
            input_np = np.fromstring(input_str, dtype="complex", sep=',')
            output_np = np.fromstring(output_str, dtype="complex", sep=',')
            targ_obj.append((input_np, output_np))
    elif json_data[SPEC_TYPE] == Spec.FULL_SYN_SPEC:
        targ_obj = []
        for row in json_data[SPEC][SPEC_OBJ]:
            targ_obj.append(np.fromstring(row, dtype="complex", sep=','))
        targ_obj = np.array(targ_obj)
    else:
        raise Exception(
            f"The \"Spec Type\" must be set in either {Spec.PART_SYN_SPEC} or {Spec.FULL_SYN_SPEC}. The Data was {json_data[SPEC_TYPE]}")

    return Spec(
        ID=json_data[ID],
        target_object=targ_obj,
        qreg_size=int(json_data[SPEC][QREG_SIZE]),
        max_inst_number=int(json_data[SPEC][MAX_INST_ALLOW]),
        components=load_gate_components(json_data[SPEC][COMPONENTS]),
        equiv_phase=json_data[EQUIV_PHASE]
    )
