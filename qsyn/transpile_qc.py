def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn


import numpy as np
from qiskit import QuantumCircuit
from cirq import Simulator
from qiskit import transpile
from qiskit import circuit
import argparse
from qsyn.set_synthesis import *
from qiskit.extensions import UnitaryGate
from qsyn.state_search import * 
from warnings import filterwarnings
# filterwarnings(action='ignore', category=DeprecationWarning, message='`np.int` is a deprecated alias')
# filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')


if __name__=="__main__":
    basis_gates = ['u3', 'cx']
    CT = cirq.ops.ControlledGate(sub_gate = cirq.T)
    CS = cirq.ops.ControlledGate(sub_gate = cirq.S)
    circuit = cirq.Circuit()
    c,b,a =cirq.LineQubit.range(3)
    circuit.append(cirq.H(a))
    circuit.append(CS(b,a))
    circuit.append(CT(c,a))
    circuit.append(cirq.H(b))
    circuit.append(CS(c,b))
    circuit.append(cirq.H(c))
    circuit.append(cirq.SWAP(a,c))
    spec_QFT = cirq.unitary(circuit)
    spec_cluster = [0.5, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, -0.5]
    # construct by cirq then get unitary 
    QFT_2_WO_REV = cirq.ops.QuantumFourierTransformGate(2, without_reverse=True)
    QFT_2_WO_REV_INV = cirq.ops.QuantumFourierTransformGate(2, without_reverse=True) ** -1
    C_R_2 = cirq.ops.ControlledGate(sub_gate=cirq.ZPowGate(exponent=2 / (2**2)))
    C_R_1 = cirq.ops.ControlledGate(sub_gate=cirq.ZPowGate(exponent=2 / (2**1)))
    d,c,b,a, = cirq.LineQubit.range(4)
    circuit = cirq.Circuit()
    circuit.append(QFT_2_WO_REV(c, d))
    circuit.append(C_R_1(a, c))
    circuit.append(C_R_1(b, d))
    circuit.append(C_R_2(b, d))
    circuit.append(QFT_2_WO_REV_INV(c, d))
    spec_draper_adder = cirq.unitary(circuit)


    # for cluster
    circuit_to_transpile = QuantumCircuit(4)
    circuit_to_transpile.initialize(np.array(spec_cluster), list(range(len(circuit_to_transpile.qubits))))
    transpiled_circuit = transpile(circuit_to_transpile, basis_gates=basis_gates, optimization_level=3)
    print("==============================")
    print("Transpiled Circuit of \'cluster\'")
    print("==============================")
    print(transpiled_circuit)
    # for QFT
    circuit_to_transpile = QuantumCircuit(3)
    circuit_to_transpile.append(UnitaryGate(spec_QFT), circuit_to_transpile.qubits)
    transpiled_circuit = transpile(circuit_to_transpile, basis_gates=basis_gates,optimization_level=3)
    print("==============================")
    print("Transpiled Circuit of \'QFT\'")
    print("==============================")

    print(transpiled_circuit)
    # for draper
    circuit_to_transpile = QuantumCircuit(4)
    circuit_to_transpile.append(UnitaryGate(spec_draper_adder), circuit_to_transpile.qubits)
    transpiled_circuit = transpile(circuit_to_transpile, basis_gates=basis_gates,optimization_level=3)
    print("==============================")
    print("Transpiled Circuit of \'draper\'")
    print("==============================")

    print(transpiled_circuit)


