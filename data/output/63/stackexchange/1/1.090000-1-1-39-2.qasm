OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
tdg qubits[1];
sdg qubits[1];
h qubits[1];
s qubits[1];
s qubits[0];
U qubits[1], qubits[0];
sdg qubits[0];
sdg qubits[0];
tdg qubits[1];
h qubits[0];

