OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
t qubits[1];
h qubits[1];
sdg qubits[0];
s qubits[1];
U qubits[0], qubits[1];
s qubits[0];
s qubits[0];
h qubits[0];
tdg qubits[1];
s qubits[0];
s qubits[0];

