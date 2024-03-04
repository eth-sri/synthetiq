OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
s qubits[0];
s qubits[0];
cx qubits[1], qubits[0];
t qubits[1];
tdg qubits[0];
h qubits[1];
h qubits[0];
cx qubits[0], qubits[1];
h qubits[0];
tdg qubits[0];
sdg qubits[0];
h qubits[1];

