OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
h qubits[1];
h qubits[0];
t qubits[1];
t qubits[0];
scz qubits[1], qubits[0];
h qubits[1];
h qubits[0];
t qubits[1];
t qubits[0];
scz qubits[1], qubits[0];
h qubits[1];
h qubits[0];
s qubits[1];
s qubits[0];

