OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
sdg qubits[1];
h qubits[1];
cx qubits[1], qubits[0];
s qubits[0];
h qubits[0];
sdg qubits[1];
t qubits[0];
h qubits[1];
cx qubits[1], qubits[0];
tdg qubits[1];
t qubits[0];
h qubits[1];
cx qubits[1], qubits[0];

