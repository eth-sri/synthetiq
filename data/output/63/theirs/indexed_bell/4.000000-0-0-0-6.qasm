OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[1];
cz qubits[1], qubits[2];
cx qubits[3], qubits[0];
cx qubits[1], qubits[0];

