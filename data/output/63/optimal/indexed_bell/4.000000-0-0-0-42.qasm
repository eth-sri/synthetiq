OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[0];
cz qubits[0], qubits[2];
cx qubits[0], qubits[1];
cx qubits[3], qubits[0];

