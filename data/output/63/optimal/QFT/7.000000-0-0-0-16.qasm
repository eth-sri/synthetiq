OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
h qubits[2];
cs qubits[1], qubits[2];
h qubits[1];
cs qubits[0], qubits[1];
swap qubits[2], qubits[0];
ct qubits[0], qubits[2];
h qubits[2];

