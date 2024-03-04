OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
h qubits[2];
swap qubits[2], qubits[0];
cs qubits[0], qubits[1];
h qubits[1];
ct qubits[0], qubits[2];
cs qubits[2], qubits[1];
h qubits[2];

