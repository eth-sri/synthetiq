OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
swap qubits[0], qubits[2];
h qubits[0];
cs qubits[0], qubits[1];
h qubits[1];
ct qubits[0], qubits[2];
cs qubits[2], qubits[1];
h qubits[2];

