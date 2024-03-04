OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
h qubits[2];
cx qubits[2], qubits[1];
h qubits[2];
cx qubits[1], qubits[2];
h qubits[1];
cx qubits[1], qubits[0];

