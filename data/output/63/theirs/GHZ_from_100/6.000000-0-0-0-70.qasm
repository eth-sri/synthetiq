OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
h qubits[0];
cx qubits[0], qubits[2];
h qubits[0];
cx qubits[0], qubits[1];
cx qubits[1], qubits[2];
h qubits[2];

