OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
h qubits[2];
cx qubits[2], qubits[0];
h qubits[0];
cx qubits[2], qubits[0];
h qubits[0];
cx qubits[0], qubits[1];

