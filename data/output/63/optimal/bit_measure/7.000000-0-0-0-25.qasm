OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
G3 qubits[0];
ch qubits[2], qubits[0];
h qubits[1];
cx qubits[2], qubits[3];
cx qubits[1], qubits[0];
cx qubits[1], qubits[3];
cx qubits[0], qubits[2];

