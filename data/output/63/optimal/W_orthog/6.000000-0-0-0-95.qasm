OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
G_13 qubits[0];
cx qubits[0], qubits[2];
cx qubits[0], qubits[1];
ch qubits[2], qubits[0];
cx qubits[2], qubits[0];
x qubits[2];

