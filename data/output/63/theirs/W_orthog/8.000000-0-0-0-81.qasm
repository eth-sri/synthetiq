OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
ch qubits[0], qubits[2];
ch qubits[1], qubits[0];
G_13 qubits[0];
cx qubits[0], qubits[1];
ch qubits[0], qubits[1];
cx qubits[0], qubits[2];
cx qubits[1], qubits[0];
x qubits[2];

