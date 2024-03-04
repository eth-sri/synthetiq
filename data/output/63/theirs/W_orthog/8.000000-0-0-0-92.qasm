OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
cx qubits[1], qubits[0];
x qubits[2];
G_13 qubits[0];
ch qubits[2], qubits[0];
cx qubits[2], qubits[1];
x qubits[0];
cx qubits[0], qubits[2];
x qubits[2];

