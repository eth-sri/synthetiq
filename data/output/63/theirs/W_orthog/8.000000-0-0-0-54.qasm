OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
cx qubits[2], qubits[0];
x qubits[0];
G_13 qubits[1];
ch qubits[0], qubits[1];
x qubits[0];
cx qubits[1], qubits[2];
cx qubits[0], qubits[2];
x qubits[2];

