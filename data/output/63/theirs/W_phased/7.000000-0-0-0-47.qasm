OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
G_13 qubits[0];
cx qubits[0], qubits[2];
ch qubits[1], qubits[2];
exp_pi23 qubits[1];
exp_pi23 qubits[2];
cx qubits[0], qubits[1];
x qubits[2];

