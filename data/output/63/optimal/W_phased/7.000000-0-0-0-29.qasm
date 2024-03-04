OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
G_13 qubits[0];
ch qubits[1], qubits[0];
exp_pi23 qubits[0];
cx qubits[1], qubits[0];
cx qubits[0], qubits[2];
x qubits[0];
exp_pi23 qubits[1];

