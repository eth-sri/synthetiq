OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
G_13 qubits[2];
cx qubits[2], qubits[0];
ch qubits[1], qubits[2];
x qubits[2];
exp_pi23 qubits[2];
exp_pi23 qubits[1];
exp_pi23 qubits[1];

