OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
G_13 qubits[1];
ch qubits[2], qubits[1];
exp_pi23 qubits[1];
cx qubits[2], qubits[1];
cx qubits[2], qubits[0];
x qubits[0];
exp_pi23 qubits[1];

