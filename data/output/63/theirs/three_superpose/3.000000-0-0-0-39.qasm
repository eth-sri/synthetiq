OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
G3 qubits[0];
ch qubits[1], qubits[0];
x qubits[0];

