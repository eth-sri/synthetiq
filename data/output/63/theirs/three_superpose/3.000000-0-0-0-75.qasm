OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
G3 qubits[1];
ch qubits[0], qubits[1];
x qubits[1];

