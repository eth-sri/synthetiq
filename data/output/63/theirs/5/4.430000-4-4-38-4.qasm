OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[3];
t qubits[3];
cx qubits[2], qubits[3];
tdg qubits[3];
cx qubits[1], qubits[3];
cx qubits[0], qubits[3];
t qubits[3];
cx qubits[2], qubits[3];
tdg qubits[3];
h qubits[3];
sdg qubits[3];

