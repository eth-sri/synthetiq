OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
cx qubits[0], qubits[3];
h qubits[3];
cx qubits[1], qubits[0];
t qubits[3];
s qubits[2];
cx qubits[0], qubits[3];
tdg qubits[0];
sdg qubits[0];
cx qubits[1], qubits[2];
cx qubits[1], qubits[0];
s qubits[1];
tdg qubits[2];
cx qubits[2], qubits[1];
tdg qubits[3];
s qubits[1];
cx qubits[2], qubits[3];
cx qubits[2], qubits[1];
t qubits[3];
cx qubits[1], qubits[3];
cx qubits[0], qubits[3];
tdg qubits[3];
h qubits[3];
cx qubits[1], qubits[3];
tdg qubits[0];
cx qubits[0], qubits[3];
sdg qubits[0];
sdg qubits[0];

