OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[3];
cx qubits[3], qubits[0];
cx qubits[2], qubits[3];
tdg qubits[3];
t qubits[0];
cx qubits[2], qubits[3];
cx qubits[2], qubits[0];
tdg qubits[3];
t qubits[0];
cx qubits[0], qubits[3];
cx qubits[2], qubits[0];
h qubits[0];
cx qubits[2], qubits[3];
cx qubits[0], qubits[3];
cx qubits[3], qubits[0];
t qubits[0];
s qubits[0];
h qubits[0];
cx qubits[0], qubits[1];
cx qubits[3], qubits[0];
tdg qubits[0];
t qubits[1];
cx qubits[3], qubits[0];
cx qubits[3], qubits[1];
tdg qubits[0];
t qubits[1];
cx qubits[1], qubits[0];
cx qubits[3], qubits[1];
h qubits[1];
cx qubits[3], qubits[0];
cx qubits[1], qubits[0];
cx qubits[0], qubits[1];
tdg qubits[0];
h qubits[0];
cx qubits[0], qubits[1];
cx qubits[3], qubits[0];
tdg qubits[0];
t qubits[1];
cx qubits[3], qubits[0];
cx qubits[3], qubits[1];
tdg qubits[0];
t qubits[1];
cx qubits[1], qubits[0];
cx qubits[3], qubits[1];
cx qubits[3], qubits[0];
h qubits[3];
h qubits[1];
cx qubits[3], qubits[2];
cx qubits[1], qubits[0];
cx qubits[0], qubits[3];
tdg qubits[3];
t qubits[2];
cx qubits[0], qubits[3];
cx qubits[0], qubits[2];
tdg qubits[3];
t qubits[2];
cx qubits[2], qubits[3];
cx qubits[0], qubits[2];
h qubits[2];
cx qubits[2], qubits[3];
cx qubits[0], qubits[3];
cx qubits[3], qubits[2];
s qubits[2];
cx qubits[0], qubits[1];
cx qubits[3], qubits[2];
sdg qubits[3];
sdg qubits[1];
sdg qubits[1];

