OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
t qubits[3];
s qubits[3];
h qubits[3];
cx qubits[2], qubits[3];
tdg qubits[3];
sdg qubits[1];
tdg qubits[2];
cx qubits[2], qubits[3];
sdg qubits[2];
cx qubits[3], qubits[1];
cx qubits[3], qubits[2];
cx qubits[1], qubits[2];
tdg qubits[1];
cx qubits[0], qubits[3];
cx qubits[2], qubits[1];
sdg qubits[2];
tdg qubits[3];
cx qubits[3], qubits[1];
cx qubits[0], qubits[2];
cx qubits[0], qubits[3];
cx qubits[2], qubits[1];
s qubits[3];
s qubits[1];
cx qubits[2], qubits[3];
cx qubits[0], qubits[1];
t qubits[3];
tdg qubits[2];
t qubits[1];
cx qubits[0], qubits[3];
cx qubits[0], qubits[2];
cx qubits[0], qubits[1];
s qubits[0];
h qubits[3];

