OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
sdg qubits[0];
h qubits[3];
cx qubits[3], qubits[0];
s qubits[0];
sdg qubits[2];
sdg qubits[2];
cx qubits[1], qubits[0];
cx qubits[1], qubits[2];
t qubits[0];
cx qubits[2], qubits[0];
t qubits[0];
tdg qubits[3];
cx qubits[3], qubits[0];
sdg qubits[3];
s qubits[0];
cx qubits[2], qubits[3];
cx qubits[2], qubits[0];
cx qubits[1], qubits[0];
tdg qubits[3];
tdg qubits[0];
cx qubits[0], qubits[3];
t qubits[1];
s qubits[2];
s qubits[1];
h qubits[3];
