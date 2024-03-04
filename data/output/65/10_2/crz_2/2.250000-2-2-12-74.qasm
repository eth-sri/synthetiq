OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
t qubits[0];
s qubits[0];
s qubits[0];
cx qubits[1], qubits[0];
t qubits[0];
s qubits[0];
cx qubits[1], qubits[0];
sdg qubits[1];
sdg qubits[1];

