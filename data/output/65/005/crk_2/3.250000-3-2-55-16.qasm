OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
tdg qubits[0];
sdg qubits[0];
cx qubits[1], qubits[0];
tdg qubits[0];
sdg qubits[0];
sdg qubits[0];
t qubits[1];
cx qubits[1], qubits[0];
s qubits[1];
s qubits[1];

