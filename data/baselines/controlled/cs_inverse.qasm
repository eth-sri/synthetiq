OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
tdg qubits[0];
tdg qubits[1];
cx qubits[0],qubits[1];
t qubits[1];
cx qubits[0],qubits[1];
