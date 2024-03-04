OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[2];
h qubits[0];
toffoli qubits[1], qubits[2], qubits[0];
toffoli qubits[0], qubits[1], qubits[2];
ncncx qubits[3], qubits[2], qubits[0];
cx qubits[1], qubits[2];

