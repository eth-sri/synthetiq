OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[3];
h qubits[2];
cx qubits[0], qubits[1];
toffoli qubits[1], qubits[2], qubits[3];
ncncx qubits[0], qubits[2], qubits[3];
cx qubits[1], qubits[3];
cx qubits[1], qubits[2];

