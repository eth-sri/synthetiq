OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
toffoli qubits[0], qubits[3], qubits[1];
h qubits[3];
h qubits[1];
toffoli qubits[2], qubits[1], qubits[3];
toffoli qubits[1], qubits[3], qubits[2];
ncncx qubits[0], qubits[1], qubits[3];
cx qubits[2], qubits[3];

