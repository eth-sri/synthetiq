OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
toffoli qubits[3], qubits[1], qubits[0];
h qubits[3];
h qubits[2];
toffoli qubits[1], qubits[2], qubits[3];
toffoli qubits[2], qubits[1], qubits[3];
ncncx qubits[0], qubits[3], qubits[2];
cx qubits[1], qubits[3];

