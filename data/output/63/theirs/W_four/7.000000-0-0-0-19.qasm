OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[3];
h qubits[2];
toffoli qubits[0], qubits[2], qubits[3];
cx qubits[0], qubits[1];
toffoli qubits[3], qubits[0], qubits[1];
toffoli qubits[2], qubits[0], qubits[1];
ncncx qubits[1], qubits[3], qubits[2];

