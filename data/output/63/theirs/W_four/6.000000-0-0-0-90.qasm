OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[2];
h qubits[0];
toffoli qubits[3], qubits[0], qubits[2];
toffoli qubits[2], qubits[0], qubits[3];
ncncx qubits[1], qubits[0], qubits[2];
cx qubits[3], qubits[0];

