OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
cx qubits[2], qubits[3];
h qubits[1];
h qubits[0];
ncncx qubits[2], qubits[1], qubits[0];
ncncx qubits[3], qubits[0], qubits[2];
ncncx qubits[1], qubits[2], qubits[0];
cx qubits[1], qubits[0];

