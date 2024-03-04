OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[2];
h qubits[0];
ncncx qubits[1], qubits[2], qubits[0];
ncncx qubits[3], qubits[0], qubits[1];
ncncx qubits[2], qubits[0], qubits[1];
cx qubits[2], qubits[0];

