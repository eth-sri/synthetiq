OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[2];
h qubits[0];
ncncx qubits[3], qubits[0], qubits[2];
ncncx qubits[0], qubits[3], qubits[1];
ncncx qubits[1], qubits[3], qubits[2];
cx qubits[0], qubits[2];

