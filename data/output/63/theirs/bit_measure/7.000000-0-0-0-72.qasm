OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[0];
G3 qubits[2];
cx qubits[0], qubits[1];
ch qubits[3], qubits[2];
toffoli qubits[3], qubits[1], qubits[0];
cx qubits[2], qubits[0];
cx qubits[3], qubits[2];

