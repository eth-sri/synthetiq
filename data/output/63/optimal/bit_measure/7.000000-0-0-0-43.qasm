OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
G3 qubits[2];
cx qubits[2], qubits[0];
h qubits[3];
cx qubits[3], qubits[2];
ch qubits[1], qubits[0];
toffoli qubits[1], qubits[2], qubits[3];
cx qubits[1], qubits[0];

