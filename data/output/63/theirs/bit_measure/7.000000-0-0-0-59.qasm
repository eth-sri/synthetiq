OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[3];
cx qubits[3], qubits[2];
G3 qubits[0];
cx qubits[0], qubits[3];
toffoli qubits[1], qubits[2], qubits[3];
ch qubits[1], qubits[0];
cx qubits[1], qubits[0];

