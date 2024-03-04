OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
G3 qubits[3];
ch qubits[2], qubits[3];
h qubits[1];
toffoli qubits[0], qubits[1], qubits[2];
cx qubits[3], qubits[0];
toffoli qubits[3], qubits[0], qubits[2];
cx qubits[1], qubits[2];

