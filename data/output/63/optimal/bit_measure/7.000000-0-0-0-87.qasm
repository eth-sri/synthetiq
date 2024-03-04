OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
G3 qubits[0];
ch qubits[2], qubits[0];
h qubits[3];
toffoli qubits[1], qubits[2], qubits[0];
cx qubits[3], qubits[1];
toffoli qubits[0], qubits[1], qubits[2];
cx qubits[0], qubits[2];

