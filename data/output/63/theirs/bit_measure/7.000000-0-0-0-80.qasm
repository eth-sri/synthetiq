OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
G3 qubits[3];
ch qubits[1], qubits[3];
toffoli qubits[0], qubits[1], qubits[3];
h qubits[2];
toffoli qubits[3], qubits[2], qubits[0];
cx qubits[3], qubits[1];
cx qubits[2], qubits[0];

