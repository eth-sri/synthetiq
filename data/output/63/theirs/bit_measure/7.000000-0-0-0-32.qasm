OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
G3 qubits[2];
cx qubits[2], qubits[1];
ch qubits[3], qubits[1];
h qubits[0];
cx qubits[0], qubits[1];
toffoli qubits[3], qubits[0], qubits[1];
cx qubits[3], qubits[2];

