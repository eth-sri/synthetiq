OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[0];
G3 qubits[1];
cx qubits[0], qubits[2];
ch qubits[3], qubits[1];
cx qubits[3], qubits[0];
toffoli qubits[1], qubits[0], qubits[3];
cx qubits[1], qubits[3];

