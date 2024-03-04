OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
h qubits[1];
G3 qubits[0];
cx qubits[1], qubits[2];
ch qubits[3], qubits[0];
toffoli qubits[1], qubits[3], qubits[0];
toffoli qubits[0], qubits[3], qubits[1];
cx qubits[0], qubits[3];

