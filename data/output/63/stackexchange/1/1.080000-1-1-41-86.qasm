OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[2];
h qubits[1];
s qubits[0];
s qubits[1];
h qubits[0];
U qubits[0], qubits[1];
h qubits[0];
sdg qubits[0];
h qubits[0];
s qubits[0];

