OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
h qubits[1];
h qubits[0];
cz qubits[1], qubits[0];
x qubits[1];
h qubits[2];
cz qubits[2], qubits[1];
h qubits[1];
cz qubits[2], qubits[1];

