OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
x qubits[2];
h qubits[2];
h qubits[0];
h qubits[1];
cz qubits[0], qubits[2];
cz qubits[0], qubits[1];
h qubits[0];
cz qubits[0], qubits[1];

