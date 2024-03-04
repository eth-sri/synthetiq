OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
qft qubits[0], qubits[1];
h qubits[2];
qft qubits[2], qubits[0];
qft qubits[1], qubits[2];
qft qubits[1], qubits[2];
qft qubits[2], qubits[0];
qft qubits[2], qubits[0];

