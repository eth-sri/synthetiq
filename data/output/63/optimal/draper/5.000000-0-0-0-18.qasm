OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[4];
qft2worev qubits[0], qubits[1];
cr2 qubits[2], qubits[1];
cr1 qubits[1], qubits[3];
cr1 qubits[0], qubits[2];
qft2worevinv qubits[0], qubits[1];

