OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
rx qubits[1];
rx qubits[0];
iswap qubits[1], qubits[0];
rx qubits[0];
ry qubits[0];
rx qubits[2];
iswap qubits[0], qubits[2];
ry qubits[0];

