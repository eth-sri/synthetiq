OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
rx qubits[2];
iswap qubits[2], qubits[1];
rx qubits[0];
iswap qubits[0], qubits[1];
rx qubits[2];
iswap qubits[0], qubits[2];
ry qubits[1];
ry qubits[0];

