OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
rx qubits[2];
rx qubits[0];
iswap qubits[2], qubits[0];
iswap qubits[2], qubits[1];
ry qubits[0];
rx qubits[2];
iswap qubits[2], qubits[0];
ry qubits[0];

