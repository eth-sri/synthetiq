OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
cx qubits[1], qubits[0];
cx12rev qubits[2], qubits[0];
cx qubits[1], qubits[0];
cx12 qubits[2], qubits[1];
cx12 qubits[2], qubits[0];

