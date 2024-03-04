OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
cx qubits[0], qubits[1];
cx12rev qubits[2], qubits[1];
cx qubits[0], qubits[1];
cx12 qubits[2], qubits[1];
cx12 qubits[2], qubits[0];

