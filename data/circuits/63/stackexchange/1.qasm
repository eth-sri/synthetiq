OPENQASM 2.0;
include "qelib1.inc";
gate unitary140489114080912 p0,p1 {
	u3(pi/2,-pi/2,pi/2) p0;
	u3(0.0013232938,-pi/2,-1.1857182e-13) p1;
	cx p0,p1;
	u3(0,-0.21152499,-1.3592713) p0;
	u3(pi/2,3.1402694,-pi) p1;
}
qreg qubits[2];
s qubits[0];
h qubits[1];
s qubits[1];
unitary140489114080912 qubits[0],qubits[1];
s qubits[0];
s qubits[0];
h qubits[0];
s qubits[1];
s qubits[1];
