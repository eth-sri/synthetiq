OPENQASM 2.0;
include "qelib1.inc";
gate unitary140489220481632 p0,p1 {
	u3(pi,0.58117859,-0.98961773) p0;
	u3(1.569473,pi/2,-pi) p1;
	cx p0,p1;
	u3(pi,0.84252689,2.4133232) p0;
	u3(1.5711839,1*pi/4,-1.5698606) p1;
	cx p0,p1;
	u3(0,-1.211223,-1.1449714) p0;
	u3(0.0013232938,-3*pi/4,-pi/2) p1;
}
qreg qubits[2];
h qubits[0];
t qubits[0];
h qubits[1];
t qubits[1];
unitary140489220481632 qubits[1],qubits[0];
s qubits[0];
s qubits[0];
h qubits[0];
tdg qubits[0];
s qubits[1];
h qubits[1];
sdg qubits[1];
h qubits[1];
t qubits[1];
unitary140489220479424 qubits[0],qubits[1];
sdg qubits[0];
h qubits[0];
sdg qubits[0];
s qubits[1];
h qubits[1];
sdg qubits[1];
