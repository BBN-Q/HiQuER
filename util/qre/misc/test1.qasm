OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];
creg c[2];

h q[0];
t q[1];
x q[1];
x q[1];
t q[1];
tdg q[1];
tdg q[0];
cx q[0], q[1];
