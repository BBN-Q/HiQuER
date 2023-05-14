OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];
s q[0];
h q[1];
cx q[1],q[0];
s q[0];
s q[1];
h q[0];
h q[1];
rz(0.0183674693107605) q[0];
rz(0.0183674693107605) q[1];
rz(0.0183674693107605) q[0];
rz(0.0183674693107605) q[1];
h q[0];
h q[1];
sdg q[0];
sdg q[1];
cx q[1],q[0];
sdg q[0];
h q[1];

