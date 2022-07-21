@testset "QASM" begin 

	c = Circuit([S[1], CNOT[1,2], T'[2]])
	@test string(to_qasm(c)) == """OPENQASM 2.0;
								include "qelib1.inc";
								qreg q[2];
								s q[0];
								cx q[0], q[1];
								tdg q[1];"""

end