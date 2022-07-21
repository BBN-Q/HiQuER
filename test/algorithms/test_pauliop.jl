using QuantumClifford

matX = Complex{Float64}[0 1; 1 0];
matZ = Complex{Float64}[1 0; 0 -1];
matY = Complex{Float64}[0 -1im; 1im 0;];
matI = Complex{Float64}[1 0; 0 1;]

@testset "Pauli Operations" begin

	@test HiQuER.complex(P"ZY") ≈ kron(matZ, matY)

	@test Set(collect(enumerate_paulis(1))) == Set([P"X", P"I", P"Z", P"Y"])

	D = [-1, 1, -1, -1]
	@test D ≈ diag(sum(HiQuER.complex(k)*v for (k,v) in mat_to_paulis(diagm(D))))

	A = randn(4,4)
	M = A'A
	part = greedy_partition(mat_to_paulis(M))
	comm(P"X", P"Y")
	@test all(QuantumClifford.comm(x,y) == 0x0 for g in part for x in g for y in g)

	tt = Tableau(Stabilizer([P"IXX", P"ZYZ", P"XXI"]))
	tg = copy(tt)
	c = Circuit()
	HiQuER.diagX!(tt, tg, c)
	HiQuER.clearX!(tt, tg, c)
	@test iszero(tt.X)
	@test isdiag(tt.Z)
	@test Stabilizer(tg) == Stabilizer([P"-IZZ", P"-ZZI", P"IZI"])

end