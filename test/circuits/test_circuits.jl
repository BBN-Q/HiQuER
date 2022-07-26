@testset "Circuits" begin 

	@test H[1] isa Pair
	@test all([x in (H[1]*CNOT[1,2]).qubits for x in (1,2,3)])
	@test_throws ArgumentError H[1]*H[2]

	slices = [Slice(H[1]), Slice(CNOT[1,3]), Slice(S[2])]
	circ   = Circuit([H[1], CNOT[1,3]])
	push!(circ, S[2], when=:last)
	@test circ.slices == slices 

	slices = [Slice(H[1]*S[2]), Slice(CNOT[1,3])]
	circ   = Circuit([H[1], CNOT[1,3]])
	push!(circ, S[2], when=:last)
	@test circ.slices == slices 

end
	

