module HiQuER

include("graphs/graph_coloring.jl")

include("circuits/angle_types.jl")
include("circuits/gate_types.jl")
include("circuits/circuit_types.jl")

include("io/qasm_util.jl")
include("io/qasm_circuit.jl")

include("algorithms/pauli_hamiltonian.jl")

end
