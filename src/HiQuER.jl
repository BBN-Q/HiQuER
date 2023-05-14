module HiQuER

using Documenter

include("graphs/graph_coloring.jl")

include("circuits/angle_types.jl")
include("circuits/gate_types.jl")
include("circuits/circuit_types.jl")

include("circuits/unitaries.jl")

include("io/qasm_util.jl")
include("io/qasm_parse.jl")
include("io/qasm_circuit.jl")
include("io/circuit_viz.jl")

include("algorithms/pauli_hamiltonian.jl")
include("algorithms/trotterize.jl")

include("applications/fermihubbard.jl")
include("applications/plasmaphysics.jl")

end
