#=
Copyright 2022 Raytheon BBN

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.

You may obtain a copy of the License at
   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#

#= 

trotterize.jl

Tools for Trotter-Suzuki approximations. 

=#

function product_formula_empirical_bounds(k, N)
    @assert k in (1, 2, 4, 6, 8) "Do not have empirical bound for $(k)!"
    if k == 1
        return 2417*N^1.964
    elseif k == 2
        return 39.47*N^1.883;
    elseif k == 4
        return 4.035*N^1.555
    elseif k == 6
        return 1.789*N^1.311
    elseif k == 8 
        return 1.144*N^1.141
    end
    return nothing
end

function Suzuki_Trotter_circuits_4thorder(pauli_dict, N, Tmax)
    q  = 2
    r  = product_formula_empirical_bounds(4, N)
    pq = 1/(4 - 4^(1/(2q-1)))
    
    outer = Dict(x => pq*y*Tmax/r for (x,y) in pauli_dict)
    inner = Dict(x => (1-4*pq)*y*Tmax/r for (x,y) in pauli_dict)
    
    outer_circuit = to_circuit(greedy_partition(outer))
    inner_circuit = to_circuit(greedy_partition(inner))

    circuit = Circuit()
    push!(circuit, outer_circuit)
    push!(circuit, inner_circuit)
    push!(circuit, inner_circuit)
    push!(circuit, outer_circuit)
    
    qasm = to_qasm(circuit)
    
    return circuit
end

export Suzuki_Trotter_circuits_4thorder
export product_formula_empirical_bounds