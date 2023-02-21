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

fermihubbard.jl

Generate Fermi-Hubbard model instances using OpenFermion.

=#

using LinearAlgebra
using QuantumClifford
using PyCall

of = pyimport("openfermion")

function term_to_pauli(N, term)
    xvec = zeros(Bool, 2*N)
    zvec = zeros(Bool, 2*N)
    for p in term
        if p[2] == "Z"
            zvec[p[1]+1] = true
        elseif p[2] == "X"
            xvec[p[1]+1] = true
        elseif p[2] == "Y"
            xvec[p[1]+1] = true
            zvec[p[1]+1] = true
        end
    end
    return PauliOperator(0x0, xvec, zvec)
end 

function terms_to_paulidict(N, terms)
    pdict = Dict()
    for (k, v) in terms
        pdict[term_to_pauli(N, k)] = real.(v) #hack!
    end
    return pdict
end

function gen_fermi_hubbard(N)
    U = 2.0
    J = -1.0
    hubbard = of.fermi_hubbard(1, N, tunneling=-J, coulomb=U, periodic=false)
    return of.jordan_wigner(hubbard).terms
end

export gen_fermi_hubbard