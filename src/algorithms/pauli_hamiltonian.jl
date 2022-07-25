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

pauli_hamiltonian.jl

Tools for converting an arbitrary Hamiltonian to a Pauli-product form and diagonalizing
   sets of commuting Pauli operators so that it can be "efficiently" represented as a circuit.

References:
   
   [1] van den Berg, E. and Temme, K. "Circuit optimization of Hamiltonian simulation by
      simultaneous diagonalization of Pauli clusters" Quantum 4, 322 (2020). 
      https://doi.org/10.22331/q-2020-09-12-322

   [2] Aaronson, S. and Gottesman, D. "Improved simulation of stabilizer circuits",
      Phys. Rev. A 70, 052328 (2004) 
      https://doi.org/10.1103/PhysRevA.70.052328

   [3] Crawford, O., van Straaten, B., Wang, D., Parks, T., Campbell, E. and Brierley, S. "Efficient
      quantum measurement of Pauli operators in the presence of finite sampling error". Quantum 5,
      385 (2021). 
      https://doi.org/10.22331/q-2021-01-20-385

=#

using QuantumClifford
using Base.Iterators: product
using Combinatorics, SparseArrays

################################################################################################
# Pauli Operator Convenience Functions

const _pauli_phase = Dict{UInt8, ComplexF64}(0x0 => 1.0, 
                                             0x1 => 1im,
                                             0x2 => -1.0,
                                             0x3 => -1im)

function convert(::Type{Matrix{Complex{T}}}, p::PauliOperator) where {T}
   matX = sparse(Complex{T}[0 1; 1 0])
   matZ = sparse(Complex{T}[1 0; 0 -1])
   return reduce(kron, [(1im)^(x*z)*(matX^x)*(matZ^z) for (x,z) in zip(xbit(p), zbit(p))]);
end

complex(p::PauliOperator) = convert(Matrix{ComplexF64}, p)

"""
Enumerate all n-qubit Pauli operators.
"""
function enumerate_paulis(n)
   pgroup = [P"I", P"X", P"Y", P"Z"]
   return (reduce(⊗, x) for x in product(repeat([pgroup], n)...))
end

"""
Represent a Matrix `A` as a sum of weighted Pauli operators.
Returns a dictionary that maps paulis to their corresponding weight.

"""
function mat_to_paulis(A)
   (N,M) = size(A)
   @assert N == M "Matrix must be square!"
   n = round(Int, log2(N))
   @assert n ≈ log2(N) "Matrix dimension must be a power of 2!"
   d = Dict{PauliOperator, Float64}()
   for p in enumerate_paulis(n)
      h = real((1.0/M)*tr(complex(p)*A))
      if !(h ≈ 0)
         d[p] = h
      end 
   end 
   return d 
end

"""
Parition Pauli strings into mutually commuting subsets. Uses the greedy algorithm of [3]
"""

function greedy_partition(d::Dict{T, Float64}) where T<:PauliOperator
    sets = [[]]
    for (pauli, θ) in sort(collect(d), by=x->x[2])
        inserted = false
        for idx in eachindex(sets)
            if all(QuantumClifford.comm(pauli, p) == 0x0 for p in sets[idx])
                push!(sets[idx], pauli)
                inserted = true
                break
            end
        end
        if !inserted
            push!(sets, [pauli])
        end
    end
    return sets
end    

################################################################################################
# Tableau Representations
#
# TODO: We should just be able to use the built-in Stabilizer tableau

mutable struct Tableau
    X::BitArray
    Z::BitArray
    s::BitVector
end

function Tableau(s::Stabilizer)
    Np = length(s)
    Nq = s.nqubits
    X  = stab_to_gf2(s)[:,1:Nq]
    Z  = stab_to_gf2(s)[:,Nq+1:end]
    s  = [ph > 0x1 for ph in s.phases]
    return Tableau(X,Z,s)
end

LinearAlgebra.rank(t::Tableau) = rank([t.X t.Z])

Base.copy(t::Tableau) = Tableau(copy(t.X), copy(t.Z), copy(t.s))

function QuantumClifford.Stabilizer(t::Tableau)
    return Stabilizer(Vector{UInt8}(t.s .<< 1),
                        Matrix{Bool}(t.X),
                        Matrix{Bool}(t.Z))
end

################################################################################################
# Tableau Operations
#
# TODO: We should just be able to use the built-in Stabilizer tableau

"""
Swap two rows `i` and `j` in the tableau.
"""
function rowswap!(t::Tableau, i, j)
    (i == j) && return
    t.s[i], t.s[j] = t.s[j], t.s[i]
    t.X[i,:], t.X[j,:] = t.X[j,:], t.X[i,:]
    t.Z[i,:], t.Z[j,:] = t.Z[j,:], t.Z[i,:]
    nothing
end

"""
Swap two column `i` and `j` in the tableau
"""
function colswap!(t::Tableau, i, j)
    (i==j) && return
    t.X[:,i], t.X[:,j] = t.X[:,j], t.X[:,i]
    t.Z[:,i], t.Z[:,j] = t.Z[:,j], t.Z[:,i]
    nothing
end

function bitg(x1::T, x2::T, z1::T, z2::T) where T<:Bool
    if !x1 & !z1
        return 0
    elseif x1 & z1
        return z2 - x2
    elseif x1 & !z1
        return z2*(2*x2-1)
    elseif !x1 & z1
        return x2*(1-2*z2)
    end
end

"""
Perform the rowsum algorithm of [2] on a tableau. This is what [1] refers to as a "sweep".      
""" 
function rowsum!(t::Tableau, h, i)
    (i==h) && return
    phase_sum = 2*t.s[h] + 2*t.s[i]
    phase_sum += sum([bitg(t.X[i,j], t.Z[i,j], t.X[h,j], t.Z[h,j]) for j=1:size(t.X,2)])
    if mod(phase_sum, 4) == 0
        t.s[h] = false
    elseif mod(phase_sum, 4) == 2
        t.s[h] = true
    else
        throw(DomainError(phase_sum, "Something went wrong!"))
    end
    @inbounds @simd for j = 1:size(t.X,2)
        t.X[h,j] = t.X[i,j] ⊻ t.X[h,j]
        t.Z[h,j] = t.Z[i,j] ⊻ t.Z[h,j]
    end
    nothing
end

################################################################################################
# Tableau Clifford Gates
#
# TODO: We should just be able to use the built-in Stabilizer tableau

"""
Apply a Hadamard gate to qubit `a`.
"""
function hadamard!(t::Tableau, a)
    t.X[:,a], t.Z[:,a] = t.Z[:,a], t.X[:,a]
    nothing
end

"""
Apply a phase (S) gate to qubit `a`.
"""
function phase!(t::Tableau, a)
    @inbounds @simd for j = 1:length(t.s)
        t.s[j] = t.s[j] ⊻ (t.X[j,a]*t.Z[j,a])
    end
    @inbounds @simd for j = 1:size(t.X, 1)
        t.Z[j,a] = t.Z[j,a] ⊻ t.X[j,a]
    end
    nothing
end

"""
Apply a CNOT(a,b) gate to qubits `a` (control) and `b` (target).
"""
function cnot!(t::Tableau, a, b)
    @inbounds @simd for j=1:length(t.s)
        t.s[j] = t.s[j] ⊻ (t.X[j,a]*t.X[j,b]*(t.X[j,b] ⊻ t.Z[j,a] ⊻ true))
    end
    @inbounds @simd for j=1:size(t.X, 1)
        t.X[j,b] = t.X[j,b] ⊻ t.X[j,a]
        t.Z[j,a] = t.Z[j,a] ⊻ t.Z[j,b]
    end
    nothing
end

"""
Apply a control-Z(a,b) gate to qubits `a` (control) and `b` (target).
"""
function cz!(t::Tableau, a, b)
    @inbounds @simd for j=1:length(t.s)
        t.s[j] = t.s[j] ⊻ (t.X[j,a]*t.X[j,b]*(t.Z[j,a] ⊻ t.Z[j,b] ⊻ true))
    end
    @inbounds @simd for j=1:size(t.X, 1)
        t.Z[j,a] = t.Z[j,a] ⊻ t.X[j,b]
        t.Z[j,b] = t.Z[j,b] ⊻ t.X[j,a]
    end
    nothing
end

################################################################################################
# Diagonalizaton Algorithms

"""
Diagonalize x block of tableau. Algorithm 1 in [1].

"""
function diagX!(t::Tableau, tg::Tableau, c::Circuit)
    m, n = size(t.X)
    k = 1
    while true
        idx = nothing
        for ij in CartesianIndices((k:m, k:n))
            if t.X[ij] == 1
                idx = ij
                break
            end
        end
        if isnothing(idx)
            break
        end
        rowswap!(t, idx[1], k)
        colswap!(t, idx[2], k)
        for i = 1:m
            if (t.X[i,k] == 1) && (i != k)
                rowsum!(t, i, k)
            end
        end
        k = k + 1
    end
    kx = k
    while true
        idx = nothing
        for ij in CartesianIndices((k:m, k:n))
            if t.Z[ij] == 1
                idx = ij
            end
        end
        if isnothing(idx)
            break
        end
        rowswap!(t, idx[1], k)
        colswap!(t, idx[2], k)
        for i = 1:m
            if (t.Z[i,k] == 1) && (i !=k)
                rowsum!(t, i, k)
            end
        end
        k = k + 1
    end
    for j=kx:k-1
        hadamard!(t, j)
        hadamard!(tg, j)
        push!(c, H[j])
    end
    for i = 1:k
        for j=k:n
            if t.X[i,j] == 1
                cnot!(t, i, j)
                cnot!(tg, i, j)
                push!(c, CNOT[i,j])
            end
        end
    end
    nothing
end

"""
Diagonalize Z block and clear X block of tableau. Corresponds to Algorithm 3 of [1].
"""
function clearX!(t::Tableau, tg::Tableau, c::Circuit)
    k = rank(t)
    for i = 1:k
        if mod(sum([t.Z[i,j] for j=1:i]), 2) == 0
            phase!(t, i)
            phase!(tg, i)
            push!(c, S[i])
        end
        for j=1:i-1
            if t.Z[i,j] == 1
                cnot!(t, i, j)
                cnot!(tg, i, j)
                rowsum!(t, i, j)
                push!(c, CNOT[i,j])
            end
        end
    end
    for i=1:k
        phase!(t, i)
        phase!(tg, i)
        push!(c, S[i])
        hadamard!(t, i)
        hadamard!(tg, i)
        push!(c, H[i])
    end
    nothing
end    

"""
Recursively partition a bit-matrix M such that the number of 0->1 and 1->0 transitions are 
minimized in each row.
"""
function partition!(M, k, range, rev=false)
    K = sortperm(M[k, range], rev=rev)
    if length(K) <= 1
        return
    end
    
    M[:,range] = M[:,range[K]]
    
    idx = findfirst(M[k,range] .== 1)
    if !isnothing(idx) && (k+1 < size(M,1))
        if idx > 1
            partition!(M, k+1, range[1]:range[idx-1])
        end
        partition!(M, k+1, range[idx]:range[end], true)
    end
end

partition!(M) = partition!(M, 1, 1:size(M,2), false)

function phase_gadgets(t::Tableau, θ::Vector{T}; opt=true) where T <: Number
   (N,M) = size(t.Z)

   @assert length(θ) == N "Must have as many angles as Pauli strings."

   if opt 
      K = copy([t.Z t.s])' #is the copy needed?? 
      partition!(K)
      t.Z = K'[:, 1:end-1]
      t.s = K'[:, end]
   end 

   signs = [s == 0 ? 1.0 : -1.0 for s in t.s]

   c = Circuit() 

   for row in eachrow(t.Z)
      rz_placed = false 
      cnots = []
      for j = 1:M-1
         k = findfirst(row[j+1:end] .== 1)
         if row[j] == 1 && !isnothing(k)
            push!(c, CNOT[j,k+j])
            push!(cnots, (j,k+j))
         end 
         if row[j] == 1 && isnothing(k)
            push!(c, Rz(FloatAngle(-2*signs[j]*θ[j]))[j])
            rz_placed = true
         end 
      end 
      if row[M] == 1 && !rz_placed
         push!(c, Rz(FloatAngle(-2*signs[M]*θ[M]))[M])
      end 
      #uncompute 
      for (j,k) in cnots 
         push!(c, CNOT[j,k])
      end

   end 
   return c 
end





export complex
export Tableau
export enumerate_paulis
export greedy_partition 
export mat_to_paulis 
export partition! 