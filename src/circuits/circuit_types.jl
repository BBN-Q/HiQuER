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

using LaTeXStrings
using Base.Iterators: flatten

####################################################################################################
# Circuits
#

mutable struct Slice
    gates::Dict{QubitId, T} where {T<:AbstractGate}
    qubits::Set{Int}
    time_idx::Int
end

function _get_qb_ids(gates::Dict{QubitId, T}) where {T<:AbstractGate}
    return Set([n for k in keys(gates) for n in k.id])
end

function Slice(gates)
    s = Slice()
    for (k,v) in gates
        push!(s, v, k)
    end
    return s
end

function Slice(p::Pair{QubitId, T}) where T <: AbstractGate 
    s = Slice()
    push!(s, p)
    return s 
end


function Base.show(io::IO, s::Slice)
    str = [string(g, k)  for (k,g) in s.gates]
    return print(io, join(str, "*"))
end


function Base.:(*)(A::Pair{QubitId,T1}, B::Pair{QubitId,T2}) where {T1<:AbstractGate, T2<:AbstractGate}
    return Slice([A, B])
end

function Base.:(*)(A::Pair{Int,T1}, B::Pair{Int,T2}) where {T1<:AbstractGate, T2<:AbstractGate}
    _check_multiqb(A[2], QubitId(A[1]))
    _check_multiqb(B[2], QubitId(B[1]))
    return Slice([A, B])
end

Slice() = Slice(Dict{QubitId, AbstractGate}(), Set(), -1)

Base.adjoint(s::Slice) = Slice(Dict(k=>v' for (k,v) in s.gates), s.qubits, s.time_idx)

#For equality, we don't compare the time index 
function Base.:(==)(x::Slice, y::Slice)
    return x.gates == y.gates && x.qubits == y.qubits 
end



####################################################################################################
# Circuits
#

mutable struct Circuit
    slices::Vector{Slice}
    wires::Vector{Int}
end

qubits(c::Circuit) = Set(flatten([s.qubits for s in c.slices]))

Circuit() = Circuit([], [])

function available(s::Slice, qb_idx::QubitId)
    return !any(id in s.qubits for id in qb_idx)
end
available(s::Slice, qb_idx::Set{Int}) = !any(id in s.qubits for id in qb_idx)

single_gates(s::Slice) = Dict(k=>v for (k,v) in s.gates if length(k) == 1)
controlled_gates(s::Slice) = Dict(k=>v for (k,v) in s.gates if length(k) > 1)

function Base.push!(s::Slice, g::T, qb_idx::QubitId) where T<:AbstractGate
    _check_multiqb(g, qb_idx)
    if available(s, qb_idx)
        s.gates[qb_idx] = g
        push!(s.qubits, qb_idx...)
    else
        throw(ArgumentError("Cannot add to this slice!"))
    end
end

function Base.push!(s::Slice, g::Pair{QubitId,T}) where {T<:AbstractGate} 
    push!(s, g[2], g[1])
end

Base.push!(s::Slice, g::T, qb_idx::Int) where T<:AbstractGate = push!(s, g, QubitId(qb_idx))

function Base.:(*)(s::Slice, g::Pair{QubitId,T}) where {T<:AbstractGate}
    push!(s, g[2], g[1])
    return s
end

Base.:(*)(g::Pair{QubitId,T}, s::Slice) where {T<:AbstractGate} = s*g

Base.getindex(c::Circuit, i::Int) = c.slices[i]
Base.setindex!(c::Circuit, s::Slice, i::Int) = setindex!(c.slices, s, i)
Base.firstindex(c::Circuit) = c.slices[1]
Base.lastindex(c::Circuit)  = c.slices[end]


function Base.push!(c::Circuit, s::Slice)
    
    idx = findfirst(!available(sl, QubitId(s.qubits)) for sl in c.slices[end:-1:1])
    if idx == 1 || isnothing(idx)
        s.time_idx = length(c.slices)+1
        push!(c.slices, s)
        return nothing
    end 

    idx = length(c.slices) - idx + 1

    s.time_idx = idx+1
    insert!(c.slices, idx+1, s)
    for j = idx+1:length(c.slices)
        c.slices[idx].time_idx += 1
    end 
    return nothing
end


function Base.push!(c::Circuit, sg::Vector)
    for x in sg
        if x isa Slice
            push!(c, x)
        else
            push!(c, x; when=:earliest)
        end
    end
end

function Circuit(sg::Vector)
    c = Circuit()
    push!(c, sg)
    return c
end

Base.push!(x::Circuit, y::Circuit) = push!(x, y.slices)

function Base.insert!(x::Circuit, index::Integer, y::Circuit)
    for (jj, slice) in y.slices
        insert!(x.slices, index+jj-1, slice)
    end 
    nothing 
end

function Base.push!(c::Circuit, g::T, qb_idx::QubitId; when=:earliest) where T<:AbstractGate
    if when == :earliest
        
        idx = findfirst(!available(sl, qb_idx) for sl in reverse(c.slices))

        #the last slice has a gate on this wire
        if idx == 1 || length(c.slices) == 0 
            push!(c.slices, Slice(Dict{QubitId, AbstractGate}(qb_idx=>g)))
            c.slices[end].time_idx = length(c.slices)
            return nothing
        end 

        if isnothing(idx)
            push!(c.slices[1], g, qb_idx)
            return nothing
        end

        #otherwise put the gate in the first available slot
        push!(c.slices[end-idx+2], g, qb_idx)
        nothing
    end
    
    if when == :last
        push!(c.slices, Slice(Dict{QubitId, AbstractGate}(qb_idx=>g)))
        c.slices[end].time_idx = length(c.slices)
    end
end

function Base.push!(c::Circuit, g::T, qb_idx::Int; when=:earliest) where T<:AbstractGate 
    push!(c, g, QubitId(qb_idx); when=when)
    nothing
end

function Base.push!(c::Circuit, g::Pair{QubitId,T}; when=:earliest) where {T<:AbstractGate}
    push!(c, g[2], g[1]; when=when)
end

function Base.push!(c::Circuit, g::Pair{Int,T}; when=:earliest) where {T<:AbstractGate}
    push!(c, g[2], g[1]; when=when)
end

function Base.adjoint(c::Circuit)
    cdg = Circuit()
    for idx in length(c.slices):-1:1
        push!(cdg, c.slices[idx]')
    end
    return cdg
end

function Base.show(io::IO, c::Circuit)
    return print(io, join([string(s.time_idx, ": ", s) for s in c.slices], "\n"))
end

export Circuit
export Slice