using LaTeXStrings
using Base.Iterators: flatten
using TikzPictures

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


function Base.:(*)(A::Pair{QubitId,T1}, B::Pair{QubitId,T2}) where {T1<:AbstractGate, T2<:AbstractGate}
    return Slice([A, B])
end

Slice() = Slice(Dict{QubitId, AbstractGate}(), Set(), -1)

mutable struct Circuit
    slices::Vector{Slice}
    wires::Vector{Int}
end

qubits(c::Circuit) = Set(flatten([s.qubits for s in c.slices]))

Circuit() = Circuit([], [])

function available(s::Slice, qb_idx::QubitId; check_overlap=false)
    if check_overlap
        all_qb = flatten([span(k) for k in keys(s.gates)])
        return !any(id in all_qb for id in span(qb_idx))
    else
        return !any(id in s.qubits for id in qb_idx)
    end
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

function remove_overlaps(c::Circuit)
    newc = Circuit()
    for s in c.slices
        for (k, v) in single_gates(s)
            push!(newc, v, k; when=:earliest, no_overlaps=true)
        end
        for (k, v) in controlled_gates(s)
            push!(newc, v, k; when=:earliest, no_overlaps=true)
        end
    end
    return newc
end

function Base.push!(c::Circuit, s::Slice; no_overlaps=false)
    @inbounds for idx in eachindex(c.slices)
        if available(c.slices[idx], QubitId(s.qubits), check_overlap=no_overlaps)
            s.time_idx = idx
            insert!(c.slices, idx, s)
            return nothing
        end
    end
    s.time_idx = length(c.slices)+1
    push!(c.slices, s)
    nothing
end    

function Base.push!(c::Circuit, sg::Vector)
    for x in sg
        if x isa Slice
            push!(c, x; no_overlaps=true)
        else
            push!(c, x; no_overlaps=true, when=:earliest)
        end
    end
end

function Circuit(sg::Vector)
    c = Circuit()
    push!(c, sg)
    return c
end

function Base.push!(c::Circuit, g::T, qb_idx::QubitId; when=:earliest, no_overlaps=true) where T<:AbstractGate
    if when == :earliest
        @inbounds for idx in eachindex(c.slices)
            if available(c.slices[idx], qb_idx, check_overlap=no_overlaps)
                push!(c.slices[idx], g, qb_idx)
                return nothing
            end
        end
        push!(c.slices, Slice(Dict{QubitId, AbstractGate}(qb_idx=>g)))
        c.slices[end].time_idx = length(c.slices)
    end
    
    if when == :last
        push!(c.slices, Slice(Dict{QubitId, AbstractGate}(qb_idx=>g)))
        c.slices[end].time_idx = length(c.slices)
    end
end

function Base.push!(c::Circuit, g::T, qb_idx::Int; when=:earliest, no_overlaps=false) where T<:AbstractGate 
    push!(c, g, QubitId(qb_idx); when=when, no_overlaps=no_overlaps)
    nothing
end

function Base.push!(c::Circuit, g::Pair{QubitId,T}; when=:earliest, no_overlaps=false) where {T<:AbstractGate}
    push!(c, g[2], g[1]; when=when, no_overlaps=no_overlaps)
end