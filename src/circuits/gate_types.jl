using LaTeXStrings

abstract type AbstractGate end

struct Pauli
    name::Symbol
    prod::ComplexF64
end

const X = Pauli(:X, 1.0)
const Y = Pauli(:Y, 1.0)
const Z = Pauli(:Z, 1.0)
const Id = Pauli(:I, 1.0)

function Base.:*(a::Pauli, b::Pauli)
    if (a.name == :X) && (b.name == :Y)
        return Pauli(:Z, 1im*a.prod*b.prod)
    elseif (a.name == :Y) && (b.name == :X)
        return Pauli(:Z, -1im*a.prod*b.prod)
    elseif (a.name == :Y) && (b.name == :Z)
        return Pauli(:X, 1im*a.prod*b.prod)
    elseif (a.name == :Z) && (b.name == :Y)
        return Pauli(:X, -1im*a.prod*b.prod)
    elseif (a.name == :Z) && (b.name == :X)
        return Pauli(:Y,1im*a.prod*b.prod)
    elseif (a.name == :X) && (b.name == :Z)
        return Pauli(:Y, -1im*a.prod*b.prod)
    elseif (a.name == :X) && (b.name == :X)
        return Pauli(:I, a.prod*b.prod)
    elseif (a.name == :Y) && (b.name == :Y)
        return Pauli(:I, a.prod*b.prod)
    elseif (a.name == :Z) && (b.name == :Z)
        return Pauli(:I, a.prod*b.prod)
    else
        println("????")
    end
end

Base.adjoint(a::Pauli) = Pauli(a.name, conj(a.prod))

Base.:*(a::Pauli, b::T) where T<:Number = Pauli(a.name, a.prod*b)
Base.:*(b::T, a::Pauli) where T<:Number = Pauli(a.name, a.prod*b)

function Base.:+(a::Pauli, b::Pauli)
    if a.name == b.name
        return Pauli(a.name, a.prod+b.prod)
    end
end

function Base.:-(a::Pauli, b::Pauli)
    if a.name == b.name
        return Pauli(a.name, a.prod-b.prod)
    end
end

comm(a::Pauli, b::Pauli) = a*b - b*a
anticomm(a::Pauli, b::Pauli) = a*b + b*a
Base.iszero(x::Pauli) = iszero(x.prod)

commute(A::Pauli, B::Pauli) = iszero(comm(A,B))


struct QubitId
    id::Tuple{Vararg{Int}}
end

Base.getindex(Q::QubitId, i) = Q.id[i]
Base.firstindex(Q::QubitId) = Q.id[1]
Base.lastindex(Q::QubitId)  = Q.id[end]
Base.iterate(Q::QubitId) = Base.iterate(Q.id)
Base.iterate(Q::QubitId, state) = Base.iterate(Q.id, state)
Base.eltype(::Type{QubitId}) = Int
Base.length(Q::QubitId) = length(Q.id)

QubitId(x::Int...) = QubitId(x)
QubitId(x::Set{Int}) = QubitId(x...)

target(id::QubitId) = id.id[end]
controls(id::QubitId) = id.id[1:end-1]
is_single(id::QubitId) = length(id.id) == 1
span(id::QubitId) = is_single(id) ? id[1] : minimum(id):maximum(id)

struct PauliGate{T<:AbstractAngle} <: AbstractGate
    pauli::Pauli
    angle::T
end

Base.adjoint(g::PauliGate) = PauliGate(adjoint(g.pauli), -1.0*g.angle)

const S = PauliGate(Z, Angle(1//4))
const T = PauliGate(Z, Angle(1//8))    

Rz(θ::T) where T<:AbstractAngle = PauliGate(Z, θ)
Rx(θ::T) where T<:AbstractAngle = PauliGate(X, θ)
Ry(θ::T) where T<:AbstractAngle = PauliGate(Y, θ)


function _to_latex_raw(p::Pauli)
    return string("\\gate{", p.name, "}")
end

function _to_latex_raw(g::PauliGate; st=true)
    if st && g == S
        return "\\gate{S}"
    end
    if st && g == T
        return "\\gate{T}"
    end
    return string("\\gate{", g.pauli.name, _to_latex_raw(g.angle), "}")
end

struct ControlledGate <: AbstractGate
    ctrl_pauli::Vector{Pauli}
    tgt_pauli::Pauli
end

ControlledGate(ctrl::Pauli, tgt::Pauli) where T<:AbstractAngle = ControlledGate([ctrl], tgt)

const CNOT = ControlledGate(Z, X)
const XZ   = ControlledGate(X, Z)

struct PauliMultiGate{T <: AbstractAngle} <:AbstractGate
    paulis::Vector{Pauli}
    angle::T
end

function _check_multiqb(g::T, id::QubitId) where T<:AbstractGate
    if (g isa ControlledGate)
        @assert length(id) > 1 "Multi-qubit gates must have more than one qubit ID."
        @assert length(id) == length(g.ctrl_pauli) + 1 "Incorrect number of qubit IDs."
    elseif (g isa PauliMultiGate)
        @assert length(id) > 1 "Multi-qubit gates must have more than one qubit ID."
        @assert length(id) == length(g.paulis) "Incorrect number of qubit IDs."
    end
end
        
    
struct Measure <: AbstractGate
    pauli::Pauli
end

Measure() = Measure(Z)

Base.getindex(A::G, i::Int...) where G<:AbstractGate = QubitId(i)=>A