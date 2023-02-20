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

using QuantumClifford
using LaTeXStrings

abstract type AbstractGate end


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

function Base.show(io::IO, Q::QubitId) 
    if is_single(Q)
        return print(io, string("(", Q.id[1], ")"))
    else
        return print(io, string(Q.id))
    end 
end

target(id::QubitId) = id.id[end]
controls(id::QubitId) = id.id[1:end-1]
is_single(id::QubitId) = length(id.id) == 1
span(id::QubitId) = is_single(id) ? id[1] : minimum(id):maximum(id)

struct Gate <: AbstractGate
    name::Symbol
    dag::Bool
    self_adjoint::Bool
end

Gate(name::Symbol) = Gate(name, false, false)

const H = Gate(:H, false, true)

Base.adjoint(g::Gate) = Gate(g.name, g.self_adjoint ? g.dag : !g.dag, g.self_adjoint)

struct PauliGate{T<:AbstractAngle} <: AbstractGate
    pauli::PauliOperator
    angle::T
end

PauliGate(pauli::PauliOperator, θ::Rational) = PauliGate(pauli, Angle(θ))
PauliGate(pauli::PauliOperator, θ::AbstractFloat) = PauliGate(pauli, FloatAngle(θ))

Base.adjoint(g::PauliGate) = PauliGate(g.pauli, -1.0*g.angle)

const S = PauliGate(Z, Angle(1//4))
const T = PauliGate(Z, Angle(1//8))    

Rz(θ::T) where T<:AbstractAngle = PauliGate(Z, θ)
Rx(θ::T) where T<:AbstractAngle = PauliGate(X, θ)
Ry(θ::T) where T<:AbstractAngle = PauliGate(Y, θ)

const X180 = Rx(Angle(1))
const Y180 = Ry(Angle(1))
const Z180 = Rz(Angle(1))

const X90 = Rx(Angle(1//2))
const Y90 = Ry(Angle(1//2))
const Z90 = Rz(Angle(1//2))

function Base.show(io::IO, g::Gate)
    return print(io, string(g.name, g.dag ? "†" : ""))
end

function Base.show(io::IO, g::PauliGate)
    g == T  && return print(io, "T")
    g == T' && return print(io, "T†")
    g == S  && return print(io, "S")
    g == S' && return print(io, "S†")
    g == X180  && return print(io, "X")
    g == Y180  && return print(io, "Y")
    g == Z180  && return print(io, "Z")

    g.pauli == X && return print(io, string("Rx(", string(g.angle), ")"))
    g.pauli == Y && return print(io, string("Ry(", string(g.angle), ")"))
    g.pauli == Z && return print(io, string("Rz(", string(g.angle), ")"))
end


const _pauli_name = Dict(P"Z" => "Z", P"X" => "X", P"Y" => "Y", P"I" => "I")

function _to_latex_raw(g::PauliGate)
    if g == S
        return "\\gate{S}"
    end
    if g == T
        return "\\gate{T}"
    end
    return string("\\gate{", _pauli_name[g.pauli], _to_latex_raw(g.angle), "}")
end

struct ControlledGate <: AbstractGate
    ctrl_pauli::Vector{PauliOperator}
    tgt_pauli::PauliOperator
end

ControlledGate(ctrl::PauliOperator, tgt::PauliOperator) where T<:AbstractAngle = ControlledGate([ctrl], tgt)

##TODO:FIXME! 
##Only true for CNOT...
Base.adjoint(cg::ControlledGate) = cg

const CNOT = ControlledGate(Z, X)
const XZ   = ControlledGate(X, Z)

function Base.show(io::IO, g::ControlledGate)
    g == CNOT && return print(io, "CNOT")
end

struct PauliMultiGate{T <: AbstractAngle} <:AbstractGate
    paulis::Vector{PauliOperator}
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
    pauli::PauliOperator
end

Measure() = Measure(Z)

const MEAS = Measure(Z)

Base.getindex(A::G, i::Int...) where G<:AbstractGate = QubitId(i)=>A

export Gate, PauliGate, ControlledGate, Measure 
export H, T, S, Rx, Rz, Ry, X180, Y180, Z180, X90, Y90, Z90, CNOT, MEAS