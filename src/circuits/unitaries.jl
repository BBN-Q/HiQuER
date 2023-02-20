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
using LinearAlgebra 
using SparseArrays
using Folds

const cxu = Matrix{ComplexF64}([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]) |>sparse

@doc raw"""
Compute the trace norm of a matrix `A` given by:

   ``\Tr \sqrt{A^\dagger A}``

From: https://github.com/pyGSTio/pyGSTi/blob/master/pygsti/tools/optools.py
"""
function tracenorm(A::AbstractMatrix; tol=sqrt(eps()))
   if norm(A - A') < tol
      #Hermitian, so just sum eigenvalue magnitudes
      return sum(abs.(eigvals(A)))
   else
      #Sum of singular values, which is positive by construction
      F = svd(A)
      return sum(F.S)
   end 
end 

function isunitary(A::AbstractMatrix; tol=sqrt(eps()))
   return tracenorm(A'A - LinearAlgebra.I) < tol
end

@doc raw"""
Compute the trace distance between matrices A and B given by:
  ``0.5 * Tr(\sqrt{(A-B)^\dagger * (A-B)})``
From: https://github.com/pyGSTio/pyGSTi/blob/master/pygsti/tools/optools.py
"""
function tracedist(A::AbstractMatrix, B::AbstractMatrix)
   return 0.5*tracenorm(A - B)
end 

function is_close(A::AbstractMatrix, B::AbstractMatrix; tol=sqrt(eps()))
   size(A) != size(B) && return false

   ϕa = angle(A[findfirst(abs.(A[:]) .> tol)])
   ϕb = angle(B[findfirst(abs.(B[:]) .> tol)])
   return (tracedist(A*exp(-1im*ϕa),B*exp(-1im*ϕb)) < tol)
end

function average_fidelity(A::AbstractMatrix, B::AbstractMatrix)
   @assert size(A) == size(B)
   d = size(A)[1]
   return (d + abs(tr(A*B'))^2)/(d*(d+1))
end

"""
Generate a Haar-random unitary matrix of dimension d.
Algorithm from:
  F. Medrazzi. "How to generate random matrices from the classical
  compact groups" arXiv: math-ph/0609050
"""
function haar_random_unitary(d::Int)
   @assert d > 1 "Dimension must be > 1"
   re_z = randn((d,d))
   im_z = randn((d,d))
   Z = @. (re_z + im_z)/sqrt(2.0)
   F = qr(Z)
   L = diagm(diag(F.R) ./ abs.(diag(F.R)))
   return F.Q * L * F.Q
end

struct UnitaryGate <: AbstractGate
   U::SparseMatrixCSC
   function UnitaryGate(U::AbstractMatrix; tol=sqrt(eps()))
      Usp = sparse(U)
      Usp_re = real.(Usp)
      Usp_im = imag.(Usp)
      droptol!(Usp_re, tol)
      droptol!(Usp_im, tol)
      return new(Usp_re + 1im*Usp_im)
   end
end 


Matrix(G::UnitaryGate) = Matrix(G.U)

Base.:(*)(A::UnitaryGate, B::UnitaryGate) = UnitaryGate(A.U*B.U)
Base.:(*)(A::UnitaryGate, B::AbstractMatrix) = UnitaryGate(A.U*B)
Base.:(*)(A::AbstractMatrix, B::UnitaryGate) = UnitaryGate(A*B.U)
Base.:(*)(A::UnitaryGate, x::T) where {T <: Number} = UnitaryGate(A.U*x)
Base.:(*)(x::T, A::UnitaryGate) where {T <: Number} = UnitaryGate(x*A.U)
Base.:(+)(A::UnitaryGate, B::UnitaryGate) = UnitaryGate(A.U + B.U)
Base.:(-)(A::UnitaryGate, B::UnitaryGate) = UnitaryGate(A.U - B.U)

function LinearAlgebra.exp(A::UnitaryGate)
   size(A.U)[1] > 128 && @warn "Taking exponential of large matrix..."
   return UnitaryGate(LinearAlgebra.exp(Matrix(A.U)))
end

LinearAlgebra.kron(A::UnitaryGate, B::UnitaryGate) = UnitaryGate(kron(A.U, B.U))
⊗(A::UnitaryGate, B::UnitaryGate) = kron(A,B)

Base.adjoint(A::UnitaryGate) = UnitaryGate(A.U')
LinearAlgebra.tr(A::UnitaryGate) = tr(UnitaryGate)

dim(A::UnitaryGate) = size(A.U)[1]

tracenorm(A::UnitaryGate) = tracenorm(A.U)
tracedist(A::UnitaryGate, B::UnitaryGate) = tracedist(A.U, B.U)
average_fidelity(A::UnitaryGate, B::UnitaryGate) = average_fidelity(A.U, B.U)

Base.:(≈)(A::UnitaryGate, B::UnitaryGate) = is_close(A.U, B.U)

function UnitaryGate(g::Gate; U::Union{Nothing, AbstractMatrix}=nothing)
   if !isnothing(U)
      return UnitaryGate(U)
   end 

   g == H && return UnitaryGate([1.0 1.0; 1.0 -1.0]./sqrt(2))

   throw("not implemented!")
end  

function UnitaryGate(g::ControlledGate)
   g == cx && return UnitaryGate(cxu)
   throw("not implemented!")
end 

Ugate(θ, ϕ, λ) = UnitaryGate([cos(θ/2) -exp(im*λ)*sin(θ/2);
                              exp(im*ϕ)*sin(θ/2) exp(1im*(ϕ+λ))*cos(θ/2)])

function UnitaryGate(g::PauliGate)

   g == S  && return Ugate(0, 0, π/2)
   g == S' && return Ugate(0, 0, -π/2)
   g == T  && return Ugate(0, 0, π/4)
   g == T' && return Ugate(0, 0, -π/4)

   g.pauli == X && return Ugate(g.angle.value*pi, -π/2, π/2)
   g.pauli == Y && return Ugate(g.angle.value*pi, 0, 0)
   g.pauli == Z && return UnitaryGate([exp(-im*g.angle.value*pi/2) 0; 
                                       0 exp(im*g.angle.value*pi/2)])
end

const Ui = UnitaryGate([1 0; 0 1])

unitary(g::AbstractGate) = UnitaryGate(g)

function build_cnot(control::Int, target::Int, qubits::Vector{Int})
   #Write a CNOT in the form
   # |0><0| ⊗ I + |1><1> ⊗ X
   function mat00(i)
      i == control && return sparse([1 0; 0 0])
      return sparse([1 0; 0 1])
   end
   function mat11(i)
      i == control && return sparse([0 0; 0 1])
      i == target  && return sparse([0 1; 1 0])  
      return sparse([1 0; 0 1])
   end
   U00 = kron([mat00(i) for i in qubits]...)
   U11 = kron([mat11(i) for i in qubits]...)
   return U00+U11
end

function UnitaryGate(s::Slice, qubits::Vector{Int})

   N = length(qubits)
   sg = single_gates(s)
   if length(sg) > 0
      U = kron([QubitId(qb) in keys(sg) ? UnitaryGate(sg[QubitId(qb)]).U : Ui.U for qb in qubits]...)
   else
      U = spdiagm(ones(2^length(qubits)))
   end

   cg = controlled_gates(s)
   if length(cg) > 0
      for (ids, g) in cg 
         U = U*build_cnot(ids[1], ids[2], qubits)
      end
   end 
   return UnitaryGate(U)
end 

function UnitaryGate(c::Circuit)
   qubits = sort(collect(HiQuER.qubits(c)))

   return Folds.reduce(*, (UnitaryGate(c.slices[idx], qubits) for idx in length(c.slices):-1:1), ThreadedEx())
end

unitary(c::Circuit) = UnitaryGate(c)

export UnitaryGate, Ui, unitary


 
   




