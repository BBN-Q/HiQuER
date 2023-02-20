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

Tools for circuit I/O into OpenQASM 2.0 format

	TODO: Write a real parser.

=#

using Quantikz: U as Uqtz
using Quantikz: Id as Idqtz
using Quantikz: CNOT as CNOTqtz
using Quantikz: Measurement as Mqtz 
using Quantikz: displaycircuit, circuit2string, savecircuit

function gate_to_tikz(g::Gate) 
	name = string(g.name, g.dag ? raw"\dagger" : "")
	return  x -> Uqtz(string(g.name), x)
end

function gate_to_tikz(g::PauliGate)
	return x -> Uqtz(string(_pauli_name[g.pauli],"(", _to_latex_raw(g.angle), ")"), x)
end 

function gate_to_tikz(g::ControlledGate)
	if g == CNOT 
		return (x,y) -> CNOTqtz(x,y)
	else
		throw("unimplemented!")
	end 
end

function gate_to_tikz(g::Measure)
	return x -> Mqtz(_pauli_name[g.pauli], x)
end


function to_qtikz(c::Circuit)
	qbs = qubits(c)
	qtc = [] 

	for sl in c.slices 
		if isempty(sl.gates)
			continue 
		end
		for qb in setdiff(qbs, sl.qubits)
			push!(qtc, Idqtz(qb))
		end 
		for (qbs, g) in sl.gates 
			push!(qtc, gate_to_tikz(g)(qbs.id...))
		end 
	end 
	return qtc 
end 

to_tex(c::Circuit) = circuit2string(to_qtikz(c))

draw(c::Circuit) = displaycircuit(to_qtikz(c))


function save(c::Circuit, filename::AbstractString)
	qtc = to_qtikz(c)
	savecircuit(qtc, filename)
end

export draw, save, to_tex