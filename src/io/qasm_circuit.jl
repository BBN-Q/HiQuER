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

function gate_name(g::Gate)
    return lowercase(string(g.name, g.dag ? "dg" : ""))
end

function gate_name(g::PauliGate)
    
    g == S && return "s"
    g == S' && return "sdg"
    g == T && return "t"
    g == T' && return "tdg"
    g == X180 && return "x"
    g == Y180 && return "y"
    g == Z180 && return "z"

    if g.pauli == X
        return string("rx(", g.angle, ")")
    elseif g.pauli == Y
        return string("ry(", g.angle, ")")
    elseif g.pauli == Z
        return string("rz(", g.angle, ")")
    end
end

function gate_name(g::ControlledGate)
    g == CNOT && return "cx"
    throw(ArgumentError("Unknown gate."))
end

function QASMGateDef(c::Circuit, name::Symbol; prefix="q")
	 sorted_qb = sort(collect(qubits(c)))
    input_names = [Symbol(string(prefix, idx)) for idx in sorted_qb]
    inputs = Dict(idx=>QASMQubit(n, -1) for (idx,n) in zip(sorted_qb, input_names))
    
    instr = Vector{QASMGate}()
    
    for sl in c.slices
        for (id, g) in sl.gates
            targets = [inputs[idx] for idx = id]
            push!(instr, QASMGate(Symbol(gate_name(g)), targets))
        end
    end
    
    return QASMGateDef(name, input_names, instr)

end


function to_qasm(c::Circuit; gatedefs=[])
    
    regs = [QASMRegister(:q, length(qubits(c)))]
    gates = Vector{QASMGate}()
    for sl in c.slices
        for (id, g) in sl.gates
            targets = [QASMQubit(regs[1].name, idx-1) for idx in id]
            push!(gates, QASMGate(Symbol(gate_name(g)), targets))
        end
    end
    return QASMListing([], gatedefs, regs, [], gates)
end  

export QASMGateDef
export to_qasm
