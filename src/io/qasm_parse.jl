#=
Copyright 2023 Raytheon BBN

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

Basic OpenQASM 2.0 parser.

=#


using OpenQASM
using Match
using Base.Iterators: flatten

function load_qasm(path)
    open(path, "r") do io
        qasm = read(io, String)
        return OpenQASM.parse(qasm)
    end 
end

function parse_qasm(reg::OpenQASM.Types.RegDecl)
    name = reg.name.str
    width = parse(Int, reg.size.str)
    return (HiQuER.QRegister(name, width),)
end

function parse_qasm(bit::OpenQASM.Types.Bit)
    return (bit.name.str, 1+parse(Int, bit.address.str))
end 

#TODO: Fixme!
function parse_carg(carg)


    if carg isa OpenQASM.Types.Neg
        return -1.0*parse(Float64, carg.val.str)
    else
        return parse(Float64, carg.str)
    end
end
    

function parse_qasm(instr::OpenQASM.Types.Instruction)
    qargs = [parse_qasm(q) for q in instr.qargs]
    cargs = instr.cargs
    gate = @match instr.name begin
        
        "s"   => S
        "sdg" => S'
        "h"   => H
        "t"   => T
        "tdg" => T'
        "x"   => X180
        "y"   => Y180
        "z"   => Z180
        "cx"  => CNOT
        "rx"  => Rx(FloatAngle(parse_carg(cargs[1])))
        "ry"  => Ry(FloatAngle(parse_carg(cargs[1])))
        "rz"  => Rz(FloatAngle(parse_carg(cargs[1])))
        "u1"  => Rz(FloatAngle(parse_carg(cargs[1])/2))
        "u2"  => begin
                    θ1 = FloatAngle(parse_carg(cargs[2]) - π/2)
                    θ3 = FloatAngle(parse_carg(cargs[1]) - π/2)
                    [Rz(θ1), X90, Rz(θ3)]
                end
        "u3"  => begin 
                    θ1 = FloatAngle(parse_carg(cargs[3]))
                    θ2 = FloatAngle(parse_carg(cargs[2]) + π)
                    θ3 = FloatAngle(parse_carg(cargs[1]) + 3π)
                    [Rz(θ1), X90, Rz(θ2), X90, Rz(θ3)]
                end
        _     => error("Unknown gate $instr.name")
    end

    
    return (gate, qargs)
end

function parse_qasm(instr::OpenQASM.Types.Measure)
    qargs = [parse_qasm(instr.qarg)]
    return (MEAS, qargs)
end

function parse_qasm(include::OpenQASM.Types.Include)
    return nothing
end

function parse_qasm_circuit(circuit::String)

    ast = OpenQASM.parse(circuit)
    prog = ast.prog
    
    c = Circuit()
    
    for (idx, instr) in enumerate(prog)
    	#println(idx, " ", instr)

        if !isnothing(local m = parse_qasm(instr))
            push!(c, m...)
        end
    end
    return c
end 

function load_qasm_circuit(path)
    open(path, "r") do io
        qasm = read(io, String)
        return parse_qasm_circuit(qasm)
    end 
end

export load_qasm_circuit, parse_qasm_circuit



