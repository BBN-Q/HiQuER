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

Tools for interacting with OpenQASM 2.0

	TODO: Write a real parser.

=#

struct ParseError <: Exception; end

struct QASMComment
    string::String
    function QASMComment(s::String)
        regex = r"^\/\/([^\n\r]*)$"
        m = match(regex, s)
        if isnothing(m)
            return nothing
        else
            return new(m.captures[1])
        end
    end
end

function Base.show(io::IO, x::QASMComment) 
    return print(io, string("// ", x.string))
end

struct QASMRegister
    name::Symbol
    width::Int
end

Base.:(==)(q1::QASMRegister, q2::Symbol) = (q1.name == q2)

function Base.show(io::IO, x::QASMRegister) 
    #if x.width > 1
        wstring = string("[", x.width, "]")
    #else
    #    wstring = ""
    #end
    return print(io, string("qreg ", x.name, wstring, ";"))
end

function QASMRegister(s::T) where T<:AbstractString
    regex = r"qreg ([_a-zA-z]\w*)\s*(?:\[(\d+)\])*\s*;"
    m = match(regex, s)
    if isnothing(m)
        return nothing
    else
        try
            return QASMRegister(Symbol(m.captures[1]), 
                    isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2]))
        catch
            throw(ParseError("Could not parse "*s));
        end
    end
end   

struct QASMCRegister
    name::Symbol
    width::Int
end

function Base.show(io::IO, x::QASMCRegister) 
    #if x.width > 1
        wstring = string("[", x.width, "]")
    #else
    #    wstring = ""
    #end
    return print(io, string("creg ", x.name, wstring, ";"))
end

function QASMCRegister(s::T) where T<:AbstractString
    regex = r"creg ([_a-zA-z]\w*)\s*(?:\[(\d+)\])*\s*;"
    m = match(regex, s)
    if isnothing(m)
        return nothing
    else
        try
            return QASMCRegister(Symbol(m.captures[1]), 
                    isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2]))
        catch
            throw(ParseError("Could not parse "*s));
        end
    end
end   

struct QASMQubit
    register::Symbol
    index::Int
end

function Base.show(io::IO, x::QASMQubit) 
    if x.index >= 0
        wstring = string("[", x.index, "]")
    else
        wstring = ""
    end
    return print(io, string(x.register, wstring))
    
end

struct QASMGate
    gate::Symbol
    targets::Vector{QASMQubit}
end


function Base.show(io::IO, x::QASMGate) 
    return print(io, string(x.gate, " ", 
                            join([string(z) for z in x.targets], ", "), ";"))
end

struct QASMGateDef
	name::Symbol
	inputs::Vector{Symbol}
	instructions::Vector{QASMGate}
end

function Base.show(io::IO, gdef::QASMGateDef)
    output = string("gate ", gdef.name, " ", join([string(x) for x in gdef.inputs], " ")..., " {\n")
    output = output * join([string(x) for x in gdef.instructions], "\n")*"\n}"
    return print(io, output)
end

function QASMGate(s::T) where T<:AbstractString
    regex = r"^([_a-zA-Z]\w*)(\s(?:[_a-zA-Z]\w*)(?:\[\d+\])*\,*)*;"
    m = match(regex, s)
    if isnothing(m)
        return nothing
    else
        gstr = m.captures[1]
        qstr = replace(s[length(m.captures[1])+1:end], " "=>"", ";"=>"")
        qm   = eachmatch(r"(?:([_a-zA-Z]\w*)(?:\[(\d+)\])*)", qstr)
        targets = [QASMQubit(Symbol(x.captures[1]), 
                   isnothing(x.captures[2]) ? 0 : parse(Int, x.captures[2]))
                   for x in qm]
        return QASMGate(Symbol(gstr), targets)
    end
end

function extras(s::String)
    regexes = [r"OPENQASM\s+\d\.\d",
               r"""include\s+\"[\w\.]+\"""" ]
    if all([isnothing(match(r, s)) for r in regexes])
        return nothing
    else
        return true
    end
end

mutable struct QASMListing
    comments::Vector{QASMComment}
    gatedefs::Vector{QASMGateDef}
    qregisters::Vector{QASMRegister}
    cregisters::Vector{QASMCRegister}
    instructions::Vector{QASMGate}
end

function QASMListing(fs::Vector{String})
    
    comments = Vector{QASMComment}()
    qregs    = Vector{QASMRegister}()
    cregs    = Vector{QASMCRegister}()
    gates    = Vector{QASMGate}()
    
    for str in fs  
        println(str)
        if !isnothing(extras(str))
            continue
        elseif !isnothing(local qr = QASMRegister(str))
            push!(qregs, qr)
        elseif !isnothing(local cr = QASMCRegister(str))
            push!(cregs,cr)
        elseif !(isnothing(local comm = QASMComment(str)))
            push!(comments, comm)
        elseif !(isnothing(local g = QASMGate(str)))
            push!(gates, g)
        else
            throw(ParseError("Not a QASM file??"))
        end
    end
    
    return QASMListing(comments, qregs, cregs, gates)
end

function Base.show(io::IO, q::QASMListing)
    strs = Vector{String}();
    
    for c in q.comments
        push!(strs, string(c))
    end
    
    push!(strs, "OPENQASM 2.0;")
    push!(strs, """include "qelib1.inc";""")

    for gdef in q.gatedefs
    	push!(strs, string(gdef))
    end
        
    for qr in q.qregisters
        push!(strs, string(qr))
    end
    
    for cr in q.cregisters
        continue
        #push!(strs, string(cr))
    end
    
    for g in q.instructions
        push!(strs, string(g))
    end
    
    return print(io, join(strs, "\n"))
end
    
function used_qubits(q::QASMListing)
    qmap = Dict{Symbol, Vector{Bool}}()
    for  qr in q.qregisters
        qmap[qr.name] = falses(qr.width)
    end
    
    for gate in q.instructions
        for qb in gate.targets
            qmap[qb.register][qb.index+1] = true
        end
    end
    return qmap
end   

function istarget(g::QASMGate, qb::QASMQubit)
    return any([x == qb for x in g.targets])
end

function istarget(g::QASMGate, qr::QASMRegister)
    return any([x.register == qr.name for x in g.targets])
end

function istarget(g::QASMGate, qr::Symbol)
    return any([x.register == qr for x in g.targets])
end

function istarget(g::QASMGate, qr::QASMRegister, idx)
    return any([qb.register == qr.name && qb.index in idx for qb in g.targets])
end

function istarget(g::QASMGate, qr::Symbol, idx)
    return any([qb.register == qr && qb.index in idx for qb in g.targets])
end

function rename_register(g::QASMGate, old_reg::Symbol, new_reg::QASMRegister)
    new_targets = [qb.register == old_reg ? QASMQubit(new_reg.name, qb.index) : qb for qb in g.targets]
    return QASMGate(g.gate, new_targets)
end

function remap_qubits(g::QASMGate, old_reg::Symbol, new_reg::QASMRegister, map)
    new_targets = []
    for qb in g.targets
        if qb.register == old_reg && qb.index in keys(map)
            push!(new_targets, QASMQubit(new_reg.name, map[qb.index]))
        else
            push!(new_targets, qb)
        end
    end
    return QASMGate(g.gate, new_targets)
end

function clear_unused!(q::QASMListing)
    qmap = used_qubits(q)
    new_regs = []
    for qr in q.qregisters
        N = sum(qmap[qr.name])
        if N == 0
            continue
        elseif N < qr.width
            nr = QASMRegister(qr.name, N)
            sub_idx = Dict(j=>(i-1) for (i,j) in enumerate([x for x=1:N if qmap[qr.name][x]]))
            rename_register!(q, qr.name, nr, mapping=sub_idx)
            
        else
            push!(new_regs, qr)
        end
    end
    q.qregisters = new_regs
end

function rename_register!(q::QASMListing, old_reg::Symbol, new_reg::QASMRegister; mapping=nothing)
    #simple case, just rename names
    
    old_width = [x.width for x in q.qregisters if x == old_reg][1]
        
    if old_width <= new_reg.width
        @inbounds for idx in eachindex(q.qregisters)
            if q.qregisters[idx] == old_reg
                q.qregisters[idx] = new_reg
            end
        end
        @inbounds for idx in eachindex(q.instructions)
            if istarget(q.instructions[idx], old_reg)
                q.instructions[idx] = rename_register(q.instructions[idx], old_reg, new_reg)
            end
        end
    else
        @assert !isnothing(mapping) "Must provide a mapping to rename registers!"
        
        push!(q.qregisters, new_reg)
        @inbounds for idx in eachindex(q.instructions)
            if istarget(q.instructions[idx], old_reg, keys(mapping))
                q.instructions[idx] = remap_qubits(q.instructions[idx], old_reg, new_reg, mapping)
            end        
        end
        
    end        
end     

export QASMListing