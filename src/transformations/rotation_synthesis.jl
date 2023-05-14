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

Use gridsynth for rotation synthesis.

=#

using LRUCache

cliffordT = Set([H, S, T, X, Y, Z, CNOT])
union!(cliffordT, Set([g' for g in cliffordT]))

function gridsynth(angle::Float64, ϵ=1e-10)
    decomp = readchomp(`gridsynth -e $ϵ -p $angle`)
    gates = Vector{HiQuER.AbstractGate}()
    for u in decomp
        if u == 'X'
            push!(gates, X180)
        else
            push!(gates, eval(Symbol(u)))
        end
    end
    return gates
end

const synth_lru = LRU{Tuple{Float64, Float64}, Vector{HiQuER.AbstractGate}}(maxsize=10000)

function cached_gridsynth(θ::FloatAngle, ϵ=1e-10)
    angle = θ.value
    old = filter(x -> abs(x[1] - angle) < ϵ, collect(keys(synth_lru)))
    if length(old) > 0
        angle = old[1][1]
    end
    get!(synth_lru, (angle, ϵ)) do
        gridsynth(angle, ϵ)
    end
end   
