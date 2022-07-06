#=
Functions for implementing the Misra-Gries edge coloring algorithm.

References:

Husfeldt, T. (2015). Graph colouring algorithms. In L. Beineke & R. Wilson (Eds.), Topics in 
Chromatic Graph Theory (Encyclopedia of Mathematics and its Applications, pp. 277-303). 
Cambridge: Cambridge University Press. doi:10.1017/CBO9781139519793.016

Misra, J, Gries, D. "A constructive proof of Vizing's Theorem". Inf. Proc. Lett. 41(3),
pp. 131-133 (1992). doi:10.1016/0020-0190(92)90041-S

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

using LinearAlgebra
using Graphs, MetaGraphs, GraphPlot
using Colors, Random


##################################################
## General Utility Functions							  

"""
Plot a graph with edge colors; assumes we have colored it with ``\\Delta + 1`` colors. 
"""
function cplot(g::T; layout=spring_layout) where T <: MetaGraph
	nc = Δ(g) + 1
    nodes = collect(1:nv(g))
    cols = distinguishable_colors(nc, [colorant"black", colorant"white"], dropseed=true)
    edgec = []
    black = colorant"black"
    for e in edges(g)
        c = color(g, e)
        push!(edgec, c == 0 ? black : cols[c])
    end
    
    gplot(cg, layout=layout, nodelabel=nodes, linetype="curve",
    	  edgestrokec=edgec, nodefillc=colorant"lightgray")
end

"""
Get the number of colors on the edges of a graph.
"""
function ncolors(g)
    colors = Set()
    for e in edges(g)
        push!(colors, color(g, e))
    end
    return length(colors)
end

"""
Check that no vertex has more than once incident edge sharing a color.
"""
function test_colors(g)
    for v in vertices(g)
        cin = [color(g, Edge(u, v)) for u in neighbors(g, v)]
        if length(cin) != length(Set(cin))
            return false
        end
    end
    return true
end 

"""
Build a test graph.
"""
function test_graph(kind=:bull)
    cg = MetaGraph(smallgraph(kind))
    nc = Δ(cg) + 1
    for e in edges(cg)
        set_color!(cg, e, 0)
    end
    return cg
end

##################################################
## Utility Functions for Graph Colors

"""
Get the color of a graph edge between vertices (`u`, `v`).
"""
function color(g::T, u::Int, v::Int) where T <: MetaGraph
    return get_prop(g, Edge(u,v), :color)
end

"""
Get the color of a graph edge `e`.
"""
function color(g::T, e) where T <: MetaGraph
    return get_prop(g, e, :color)
end

"""
Set the color of a graph edge between vertices (`u`, `v`).
"""
function set_color!(g::T, u::Int, v::Int, c) where T <: MetaGraph
    set_prop!(g, Edge(u,v), :color, c)
end

"""
Set the color of a graph edge `e`.
"""
function set_color!(g::T, e, c) where T <: MetaGraph
    set_prop!(g, e, :color, c)
end

"""
Check that an edge has a color: the color value of '0' means no color has been assigned.  
"""
function has_color(g::T, e) where T <: MetaGraph
    return color(g, e) > 0
end

"""
Check that color `c` on node `v` is free: no edges incident to `v` carry that color.
"""
function color_free(g, v, c)
    for u in neighbors(g, v)
        if color(g, Edge(v, u)) == c
            return false
        end
    end
    return true
end

##################################################
## Misra-Gries Algorithm	

"""
Build the maximum fan between nodes `x` and `f`.
"""
function max_fan(g, x, f)
    x_fan = [f]
    x_n   = neighbors(g, x)
    is_max = false
    prev_node = f
    while !(is_max)
        is_max = true
        for v in x_n
            if !(v in x_fan) &&
                has_color(g, Edge(x, v)) &&
                color_free(g, prev_node, color(g, Edge(x, v)))
                push!(x_fan, v)
                prev_node = v
                is_max = false
            end
        end
    end
    return x_fan
end 

"""
Find colors `c`, `d` that a free between node `x` and the terminal node of `fan`.
"""
function cd_colors(g, x, fan)
    c = 1
    d = 1
    while !color_free(g, x, c)
        c += 1
    end
    while !color_free(g, fan[end], d)
        d += 1
    end
    return c,d
end

"""
Invert the cd_u path. Returns the length of the path, and modifies the graph.
"""
function invert_cd!(g, u, c, d)
    is_max = false
    visited = Set(u)
    while !is_max
        is_max = true
        for v in neighbors(g, u)
            if color(g, u, v) == d && !(v in visited)
                set_color!(g, u, v, c)
                u = v
                c, d = d, c
                is_max = false
                push!(visited, v)
                break
            end
        end
    end
    return length(visited) - 1
end     


function rotate_fan!(g, x, fan)
    for (a, b) in zip(fan, fan[2:end])
        c = color(g, x, b)
        set_color!(g, x, a, c)
    end
end

"""
Find the node in a fan that has color `d` free.
"""
function find_in_fan(g, d, fan)
    for (idx, u) in enumerate(fan)
        if color_free(g, u, d)
            return (idx, u)
        end
    end
    return (0, nothing)
end

"""
Implement Misra-Gries algorithm for edge coloring of graph 'g'
"""
function color_graph_edges!(g)
    idx = 0
    for e in edges(g)
        u, v = src(e), dst(e)
        fan = max_fan(g, u, v)
        c, d = cd_colors(g, u, fan)
        plen = invert_cd!(g, u, c, d)
        if plen > 0
            wi, w = find_in_fan(g, d, fan)
        else
            wi, w = length(fan), fan[end]
        end
        rotate_fan!(g, u, fan[1:wi])
        set_color!(g, u, w, d) 
    end
    @assert test_colors(g) "Graph coloring failed!"
end