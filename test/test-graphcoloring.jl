using Graphs, MetaGraphs


@testset "Graph Coloring" begin
	
	graph = HiQuER.test_graph(:bull)
	HiQuER.color_graph_edges!(graph)
	@test HiQuER.test_colors(graph)

	graph = HiQuER.test_graph(:pappus)
	HiQuER.color_graph_edges!(graph)
	@test HiQuER.test_colors(graph)

	graph = HiQuER.test_graph(:karate)
	HiQuER.color_graph_edges!(graph)
	@test HiQuER.test_colors(graph)

	#now a random graph, including self-loops
	M = rand(8,8)
	graph = MetaGraph(Graph(M'*M))
 	for e in edges(graph)
    	HiQuER.set_color!(graph, e, 0)
	end
	HiQuER.color_graph_edges!(graph)
	@test HiQuER.test_colors(graph)

end