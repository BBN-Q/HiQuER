
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

end