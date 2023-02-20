using Graphs, MetaGraphs, LinearAlgebra


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
	graph = HiQuER.graph_from_matrix(M'*M)
	HiQuER.color_graph_edges!(graph)
	@test HiQuER.test_colors(graph)

	M = rand(12,12)
	H = Tridiagonal(M'*M)
	Hc = HiQuER.split_by_color(H)
	@test all(H .≈ sum(Hc))

end