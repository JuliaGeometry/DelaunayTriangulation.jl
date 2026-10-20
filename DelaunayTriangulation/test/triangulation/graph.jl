using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures

@testset "Constructing empty graphs" begin
    g1 = DT.Graph{Int64}()
    g2 = DT.Graph{Int32}()
    @test g1.vertices == Set{Int64}()
    @test g1.edges == Set{NTuple{2, Int64}}()
    @test g1.neighbours == Dict{Int64, Set{Int64}}()
    @test g2.vertices == Set{Int32}()
    @test g2.edges == Set{NTuple{2, Int32}}()
    @test g2.neighbours == Dict{Int32, Set{Int32}}()
end

global sg = [
    0 0 0 1 0 1
    0 0 1 1 0 0
    0 1 0 1 0 1
    1 1 1 0 1 1
    0 0 0 1 0 0
    1 0 1 1 0 0
]
global g = _make_graph_from_adjacency(sg)

@testset "Getting graph and getting statistics" begin
    @test DT.get_edges(g) == g.edges
    @test DT.get_neighbours(g) == g.neighbours
    @test DT.get_vertices(g) == g.vertices
    @test DT.get_neighbours(g, 1) == g.neighbours[1]
    @test DT.get_neighbours(g, 3) == Set{Int}((2, 4, 6))
    @test DT.num_neighbours(g, 5) == 1
    @test DT.num_neighbours(g, 2) == 2
    @test DT.num_neighbours(g, 4) == 5
    @test DT.num_edges(g) == 8
end

@testset "Adding vertices" begin
    DT.add_vertex!(g, 7)
    @test 7 ∈ g.vertices
    DT.add_vertex!(g, 13, 20)
    @test 13 ∈ g.vertices && 20 ∈ g.vertices
end

@testset "Adding neighbours" begin
    DT.add_neighbour!(g, 3, 1)
    @test 1 ∈ g.neighbours[3] && 3 ∈ g.neighbours[1]
    DT.add_neighbour!(g, 10, 15, 21)
    @test 10 ∈ g.vertices && 15 ∈ g.vertices && 21 ∈ g.vertices
    @test 15 ∈ g.neighbours[10] && 21 ∈ g.neighbours[10]
    @test (3, 1) ∈ g.edges
end

global T1 = (25, 26, 27)
global T2 = (28, 30, 50)
global T3 = (51, 52, 53)
global T4 = (60, 65, 0)
global T5 = (-1, 7, 5)

@testset "Adding triangles" begin
    DT.add_triangle!(g, T1...)
    DT.add_triangle!(g, T2)
    DT.add_triangle!.(Ref(g), (T3, T4, T5))
    for ((i, j, k)) in (T1, T2, T3, T4, T5)
        @test i ∈ g.vertices && j ∈ g.vertices && k ∈ g.vertices
        @test i ∈ g.neighbours[j] && i ∈ g.neighbours[k]
        @test j ∈ g.neighbours[i] && j ∈ g.neighbours[k]
        @test k ∈ get_neighbours(g, i) && k ∈ get_neighbours(g, j)
    end
end

@testset "Deleting triangles" begin
    DT.delete_triangle!(g, T1...)
    DT.delete_triangle!(g, T2)
    DT.delete_triangle!.(Ref(g), (T3, T4, T5))
    for ((i, j, k)) in (T1, T2, T3, T4, T5)
        i_in_j = haskey(g.neighbours, j) && i ∈ g.neighbours[j] && i ∈ get_neighbours(g, j)
        i_in_k = haskey(g.neighbours, k) && i ∈ g.neighbours[k] && i ∈ get_neighbours(g, k)
        j_in_i = haskey(g.neighbours, i) && j ∈ g.neighbours[i] && j ∈ get_neighbours(g, i)
        j_in_k = haskey(g.neighbours, k) && j ∈ g.neighbours[k] && j ∈ get_neighbours(g, k)
        k_in_i = haskey(g.neighbours, i) && k ∈ g.neighbours[i] && k ∈ get_neighbours(g, i)
        k_in_j = haskey(g.neighbours, j) && k ∈ g.neighbours[j] && k ∈ get_neighbours(g, j)
        @test !(i_in_j || i_in_k)
        @test !(j_in_i || j_in_k)
        @test !(k_in_i || k_in_j)
    end
end

@testset "Deleting neighbours" begin
    DT.delete_neighbour!(g, 3, 6)
    @test 3 ∉ get_neighbours(g, 6) && 6 ∉ get_neighbours(g, 3)
    DT.delete_neighbour!(g, 4, 6, 2, 3, 1)
    @test 4 ∉ get_neighbours(g, 6) && 6 ∉ get_neighbours(g, 4)
    @test 4 ∉ get_neighbours(g, 2) && 2 ∉ get_neighbours(g, 4)
    @test 4 ∉ get_neighbours(g, 3) && 3 ∉ get_neighbours(g, 4)
    @test 4 ∉ get_neighbours(g, 1) && 1 ∉ get_neighbours(g, 4)
    DT.delete_neighbour!(g, 10, 15, 21)
    @test !haskey(g.neighbours, 10) && !haskey(g.neighbours, 15) && !haskey(g.neighbours, 21)
    DT.delete_neighbour!(g, 1, 6)
    @test haskey(g.neighbours, 1) && !haskey(g.neighbours, 6) && 6 ∉ get_neighbours(g, 1)
end

@testset "Deleting vertices" begin
    DT.delete_vertex!(g, 3)
    @test 3 ∉ g.vertices
    DT.delete_vertex!(g, 60, 65, 0)
    @test 60 ∉ g.vertices && 65 ∉ g.vertices && 0 ∉ g.vertices
    @test 60 ∉ keys(g.neighbours) && 65 ∉ keys(g.neighbours) && 0 ∉ keys(g.neighbours)
    @test 3 ∉ keys(g.neighbours)
end

@testset "Deleting ghost vertices" begin
    bidx = DT.𝒢
    DT.add_vertex!(g, bidx, bidx - 1, bidx - 2, bidx - 3, bidx - 4)
    DT.delete_ghost_vertices_from_graph!(g)
    @test bidx ∉ g.vertices && bidx - 1 ∉ g.vertices && bidx - 2 ∉ g.vertices && bidx - 3 ∉ g.vertices &&
        bidx - 4 ∉ g.vertices
end

global sg = [
    0 0 0 1 0 1
    0 0 1 1 0 0
    0 1 0 1 0 1
    1 1 1 0 1 1
    0 0 0 1 0 0
    1 0 1 1 0 0
]
global g = _make_graph_from_adjacency(sg)

@testset "Removing empty parts" begin
    clean_dg = deepcopy(g)
    DT.add_vertex!(g, 13, 5)
    @test g ≠ clean_dg
    empty!(g.neighbours[13])
    delete!(g.neighbours[5], 13)
    @test g ≠ clean_dg
    DT.clear_empty_vertices!(g)
    @test g == clean_dg
end

@testset "Number of vertices" begin
    @test num_vertices(g) == 6
end

@testset "has_vertex and has_ghost_vertices" begin
    bidx = DT.𝒢
    DT.add_vertex!(g, bidx, bidx - 1, bidx - 2, bidx - 3, bidx - 4)
    @test DT.has_vertex(g, bidx)
    @test DT.has_vertex(g, bidx - 4)
    @test DT.has_vertex(g, 5)
    @test DT.has_ghost_vertices(g)
    DT.delete_ghost_vertices_from_graph!(g)
    @test !DT.has_ghost_vertices(g)
end

@testset "empty!" begin
    tri = triangulate(rand(2, 50))
    graph = get_graph(tri)
    @test !isempty(get_neighbours(graph))
    @test !isempty(DT.get_edges(graph))
    @test !isempty(DT.get_vertices(graph))
    empty!(graph)
    @test isempty(get_neighbours(graph))
    @test isempty(DT.get_edges(graph))
    @test isempty(DT.get_vertices(graph))
end

@testset "copy/deepcopy" begin
    tri = triangulate(rand(2, 50))
    graph = get_graph(tri)
    graph2 = copy(graph)
    @test graph == graph2 && !(graph === graph2)
    S = get_neighbours(tri, 1)
    S2 = get_neighbours(graph2, 1)
    @test S === S2
    graph2 = deepcopy(graph)
    @test graph == graph2 
    S = get_neighbours(tri, 1)
    S2 = get_neighbours(graph2, 1)
    @test S == S2 && !(S === S2)
end

@testset "Deleting a vertex deletes its edges" begin
    g = DT.Graph{Int}()
    DT.add_triangle!(g, 1, 2, 3)
    DT.add_triangle!(g, 1, 3, 4)
    DT.delete_vertex!(g, 3)
    @test 3 ∉ g.vertices && !haskey(g.neighbours, 3)
    @test all(e -> 3 ∉ e, g.edges)
    @test DT.num_edges(g) == 2 # (1, 2) and (1, 4)
    @test all(v -> 3 ∉ get_neighbours(g, v), (1, 2, 4))
    DT.clear_empty_vertices!(g)
    @test g.vertices == Set((1, 2, 4))
end
