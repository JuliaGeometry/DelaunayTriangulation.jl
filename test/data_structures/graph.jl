using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using SimpleGraphs

@testset "Constructing empty graphs" begin
    sg1 = UndirectedGraph{Int64}()
    sg2 = UndirectedGraph{Int32}()
    g1 = DT.Graph{Int64}()
    g2 = DT.Graph{Int32}()
    @test g1.graph == sg1
    @test g2.graph == sg2
    @test g1 == DT.Graph()
    @test g1 == DT.Graph(sg1)
    @test g2 == DT.Graph(sg2)
end

global sg = UndirectedGraph([0 0 0 1 0 1
    0 0 1 1 0 0
    0 1 0 1 0 1
    1 1 1 0 1 1
    0 0 0 1 0 0
    1 0 1 1 0 0])
global g = DT.Graph(sg)

@testset "Getting graph and getting statistics" begin
    @test DT.get_graph(g) == sg
    @test DT.get_edges(g) == sg.E
    @test DT.get_neighbours(g) == sg.N
    @test DT.get_vertices(g) == sg.V
    @test DT.get_neighbours(g, 1) == sg.N[1]
    @test DT.get_neighbours(g, 3) == Set{Int64}((2, 4, 6))
    @test DT.num_neighbours(g, 5) == 1
    @test DT.num_neighbours(g, 2) == 2
    @test DT.num_neighbours(g, 4) == 5
    @test DT.num_edges(g) == 8
end

@testset "Adding vertices" begin
    DT.add_vertex!(g, 7)
    @test 7 ∈ sg.V
    DT.add_vertex!(g, 13, 20)
    @test 13 ∈ sg.V && 20 ∈ sg.V
end

@testset "Adding neighbours" begin
    DT.add_neighbour!(g, 3, 1)
    @test 1 ∈ sg.N[3] && 3 ∈ sg.N[1]
    DT.add_neighbour!(g, 10, 15, 21)
    @test 10 ∈ sg.V && 15 ∈ sg.V && 21 ∈ sg.V
    @test 15 ∈ sg.N[10] && 21 ∈ sg.N[10]
end

global T1 = (25, 26, 27)
global T2 = (28, 30, 50)
global T3 = (51, 52, 53)
global T4 = (60, 65, 0)
global T5 = (-1, 7, 5)

@testset "Adding triangles" begin
    DT.add_triangle!(g, T1...)
    DT.add_triangle!(g, T2)
    DT.add_triangle!(g, T3, T4, T5)
    for ((i, j, k)) in (T1, T2, T3, T4, T5)
        @test i ∈ sg.V && j ∈ sg.V && k ∈ sg.V
        @test i ∈ sg.N[j] && i ∈ sg.N[k]
        @test j ∈ sg.N[i] && j ∈ sg.N[k]
        @test k ∈ get_neighbours(g, i) && k ∈ get_neighbours(g, j)
    end
end

@testset "Deleting triangles" begin
    DT.delete_triangle!(g, T1...)
    DT.delete_triangle!(g, T2)
    DT.delete_triangle!(g, T3, T4, T5)
    for ((i, j, k)) in (T1, T2, T3, T4, T5)
        @test !(i ∈ sg.N[j] && i ∈ sg.N[k])
        @test !(j ∈ sg.N[i] && j ∈ sg.N[k])
        @test !(k ∈ get_neighbours(g, i) && k ∈ get_neighbours(g, j))
    end
end

@testset "Deleting neighbours" begin
    DT.delete_neighbour!(g, 30, 50)
    @test 50 ∉ get_neighbours(g, 30) && 30 ∉ get_neighbours(g, 50)
    DT.delete_neighbour!(g, 51, 52, 53)
    @test 52 ∉ get_neighbours(g, 51) && 53 ∉ get_neighbours(g, 51)
end

@testset "Deleting vertices" begin
    DT.delete_vertex!(g, 3)
    @test 3 ∉ sg.V
    DT.delete_vertex!(g, 60, 65, 0)
    @test 60 ∉ sg.V && 65 ∉ sg.V && 0 ∉ sg.V
end

@testset "Deleting boundary indices" begin
    bidx = DT.BoundaryIndex
    DT.add_vertex!(g, bidx, bidx - 1, bidx - 2, bidx - 3, bidx - 4)
    DT.delete_boundary_vertices_from_graph!(g)
    @test bidx ∉ sg.V && bidx - 1 ∉ sg.V && bidx - 2 ∉ sg.V && bidx - 3 ∉ sg.V &&
          bidx - 4 ∉ sg.V
end

global sg = UndirectedGraph([0 0 0 1 0 1
    0 0 1 1 0 0
    0 1 0 1 0 1
    1 1 1 0 1 1
    0 0 0 1 0 0
    1 0 1 1 0 0])
global g = DT.Graph(sg)

@testset "Removing empty parts" begin
    clean_dg = deepcopy(g)
    DT.add_vertex!(g, 13, 5)
    @test g ≠ clean_dg
    DT.delete_neighbour!(g, 13, 5)
    @test g ≠ clean_dg
    DT.clear_empty_points!(g)
    @test g == clean_dg
end

@testset "Number of vertices and each_vertex" begin
    @test num_vertices(g) == 6
    @test each_vertex(g) == get_vertices(g)
end