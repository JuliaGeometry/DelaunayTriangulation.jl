using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StableRNGs

# These tests reproduce the examples shown in the docstrings, which are not run as doctests.
# Sets and Dicts are compared by value, so that the tests do not depend on hash iteration order.

@testset "constrained_triangulation.jl" begin
    @testset "fix_segments!" begin
        segments = [(2, 15), (2, 28), (2, 41)]
        bad_indices = [1, 2, 3]
        @test DT.fix_segments!(segments, bad_indices) == [(2, 15), (15, 28), (28, 41)]
        segments = [(2, 7), (2, 12), (12, 17), (2, 22), (2, 27), (2, 32), (32, 37), (2, 42), (42, 47)]
        bad_indices = [2, 4, 5, 6, 8]
        @test DT.fix_segments!(segments, bad_indices) == [(2, 7), (7, 12), (12, 17), (17, 22), (22, 27), (27, 32), (32, 37), (37, 42), (42, 47)]
    end

    @testset "connect_segments!" begin
        segments = [(7, 12), (12, 17), (17, 22), (32, 37), (37, 42), (42, 47)]
        @test DT.connect_segments!(segments) == [(7, 12), (12, 17), (17, 22), (22, 32), (32, 37), (37, 42), (42, 47)]
    end

    @testset "extend_segments!" begin
        segments = [(2, 7), (7, 12), (12, 49)]
        segment = (1, 68)
        @test DT.extend_segments!(segments, segment) == [(1, 2), (2, 7), (7, 12), (12, 49), (49, 68)]
    end

    @testset "split_segment!" begin
        segments = Set(((2, 3), (3, 5), (10, 12)))
        collinear_segments = [(2, 10), (11, 15), (2, 3)]
        segment = (3, 5)
        @test DT.split_segment!(segments, segment, collinear_segments) == Set(((2, 10), (2, 3), (11, 15), (10, 12)))
    end
end

@testset "clipped.jl" begin
    @testset "get_shared_vertex" begin
        @test DT.get_shared_vertex((1, 3), (5, 7)) == 0
        @test DT.get_shared_vertex((1, 3), (3, 7)) == 3
        @test DT.get_shared_vertex((10, 3), (10, 5)) == 10
        @test DT.get_shared_vertex((9, 4), (9, 5)) == 9
    end
end

@testset "adjacent.jl" begin
    @testset "get_adjacent(adj)" begin
        d = Dict((1, 2) => 3, (2, 3) => 1, (3, 1) => 2)
        adj = DT.Adjacent(d)
        @test adj isa DT.Adjacent{Int64, Tuple{Int64, Int64}}
        @test get_adjacent(adj) == Dict((1, 2) => 3, (3, 1) => 2, (2, 3) => 1)
        @test get_adjacent(adj) == d
    end

    @testset "get_adjacent(adj, uv)" begin
        adj = DT.Adjacent(Dict((1, 2) => 3, (2, 3) => 1, (3, 1) => 2, (4, 5) => -1))
        @test get_adjacent(adj, 4, 5) == -1
        @test get_adjacent(adj, (3, 1)) == 2
        @test get_adjacent(adj, (1, 2)) == 3
        @test get_adjacent(adj, 17, 5) == 0
        @test get_adjacent(adj, (1, 6)) == 0
    end

    @testset "add_adjacent!" begin
        adj = DT.Adjacent{Int64, NTuple{2, Int64}}()
        @test get_adjacent(DT.add_adjacent!(adj, 1, 2, 3)) == Dict((1, 2) => 3)
        @test get_adjacent(DT.add_adjacent!(adj, (2, 3), 1)) == Dict((1, 2) => 3, (2, 3) => 1)
        @test get_adjacent(DT.add_adjacent!(adj, 3, 1, 2)) == Dict((1, 2) => 3, (3, 1) => 2, (2, 3) => 1)
    end

    @testset "delete_adjacent!" begin
        adj = DT.Adjacent(Dict((2, 7) => 6, (7, 6) => 2, (6, 2) => 2, (17, 3) => -1, (-1, 5) => 17, (5, 17) => -1))
        @test get_adjacent(DT.delete_adjacent!(adj, 2, 7)) == Dict((-1, 5) => 17, (17, 3) => -1, (6, 2) => 2, (5, 17) => -1, (7, 6) => 2)
        @test get_adjacent(DT.delete_adjacent!(adj, (6, 2))) == Dict((-1, 5) => 17, (17, 3) => -1, (5, 17) => -1, (7, 6) => 2)
        @test get_adjacent(DT.delete_adjacent!(adj, 5, 17)) == Dict((-1, 5) => 17, (17, 3) => -1, (7, 6) => 2)
    end

    @testset "add_triangle!(adj, u, v, w)" begin
        adj = DT.Adjacent{Int32, NTuple{2, Int32}}()
        res = add_triangle!(adj, 1, 2, 3)
        @test res isa DT.Adjacent{Int32, Tuple{Int32, Int32}}
        @test get_adjacent(res) == Dict((1, 2) => 3, (3, 1) => 2, (2, 3) => 1)
        @test get_adjacent(add_triangle!(adj, 6, -1, 7)) == Dict((1, 2) => 3, (3, 1) => 2, (6, -1) => 7, (-1, 7) => 6, (2, 3) => 1, (7, 6) => -1)
    end

    @testset "delete_triangle!(adj, u, v, w)" begin
        adj = DT.Adjacent{Int32, NTuple{2, Int32}}()
        add_triangle!(adj, 1, 6, 7)
        add_triangle!(adj, 17, 3, 5)
        @test get_adjacent(adj) == Dict((17, 3) => 5, (1, 6) => 7, (6, 7) => 1, (7, 1) => 6, (5, 17) => 3, (3, 5) => 17)
        @test get_adjacent(delete_triangle!(adj, 3, 5, 17)) == Dict((1, 6) => 7, (6, 7) => 1, (7, 1) => 6)
        @test isempty(get_adjacent(delete_triangle!(adj, 7, 1, 6)))
    end
end

@testset "adjacent2vertex.jl" begin
    @testset "get_adjacent2vertex(adj2v)" begin
        e1 = Set(((1, 2), (5, 3), (7, 8)))
        e2 = Set(((2, 3), (13, 5), (-1, 7)))
        d = Dict(9 => e1, 6 => e2)
        adj2v = DT.Adjacent2Vertex(d)
        @test adj2v isa DT.Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}}
        @test get_adjacent2vertex(adj2v) == Dict(6 => Set([(13, 5), (-1, 7), (2, 3)]), 9 => Set([(1, 2), (7, 8), (5, 3)]))
        @test get_adjacent2vertex(adj2v) == d
    end

    @testset "get_adjacent2vertex(adj2v, w)" begin
        adj2v = DT.Adjacent2Vertex(Dict(1 => Set(((2, 3), (5, 7), (8, 9))), 5 => Set(((1, 2), (7, 9), (8, 3)))))
        @test get_adjacent2vertex(adj2v, 1) == Set(((8, 9), (5, 7), (2, 3)))
        @test get_adjacent2vertex(adj2v, 5) == Set(((1, 2), (8, 3), (7, 9)))
    end

    @testset "add_adjacent2vertex!" begin
        adj2v = DT.Adjacent2Vertex{Int64, Set{NTuple{2, Int64}}}()
        @test isempty(get_adjacent2vertex(adj2v))
        @test get_adjacent2vertex(DT.add_adjacent2vertex!(adj2v, 1, (2, 3))) == Dict(1 => Set([(2, 3)]))
        @test get_adjacent2vertex(DT.add_adjacent2vertex!(adj2v, 1, 5, 7)) == Dict(1 => Set([(5, 7), (2, 3)]))
        @test get_adjacent2vertex(DT.add_adjacent2vertex!(adj2v, 17, (5, -1))) == Dict(17 => Set([(5, -1)]), 1 => Set([(5, 7), (2, 3)]))
    end

    @testset "delete_adjacent2vertex!(adj2v, w, uv)" begin
        adj2v = DT.Adjacent2Vertex(Dict(1 => Set(((2, 3), (5, 7), (8, 9))), 5 => Set(((1, 2), (7, 9), (8, 3)))))
        @test get_adjacent2vertex(DT.delete_adjacent2vertex!(adj2v, 5, 8, 3)) == Dict(5 => Set([(1, 2), (7, 9)]), 1 => Set([(8, 9), (5, 7), (2, 3)]))
        @test get_adjacent2vertex(DT.delete_adjacent2vertex!(adj2v, 1, (2, 3))) == Dict(5 => Set([(1, 2), (7, 9)]), 1 => Set([(8, 9), (5, 7)]))
    end

    @testset "delete_adjacent2vertex!(adj2v, w)" begin
        adj2v = DT.Adjacent2Vertex(Dict(1 => Set(((2, 3), (5, 7))), 5 => Set(((-1, 2), (2, 3)))))
        @test get_adjacent2vertex(DT.delete_adjacent2vertex!(adj2v, 1)) == Dict(5 => Set([(-1, 2), (2, 3)]))
        @test isempty(get_adjacent2vertex(DT.delete_adjacent2vertex!(adj2v, 5)))
    end

    @testset "add_triangle!(adj2v, u, v, w)" begin
        adj2v = DT.Adjacent2Vertex{Int32, Set{NTuple{2, Int32}}}()
        res = add_triangle!(adj2v, 17, 5, 8)
        @test res isa DT.Adjacent2Vertex{Int32, Set{Tuple{Int32, Int32}}}
        @test get_adjacent2vertex(res) == Dict(5 => Set([(8, 17)]), 8 => Set([(17, 5)]), 17 => Set([(5, 8)]))
        @test get_adjacent2vertex(add_triangle!(adj2v, 1, 5, 13)) == Dict(5 => Set([(8, 17), (13, 1)]), 13 => Set([(1, 5)]), 8 => Set([(17, 5)]), 17 => Set([(5, 8)]), 1 => Set([(5, 13)]))
    end

    @testset "delete_triangle!(adj2v, u, v, w)" begin
        adj2v = DT.Adjacent2Vertex{Int32, Set{NTuple{2, Int32}}}()
        @test get_adjacent2vertex(add_triangle!(adj2v, 1, 2, 3)) == Dict(2 => Set([(3, 1)]), 3 => Set([(1, 2)]), 1 => Set([(2, 3)]))
        @test get_adjacent2vertex(add_triangle!(adj2v, 17, 5, 2)) == Dict(5 => Set([(2, 17)]), 2 => Set([(3, 1), (17, 5)]), 17 => Set([(5, 2)]), 3 => Set([(1, 2)]), 1 => Set([(2, 3)]))
        @test get_adjacent2vertex(delete_triangle!(adj2v, 5, 2, 17)) == Dict(5 => Set(), 2 => Set([(3, 1)]), 17 => Set(), 3 => Set([(1, 2)]), 1 => Set([(2, 3)]))
        @test get_adjacent2vertex(delete_triangle!(adj2v, 2, 3, 1)) == Dict(5 => Set(), 2 => Set(), 17 => Set(), 3 => Set(), 1 => Set())
    end

    @testset "clear_empty_keys!(adj2v)" begin
        adj2v = DT.Adjacent2Vertex{Int64, Set{NTuple{2, Int64}}}()
        @test get_adjacent2vertex(add_triangle!(adj2v, 1, 2, 3)) == Dict(2 => Set([(3, 1)]), 3 => Set([(1, 2)]), 1 => Set([(2, 3)]))
        @test get_adjacent2vertex(delete_triangle!(adj2v, 2, 3, 1)) == Dict(2 => Set(), 3 => Set(), 1 => Set())
        @test isempty(get_adjacent2vertex(DT.clear_empty_keys!(adj2v)))
    end
end

@testset "boundary_nodes.jl" begin
    @testset "has_multiple_curves" begin
        @test !DT.has_multiple_curves([1, 2, 3, 1])
        @test !DT.has_multiple_curves([[1, 2, 3], [3, 4, 1]])
        @test DT.has_multiple_curves([[[1, 2, 3], [3, 4, 1]], [[5, 6, 7, 8, 5]]])
    end

    @testset "has_multiple_sections" begin
        @test !DT.has_multiple_sections([1, 2, 3, 1])
        @test DT.has_multiple_sections([[1, 2, 3], [3, 4, 1]])
        @test DT.has_multiple_sections([[[1, 2, 3], [3, 4, 1]], [[5, 6, 7, 8, 5]]])
    end

    @testset "num_curves" begin
        @test DT.num_curves([1, 2, 3, 1]) == 1
        @test DT.num_curves([[1, 2, 3], [3, 4, 1]]) == 1
        @test DT.num_curves([[[1, 2, 3], [3, 4, 1]], [[5, 6, 7, 8, 5]]]) == 2
    end

    @testset "num_sections" begin
        @test DT.num_sections([1, 2, 3, 4, 5, 1]) == 1
        @test DT.num_sections([[1, 2, 3, 4], [4, 5, 1]]) == 2
        @test DT.num_sections([[1, 2, 3], [3, 4, 5, 6, 7, 8], [8, 9], [9, 1]]) == 4
    end

    @testset "get_boundary_nodes" begin
        @test get_boundary_nodes([[[1, 2, 3, 4], [4, 5, 1]], [[6, 7, 8, 9], [9, 10, 6]]], 2) == [[6, 7, 8, 9], [9, 10, 6]]
        @test get_boundary_nodes([[1, 2, 3, 4], [4, 5, 1]], 1) == [1, 2, 3, 4]
        @test get_boundary_nodes([1, 2, 3, 4, 5, 6, 1], 4) == 4
        @test get_boundary_nodes([[[1, 2, 3, 4], [4, 5, 1]], [[6, 7, 8, 9], [9, 10, 6]]], 1, 2) == [4, 5, 1]
        @test get_boundary_nodes([[1, 2, 3, 4], [4, 5, 6, 1]], 2, 3) == 6
        @test get_boundary_nodes([1, 2, 3, 4, 5, 1], [1, 2, 3, 4, 5, 1]) == [1, 2, 3, 4, 5, 1]
    end

    @testset "each_boundary_node" begin
        @test DT.each_boundary_node([7, 8, 19, 2, 17]) == [7, 8, 19, 2, 17]
        @test DT.each_boundary_node([7, 8, 19, 2, 17, 7]) == [7, 8, 19, 2, 17, 7]
    end

    @testset "construct_ghost_vertex_map" begin
        gv_map = DT.construct_ghost_vertex_map([1, 2, 3, 4, 5, 1])
        @test gv_map isa Dict{Int64, Vector{Int64}}
        @test gv_map == Dict(-1 => [1, 2, 3, 4, 5, 1])
        gv_map = DT.construct_ghost_vertex_map([[17, 29, 23, 5, 2, 1], [1, 50, 51, 52], [52, 1]])
        @test gv_map isa Dict{Int64, Int64}
        @test gv_map == Dict(-1 => 1, -3 => 3, -2 => 2)
        gv_map = DT.construct_ghost_vertex_map([[[1, 5, 17, 18, 1]], [[23, 29, 31, 33], [33, 107, 101], [101, 99, 85, 23]]])
        @test gv_map isa Dict{Int64, Tuple{Int64, Int64}}
        @test gv_map == Dict(-1 => (1, 1), -3 => (2, 2), -2 => (2, 1), -4 => (2, 3))
    end

    @testset "construct_boundary_edge_map" begin
        bn = [17, 18, 15, 4, 3, 17]
        @test DT.construct_boundary_edge_map(bn) == Dict((18, 15) => (bn, 2), (3, 17) => (bn, 5), (17, 18) => (bn, 1), (4, 3) => (bn, 4), (15, 4) => (bn, 3))
        @test DT.construct_boundary_edge_map([[5, 17, 3, 9], [9, 18, 13, 1], [1, 93, 57, 5]]) == Dict(
            (18, 13) => (2, 2), (17, 3) => (1, 2), (9, 18) => (2, 1), (13, 1) => (2, 3), (3, 9) => (1, 3),
            (93, 57) => (3, 2), (5, 17) => (1, 1), (57, 5) => (3, 3), (1, 93) => (3, 1),
        )
        @test DT.construct_boundary_edge_map([[[2, 5, 10], [10, 11, 2]], [[27, 28, 29, 30], [30, 31, 85, 91], [91, 92, 27]]]) == Dict(
            (92, 27) => ((2, 3), 2), (2, 5) => ((1, 1), 1), (11, 2) => ((1, 2), 2), (10, 11) => ((1, 2), 1),
            (30, 31) => ((2, 2), 1), (91, 92) => ((2, 3), 1), (29, 30) => ((2, 1), 3), (31, 85) => ((2, 2), 2),
            (27, 28) => ((2, 1), 1), (5, 10) => ((1, 1), 2), (28, 29) => ((2, 1), 2), (85, 91) => ((2, 2), 3),
        )
    end

    @testset "insert_boundary_node!" begin
        boundary_nodes = [1, 2, 3, 4, 5, 1]
        @test DT.insert_boundary_node!(boundary_nodes, (boundary_nodes, 4), 23) == [1, 2, 3, 23, 4, 5, 1]
        boundary_nodes = [[7, 13, 9, 25], [25, 26, 29, 7]]
        @test DT.insert_boundary_node!(boundary_nodes, (2, 1), 57) == [[7, 13, 9, 25], [57, 25, 26, 29, 7]]
        boundary_nodes = [[[17, 23, 18, 25], [25, 26, 81, 91], [91, 101, 17]], [[1, 5, 9, 13], [13, 15, 1]]]
        @test DT.insert_boundary_node!(boundary_nodes, ((1, 3), 3), 1001) == [[[17, 23, 18, 25], [25, 26, 81, 91], [91, 101, 1001, 17]], [[1, 5, 9, 13], [13, 15, 1]]]
    end

    @testset "delete_boundary_node!" begin
        boundary_nodes = [71, 25, 33, 44, 55, 10]
        @test DT.delete_boundary_node!(boundary_nodes, (boundary_nodes, 4)) == [71, 25, 33, 55, 10]
        boundary_nodes = [[7, 13, 9, 25], [25, 26, 29, 7]]
        @test DT.delete_boundary_node!(boundary_nodes, (2, 3)) == [[7, 13, 9, 25], [25, 26, 7]]
        boundary_nodes = [[[17, 23, 18, 25], [25, 26, 81, 91], [91, 101, 17]], [[1, 5, 9, 13], [13, 15, 1]]]
        @test DT.delete_boundary_node!(boundary_nodes, ((2, 2), 2)) == [[[17, 23, 18, 25], [25, 26, 81, 91], [91, 101, 17]], [[1, 5, 9, 13], [13, 1]]]
    end

    @testset "get_curve_index" begin
        @test DT.get_curve_index(-1) == 1
        @test DT.get_curve_index((5, 3)) == 5
        gv_map = DT.construct_ghost_vertex_map([[[1, 5, 17, 18, 1]], [[23, 29, 31, 33], [33, 107, 101], [101, 99, 85, 23]]])
        @test DT.get_curve_index(gv_map, -1) == 1
        @test DT.get_curve_index(gv_map, -2) == 2
        @test DT.get_curve_index(gv_map, -3) == 2
        @test DT.get_curve_index(gv_map, -4) == 2
    end

    @testset "get_section_index" begin
        @test DT.get_section_index((2, 3)) == 3
        @test DT.get_section_index(4) == 4
        @test DT.get_section_index([1, 2, 3, 4, 5, 1]) == 1
        gv_map = DT.construct_ghost_vertex_map([[[1, 5, 17, 18, 1]], [[23, 29, 31, 33], [33, 107, 101], [101, 99, 85, 23]]])
        @test DT.get_section_index(gv_map, -1) == 1
        @test DT.get_section_index(gv_map, -2) == 1
        @test DT.get_section_index(gv_map, -3) == 2
        @test DT.get_section_index(gv_map, -4) == 3
    end

    @testset "construct_ghost_vertex_ranges" begin
        boundary_nodes = [
            [
                [1, 2, 3, 4], [4, 5, 6, 1],
            ],
            [
                [18, 19, 20, 25, 26, 30],
            ],
            [
                [50, 51, 52, 53, 54, 55], [55, 56, 57, 58], [58, 101, 103, 105, 107, 120], [120, 121, 122, 50],
            ],
        ]
        @test DT.construct_ghost_vertex_ranges(boundary_nodes) == Dict(
            -1 => -2:-1, -2 => -2:-1, -3 => -3:-3,
            -4 => -7:-4, -5 => -7:-4, -6 => -7:-4, -7 => -7:-4,
        )
    end
end

@testset "edges.jl" begin
    @testset "construct_edge" begin
        @test DT.construct_edge(NTuple{2, Int}, 2, 5) == (2, 5)
        e = DT.construct_edge(Vector{Int32}, 5, 15)
        @test e isa Vector{Int32}
        @test e == [5, 15]
    end

    @testset "initial" begin
        @test DT.initial((1, 3)) == 1
        @test DT.initial([2, 5]) == 2
    end

    @testset "terminal" begin
        @test DT.terminal((1, 7)) == 7
        @test DT.terminal([2, 13]) == 13
    end

    @testset "edge_vertices" begin
        @test edge_vertices((1, 5)) == (1, 5)
        @test edge_vertices([23, 50]) == (23, 50)
    end

    @testset "reverse_edge" begin
        @test DT.reverse_edge((17, 3)) == (3, 17)
        @test DT.reverse_edge([1, 2]) == [2, 1]
    end

    @testset "compare_unoriented_edges" begin
        u = (1, 3)
        @test !DT.compare_unoriented_edges(u, (5, 3))
        @test DT.compare_unoriented_edges(u, (1, 3))
        @test DT.compare_unoriented_edges(u, (3, 1))
    end

    @testset "num_edges" begin
        @test num_edges([(1, 2), (3, 4), (1, 5)]) == 3
    end

    @testset "edge_type" begin
        @test DT.edge_type(Set(((1, 2), (2, 3), (17, 5)))) == Tuple{Int64, Int64}
        @test DT.edge_type([[1, 2], [3, 4], [17, 3]]) == Vector{Int64}
    end

    @testset "contains_edge" begin
        E = Set(((1, 3), (17, 3), (1, -1)))
        @test !DT.contains_edge((1, 2), E)
        @test DT.contains_edge((17, 3), E)
        @test !DT.contains_edge(3, 17, E)
        E = [[1, 2], [5, 13], [-1, 1]]
        @test DT.contains_edge(1, 2, E)
    end

    @testset "add_to_edges!" begin
        E = Set(((1, 2), (3, 5)))
        @test DT.add_to_edges!(E, (1, 5)) == Set(((1, 2), (3, 5), (1, 5)))
    end

    @testset "add_edge!" begin
        E = Set(((1, 5), (17, 10), (5, 3)))
        @test DT.add_edge!(E, (3, 2)) === nothing
        @test E == Set(((3, 2), (5, 3), (17, 10), (1, 5)))
        DT.add_edge!(E, (1, -3), (5, 10), (1, -1))
        @test E == Set(((3, 2), (5, 10), (1, -3), (1, -1), (5, 3), (17, 10), (1, 5)))
    end

    @testset "delete_from_edges!" begin
        E = Set(([1, 2], [5, 15], [17, 10], [5, -1]))
        @test DT.delete_from_edges!(E, [5, 15]) == Set(([17, 10], [5, -1], [1, 2]))
    end

    @testset "delete_edge!" begin
        E = Set(([1, 2], [10, 15], [1, -1], [13, 23], [1, 5]))
        @test DT.delete_edge!(E, [10, 15]) === nothing
        @test E == Set(([1, 5], [1, 2], [1, -1], [13, 23]))
        DT.delete_edge!(E, [1, 5], [1, -1])
        @test E == Set(([1, 2], [13, 23]))
    end

    @testset "each_edge" begin
        E = Set(((1, 2), (1, 3), (2, -1)))
        @test each_edge(E) == Set(((1, 2), (1, 3), (2, -1)))
    end

    @testset "random_edge" begin
        E = Set(((1, 2), (10, 15), (23, 20)))
        rng = StableRNG(123)
        for _ in 1:3
            @test DT.random_edge(rng, E) ∈ E # which edge is returned depends on the Set's iteration order
        end
    end
end

@testset "points.jl" begin
    @testset "getx" begin
        @test getx((0.3, 0.7)) == 0.3
    end

    @testset "gety" begin
        @test gety((0.9, 1.3)) == 1.3
    end

    @testset "getxy" begin
        @test getxy([0.9, 23.8]) == (0.9, 23.8)
    end

    @testset "_getx" begin
        @test DT._getx((0.37, 0.7)) === 0.37
        @test DT._getx((0.37f0, 0.7f0)) === 0.3700000047683716
    end

    @testset "_gety" begin
        @test DT._gety((0.5, 0.5)) === 0.5
        @test DT._gety((0.5f0, 0.5f0)) === 0.5
    end

    @testset "_getxy" begin
        @test DT._getxy([0.3, 0.5]) === (0.3, 0.5)
        @test DT._getxy([0.3f0, 0.5f0]) === (0.30000001192092896, 0.5)
    end

    @testset "getpoint" begin
        points = [(0.3, 0.7), (1.3, 5.0), (5.0, 17.0)]
        @test DT.getpoint(points, 2) == (1.3, 5.0)
        points = [0.3 1.3 5.0; 0.7 5.0 17.0]
        @test DT.getpoint(points, 2) == (1.3, 5.0)
        @test DT.getpoint(points, (17.3, 33.0)) == (17.3, 33.0)
    end

    @testset "get_point" begin
        points = [(1.0, 2.0), (3.0, 5.5), (1.7, 10.3), (-5.0, 0.0)]
        @test get_point(points, 1) == (1.0, 2.0)
        @test get_point(points, 1, 2, 3, 4) == ((1.0, 2.0), (3.0, 5.5), (1.7, 10.3), (-5.0, 0.0))
        points = [1.0 3.0 1.7 -5.0; 2.0 5.5 10.3 0.0]
        @test get_point(points, 1) == (1.0, 2.0)
        pts = get_point(points, 1, 2, 3, 4)
        @test pts == ((1.0, 2.0), (3.0, 5.5), (1.7, 10.3), (-5.0, 0.0))
        @test pts isa NTuple{4, Tuple{Float64, Float64}}
    end

    @testset "each_point_index" begin
        @test DT.each_point_index([(1.0, 2.0), (-5.0, 2.0), (2.3, 2.3)]) == Base.OneTo(3)
        @test DT.each_point_index([1.0 -5.0 2.3; 2.0 2.0 2.3]) == Base.OneTo(3)
    end

    @testset "each_point" begin
        @test DT.each_point([(1.0, 2.0), (5.0, 13.0)]) == [(1.0, 2.0), (5.0, 13.0)]
        @test collect(DT.each_point([1.0 5.0 17.7; 5.5 17.7 0.0])) == [[1.0, 5.5], [5.0, 17.7], [17.7, 0.0]]
    end

    @testset "num_points" begin
        @test DT.num_points([(1.0, 1.0), (2.3, 1.5), (0.0, -5.0)]) == 3
        @test DT.num_points([1.0 5.5 10.0 -5.0; 5.0 2.0 0.0 0.0]) == 4
    end

    @testset "points_are_unique" begin
        points = [1.0 2.0 3.0 4.0 5.0; 0.0 5.5 2.0 1.3 17.0]
        @test DT.points_are_unique(points)
        points[:, 4] .= points[:, 1]
        @test !DT.points_are_unique(points)
    end

    @testset "lexicographic_order" begin
        points = [(1.0, 5.0), (0.0, 17.0), (0.0, 13.0), (5.0, 17.3), (3.0, 1.0), (5.0, -2.0)]
        order = DT.lexicographic_order(points)
        @test order == [3, 2, 1, 5, 6, 4]
        @test points[order] == [(0.0, 13.0), (0.0, 17.0), (1.0, 5.0), (3.0, 1.0), (5.0, -2.0), (5.0, 17.3)]
    end

    @testset "push_point!" begin
        points = [(1.0, 3.0), (5.0, 1.0)]
        @test DT.push_point!(points, 2.3, 5.3) == [(1.0, 3.0), (5.0, 1.0), (2.3, 5.3)]
        @test DT.push_point!(points, (17.3, 5.0)) == [(1.0, 3.0), (5.0, 1.0), (2.3, 5.3), (17.3, 5.0)]
    end

    @testset "pop_point!" begin
        points = [(1.0, 2.0), (1.3, 5.3)]
        @test DT.pop_point!(points) == (1.3, 5.3)
        @test points == [(1.0, 2.0)]
    end

    @testset "mean_points" begin
        points = [(1.0, 2.0), (2.3, 5.0), (17.3, 5.3)]
        @test collect(DT.mean_points(points)) ≈ [(1.0 + 2.3 + 17.3) / 3, (2.0 + 5.0 + 5.3) / 3]
        points = [1.0 2.3 17.3; 2.0 5.0 5.3]
        @test collect(DT.mean_points(points)) ≈ [(1.0 + 2.3 + 17.3) / 3, (2.0 + 5.0 + 5.3) / 3]
        @test collect(DT.mean_points(points, (1, 3))) ≈ [(1.0 + 17.3) / 2, (2.0 + 5.3) / 2]
    end

    @testset "set_point!" begin
        points = [(1.0, 3.0), (5.0, 17.0)]
        @test DT.set_point!(points, 1, 0.0, 0.0) == (0.0, 0.0)
        @test points == [(0.0, 0.0), (5.0, 17.0)]
        points = [1.0 2.0 3.0; 4.0 5.0 6.0]
        @test DT.set_point!(points, 2, (17.3, 0.0)) == [17.3, 0.0]
        @test points == [1.0 17.3 3.0; 4.0 0.0 6.0]
    end
end

@testset "triangles.jl" begin
    @testset "construct_triangle" begin
        @test DT.construct_triangle(NTuple{3, Int}, 1, 2, 3) == (1, 2, 3)
        T = DT.construct_triangle(Vector{Int32}, 1, 2, 3)
        @test T isa Vector{Int32}
        @test T == [1, 2, 3]
    end

    @testset "geti" begin
        @test DT.geti((1, 2, 3)) == 1
        @test DT.geti([2, 5, 1]) == 2
    end

    @testset "getj" begin
        @test DT.getj((5, 6, 13)) == 6
        @test DT.getj([10, 19, 21]) == 19
    end

    @testset "getk" begin
        @test DT.getk((1, 2, 3)) == 3
        @test DT.getk([1, 2, 3]) == 3
    end

    @testset "triangle_vertices" begin
        @test triangle_vertices((1, 5, 17)) == (1, 5, 17)
        @test triangle_vertices([5, 18, 23]) == (5, 18, 23)
    end

    @testset "triangle_type" begin
        @test DT.triangle_type(Set{NTuple{3, Int64}}) == Tuple{Int64, Int64, Int64}
        @test DT.triangle_type(Vector{NTuple{3, Int32}}) == Tuple{Int32, Int32, Int32}
        @test DT.triangle_type(Vector{Vector{Int64}}) == Vector{Int64}
    end

    @testset "num_triangles" begin
        T1, T2, T3 = (1, 5, 10), (17, 23, 10), (-1, 10, 5)
        @test num_triangles(Set((T1, T2, T3))) == 3
    end

    @testset "triangle_edges" begin
        @test DT.triangle_edges((1, 2, 3)) == ((1, 2), (2, 3), (3, 1))
        @test DT.triangle_edges(1, 2, 3) == ((1, 2), (2, 3), (3, 1))
    end

    @testset "rotate_triangle" begin
        T = (1, 2, 3)
        @test DT.rotate_triangle(T, 0) == (1, 2, 3)
        @test DT.rotate_triangle(T, 1) == (2, 3, 1)
        @test DT.rotate_triangle(T, 2) == (3, 1, 2)
        @test DT.rotate_triangle(T, 3) == (1, 2, 3)
    end

    @testset "construct_positively_oriented_triangle" begin
        points = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)]
        @test DT.construct_positively_oriented_triangle(NTuple{3, Int}, 1, 2, 3, points) == (2, 1, 3)
        @test DT.construct_positively_oriented_triangle(NTuple{3, Int}, 2, 3, 1, points) == (3, 2, 1)
        @test DT.construct_positively_oriented_triangle(NTuple{3, Int}, 2, 1, 3, points) == (2, 1, 3)
        @test DT.construct_positively_oriented_triangle(NTuple{3, Int}, 3, 2, 1, points) == (3, 2, 1)
        points = [(1.0, 1.0), (2.5, 2.3), (17.5, 23.0), (50.3, 0.0), (-1.0, 2.0), (0.0, 0.0), (5.0, 13.33)]
        @test DT.construct_positively_oriented_triangle(Vector{Int}, 5, 3, 2, points) == [3, 5, 2]
        @test DT.construct_positively_oriented_triangle(Vector{Int}, 7, 1, 2, points) == [7, 1, 2]
        @test DT.construct_positively_oriented_triangle(Vector{Int}, 7, 2, 1, points) == [2, 7, 1]
        @test DT.construct_positively_oriented_triangle(Vector{Int}, 5, 4, 3, points) == [5, 4, 3]
    end

    @testset "compare_triangles" begin
        T1 = (1, 5, 10)
        @test !DT.compare_triangles(T1, (17, 23, 20))
        @test DT.compare_triangles(T1, (5, 10, 1))
        @test DT.compare_triangles(T1, (10, 1, 5))
        @test !DT.compare_triangles(T1, (10, 5, 1))
    end

    @testset "contains_triangle" begin
        V = Set(((1, 2, 3), (4, 5, 6), (7, 8, 9)))
        @test DT.contains_triangle((1, 2, 3), V) == ((1, 2, 3), true)
        @test DT.contains_triangle((2, 3, 1), V) == ((1, 2, 3), true)
        @test DT.contains_triangle((10, 18, 9), V) == ((10, 18, 9), false)
        @test DT.contains_triangle(9, 7, 8, V) == ((7, 8, 9), true)
    end

    @testset "sort_triangle" begin
        @test DT.sort_triangle((1, 5, 3)) == (5, 3, 1)
        @test DT.sort_triangle((1, -1, 2)) == (2, 1, -1)
        @test DT.sort_triangle((3, 2, 1)) == (3, 2, 1)
    end

    @testset "add_to_triangles!" begin
        T = Set(((1, 2, 3), (17, 8, 9)))
        @test DT.add_to_triangles!(T, (1, 5, 12)) == Set(((1, 5, 12), (1, 2, 3), (17, 8, 9)))
        @test DT.add_to_triangles!(T, (-1, 3, 6)) == Set(((1, 5, 12), (1, 2, 3), (17, 8, 9), (-1, 3, 6)))
    end

    @testset "add_triangle!(T, V...)" begin
        T = Set(((1, 2, 3), (4, 5, 6)))
        add_triangle!(T, (7, 8, 9))
        add_triangle!(T, (10, 11, 12), (13, 14, 15))
        add_triangle!(T, 16, 17, 18)
        @test T == Set(((7, 8, 9), (10, 11, 12), (4, 5, 6), (13, 14, 15), (16, 17, 18), (1, 2, 3)))
    end

    @testset "delete_from_triangles!" begin
        V = Set(((1, 2, 3), (4, 5, 6), (7, 8, 9)))
        @test DT.delete_from_triangles!(V, (4, 5, 6)) == Set(((7, 8, 9), (1, 2, 3)))
        @test DT.delete_from_triangles!(V, (9, 7, 8)) == Set(((1, 2, 3),))
    end

    @testset "delete_triangle!(V, T...)" begin
        V = Set(((1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12), (13, 14, 15)))
        @test delete_triangle!(V, (6, 4, 5)) == Set(((7, 8, 9), (10, 11, 12), (13, 14, 15), (1, 2, 3)))
        @test delete_triangle!(V, (10, 11, 12), (1, 2, 3)) == Set(((7, 8, 9), (13, 14, 15)))
        @test delete_triangle!(V, 8, 9, 7) == Set(((13, 14, 15),))
    end

    @testset "each_triangle" begin
        T = Set(((1, 2, 3), (-1, 5, 10), (17, 13, 18)))
        @test each_triangle(T) == Set(((-1, 5, 10), (1, 2, 3), (17, 13, 18)))
        T = [[1, 2, 3], [10, 15, 18], [1, 5, 6]]
        @test each_triangle(T) == [[1, 2, 3], [10, 15, 18], [1, 5, 6]]
    end

    @testset "compare_triangle_collections" begin
        T = Set(((1, 2, 3), (4, 5, 6), (7, 8, 9)))
        V = [[2, 3, 1], [4, 5, 6], [9, 7, 8]]
        @test DT.compare_triangle_collections(T, V)
        V[1] = [17, 19, 20]
        @test !DT.compare_triangle_collections(T, V)
        V = [[1, 2, 3], [8, 9, 7]]
        @test !DT.compare_triangle_collections(T, V)
    end

    @testset "sort_triangles" begin
        T = Set(((1, 3, 2), (5, 2, 3), (10, 1, 13), (-1, 10, 12), (10, 1, 17), (5, 8, 2)))
        @test DT.sort_triangles(T) == Set(((13, 10, 1), (3, 5, 2), (10, 12, -1), (5, 8, 2), (17, 10, 1), (3, 2, 1)))
    end
end

@testset "utils.jl" begin
    @testset "number_type" begin
        @test DT.number_type([1, 2, 3]) == Int64
        @test DT.number_type((1, 2, 3)) == Int64
        @test DT.number_type([1.0 2.0 3.0; 4.0 5.0 6.0]) == Float64
        @test DT.number_type([[[1, 2, 3, 4, 5, 1]], [[6, 8, 9], [9, 10, 11], [11, 12, 6]]]) == Int64
        @test DT.number_type((1.0f0, 2.0f0)) == Float32
        @test DT.number_type(Vector{Float64}) == Float64
        @test DT.number_type(Vector{Vector{Float64}}) == Float64
        @test DT.number_type(NTuple{2, Float64}) == Float64
    end

    @testset "get_ghost_vertex" begin
        @test DT.get_ghost_vertex(1, 7, -2) == -2
        @test DT.get_ghost_vertex(-1, 2, 3) == -1
        @test DT.get_ghost_vertex(1, 5, 10) == 10
        @test DT.get_ghost_vertex(1, -1) == -1
        @test DT.get_ghost_vertex(-5, 2) == -5
    end

    @testset "is_true" begin
        @test DT.is_true(true)
        @test !DT.is_true(false)
        @test DT.is_true(Val(true))
        @test !DT.is_true(Val(false))
    end

    @testset "get_ordinal_suffix" begin
        @test DT.get_ordinal_suffix(1) == "st"
        @test DT.get_ordinal_suffix(2) == "nd"
        @test DT.get_ordinal_suffix(3) == "rd"
        @test DT.get_ordinal_suffix(4) == "th"
        @test DT.get_ordinal_suffix(5) == "th"
        @test DT.get_ordinal_suffix(6) == "th"
        @test DT.get_ordinal_suffix(11) == "th"
        @test DT.get_ordinal_suffix(15) == "th"
        @test DT.get_ordinal_suffix(100) == "th"
    end
end
