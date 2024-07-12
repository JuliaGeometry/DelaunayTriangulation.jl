using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using ..DelaunayTriangulation: Triangulation
using StableRNGs
using StaticArrays


_x, _y = complicated_geometry()
global x = _x
global y = _y
boundary_nodes, points = convert_boundary_points_to_indices(x, y)
rng = StableRNG(9881)
global tri = triangulate(points; boundary_nodes, rng, delete_ghosts = false)
A = get_area(tri)
refine!(tri; max_area = 1.0e-3A, rng, use_circumcenter = true)
boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
rng = StableRNG(9881)
global tri_2 = triangulate(points; boundary_nodes, rng, delete_ghosts = false)
A = get_area(tri_2)
refine!(tri_2; max_area = 1.0e-3A, rng, use_circumcenter = true)
global i, j, k = 451, 307, 227
global u, v, w = 420, 417, 418
global T = (u, v, w)


@testset "Simple add, no connections" begin
    DT.add_triangle!(tri, i, j, k)
    DT.add_triangle!(tri, T)
    for ((a, b, c)) in ((i, j, k), (u, v, w))
        @test DT.get_adjacent(tri, a, b) == c
        @test DT.get_adjacent(tri, b, c) == a
        @test DT.get_adjacent(tri, c, a) == b
        @test (a, b) ∈ DT.get_adjacent2vertex(tri, c)
        @test (b, c) ∈ DT.get_adjacent2vertex(tri, a)
        @test (c, a) ∈ DT.get_adjacent2vertex(tri, b)
        @test a ∈ DT.get_neighbours(tri, b) && a ∈ DT.get_neighbours(tri, c)
        @test b ∈ DT.get_neighbours(tri, a) && b ∈ DT.get_neighbours(tri, c)
        @test c ∈ DT.get_neighbours(tri, a) && c ∈ DT.get_neighbours(tri, b)
        @test (a, b, c) ∈ tri.triangles
    end
end


@testset "Adding to empty triangulation, ghost edges included" begin
    pts = rand(2, 500)
    u, v, w = 37, 38, 73
    pts[:, u] .= 0.0
    pts[:, v] .= [1.0, 0.0]
    pts[:, w] .= [0.0, 1.0]
    tri = Triangulation(pts)
    DT.add_triangle!(tri, (u, v, w); update_ghost_edges = true)
    BI = DT.𝒢
    @test get_triangles(tri) == Set(((u, v, w), (w, v, BI), (v, u, BI), (u, w, BI)))
    @test get_adjacent(get_adjacent(tri)) == Dict(
        (u, v) => w,
        (v, w) => u,
        (w, u) => v,
        (w, v) => BI,
        (v, BI) => w,
        (BI, w) => v,
        (v, u) => BI,
        (u, BI) => v,
        (BI, v) => u,
        (u, w) => BI,
        (w, BI) => u,
        (BI, u) => w,
    )
    @test get_adjacent2vertex(get_adjacent2vertex(tri)) ==
        Dict(
        BI => Set{NTuple{2, Int}}([(w, v), (v, u), (u, w)]),
        u => Set{NTuple{2, Int}}([(v, w), (BI, v), (w, BI)]),
        v => Set{NTuple{2, Int}}([(w, u), (BI, w), (u, BI)]),
        w => Set{NTuple{2, Int}}([(u, v), (v, BI), (BI, u)]),
    )
    @test get_graph(tri).edges ==
        Set([(BI, v), (u, w), (u, v), (BI, w), (BI, u), (v, w)])
    @test get_graph(tri).vertices == Set([u, v, w, BI])
end


@testset "Smaller example" begin
    tri = example_triangulation()
    i, j, k = 1, 3, 7
    DT.add_triangle!(tri, i, j, k)


    @testset "Adding an interior triangle" begin
        true_T = Set{NTuple{3, Int}}(
            [
                (3, 2, 5),
                (1, 3, 7),
                (4, 1, 5),
                (6, 3, 1),
                (4, 6, 1),
                (5, 1, 3),
            ],
        )
        true_adj =
            Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (3, 5) => 1,
            (4, 5) => DT.𝒢, (5, 2) => DT.𝒢,
            (2, 3) => DT.𝒢, (3, 6) => DT.𝒢,
            (6, 4) => DT.𝒢,
        )
        true_adj2v = Dict(
            DT.𝒢 => Set{NTuple{2, Int}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
            1 => Set{NTuple{2, Int}}([(3, 7), (3, 5), (6, 3), (5, 4), (4, 6)]),
            2 => Set{NTuple{2, Int}}([(5, 3)]),
            3 => Set{NTuple{2, Int}}([(2, 5), (5, 1), (7, 1), (1, 6)]),
            4 => Set{NTuple{2, Int}}([(1, 5), (6, 1)]),
            5 => Set{NTuple{2, Int}}([(4, 1), (1, 3), (3, 2)]),
            6 => Set{NTuple{2, Int}}([(3, 1), (1, 4)]),
            7 => Set{NTuple{2, Int}}([(1, 3)]),
        )
        true_DG = [
            0 0 1 1 1 1 1 0
            0 0 0 1 1 1 1 1
            1 0 0 1 0 1 0 0
            1 1 1 0 0 1 1 1
            1 1 0 0 0 1 1 0
            1 1 1 1 1 0 0 0
            1 1 0 1 1 0 0 0
            0 1 0 1 0 0 0 0
        ]
        true_DG = _make_graph_from_adjacency(true_DG, Dict(1:8 .=> [-1, (1:7)...]))
        @test get_triangles(tri) == true_T
        @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
        @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
        @test get_graph(tri) == true_DG
    end


    @testset "Adding a boundary triangle with one boundary edge" begin
        for (i, j, k) in ((5, 2, 8), (2, 8, 5), (8, 5, 2))
            local true_T, true_adj, true_adj2v, true_DG
            _tri = deepcopy(tri)
            DT.add_triangle!(_tri, i, j, k)
            true_T = Set{NTuple{3, Int}}(
                [
                    (3, 2, 5),
                    (1, 3, 7),
                    (4, 1, 5),
                    (6, 3, 1),
                    (4, 6, 1),
                    (5, 1, 3),
                    (i, j, k),
                ],
            )
            true_adj =
                Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
                (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
                (5, 1) => 3, (3, 5) => 1,
                (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
                (4, 5) => DT.𝒢, (5, 8) => DT.𝒢,
                (8, 2) => DT.𝒢, (2, 3) => DT.𝒢,
                (3, 6) => DT.𝒢, (6, 4) => DT.𝒢,
            )
            true_adj2v = Dict(
                DT.𝒢 => Set{NTuple{2, Int}}([(4, 5), (5, 8), (8, 2), (2, 3), (3, 6), (6, 4)]),
                1 => Set{NTuple{2, Int}}([(3, 7), (3, 5), (6, 3), (5, 4), (4, 6)]),
                2 => Set{NTuple{2, Int}}([(5, 3), (8, 5)]),
                3 => Set{NTuple{2, Int}}([(2, 5), (5, 1), (7, 1), (1, 6)]),
                4 => Set{NTuple{2, Int}}([(1, 5), (6, 1)]),
                5 => Set{NTuple{2, Int}}([(4, 1), (1, 3), (3, 2), (2, 8)]),
                6 => Set{NTuple{2, Int}}([(3, 1), (1, 4)]),
                7 => Set{NTuple{2, Int}}([(1, 3)]),
                8 => Set{NTuple{2, Int}}([(5, 2)]),
            )
            true_DG = _make_graph_from_adjacency(
                [
                    0 0 1 1 1 1 1 0 1
                    0 0 0 1 1 1 1 1 0
                    1 0 0 1 0 1 0 0 1
                    1 1 1 0 0 1 1 1 0
                    1 1 0 0 0 1 1 0 0
                    1 1 1 1 1 0 0 0 1
                    1 1 0 1 1 0 0 0 0
                    0 1 0 1 0 0 0 0 0
                    1 0 1 0 0 1 0 0 0
                ], Dict(1:9 .=> [-1, (1:8)...]),
            )
            DT.clear_empty_features!(_tri)
            @test get_triangles(_tri) == true_T
            @test (get_adjacent ∘ get_adjacent)(_tri) == true_adj
            @test (get_adjacent2vertex ∘ get_adjacent2vertex)(_tri) == true_adj2v
            @test get_graph(_tri) == true_DG
        end
        DT.add_triangle!(tri, 5, 2, 8) # Get an actual copy for later, test is above
    end


    @testset "Adding a boundary triangle with two boundary edges" begin
        for (i, j, k) in ((3, 6, 2), (6, 2, 3), (2, 3, 6))
            local true_T, true_adj, true_adj2v, true_DG
            _tri = deepcopy(tri)
            DT.add_triangle!(_tri, i, j, k)
            true_T = Set{NTuple{3, Int}}(
                [
                    (3, 2, 5),
                    (1, 3, 7),
                    (4, 1, 5),
                    (6, 3, 1),
                    (4, 6, 1),
                    (5, 1, 3),
                    (5, 2, 8),
                    (i, j, k),
                ],
            )
            true_adj =
                Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
                (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
                (5, 1) => 3, (3, 5) => 1,
                (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
                (2, 3) => 6, (3, 6) => 2, (6, 2) => 3,
                (4, 5) => DT.𝒢, (5, 8) => DT.𝒢,
                (8, 2) => DT.𝒢, (2, 6) => DT.𝒢,
                (6, 4) => DT.𝒢,
            )
            true_adj2v = Dict(
                DT.𝒢 => Set{NTuple{2, Int}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
                1 => Set{NTuple{2, Int}}([(3, 7), (3, 5), (6, 3), (5, 4), (4, 6)]),
                2 => Set{NTuple{2, Int}}([(5, 3), (8, 5), (3, 6)]),
                3 => Set{NTuple{2, Int}}([(2, 5), (5, 1), (7, 1), (1, 6), (6, 2)]),
                4 => Set{NTuple{2, Int}}([(1, 5), (6, 1)]),
                5 => Set{NTuple{2, Int}}([(4, 1), (1, 3), (3, 2), (2, 8)]),
                6 => Set{NTuple{2, Int}}([(3, 1), (1, 4), (2, 3)]),
                7 => Set{NTuple{2, Int}}([(1, 3)]),
                8 => Set{NTuple{2, Int}}([(5, 2)]),
            )
            true_DG = _make_graph_from_adjacency(
                [
                    0 0 1 0 1 1 1 0 1
                    0 0 0 1 1 1 1 1 0
                    1 0 0 1 0 1 1 0 1
                    0 1 1 0 0 1 1 1 0
                    1 1 0 0 0 1 1 0 0
                    1 1 1 1 1 0 0 0 1
                    1 1 1 1 1 0 0 0 0
                    0 1 0 1 0 0 0 0 0
                    1 0 1 0 0 1 0 0 0
                ], Dict(1:9 .=> [-1, (1:8)...]),
            )
            DT.clear_empty_features!(_tri)
            @test get_triangles(_tri) == true_T
            @test (get_adjacent ∘ get_adjacent)(_tri) == true_adj
            @test (get_adjacent2vertex ∘ get_adjacent2vertex)(_tri) == true_adj2v
            @test get_graph(_tri) == true_DG
        end
    end
end


@testset "Empty triangulation" begin
    tri = example_empty_triangulation()
    true_T = Set{NTuple{3, Int}}([(1, 2, 3)])
    true_adj =
        Dict(
        (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
        (3, 2) => DT.𝒢, (2, 1) => DT.𝒢, (1, 3) => DT.𝒢,
    )
    true_adj2v = Dict(
        DT.𝒢 => Set{NTuple{2, Int}}([(3, 2), (2, 1), (1, 3)]),
        1 => Set{NTuple{2, Int}}([(2, 3)]),
        2 => Set{NTuple{2, Int}}([(3, 1)]),
        3 => Set{NTuple{2, Int}}([(1, 2)]),
    )
    true_DG = _make_graph_from_adjacency(
        [
            0 1 1 1
            1 0 1 1
            1 1 0 1
            1 1 1 0
        ], Dict(1:4 .=> [-1, 1, 2, 3]),
    )
    DT.add_triangle!(tri, 1, 2, 3)
    @test get_triangles(tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
    @test get_graph(tri) == true_DG
end


@testset "Testing the boundary addition cases for updating ghost edges" begin
    p1 = @SVector[0.0, 0.0]
    p2 = @SVector[1.0, 0.0]
    p3 = @SVector[0.0, 1.0]
    pts = [p1, p2, p3]
    tri = Triangulation(pts)
    DT.add_triangle!(tri, 1, 2, 3; update_ghost_edges = true)
    true_T = Set{NTuple{3, Int}}(
        [
            (1, 2, 3),
            (2, 1, DT.𝒢),
            (1, 3, DT.𝒢),
            (3, 2, DT.𝒢),
        ],
    )
    true_adj = DT.Adjacent(
        Dict(
            (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (2, 1) => DT.𝒢, (1, DT.𝒢) => 2, (DT.𝒢, 2) => 1,
            (1, 3) => DT.𝒢, (3, DT.𝒢) => 1, (DT.𝒢, 1) => 3,
            (3, 2) => DT.𝒢, (2, DT.𝒢) => 3, (DT.𝒢, 3) => 2,
        ),
    )
    true_adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.𝒢 => Set{NTuple{2, Int}}([(2, 1), (1, 3), (3, 2)]),
            1 => Set{NTuple{2, Int}}([(2, 3), (DT.𝒢, 2), (3, DT.𝒢)]),
            2 => Set{NTuple{2, Int}}([(3, 1), (1, DT.𝒢), (DT.𝒢, 3)]),
            3 => Set{NTuple{2, Int}}([(1, 2), (DT.𝒢, 1), (2, DT.𝒢)]),
        ),
    )
    true_DG = DT.Graph{Int}()
    DT.add_neighbour!(true_DG, DT.𝒢, [1, 2, 3]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3]...)
    DT.add_neighbour!(true_DG, 3, [1, 2]...)
    @test get_triangles(tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
    p4 = @SVector[1.7, 1.7]
    push!(pts, p4)
    DT.add_triangle!(tri, 3, 2, 4; update_ghost_edges = true)
    true_T = Set{NTuple{3, Int}}(
        [
            (1, 2, 3),
            (2, 1, DT.𝒢),
            (1, 3, DT.𝒢),
            (3, 4, DT.𝒢),
            (4, 2, DT.𝒢),
            (3, 2, 4),
        ],
    )
    true_adj = DT.Adjacent(
        Dict(
            (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (2, 1) => DT.𝒢, (1, DT.𝒢) => 2, (DT.𝒢, 2) => 1,
            (1, 3) => DT.𝒢, (3, DT.𝒢) => 1, (DT.𝒢, 1) => 3,
            (3, 4) => DT.𝒢, (4, DT.𝒢) => 3, (DT.𝒢, 3) => 4,
            (4, 2) => DT.𝒢, (2, DT.𝒢) => 4, (DT.𝒢, 4) => 2,
            (3, 2) => 4, (2, 4) => 3, (4, 3) => 2,
        ),
    )
    true_adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.𝒢 => Set{NTuple{2, Int}}([(2, 1), (1, 3), (3, 4), (4, 2)]),
            1 => Set{NTuple{2, Int}}([(2, 3), (DT.𝒢, 2), (3, DT.𝒢)]),
            2 => Set{NTuple{2, Int}}([(3, 1), (1, DT.𝒢), (DT.𝒢, 4), (4, 3)]),
            3 => Set{NTuple{2, Int}}([(1, 2), (DT.𝒢, 1), (4, DT.𝒢), (2, 4)]),
            4 => Set{NTuple{2, Int}}([(DT.𝒢, 3), (2, DT.𝒢), (3, 2)]),
        ),
    )
    true_DG = DT.Graph{Int}()
    DT.add_neighbour!(true_DG, DT.𝒢, [1, 3, 4, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4]...)
    DT.add_neighbour!(true_DG, 4, [2, 3]...)
    @test get_triangles(tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
    p5 = @SVector[1.0, 3.0]
    p6 = @SVector[3.0, 1.0]
    push!(pts, p5, p6)
    DT.add_triangle!(tri, 3, 4, 5; update_ghost_edges = true)
    DT.add_triangle!(tri, 4, 2, 6; update_ghost_edges = true)
    true_T = Set{NTuple{3, Int}}(
        [
            (1, 2, 3),
            (3, 2, 4),
            (3, 4, 5),
            (4, 2, 6),
            (2, 1, DT.𝒢),
            (1, 3, DT.𝒢),
            (3, 5, DT.𝒢),
            (5, 4, DT.𝒢),
            (4, 6, DT.𝒢),
            (6, 2, DT.𝒢),
        ],
    )
    true_adj = DT.Adjacent(
        Dict(
            (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (3, 2) => 4, (2, 4) => 3, (4, 3) => 2,
            (3, 4) => 5, (4, 5) => 3, (5, 3) => 4,
            (4, 2) => 6, (2, 6) => 4, (6, 4) => 2,
            (2, 1) => DT.𝒢, (1, DT.𝒢) => 2, (DT.𝒢, 2) => 1,
            (1, 3) => DT.𝒢, (3, DT.𝒢) => 1, (DT.𝒢, 1) => 3,
            (3, 5) => DT.𝒢, (5, DT.𝒢) => 3, (DT.𝒢, 3) => 5,
            (5, 4) => DT.𝒢, (4, DT.𝒢) => 5, (DT.𝒢, 5) => 4,
            (4, 6) => DT.𝒢, (6, DT.𝒢) => 4, (DT.𝒢, 4) => 6,
            (6, 2) => DT.𝒢, (2, DT.𝒢) => 6, (DT.𝒢, 6) => 2,
        ),
    )
    true_adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.𝒢 => Set{NTuple{2, Int}}([(1, 3), (3, 5), (5, 4), (4, 6), (6, 2), (2, 1)]),
            1 => Set{NTuple{2, Int}}([(2, 3), (3, DT.𝒢), (DT.𝒢, 2)]),
            2 => Set{NTuple{2, Int}}([(3, 1), (4, 3), (6, 4), (1, DT.𝒢), (DT.𝒢, 6)]),
            3 => Set{NTuple{2, Int}}([(1, 2), (2, 4), (4, 5), (DT.𝒢, 1), (5, DT.𝒢)]),
            4 => Set{NTuple{2, Int}}([(3, 2), (5, 3), (2, 6), (DT.𝒢, 5), (6, DT.𝒢)]),
            5 => Set{NTuple{2, Int}}([(3, 4), (DT.𝒢, 3), (4, DT.𝒢)]),
            6 => Set{NTuple{2, Int}}([(4, 2), (DT.𝒢, 4), (2, DT.𝒢)]),
        ),
    )
    true_DG = DT.Graph{Int}()
    DT.add_neighbour!(true_DG, DT.𝒢, [1, 3, 5, 4, 6, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4, 6]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4, 5]...)
    DT.add_neighbour!(true_DG, 4, [2, 3, 5, 6]...)
    DT.add_neighbour!(true_DG, 5, [3, 4]...)
    DT.add_neighbour!(true_DG, 6, [2, 4]...)
    @test get_triangles(tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
    DT.add_triangle!(tri, 5, 4, 6; update_ghost_edges = true)
    true_T = Set{NTuple{3, Int}}(
        [
            (1, 2, 3),
            (3, 2, 4),
            (3, 4, 5),
            (4, 2, 6),
            (5, 4, 6),
            (2, 1, DT.𝒢),
            (1, 3, DT.𝒢),
            (3, 5, DT.𝒢),
            (5, 6, DT.𝒢),
            (6, 2, DT.𝒢),
        ],
    )
    true_adj = DT.Adjacent(
        Dict(
            (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (3, 2) => 4, (2, 4) => 3, (4, 3) => 2,
            (3, 4) => 5, (4, 5) => 3, (5, 3) => 4,
            (4, 2) => 6, (2, 6) => 4, (6, 4) => 2,
            (5, 4) => 6, (4, 6) => 5, (6, 5) => 4,
            (2, 1) => DT.𝒢, (1, DT.𝒢) => 2, (DT.𝒢, 2) => 1,
            (1, 3) => DT.𝒢, (3, DT.𝒢) => 1, (DT.𝒢, 1) => 3,
            (3, 5) => DT.𝒢, (5, DT.𝒢) => 3, (DT.𝒢, 3) => 5,
            (5, 6) => DT.𝒢, (6, DT.𝒢) => 5, (DT.𝒢, 5) => 6,
            (6, 2) => DT.𝒢, (2, DT.𝒢) => 6, (DT.𝒢, 6) => 2,
        ),
    )
    true_adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.𝒢 => Set{NTuple{2, Int}}([(1, 3), (3, 5), (5, 6), (6, 2), (2, 1)]),
            1 => Set{NTuple{2, Int}}([(2, 3), (3, DT.𝒢), (DT.𝒢, 2)]),
            2 => Set{NTuple{2, Int}}([(3, 1), (4, 3), (6, 4), (1, DT.𝒢), (DT.𝒢, 6)]),
            3 => Set{NTuple{2, Int}}([(1, 2), (2, 4), (4, 5), (DT.𝒢, 1), (5, DT.𝒢)]),
            4 => Set{NTuple{2, Int}}([(3, 2), (2, 6), (6, 5), (5, 3)]),
            5 => Set{NTuple{2, Int}}([(3, 4), (4, 6), (6, DT.𝒢), (DT.𝒢, 3)]),
            6 => Set{NTuple{2, Int}}([(4, 2), (5, 4), (DT.𝒢, 5), (2, DT.𝒢)]),
        ),
    )
    true_DG = DT.Graph{Int}()
    DT.add_neighbour!(true_DG, DT.𝒢, [1, 3, 5, 6, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4, 6]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4, 5]...)
    DT.add_neighbour!(true_DG, 4, [2, 3, 5, 6]...)
    DT.add_neighbour!(true_DG, 5, [3, 4, 6]...)
    DT.add_neighbour!(true_DG, 6, [2, 4, 5]...)
    @test get_triangles(tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
    DT.convex_hull!(tri; reconstruct = false)
    DT.compute_representative_points!(tri)
    @test validate_triangulation(tri)
end


@testset "Larger example" begin
    p1 = @SVector[-3.32, 3.53]
    p2 = @SVector[-5.98, 2.17]
    p3 = @SVector[-6.36, -1.55]
    p4 = @SVector[-2.26, -4.31]
    p5 = @SVector[6.34, -3.23]
    p6 = @SVector[-3.24, 1.01]
    p7 = @SVector[0.14, -1.51]
    p8 = @SVector[0.2, 1.25]
    p9 = @SVector[1.0, 4.0]
    p10 = @SVector[4.74, 2.21]
    p11 = @SVector[2.32, -0.27]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
    tri = Triangulation(pts)
    for (i, j, k) in (
            (1, 2, 6),
            (1, 6, 8),
            (9, 1, 8),
            (9, 8, 10),
            (10, 8, 11),
            (8, 7, 11),
            (8, 6, 7),
            (6, 2, 3),
            (6, 3, 4),
            (6, 4, 7),
            (7, 4, 5),
            (11, 7, 5),
            (10, 11, 5),
        )
        DT.add_triangle!(tri, i, j, k; update_ghost_edges = true)
    end
    true_T = Set{NTuple{3, Int}}(
        [
            (1, 6, 8),
            (1, 8, 9),
            (9, 8, 10),
            (10, 8, 11),
            (10, 11, 5),
            (11, 7, 5),
            (7, 4, 5),
            (7, 6, 4),
            (6, 3, 4),
            (6, 2, 3),
            (1, 2, 6),
            (8, 6, 7),
            (8, 7, 11),
            (2, 1, DT.𝒢),
            (1, 9, DT.𝒢),
            (9, 10, DT.𝒢),
            (10, 5, DT.𝒢),
            (5, 4, DT.𝒢),
            (4, 3, DT.𝒢),
            (3, 2, DT.𝒢),
        ],
    )
    true_adj = DT.Adjacent(
        Dict(
            (1, 6) => 8, (6, 8) => 1, (8, 1) => 6,
            (1, 8) => 9, (8, 9) => 1, (9, 1) => 8,
            (9, 8) => 10, (8, 10) => 9, (10, 9) => 8,
            (10, 8) => 11, (8, 11) => 10, (11, 10) => 8,
            (10, 11) => 5, (11, 5) => 10, (5, 10) => 11,
            (11, 7) => 5, (7, 5) => 11, (5, 11) => 7,
            (7, 4) => 5, (4, 5) => 7, (5, 7) => 4,
            (7, 6) => 4, (6, 4) => 7, (4, 7) => 6,
            (6, 3) => 4, (3, 4) => 6, (4, 6) => 3,
            (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
            (1, 2) => 6, (2, 6) => 1, (6, 1) => 2,
            (8, 6) => 7, (6, 7) => 8, (7, 8) => 6,
            (8, 7) => 11, (7, 11) => 8, (11, 8) => 7,
            (2, 1) => DT.𝒢, (1, DT.𝒢) => 2, (DT.𝒢, 2) => 1,
            (1, 9) => DT.𝒢, (9, DT.𝒢) => 1, (DT.𝒢, 1) => 9,
            (9, 10) => DT.𝒢, (10, DT.𝒢) => 9, (DT.𝒢, 9) => 10,
            (10, 5) => DT.𝒢, (5, DT.𝒢) => 10, (DT.𝒢, 10) => 5,
            (5, 4) => DT.𝒢, (4, DT.𝒢) => 5, (DT.𝒢, 5) => 4,
            (4, 3) => DT.𝒢, (3, DT.𝒢) => 4, (DT.𝒢, 4) => 3,
            (3, 2) => DT.𝒢, (2, DT.𝒢) => 3, (DT.𝒢, 3) => 2,
        ),
    )
    true_adj2v = DT.Adjacent2Vertex{Int, Set{NTuple{2, Int}}}()
    for (ij, k) in true_adj.adjacent
        DT.add_adjacent2vertex!(true_adj2v, k, ij)
    end
    true_DG = DT.Graph{Int}()
    DT.add_neighbour!(true_DG, DT.𝒢, 1, 9, 10, 5, 4, 3, 2)
    DT.add_neighbour!(true_DG, 1, 2, 6, 8, 9)
    DT.add_neighbour!(true_DG, 2, 1, 6, 3)
    DT.add_neighbour!(true_DG, 3, 2, 6, 4)
    DT.add_neighbour!(true_DG, 4, 6, 3, 7, 5)
    DT.add_neighbour!(true_DG, 5, 10, 11, 7, 4)
    DT.add_neighbour!(true_DG, 6, 1, 2, 3, 4, 7, 8)
    DT.add_neighbour!(true_DG, 7, 8, 6, 4, 5, 11)
    DT.add_neighbour!(true_DG, 8, 1, 6, 7, 11, 10, 9)
    DT.add_neighbour!(true_DG, 9, 1, 8, 10)
    DT.add_neighbour!(true_DG, 10, 9, 8, 11, 5)
    DT.add_neighbour!(true_DG, 11, 10, 8, 7, 5)
    @test DT.compare_triangle_collections(get_triangles(tri), true_T)
    @test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
    @test get_graph(tri) == true_DG
    DT.convex_hull!(tri; reconstruct = false)
    DT.compute_representative_points!(tri)
    @test validate_triangulation(tri)
end