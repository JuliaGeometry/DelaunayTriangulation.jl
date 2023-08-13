using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Random
using StableRNGs
using StatsBase

include("../helper_functions.jl")

tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)
DT.compute_representative_points!(tri)
rep = DT.get_representative_point_list(tri)
pts = get_points(tri)
adj = get_adjacent(tri)
adj2v = get_adjacent2vertex(tri)
boundary_index_ranges = get_boundary_index_ranges(tri)
boundary_map = get_boundary_map(tri)
graph = get_graph(tri)
boundary_nodes = get_boundary_nodes(tri)

@testset "Specific point and triangle examples" begin
    q1 = (2.0, 4.0)
    V1 = [(index_map["a"], index_map["s"], index_map["h"])]
    q2 = (1.1803475764382, 17.4209559058887)
    V2 = [(index_map["g"], index_map["b1"], index_map["i"])]
    q3 = (14.2726546767168, 18.7727571005127)
    V3 = [(index_map["f"], index_map["v"], index_map["e"])]
    q4 = (12.0, 2.0)
    V4 = [(index_map["t"], index_map["p"], index_map["o"]),
        (index_map["t"], index_map["b"], index_map["p"])]
    q5 = (13.0, 0.0)
    V5 = [(index_map["b"], index_map["c"], index_map["p"]),
        (index_map["c"], index_map["b"], DT.BoundaryIndex)]
    q6 = [6.0819947177817, 8.90410894457]
    V6 = [(index_map["j"], index_map["k"], DT.BoundaryIndex - 1)]
    q7 = (15.8, 9.2)
    V7 = [(index_map["m"], DT.BoundaryIndex - 2, index_map["r"]),
        (index_map["m"], DT.BoundaryIndex - 3, index_map["r"])]
    q8 = (14.5391680629781, 8.5135910023477)
    V8 = [(index_map["m"], index_map["n"], DT.BoundaryIndex - 2),
        (index_map["m"], index_map["n"], DT.BoundaryIndex - 3)]
    q9 = (22.0, 8.0)
    V9 = [(index_map["d"], index_map["c"], DT.BoundaryIndex)]
    q10 = (6.0, 11.0)
    V10 = [(index_map["j"], index_map["k"], DT.BoundaryIndex - 1),
        (index_map["k"], index_map["ℓ"], DT.BoundaryIndex - 1),
        (index_map["ℓ"], index_map["i"], DT.BoundaryIndex - 1),
        (index_map["i"], index_map["j"], DT.BoundaryIndex - 1)]
    q11 = (6.0819947177817, 8.90410894457)
    V11 = [(index_map["j"], index_map["k"], DT.BoundaryIndex - 1)]
    q12 = (-6.6465702688004, 19.4798146355492)
    V12 = [(index_map["h"], index_map["g"], DT.BoundaryIndex)]
    q13 = (-20.0, 10.0)
    V13 = [(index_map["h"], index_map["g"], DT.BoundaryIndex),
        (index_map["a"], index_map["h"], DT.BoundaryIndex)]
    q14 = (20.0, 6.0)
    V14 = [(index_map["c"], index_map["d"], index_map["q"]),
        (index_map["c"], index_map["d"], DT.BoundaryIndex)]
    allq = [q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14]
    allV = [V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14]
    for _ in 1:36
        for k in each_point_index(pts)
            for i in eachindex(allq)
                @show i, k
                T1 = jump_and_march(tri, allq[i]; k)
                T2 = jump_and_march(tri, allq[i])
                @test DT.is_positively_oriented(DT.triangle_orientation(tri, T1))
                @test DT.is_positively_oriented(DT.triangle_orientation(tri, T2))
                @inferred jump_and_march(tri, allq[i]; k)
                @inferred jump_and_march(tri, allq[i])
                for V in [T1, T2]
                    if length(allV[i]) == 1
                        @test DT.compare_triangles(V, allV[i][1]) && DT.is_inside(DT.point_position_relative_to_triangle(tri, V, allq[i]))
                    elseif length(allV[i]) == 2
                        @test (DT.compare_triangles(V, allV[i][1]) || DT.compare_triangles(V, allV[i][2])) && (DT.is_on(DT.point_position_relative_to_triangle(tri, V, allq[i])) || (DT.is_ghost_triangle(V) && DT.is_inside(DT.point_position_relative_to_triangle(tri, V, allq[i]))))
                    else
                        bool1 = any(j -> DT.compare_triangles(V, allV[i][j]), eachindex(allV[i]))
                        bool2 = (DT.is_on(DT.point_position_relative_to_triangle(tri, V, allq[i])) || (DT.is_ghost_triangle(V) && DT.is_inside(DT.point_position_relative_to_triangle(tri, V, allq[i]))))
                        @test bool1 && bool2
                    end
                end
            end
        end
    end
end

rep[1].x = 10.0
rep[1].y = 10.0
_pts = tri.points[[12, 11, 10, 9]]
rep[2].x = mean([8.0, 8.0, 4.0, 4.0])
rep[2].y = mean([16.0, 6.0, 6.0, 16.0])
_pts = tri.points[[18, 17, 16, 15, 14, 13]]
rep[3].x = mean([18.0, 18.0, 14.0, 12.0, 14.0, 14.0])
rep[3].y = mean([12.0, 6.0, 2.0, 4.0, 6.0, 10.0])

@testset "Tests with different types of triangulations" begin
    x, y = complicated_geometry()
    rng = StableRNG(919191919)
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri2 = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
    A = get_total_area(tri2)
    refine!(tri2; max_area=1e-2A)

    a, b, c, d = 2.0, 10.0, -5.0, 7.5
    nx = 20
    ny = 10
    tri3 = DT.triangulate_rectangle(a, b, c, d, nx, ny)

    for tri in (tri, tri2, tri3)
        DT.compute_representative_points!(tri)
        rep = DT.get_representative_point_list(tri)
        if !(tri === tri2 || tri === tri3)
            local _pts
            rep[1].x = 10.0
            rep[1].y = 10.0
            _pts = tri.points[[12, 11, 10, 9]]
            rep[2].x = mean([8.0, 8.0, 4.0, 4.0])
            rep[2].y = mean([16.0, 6.0, 6.0, 16.0])
            _pts = tri.points[[18, 17, 16, 15, 14, 13]]
            rep[3].x = mean([18.0, 18.0, 14.0, 12.0, 14.0, 14.0])
            rep[3].y = mean([12.0, 6.0, 2.0, 4.0, 6.0, 10.0])
        end

        @testset "Test that we can find a point in every triangle" begin
            for _ in 1:36
                for V in each_triangle(tri.triangles)
                    if !DT.is_outer_ghost_triangle(indices(V)..., tri.boundary_map)
                        i, j, k = indices(V)
                        p, q, r = get_point(tri.points, tri.representative_point_list, tri.boundary_map, i, j, k)
                        local c
                        c = (p .+ q .+ r) ./ 3
                        for k in each_point_index(tri.points)
                            _V = jump_and_march(tri, c; k)
                            @test DT.is_positively_oriented(DT.triangle_orientation(tri, _V))
                            if !DT.is_ghost_triangle(_V...)
                                @test DT.compare_triangles(_V, V) &&
                                      DT.is_inside(DT.point_position_relative_to_triangle(tri,
                                    _V,
                                    c))
                            else
                                local V1, V2
                                V1 = DT.rotate_ghost_triangle_to_standard_form(V)
                                V2 = DT.rotate_ghost_triangle_to_standard_form(_V)
                                i1 = geti(V1)
                                i2 = geti(V2)
                                if i1 ≠ i2
                                    i1 = i1 - 1
                                end
                                if i1 ≠ i2
                                    i1 = i1 + 1
                                end
                                if i1 ≠ i2
                                    @test false
                                end
                                _V = DT.construct_triangle(typeof(V), i1, getj(V1), getk(V1))
                                @test DT.compare_triangles(_V, V) &&
                                      DT.is_inside(DT.point_position_relative_to_triangle(tri,
                                    _V,
                                    c))
                            end
                        end
                    end
                end
            end

            @testset "Test that we don't break for points already in the triangulation" begin
                for _ in 1:36
                    for k in each_point_index(tri.points)
                        for j in each_point_index(tri.points)
                            _V = jump_and_march(tri, get_point(tri, k); k=j)
                            @test k ∈ indices(_V)
                            @test DT.is_positively_oriented(DT.triangle_orientation(tri, _V))
                        end
                    end
                end
            end

            @testset "Finding points in ghost triangles" begin
                # Technically this will also find points in solid triangles, but by doing random testing with large random points, 
                # we ensure that we primarily find ghost triangles
                for _ in 1:36
                    q = (50randn(), 50rand())
                    for k in each_point_index(tri.points)
                        _V1 = jump_and_march(tri, q; k)
                        @test DT.is_inside(DT.point_position_relative_to_triangle(tri, _V1, q))
                    end
                end
            end
        end
    end
end

@testset "Small geometry with some collinearities" begin
    a = [0.0, 0.0]
    b = [0.0, 1.0]
    c = [0.0, 4.0]
    d = [2.0, 0.0]
    e = [6.0, 0.0]
    f = [8.0, 0.0]
    g = [8.0, 0.5]
    h = [7.5, 1.0]
    i = [4.0, 0.5]
    j = [4.0, 4.0]
    k = [8.0, 4.0]
    pts = [a, b, c, d, e, f, g, h, i, j, k]
    rng = StableRNG(213)
    tri = triangulate(pts; rng, delete_ghosts=false, randomise=false)
    for qi in each_solid_vertex(tri)
        for k in each_solid_vertex(tri)
            q = get_point(tri, qi)
            T = jump_and_march(tri, q; k)
            @test DT.is_positively_oriented(DT.triangle_orientation(tri, T))
            @test DT.is_on(DT.point_position_relative_to_triangle(tri, T, q))
        end
    end
end

@testset "History structure" begin
    history = DT.PointLocationHistory{NTuple{3,Int},NTuple{2,Int},Int}()
    @test history.triangles == NTuple{3,Int}[]
    @test history.collinear_segments == NTuple{2,Int}[]
    @test history.left_vertices == Int[]
    @test history.right_vertices == Int[]
    @test history.collinear_point_indices == Int[]
    DT.add_triangle!(history, 2, 3, 4)
    DT.add_edge!(history, 7, 14)
    DT.add_left_vertex!(history, 10)
    DT.add_right_vertex!(history, 20)
    @test history.triangles == [(2, 3, 4)]
    @test history.collinear_segments == [(7, 14)]
    @test num_edges(history) == 1
    @test history.left_vertices == [10]
    @test history.right_vertices == [20]
    DT.add_left_vertex!(history, 17)
    DT.add_right_vertex!(history, 12)
    DT.add_left_vertex!(history, 19)
    DT.add_right_vertex!(history, 29)
    @test history.left_vertices == [10, 17, 19]
    @test history.right_vertices == [20, 12, 29]
    add_edge!(history, 37, 23)
    add_edge!(history, 50, 101)
    @test num_edges(history) == 3
    DT.add_index!(history, 2)
    DT.add_index!(history, 7)
    DT.add_index!(history, 13)
    @test history.collinear_point_indices == [2, 7, 13]
end

@testset "Reinitialising non-convex triangulations: Used to run into an infinite loop" begin
    a = (0.0, 0.0)
    b = (5.0, 0.0)
    c = (5.0, 5.0)
    d = (0.0, 5.0)
    e = (1.0, 4.0)
    f = (0.0, 3.0)
    g = (1.5, 2.0)
    h = (4.0, 2.0)
    i = (4.0, 1.0)
    j = (1.5, 1.0)
    k = (2.0, 1.5)
    ℓ = (2.4, 1.0)
    m = (2.8, 3.6)
    n = (2.2, 2.8)
    o = (7.2, 2.4)
    p = (7.3, 3.9)
    q = (2.0, 1.4)
    r = (3.0, 1.6)
    s = (3.4, 1.4)
    t = (3.6, 2.2)
    u = (1.0, 2.8)
    v = (0.8, 1.8)
    w = (0.8, 4.0)
    z = (0.4, 4.4)
    c1 = (5.0, 4.6)
    d1 = (4.9, 4.8)
    a1 = (5.6, 1.8)
    b1 = (6.0, 0.6)
    e1 = (5.1, 4.4)
    f1 = (5.0, 4.3)
    pts = [[[d, e, f, a, b, f1, e1, c1, d1, c, d]], [[g, h, i, ℓ, k, j, g]]]
    existing_points = [m, n, o, p, q, r, s, t, u, v, w, z, a1, b1]
    nodes, points = convert_boundary_points_to_indices(pts; existing_points=existing_points)
    tri = triangulate(points; boundary_nodes=nodes, delete_ghosts=false)
    i, j, k = 2, 25, 8
    T = (i, j, k)
    r = 5
    for i in 1:100
        @show i
        V = jump_and_march(tri, get_point(tri, k); point_indices=nothing, m=nothing, try_points=nothing, k=r)
        @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, get_point(tri, k)))
    end
end

@testset "Finding points in a triangulation with concave boundaries" begin
    a, b, c = (0.0, 8.0), (0.0, 6.0), (0.0, 4.0)
    d, e, f = (0.0, 2.0), (0.0, 0.0), (2.0, 0.0)
    g, h, i = (4.0, 0.0), (6.0, 0.0), (8.0, 0.0)
    j, k, ℓ = (8.0, 1.0), (7.0, 2.0), (5.0, 2.0)
    m, n, o = (3.0, 2.0), (2.0, 3.0), (2.0, 5.0)
    p, q, r = (2.0, 7.0), (1.0, 8.0), (1.0, 2.2)
    s, t, u = (0.4, 1.4), (1.2, 1.8), (2.8, 0.6)
    v, w, z = (3.4, 1.2), (1.6, 1.4), (1.6, 2.2)
    outer = [[a, b, c, d, e], [e, f, g, h, i, j, k, ℓ], [ℓ, m, n, o, p, q, a]]
    inner = [[r, z, v, u, w, t, s, r]]
    boundary_nodes, points = convert_boundary_points_to_indices([outer, inner])
    rng = StableRNG(123)
    tri = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
    refine!(tri; max_area=0.01get_total_area(tri), rng)
    qs = [
        (4.0, 5.0), (1.0, 5.6), (0.2, 5.0),
        (0.0, -1.0), (0.5, 3.5), (2.5, 1.5),
        (1.0, 2.0), (4.5, 1.0), (6.0, 1.5),
        (0.5, 8.5), (1.0, 7.5), (1.2, 1.6)
    ]
    δs = [DelaunayTriangulation.distance_to_polygon(q, get_points(tri), get_boundary_nodes(tri)) for q in qs]
    for i in 1:1000
        @show i
        Vs = [jump_and_march(tri, q, concavity_protection=true, rng=StableRNG(i + j)) for (j, q) in enumerate(qs)]
        for (q, δ, V) in zip(qs, δs, Vs)
            cert = DelaunayTriangulation.point_position_relative_to_triangle(tri, V, q)
            if δ ≥ 0.0
                @test !DelaunayTriangulation.is_outside(cert)
                @test !DelaunayTriangulation.is_ghost_triangle(V)
            else
                @test !DelaunayTriangulation.is_outside(cert)
                @test DelaunayTriangulation.is_ghost_triangle(V)
            end
        end
    end
    a, b, c, d = -1, 10, -1, 10
    for i in 1:1000
        rng = StableRNG(i)
        qs = [((b - a) * rand(rng) + a, (d - c) * rand(rng) + c) for _ in 1:1000]
        δs = [DelaunayTriangulation.distance_to_polygon(q, get_points(tri), get_boundary_nodes(tri)) for q in qs]
        Vs = [jump_and_march(tri, q, concavity_protection=true, rng=StableRNG(i + j)) for (j, q) in enumerate(qs)]
        for (q, δ, V) in zip(qs, δs, Vs)
            cert = DelaunayTriangulation.point_position_relative_to_triangle(tri, V, q)
            if δ ≥ 0.0
                @test !DelaunayTriangulation.is_outside(cert)
                @test !DelaunayTriangulation.is_ghost_triangle(V)
            else
                @test !DelaunayTriangulation.is_outside(cert)
                @test DelaunayTriangulation.is_ghost_triangle(V)
            end
        end
    end
end