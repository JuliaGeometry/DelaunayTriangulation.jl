using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Random
using StableRNGs
using StatsBase
using Preferences

@test find_triangle === jump_and_march

tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)
DT.compute_representative_points!(tri)
rep = DT.get_representative_point_list(tri)
pts = get_points(tri)
adj = get_adjacent(tri)
adj2v = get_adjacent2vertex(tri)
ghost_vertex_ranges = get_ghost_vertex_ranges(tri)
ghost_vertex_map = get_ghost_vertex_map(tri)
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
    V4 = [
        (index_map["t"], index_map["p"], index_map["o"]),
        (index_map["t"], index_map["b"], index_map["p"]),
    ]
    q5 = (13.0, 0.0)
    V5 = [
        (index_map["b"], index_map["c"], index_map["p"]),
        (index_map["c"], index_map["b"], DT.𝒢),
    ]
    q6 = (6.0819947177817, 8.90410894457)
    V6 = [(index_map["j"], index_map["k"], DT.𝒢 - 1)]
    q7 = (15.8, 9.2)
    V7 = [
        (index_map["m"], DT.𝒢 - 2, index_map["r"]),
        (index_map["m"], DT.𝒢 - 3, index_map["r"]),
    ]
    q8 = (14.5391680629781, 8.5135910023477)
    V8 = [
        (index_map["m"], index_map["n"], DT.𝒢 - 2),
        (index_map["m"], index_map["n"], DT.𝒢 - 3),
    ]
    q9 = (22.0, 8.0)
    V9 = [(index_map["d"], index_map["c"], DT.𝒢)]
    q10 = (6.0, 11.0)
    V10 = [
        (index_map["j"], index_map["k"], DT.𝒢 - 1),
        (index_map["k"], index_map["ℓ"], DT.𝒢 - 1),
        (index_map["ℓ"], index_map["i"], DT.𝒢 - 1),
        (index_map["i"], index_map["j"], DT.𝒢 - 1),
    ]
    q11 = (6.0819947177817, 8.90410894457)
    V11 = [(index_map["j"], index_map["k"], DT.𝒢 - 1)]
    q12 = (-6.6465702688004, 19.4798146355492)
    V12 = [(index_map["h"], index_map["g"], DT.𝒢)]
    q13 = (-20.0, 10.0)
    V13 = [
        (index_map["h"], index_map["g"], DT.𝒢),
        (index_map["a"], index_map["h"], DT.𝒢),
    ]
    q14 = (20.0, 6.0)
    V14 = [
        (index_map["c"], index_map["d"], index_map["q"]),
        (index_map["c"], index_map["d"], DT.𝒢),
    ]
    allq = [q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14]
    allV = [V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14]
    for _ in 1:36
        for k in DT.each_point_index(pts)
            for i in eachindex(allq)
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
    tri2 = triangulate(points; rng, boundary_nodes, delete_ghosts = false)
    A = get_area(tri2)
    refine!(tri2; max_area = 1.0e-2A, use_circumcenter = true)

    a, b, c, d = 2.0, 10.0, -5.0, 7.5
    nx = 20
    ny = 10
    tri3 = DT.triangulate_rectangle(a, b, c, d, nx, ny, delete_ghosts = false)

    for (tri_idx, tri) in enumerate((tri, tri2, tri3))
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

        if !USE_INEXACTPREDICATES
            @testset "Test that we can find a point in every triangle" begin
                for run in 1:6
                    for V in each_triangle(tri.triangles)
                        rand() < 0.5 && continue # skip 50%
                        if !DT.is_exterior_ghost_triangle(tri, triangle_vertices(V)...)
                            i, j, k = triangle_vertices(V)
                            p, q, r = get_point(tri, i, j, k)
                            local c
                            c = (p .+ q .+ r) ./ 3
                            for k in DT.each_solid_vertex(tri)
                                _V = find_triangle(tri, c; k, concavity_protection = true)
                                @test DT.is_positively_oriented(DT.triangle_orientation(tri, _V))
                                if !DT.is_ghost_triangle(_V...)
                                    @test DT.compare_triangles(_V, V) &&
                                        DT.is_inside(
                                        DT.point_position_relative_to_triangle(
                                            tri,
                                            _V,
                                            c,
                                        ),
                                    )
                                else
                                    local V1, V2
                                    V1 = DT.sort_triangle(V)
                                    V2 = DT.sort_triangle(_V)
                                    i1 = DT.geti(V1)
                                    i2 = DT.geti(V2)
                                    if i1 ≠ i2
                                        i1 = i1 - 1
                                    end
                                    if i1 ≠ i2
                                        i1 = i1 + 1
                                    end
                                    _V = DT.construct_triangle(typeof(V), i1, DT.getj(V1), DT.getk(V1))
                                    @test DT.compare_triangles(_V, V) &&
                                        DT.is_inside(
                                        DT.point_position_relative_to_triangle(
                                            tri,
                                            _V,
                                            c,
                                        ),
                                    )
                                end
                            end
                        end
                    end
                end
            end
        end

        if !!USE_INEXACTPREDICATES && tri_idx ≠ 3
            @testset "Test that we don't break for points already in the triangulation" begin
                for _ in 1:6
                    for k in DT.each_solid_vertex(tri)
                        rand() < 1 / 2 && continue
                        for j in DT.each_solid_vertex(tri)
                            _V = find_triangle(tri, get_point(tri, k); k = j)
                            @test k ∈ triangle_vertices(_V)
                            @test DT.is_positively_oriented(DT.triangle_orientation(tri, _V))
                        end
                    end
                end
            end
        end

        @testset "Finding points in ghost triangles" begin
            # Technically this will also find points in solid triangles, but by doing random testing with large random points, 
            # we ensure that we primarily find ghost triangles
            for _ in 1:6
                q = (50randn(), 50rand())
                for k in DT.each_solid_vertex(tri)
                    rand() < 1 / 2 && continue
                    _V1 = find_triangle(tri, q; k)
                    @test DT.is_inside(DT.point_position_relative_to_triangle(tri, _V1, q))
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
    tri = triangulate(pts; rng, delete_ghosts = false, randomise = false)
    for qi in each_solid_vertex(tri)
        for k in each_solid_vertex(tri)
            q = get_point(tri, qi)
            T = find_triangle(tri, q; k)
            @test DT.is_positively_oriented(DT.triangle_orientation(tri, T))
            @test DT.is_on(DT.point_position_relative_to_triangle(tri, T, q))
        end
    end
end

@testset "History structure" begin
    history = DT.PointLocationHistory{NTuple{3, Int}, NTuple{2, Int}, Int}()
    @test history.triangles == NTuple{3, Int}[]
    @test history.collinear_segments == NTuple{2, Int}[]
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
    DT.add_edge!(history, 37, 23)
    DT.add_edge!(history, 50, 101)
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
    nodes, points = convert_boundary_points_to_indices(pts; existing_points = existing_points)
    tri = triangulate(points; boundary_nodes = nodes, delete_ghosts = false)
    i, j, k = 2, 25, 8
    T = (i, j, k)
    r = 5
    for i in 1:10000
        V = find_triangle(tri, get_point(tri, k); point_indices = nothing, m = nothing, try_points = nothing, k = r, concavity_protection = true)
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
    tri = triangulate(points; rng, boundary_nodes, delete_ghosts = false)
    refine!(tri; max_area = 0.01get_area(tri), rng, use_circumcenter = true)
    qs = [
        (4.0, 5.0), (1.0, 5.6), (0.2, 5.0),
        (0.0, -1.0), (0.5, 3.5), (2.5, 1.5),
        (1.0, 2.0), (4.5, 1.0), (6.0, 1.5),
        (0.5, 8.5), (1.0, 7.5), (1.2, 1.6),
    ]
    δs = [DelaunayTriangulation.distance_to_polygon(q, get_points(tri), get_boundary_nodes(tri)) for q in qs]
    for i in 1:1000
        Vs = [find_triangle(tri, q, concavity_protection = true, rng = StableRNG(i + j)) for (j, q) in enumerate(qs)]
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
        Vs = [find_triangle(tri, q, concavity_protection = true, rng = StableRNG(i + j)) for (j, q) in enumerate(qs)]
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

@testset "Finding points in a triangulation with concave boundaries and disjoint domains" begin
    a, b, c = (0.0, 8.0), (0.0, 6.0), (0.0, 4.0)
    d, e, f = (0.0, 2.0), (0.0, 0.0), (2.0, 0.0)
    g, h, i = (4.0, 0.0), (6.0, 0.0), (8.0, 0.0)
    j, k, ℓ = (8.0, 1.0), (7.0, 1.0), (5.0, 2.0)
    m, n, o = (3.0, 2.0), (2.0, 3.0), (2.0, 5.0)
    p, q, r = (2.0, 7.0), (1.0, 8.0), (1.0, 2.2)
    s, t, u = (0.4, 1.4), (1.2, 1.8), (2.8, 0.6)
    v, w, z = (3.4, 1.2), (1.6, 1.4), (1.6, 2.2)
    outer = [[a, b, c, d, e], [e, f, g, h, i, j, k, ℓ], [ℓ, m, n, o, p, q, a]]
    inner = [[r, z, v, u, w, t, s, r]]
    m₁, n₁, o₁ = (6.0, 8.0), (8.0, 8.0), (8.0, 4.0)
    p₁, q₁, r₁ = (10.0, 4.0), (6.0, 6.0), (8.0, 6.0)
    s₁, t₁, u₁ = (9.0, 7.0), (4.0, 4.0), (5.0, 4.0)
    v₁, w₁ = (5.0, 3.0), (4.0, 3.0)
    new_domain₁ = [[m₁, q₁, o₁, p₁, r₁, s₁, n₁, m₁]]
    new_domain₂ = [[t₁, w₁, v₁, u₁, t₁]]
    boundary_nodes, points = convert_boundary_points_to_indices([outer, inner, new_domain₁, new_domain₂])
    rng = StableRNG(125123)
    tri = triangulate(points; rng, boundary_nodes, delete_ghosts = false)
    @test DT.is_disjoint(tri)
    refine!(tri; max_area = 0.001get_area(tri), rng, use_circumcenter = true)
    qs = [
        (0.6, 6.4), (1.4, 0.8), (3.1, 2.9),
        (6.3, 4.9), (4.6, 3.5), (7.0, 7.0),
        (8.9, 5.1), (5.8, 0.8), (1.0, 1.5),
        (1.5, 2.0), (8.15, 6.0),
    ]
    q = (7.0, 7.0) # When you have a point that is exactly the same as a representative point, you need to be careful of (1) the point being exactly collinear with a ghost edge, and (2) the algorithm mistakenly regarding q as if it were equal to a triangle's vertices (since this is where the ghost vertices map). This is now fixed, but we need this isolated to avoid regressions.
    V = find_triangle(tri, q, rng = StableRNG(268), k = 31, concavity_protection = true)
    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, q))
    δs = [DelaunayTriangulation.distance_to_polygon(q, get_points(tri), get_boundary_nodes(tri)) for q in qs]
    for i in 1:1000
        Vs = [find_triangle(tri, q; concavity_protection = true, rng = StableRNG(i)) for q in qs]
        for (q, δ, V) in zip(qs, δs, Vs)
            cert = DelaunayTriangulation.point_position_relative_to_triangle(tri, V, q)
            if δ > 0.0
                @test !DelaunayTriangulation.is_outside(cert)
                @test !DelaunayTriangulation.is_ghost_triangle(V)
            elseif δ < 0.0
                @test !DelaunayTriangulation.is_outside(cert)
                @test DelaunayTriangulation.is_ghost_triangle(V)
            else
                @test !DelaunayTriangulation.is_outside(cert)
            end
        end
    end
    a, b, c, d = -1, 10, -1, 10
    for i in 1:100
        rng = StableRNG(i)
        qs = [((b - a) * rand(rng) + a, (d - c) * rand(rng) + c) for _ in 1:100]
        δs = [DelaunayTriangulation.distance_to_polygon(q, get_points(tri), get_boundary_nodes(tri)) for q in qs]
        Vs = [find_triangle(tri, q, concavity_protection = true, rng = StableRNG(i + j)) for (j, q) in enumerate(qs)]
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

@testset "Stopping at barriers" begin
    @testset "A simple boundary" begin
        for _ in 1:10
            points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            tri = triangulate(points; boundary_nodes = [1, 2, 3, 4, 1], randomise = false)
            V, invisible_flag = find_triangle(tri, (1 / 2, -1.0), use_barriers = Val(true), k = 4)
            @test invisible_flag && DT.is_invisible(DT.test_visibility(tri, (1 / 2, -1.0), 4))
            @inferred find_triangle(tri, (1 / 2, -1.0), use_barriers = Val(true), k = 4)
            @test V == (1, 2, 3)
            @test DT.is_positively_oriented(DT.triangle_orientation(tri, V))
            V, invisible_flag = find_triangle(tri, (1 / 2, -1.0), use_barriers = Val(true), k = 3)
            @test invisible_flag && DT.is_invisible(DT.test_visibility(tri, (1 / 2, -1.0), 3))
            @test V == (1, 2, 3)
            @test DT.is_positively_oriented(DT.triangle_orientation(tri, V))
            V, invisible_flag = jump_and_march(tri, (1 / 2, -1.0), use_barriers = Val(true), k = 1)
            @test invisible_flag && DT.is_invisible(DT.test_visibility(tri, (1 / 2, -1.0), 1))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V, (1 / 2, -1, 0))) # starting at a boundary edge right next to the query point
            V, invisible_flag = find_triangle(tri, (1 / 2, -1.0), use_barriers = Val(true), k = 2)
            @test invisible_flag && DT.is_invisible(DT.test_visibility(tri, (1 / 2, -1.0), 2))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V, (1 / 2, -1.0)))
        end
    end

    @testset "An unconstrained triangulation with a segment inside it" begin
        for _ in 1:1000
            points = [
                (-8.0, 6.0), (-10.0, 1.0), (-2.0, -2.0), (1.38, 2.05),
                (2.0, 1.0), (4.78, 3.23), (0.0, 6.0), (-3.0, 4.0), (-3.0, 2.0),
            ]
            segments = Set([(8, 9)])
            tri = triangulate(points; segments, randomise = false)
            q1 = (-1.0, 1.0)
            q2 = (-4.0, 1.0)
            q3 = (-8.0, 2.0)
            V1, invisible_flag1 = find_triangle(tri, q1, use_barriers = Val(true), k = 1)
            V2, invisible_flag2 = find_triangle(tri, q2, use_barriers = Val(true), k = 1)
            @inferred find_triangle(tri, q1, use_barriers = Val(true), k = 1)
            V3, invisible_flag3 = find_triangle(tri, q3, use_barriers = Val(true), k = 1)
            @test DT.compare_triangles(V1, (9, 8, 1))
            @test DT.is_positively_oriented(DT.triangle_orientation(tri, V1))
            @test DT.compare_triangles(V2, (2, 3, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V2, q2))
            @test DT.compare_triangles(V3, (1, 2, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V3, q3))
            @test invisible_flag1 && DT.is_invisible(DT.test_visibility(tri, q1, 1))
            @test !invisible_flag2 && DT.is_visible(DT.test_visibility(tri, q2, 1))
            @test !invisible_flag3 && DT.is_visible(DT.test_visibility(tri, q3, 1))
            V1, invisible_flag1 = find_triangle(tri, q1, use_barriers = Val(true), k = 2)
            V2, invisible_flag2 = find_triangle(tri, q2, use_barriers = Val(true), k = 2)
            V3, invisible_flag3 = find_triangle(tri, q3, use_barriers = Val(true), k = 2)
            @test DT.compare_triangles(V1, (9, 3, 4))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V1, q1))
            @test DT.compare_triangles(V2, (2, 3, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V2, q2))
            @test DT.compare_triangles(V3, (1, 2, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V3, q3))
            @test !invisible_flag1 && DT.is_visible(DT.test_visibility(tri, q1, 2))
            @test !invisible_flag2 && DT.is_visible(DT.test_visibility(tri, q2, 2))
            @test !invisible_flag3 && DT.is_visible(DT.test_visibility(tri, q3, 2))
            V1, invisible_flag1 = find_triangle(tri, q1, use_barriers = Val(true), k = 3)
            V2, invisible_flag2 = find_triangle(tri, q2, use_barriers = Val(true), k = 3)
            V3, invisible_flag3 = find_triangle(tri, q3, use_barriers = Val(true), k = 3)
            @test DT.compare_triangles(V1, (9, 3, 4))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V1, q1))
            @test DT.compare_triangles(V2, (2, 3, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V2, q2))
            @test DT.compare_triangles(V3, (1, 2, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V3, q3))
            @test !invisible_flag1 && DT.is_visible(DT.test_visibility(tri, q1, 3))
            @test !invisible_flag2 && DT.is_visible(DT.test_visibility(tri, q2, 3))
            @test !invisible_flag3 && DT.is_visible(DT.test_visibility(tri, q3, 3))
            V1, invisible_flag1 = find_triangle(tri, q1, use_barriers = Val(true), k = 4)
            V2, invisible_flag2 = find_triangle(tri, q2, use_barriers = Val(true), k = 4)
            V3, invisible_flag3 = find_triangle(tri, q3, use_barriers = Val(true), k = 4)
            @test DT.compare_triangles(V1, (9, 3, 4))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V1, q1))
            @test DT.compare_triangles(V2, (2, 3, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V2, q2))
            @test DT.compare_triangles(V3, (8, 9, 4))
            @test DT.is_positively_oriented(DT.triangle_orientation(tri, V3))
            @test !invisible_flag1 && DT.is_visible(DT.test_visibility(tri, q1, 4))
            @test !invisible_flag2 && DT.is_visible(DT.test_visibility(tri, q2, 4))
            @test invisible_flag3 && DT.is_invisible(DT.test_visibility(tri, q3, 4))
            V1, invisible_flag1 = find_triangle(tri, q1, use_barriers = Val(true), k = 5)
            V2, invisible_flag2 = find_triangle(tri, q2, use_barriers = Val(true), k = 5)
            V3, invisible_flag3 = find_triangle(tri, q3, use_barriers = Val(true), k = 5)
            @test DT.compare_triangles(V1, (9, 3, 4))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V1, q1))
            @test DT.compare_triangles(V2, (2, 3, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V2, q2))
            @test DT.compare_triangles(V3, (1, 2, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V3, q3))
            @test !invisible_flag1 && DT.is_visible(DT.test_visibility(tri, q1, 5))
            @test !invisible_flag2 && DT.is_visible(DT.test_visibility(tri, q2, 5))
            @test !invisible_flag3 && DT.is_visible(DT.test_visibility(tri, q3, 5))
            V1, invisible_flag1 = find_triangle(tri, q1, use_barriers = Val(true), k = 6)
            V2, invisible_flag2 = find_triangle(tri, q2, use_barriers = Val(true), k = 6)
            V3, invisible_flag3 = find_triangle(tri, q3, use_barriers = Val(true), k = 6)
            @test DT.compare_triangles(V1, (9, 3, 4))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V1, q1))
            @test DT.compare_triangles(V2, (2, 3, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V2, q2))
            @test DT.compare_triangles(V3, (8, 9, 4))
            @test DT.is_positively_oriented(DT.triangle_orientation(tri, V3))
            @test !invisible_flag1 && DT.is_visible(DT.test_visibility(tri, q1, 6))
            @test !invisible_flag2 && DT.is_visible(DT.test_visibility(tri, q2, 6))
            @test invisible_flag3 && DT.is_invisible(DT.test_visibility(tri, q3, 6))
            V1, invisible_flag1 = find_triangle(tri, q1, use_barriers = Val(true), k = 7)
            V2, invisible_flag2 = find_triangle(tri, q2, use_barriers = Val(true), k = 7)
            V3, invisible_flag3 = find_triangle(tri, q3, use_barriers = Val(true), k = 7)
            @test DT.compare_triangles(V1, (9, 3, 4))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V1, q1))
            @test DT.compare_triangles(V2, (8, 9, 4))
            @test DT.is_positively_oriented(DT.triangle_orientation(tri, V2))
            @test DT.compare_triangles(V3, (1, 2, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V3, q3))
            @test !invisible_flag1 && DT.is_visible(DT.test_visibility(tri, q1, 7))
            @test invisible_flag2 && DT.is_invisible(DT.test_visibility(tri, q2, 7))
            @test !invisible_flag3 && DT.is_visible(DT.test_visibility(tri, q3, 7))
            V1, invisible_flag1 = find_triangle(tri, q1, use_barriers = Val(true), k = 8)
            V2, invisible_flag2 = find_triangle(tri, q2, use_barriers = Val(true), k = 8)
            V3, invisible_flag3 = find_triangle(tri, q3, use_barriers = Val(true), k = 8)
            @test DT.compare_triangles(V1, (9, 3, 4))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V1, q1))
            @test DT.compare_triangles(V2, (2, 3, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V2, q2))
            @test DT.compare_triangles(V3, (1, 2, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V3, q3))
            @test !invisible_flag1 && DT.is_visible(DT.test_visibility(tri, q1, 8))
            @test !invisible_flag2 && DT.is_visible(DT.test_visibility(tri, q2, 8))
            @test !invisible_flag3 && DT.is_visible(DT.test_visibility(tri, q3, 8))
            V1, invisible_flag1 = find_triangle(tri, q1, use_barriers = Val(true), k = 9)
            V2, invisible_flag2 = find_triangle(tri, q2, use_barriers = Val(true), k = 9)
            V3, invisible_flag3 = find_triangle(tri, q3, use_barriers = Val(true), k = 9)
            @test DT.compare_triangles(V1, (9, 3, 4))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V1, q1))
            @test DT.compare_triangles(V2, (2, 3, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V2, q2))
            @test DT.compare_triangles(V3, (1, 2, 9))
            @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V3, q3))
            @test !invisible_flag1 && DT.is_visible(DT.test_visibility(tri, q1, 9))
            @test !invisible_flag2 && DT.is_visible(DT.test_visibility(tri, q2, 9))
            @test !invisible_flag3 && DT.is_visible(DT.test_visibility(tri, q3, 9))
        end
    end

    @testset "Triangulation with many holes" begin
        for _ in 1:100
            a, b, c = (0.0, 0.0), (10.0, 0.0), (10.0, 10.0)
            d, e, f = (0.0, 10.0), (0.0, 8.0), (0.0, 6.0)
            g, h, i = (0.0, 4.0), (0.0, 2.0), (2.0, 0.0)
            j, k, ℓ = (4.0, 0.0), (6.0, 0.0), (8.0, 0.0)
            m, n, o = (10.0, 2.0), (10.0, 4.0), (10.0, 6.0)
            p, q, r = (10.0, 8.0), (2.0, 10.0), (4.0, 10.0)
            s, t, u = (6.0, 10.0), (8.0, 10.0), (1.0, 9.0)
            v, w, z = (1.0, 1.0), (2.0, 1.0), (2.0, 9.0)
            a1, b1, c1 = (4.0, 7.0), (7.0, 4.0), (8.0, 7.0)
            d1, e1, f1 = (5.0, 3.0), (6.0, 3.0), (6.0, 2.0)
            g1, h1, i1, j1 = (5.0, 2.0), (4.0, 9.0), (9.0, 9.0), (6.0, 8.0)
            hole1 = [d, e, f, g, h, a, i, j, k, ℓ, b, m, n, o, p, c, t, s, r, q, d]
            hole2 = [u, z, a1, w, v, u]
            hole3 = [h1, i1, j1, h1]
            hole4 = [c1, b1, e1, c1]
            hole5 = [d1, f1, g1, d1]
            boundary = [[hole1], [hole2], [hole3], [hole4], [hole5]]
            boundary_nodes, points = convert_boundary_points_to_indices(boundary)
            q1, r1 = (5.0, 6.0), (3.0, 1.0)
            push!(points, q1, r1)
            tri = triangulate(points; boundary_nodes, rng = StableRNG(123))
            K1, L1, M1, N1, O1, P1, S1 = (2.0, 6.0), (6.0, 1.0), (5.2, 2.4), (9.0, 4.0), (5.0, -2.0), (3.0, 13.0), (0.5, 4.0)
            all_q = [K1, L1, M1, N1, O1, P1, S1]
            blocked_test(Vfound, Vtrue) = begin
                _p, _q, _r = Vtrue
                _i, _j, __k = findfirst(==(_p), points), findfirst(==(_q), points), findfirst(==(_r), points)
                _Vtrue = (_i, _j, __k)
                return DT.compare_triangles(Vfound, _Vtrue) && DT.is_positively_oriented(DT.triangle_orientation(tri, Vfound))
            end
            visible_test(Vfound, Vtrue, q) = begin
                _p, _q, _r = Vtrue
                _i, _j, __k = !(_p isa Integer) ? findfirst(==(_p), points) : _p, # check for ghost vertices 
                    !(_q isa Integer) ? findfirst(==(_q), points) : _q,
                    !(_r isa Integer) ? findfirst(==(_r), points) : _r
                _Vtrue = (_i, _j, __k)
                return DT.compare_triangles(Vfound, _Vtrue) && !DT.is_outside(DT.point_position_relative_to_triangle(tri, Vfound, q))
            end

            _k = findfirst(==(a), points)
            (VK, flagVK), (VL, flagVL), (VM, flagVM), (VN, flagVN), (VO, flagVO), (VP, flagVP), (VS, flagVS) = find_triangle.(Ref(tri), all_q; k = _k, use_barriers = Val(true))
            #@test blocked_test(VK, (u, f, v))
            #@test visible_test(VL, (f1, g1, k), L1) || visible_test(VL, (f1, k, ℓ), L1) # on (f1, k)
            #@test blocked_test(VM, (d1, r1, g1))
            #@test blocked_test(VN, (d1, r1, g1))
            #@test visible_test(VO, (ℓ, k, -1), O1)
            #@test blocked_test(VP, (u, f, v))
            #@test visible_test(VS, (u, f, v), S1)
            @test flagVK && DT.is_invisible(DT.test_visibility(tri, K1, _k))
            @test !flagVL && DT.is_visible(DT.test_visibility(tri, L1, _k))
            @test flagVM && DT.is_invisible(DT.test_visibility(tri, M1, _k))
            @test flagVN && DT.is_invisible(DT.test_visibility(tri, N1, _k))
            @test flagVO && DT.is_invisible(DT.test_visibility(tri, O1, _k))
            @test flagVP && DT.is_invisible(DT.test_visibility(tri, P1, _k))
            @test !flagVS && DT.is_visible(DT.test_visibility(tri, S1, _k))
            _k = findfirst(==(q), points)
            (VK, flagVK), (VL, flagVL), (VM, flagVM), (VN, flagVN), (VO, flagVO), (VP, flagVP), (VS, flagVS) = find_triangle.(Ref(tri), all_q; k = _k, use_barriers = Val(true))
            #@test blocked_test(VK, (q, u, z))
            #@test blocked_test(VL, (h1, z, a1))
            #@test blocked_test(VM, (z, a1, h1))
            #@test blocked_test(VN, (c1, q1, e1))
            #@test blocked_test(VO, (a1, h1, z))
            #@test visible_test(VP, (r, s, -1), P1)
            #@test blocked_test(VS, (q, u, z))
            @test flagVK && DT.is_invisible(DT.test_visibility(tri, K1, _k))
            @test flagVL && DT.is_invisible(DT.test_visibility(tri, L1, _k))
            @test flagVM && DT.is_invisible(DT.test_visibility(tri, M1, _k))
            @test flagVN && DT.is_invisible(DT.test_visibility(tri, N1, _k))
            @test flagVO && DT.is_invisible(DT.test_visibility(tri, O1, _k))
            @test flagVP && DT.is_invisible(DT.test_visibility(tri, P1, _k))
            @test flagVS && DT.is_invisible(DT.test_visibility(tri, S1, _k))
            _k = findfirst(==(q1), points)
            (VK, flagVK), (VL, flagVL), (VM, flagVM), (VN, flagVN), (VO, flagVO), (VP, flagVP), (VS, flagVS) = find_triangle.(Ref(tri), all_q; k = _k, use_barriers = Val(true))
            #@test blocked_test(VK, (a1, w, d1))
            #@test blocked_test(VL, (e1, d1, f1))
            #@test blocked_test(VM, (e1, d1, f1))
            #@test blocked_test(VN, (c1, q1, e1))
            #@test blocked_test(VO, (c1, q1, e1)) || blocked_test(VO, (g1, j, k)) || blocked_test(VO, (d1, f1, e1)) # search line goes right through d1g1f1 
            #@test blocked_test(VP, (h1, a1, j1))
            #@test blocked_test(VS, (a1, w, d1))
            @test flagVK && DT.is_invisible(DT.test_visibility(tri, K1, _k))
            @test flagVL && DT.is_invisible(DT.test_visibility(tri, L1, _k))
            @test flagVM && DT.is_invisible(DT.test_visibility(tri, M1, _k))
            @test flagVN && DT.is_invisible(DT.test_visibility(tri, N1, _k))
            @test flagVO && DT.is_invisible(DT.test_visibility(tri, O1, _k))
            @test flagVP && DT.is_invisible(DT.test_visibility(tri, P1, _k))
            @test flagVS && DT.is_invisible(DT.test_visibility(tri, S1, _k))
            _k = findfirst(==(p), points)
            (VK, flagVK), (VL, flagVL), (VM, flagVM), (VN, flagVN), (VO, flagVO), (VP, flagVP), (VS, flagVS) = find_triangle.(Ref(tri), all_q; k = _k, use_barriers = Val(true))
            #@test blocked_test(VK, (a1, w, d1))
            #@test visible_test(VL, (f1, g1, k), L1) || visible_test(VL, (f1, k, ℓ), L1) # on (f1, k)
            #@test blocked_test(VM, (c1, b1, o))
            #@test visible_test(VN, (o, b1, n), N1) || visible_test(VN, (n, b1, m), N1)
            #@test blocked_test(VO, (g1, j, k))
            #@test blocked_test(VP, (i1, j1, c1))
            #@test blocked_test(VS, (a1, w, d1))
            @test flagVK && DT.is_invisible(DT.test_visibility(tri, K1, _k))
            @test !flagVL && DT.is_visible(DT.test_visibility(tri, L1, _k))
            @test flagVM && DT.is_invisible(DT.test_visibility(tri, M1, _k))
            @test !flagVN && DT.is_visible(DT.test_visibility(tri, N1, _k))
            @test flagVO && DT.is_invisible(DT.test_visibility(tri, O1, _k))
            @test flagVP && DT.is_invisible(DT.test_visibility(tri, P1, _k))
            @test flagVS && DT.is_invisible(DT.test_visibility(tri, S1, _k))
            _k = findfirst(==(b), points)
            (VK, flagVK), (VL, flagVL), (VM, flagVM), (VN, flagVN), (VO, flagVO), (VP, flagVP), (VS, flagVS) = find_triangle.(Ref(tri), all_q; k = _k, use_barriers = Val(true))
            @test blocked_test(VK, (b1, e1, f1)) || blocked_test(VK, (a1, w, d1))
            @test visible_test(VL, (f1, g1, k), L1) || visible_test(VL, (f1, k, ℓ), L1) # on (f1, k)
            @test blocked_test(VM, (f1, g1, k))
            @test visible_test(VN, (o, b1, n), N1) || visible_test(VN, (n, b1, m), N1)
            @test visible_test(VO, (ℓ, k, -1), O1)
            @test blocked_test(VP, (c1, b1, o))
            @test blocked_test(VS, (f1, g1, k))
            @test flagVK && DT.is_invisible(DT.test_visibility(tri, K1, _k))
            @test !flagVL && DT.is_visible(DT.test_visibility(tri, L1, _k))
            @test flagVM && DT.is_invisible(DT.test_visibility(tri, M1, _k))
            @test !flagVN && DT.is_visible(DT.test_visibility(tri, N1, _k))
            @test flagVO && DT.is_invisible(DT.test_visibility(tri, O1, _k))
            @test flagVP && DT.is_invisible(DT.test_visibility(tri, P1, _k))
            @test flagVS && DT.is_invisible(DT.test_visibility(tri, S1, _k))
            add_segment!(tri, findfirst(==(m), points), findfirst(==(ℓ), points))
            push!(all_q, (9.0, 1.0), (9.5, 0.5), (8.9, 1.1))
            (VK, flagVK), (VL, flagVL), (VM, flagVM), (VN, flagVN), (VO, flagVO), (VP, flagVP), (VS, flagVS), (VH4, flagVH4), (VK9, flagVK9), (VH5, flagVH5) = find_triangle.(Ref(tri), all_q; k = _k, use_barriers = Val(true))
            @test blocked_test(VK, (m, ℓ, b))
            @test blocked_test(VL, (m, ℓ, b))
            @test blocked_test(VM, (m, ℓ, b))
            @test blocked_test(VN, (m, ℓ, b))
            @test visible_test(VO, (ℓ, k, -1), O1)
            @test blocked_test(VP, (m, ℓ, b))
            @test blocked_test(VS, (m, ℓ, b))
            @test blocked_test(VH4, (m, ℓ, b)) && visible_test(VH4, (m, ℓ, b), (9.0, 1.0))
            @test visible_test(VK9, (m, ℓ, b), (9.5, 0.5))
            @test blocked_test(VH5, (m, ℓ, b))
            @test flagVK && DT.is_invisible(DT.test_visibility(tri, K1, _k))
            @test flagVL && DT.is_invisible(DT.test_visibility(tri, L1, _k))
            @test flagVM && DT.is_invisible(DT.test_visibility(tri, M1, _k))
            @test flagVN && DT.is_invisible(DT.test_visibility(tri, N1, _k))
            @test flagVO && DT.is_invisible(DT.test_visibility(tri, O1, _k))
            @test flagVP && DT.is_invisible(DT.test_visibility(tri, P1, _k))
            @test flagVS && DT.is_invisible(DT.test_visibility(tri, S1, _k))
            @test !flagVH4 && DT.is_visible(DT.test_visibility(tri, (9.0, 1.0), _k))
            @test !flagVK9 && DT.is_visible(DT.test_visibility(tri, (9.5, 0.5), _k))
            @test flagVH5 && DT.is_invisible(DT.test_visibility(tri, (8.9, 1.1), _k))
        end
    end

    @testset "Unconstrained triangulation - still get correct results even with use_barriers" begin
        for _ in 1:50
            points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            append!(points, tuple.(rand(1000), rand(1000)))
            tri = triangulate(points)
            for _ in 1:20000
                q = rand(2)
                V, flag = find_triangle(tri, q; use_barriers = Val(true))
                @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, q))
                @test !flag
            end
        end
    end

    @testset "Triangulation: Never get a ghost triangle when using barriers" begin
        points = [(-0.2, -0.2), (0.2, -0.2), (0.2, 0.2), (-0.2, 0.2)]
        for _ in 1:1000
            q = (0.19 + 0.19) * rand(2) .- 0.19
            push!(points, Tuple(q))
        end
        tri = triangulate(points; boundary_nodes = [1, 2, 3, 4, 1])
        for _ in 1:20000
            q = 2randn(2)
            k = rand(each_solid_vertex(tri))
            while DT.is_boundary_node(tri, k)[1]
                k = rand(each_solid_vertex(tri))
            end
            V, flag = find_triangle(tri, q; k, use_barriers = Val(true))
            @test !DT.is_ghost_triangle(V)
        end
    end
end
