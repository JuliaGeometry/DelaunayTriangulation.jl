module HelperFunctions
using StaticArrays
using StableRNGs
using LinearAlgebra
using ElasticArrays
using Random
using DataStructures
using DelimitedFiles
using OrderedCollections
using Test
using DelaunayTriangulation
import SpatialIndexing as SI
getxy((0.0, 0.0)) # avoid shadow

const DT = DelaunayTriangulation

using DelaunayTriangulation: validate_triangulation

function complicated_geometry()
    x1 = [collect(LinRange(0, 2, 4)),
        collect(LinRange(2, 2, 4)),
        collect(LinRange(2, 0, 4)),
        collect(LinRange(0, 0, 4))]
    y1 = [collect(LinRange(0, 0, 4)),
        collect(LinRange(0, 6, 4)),
        collect(LinRange(6, 6, 4)),
        collect(LinRange(6, 0, 4))]
    r = 0.5
    h = k = 0.6
    θ = LinRange(2π, 0, 50)
    x2 = [h .+ r .* cos.(θ)]
    y2 = [k .+ r .* sin.(θ)]
    r = 0.2
    h = 1.5
    k = 0.5
    x3 = [h .+ r .* cos.(θ)]
    y3 = [k .+ r .* sin.(θ)]
    x4 = reverse(reverse.([collect(LinRange(1, 1.5, 4)),
        collect(LinRange(1.5, 1.5, 4)),
        collect(LinRange(1.5, 1, 4)),
        collect(LinRange(1, 1, 4))]))
    y4 = reverse(reverse.([collect(LinRange(2, 2, 4)),
        collect(LinRange(2, 5, 4)),
        collect(LinRange(5, 5, 4)),
        collect(LinRange(5, 2, 4))]))
    x5 = [reverse([0.2, 0.5, 0.75, 0.75, 0.2, 0.2])]
    y5 = [reverse([2.0, 2.0, 3.0, 4.0, 5.0, 2.0])]
    x = [x1, x2, x3, x4, x5]
    y = [y1, y2, y3, y4, y5]
    return x, y, x1, x2, x3, x4, x5, y1, y2, y3, y4, y5
end

function simple_geometry()
    a = (0.0, 0.0)
    b = (10.0, 0.0)
    c = (20.0, 0.0)
    d = (20.0, 10.0)
    e = (20.0, 20.0)
    f = (10.0, 20.0)
    g = (0.0, 20.0)
    h = (0.0, 10.0)
    i = (4.0, 16.0)
    j = (4.0, 6.0)
    k = (8.0, 6.0)
    ℓ = (8.0, 16.0)
    m = (14.0, 10.0)
    n = (14.0, 6.0)
    o = (12.0, 4.0)
    p = (14.0, 2.0)
    q = (18.0, 6.0)
    r = (18.0, 12.0)
    s = (4.0, 2.0)
    t = (10.0, 2.0)
    u = (10.0, 4.0)
    v = (14.0, 18.0)
    w = (12.0, 14.0)
    z = (18.0, 16.0)
    a1 = (4.0, 18.0)
    b1 = (2.0, 12.0)
    pts = [a, b, c, d, e, f, g, h, i, j, k, ℓ,
        m, n, o, p, q, r, s, t, u, v, w, z, a1, b1]
    T = [h a s
        h s j
        b1 h j
        i b1 j
        g b1 i
        g h b1
        g i a1
        g a1 f
        a1 i f
        f i ℓ
        f ℓ w
        j s t
        j t u
        s a t
        t a b
        k j u
        w ℓ k
        w k u
        f w v
        e f v
        e v z
        e z r
        e r d
        v m z
        z m r
        v w m
        w u m
        m u n
        n u o
        o u t
        o t p
        p t b
        p b c
        q p c
        d q c
        r q d]
    T = indexin(T, pts)
    T = Set{NTuple{3,Int}}((Tuple(T) for T in eachrow(T)))
    outer = [[indexin([a, b, c, d, e, f, g, h, a], pts)...]]
    inner1 = [[indexin([ℓ, k, j, i, ℓ], pts)...]]
    inner2 = [[indexin([r, q, p], pts)...], [indexin([p, o, n, m, r], pts)...]]
    boundary_nodes = [outer, inner1, inner2]
    label_map = Dict(["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "ℓ",
        "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "z", "a1",
        "b1"] .=> pts)
    index_map = Dict(["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "ℓ",
        "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "z", "a1",
        "b1"] .=> DT.each_point_index(pts))
    return DT.Triangulation(pts, T, boundary_nodes, delete_ghosts=true), label_map, index_map
end

macro _adj(i, j, k)
    return :(($i, $j) => $k, ($j, $k) => $i, ($k, $i) => $j)
end

import SimpleGraphs: SimpleGraphs
function example_triangulation()
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    p4 = @SVector[-1.0, 2.0]
    p5 = @SVector[4.0, 2.0]
    p6 = @SVector[-2.0, -1.0]
    p7 = @SVector[2.0, 1.0]
    p8 = @SVector[5.0, 1.0]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    T = Set{NTuple{3,Int}}([
        (6, 3, 1),
        (3, 2, 5),
        (4, 1, 5),
        (4, 6, 1),
        (5, 1, 3)
    ])
    A = [
        0 0 1 1 1 1 1
        0 0 0 1 1 1 1
        1 0 0 1 0 1 0
        1 1 1 0 0 1 1
        1 1 0 0 0 1 1
        1 1 1 1 1 0 0
        1 1 0 1 1 0 0
    ]
    DG = SimpleGraphs.relabel(SimpleGraphs.UndirectedGraph(A), Dict(1:7 .=> [-1, (1:6)...]))
    DG = DT.Graph(DG.V, DG.E, DG.N)
    adj = DT.Adjacent(
        Dict(
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (1, 3) => 5, (3, 5) => 1,
            (4, 5) => DT.𝒢, (5, 2) => DT.𝒢,
            (2, 3) => DT.𝒢, (3, 6) => DT.𝒢,
            (6, 4) => DT.𝒢
        )
    )
    adj2v = DT.Adjacent2Vertex(Dict(
        DT.𝒢 => Set{NTuple{2,Int}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
        1 => Set{NTuple{2,Int}}([(5, 4), (3, 5), (6, 3), (4, 6)]),
        2 => Set{NTuple{2,Int}}([(5, 3)]),
        3 => Set{NTuple{2,Int}}([(1, 6), (5, 1), (2, 5)]),
        4 => Set{NTuple{2,Int}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int}}([(1, 4), (3, 1)])
    ))
    rep = Dict{Int,DT.RepresentativeCoordinates{Int,Float64}}()
    tri = DT.Triangulation(pts, T, Int[], Set{NTuple{2,Int}}(), Set{NTuple{2,Int}}(),
        DT.ZeroWeight(), adj, adj2v, DG, (), DT.construct_boundary_edge_map(Int[]),
        Dict{Int,Vector{Int}}(), Dict{Int,UnitRange{Int}}(), DT.ConvexHull(pts, [2, 6, 4, 5, 2]), rep, DT.construct_polygon_hierarchy(pts), nothing,
        DT._build_cache(pts, Int, NTuple{2,Int}, NTuple{3,Int},
            Set{NTuple{2,Int}}, Set{NTuple{3,Int}}, DT.ZeroWeight(), Val(true)))
    DT.compute_representative_points!(tri)
    return tri
end

function example_empty_triangulation()
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    pts = [p1, p2, p3]
    T = Set{NTuple{3,Int}}()
    A = zeros(Int, 0, 0)
    DG = SimpleGraphs.relabel(SimpleGraphs.UndirectedGraph(A), Dict(1:7 .=> [-1, (1:6)...]))
    DG = DT.Graph(DG.V, DG.E, DG.N)
    adj = DT.Adjacent(Dict{NTuple{2,Int},Int}())
    adj2v = DT.Adjacent2Vertex(Dict(DT.𝒢 => Set{NTuple{2,Int}}()))
    rep = Dict{Int,DT.RepresentativeCoordinates{Int,Float64}}()
    tri = DT.Triangulation(pts, T, Int[], Set{NTuple{2,Int}}(), Set{NTuple{2,Int}}(),
        DT.ZeroWeight(), adj, adj2v, DG, (), DT.construct_boundary_edge_map(Int[]),
        Dict{Int,Vector{Int}}(), Dict{Int,UnitRange{Int}}(), DT.ConvexHull(pts, [1, 2, 3, 1]), rep, DT.construct_polygon_hierarchy(pts), nothing,
        DT._build_cache(pts, Int, NTuple{2,Int}, NTuple{3,Int},
            Set{NTuple{2,Int}}, Set{NTuple{3,Int}}, DT.ZeroWeight(), Val(true)))
    DT.compute_representative_points!(tri)
    return tri
end

function example_with_special_corners()
    a = [8.0, 3.0]
    b = [8.0, 1.0]
    c = [6.0, 1.0]
    d = [7.0, 5.0]
    e = [6.0, 8.0]
    f = [3.0, 3.0]
    g = [-1.0, 2.0]
    h = [-5.0, 5.0]
    i = [-5.0, 9.0]
    j = [-5.0, 11.0]
    k = [-4.0, 10.0]
    ℓ = [0.0, 9.87]
    m = [3.0, 11.0]
    n = [1.0, 7.0]
    o = [-1.0, 6.0]
    p = [5.0, 5.0]
    q = [2.0, 4.0]
    r = [3.0, 8.0]
    pts = [a, b, c, d, e, f, g, h, i, j, k, ℓ, m, n, o, p, q, r]
    rng = StableRNG(29292929292)
    tri = triangulate(pts; rng, delete_ghosts=false, randomise=false)
    return tri
end

function shewchuk_example_constrained()
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
    return tri
end

function fixed_shewchuk_example_constrained()
    a = [0.0, 0.0]
    b = [0.0, 1.0]
    c = [0.0, 2.5]
    d = [2.0, 0.0]
    e = [6.0, 0.0]
    f = [8.0, 0.0]
    g = [8.0, 0.5]
    h = [7.5, 1.0]
    i = [4.0, 1.0]
    j = [4.0, 2.5]
    k = [8.0, 2.5]
    pts = [a, b, c, d, e, f, g, h, i, j, k]
    rng = StableRNG(213)
    tri = triangulate(pts; rng, delete_ghosts=false, randomise=false)
    return tri
end

function test_intersections(tri, e, allT, constrained_edges)
    for e in ((e[1], e[2]), (e[2], e[1]))
        intersecting_triangles, collinear_segments, left, right = DT.locate_intersecting_triangles(tri, e)
        @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), intersecting_triangles)
        @test all(!DT.is_none, [DT.triangle_line_segment_intersection(tri, T..., e...) for T in intersecting_triangles])
        @test DT.compare_triangle_collections(allT, intersecting_triangles)
        @test allunique(intersecting_triangles)
        if typeof(constrained_edges) <: AbstractVector
            @test collinear_segments == constrained_edges
        else # Tuple of possibilities, in case the edge's endpoints have equal degree so that we could start at any point
            @test any(==(collinear_segments), constrained_edges)
        end
    end
end

function test_split_edges(tri, edge, current_constrained_edges)
    constrained_edges = get_interior_segments(tri)
    DT.add_edge!(constrained_edges, edge)
    _, collinear_segments = DT.locate_intersecting_triangles(tri, edge)
    DT.split_segment!(tri, edge, collinear_segments)
    @test any(==(constrained_edges), current_constrained_edges)
    return nothing
end

function sort_edge_vector(E)
    sorted_E = similar(E)
    for i in eachindex(E)
        u, v = E[i]
        e = (min(u, v), max(u, v))
        sorted_E[i] = e
    end
    return sort(sorted_E)
end

function compare_edge_vectors(E1, E2)
    E1s = sort_edge_vector(collect(E1))
    E2s = sort_edge_vector(collect(E2))
    return E1s == E2s
end

function test_segment_triangle_intersections(tri, edge, true_triangles, true_collinear_segments, current_constrained_edges)
    constrained_edges = get_interior_segments(tri)
    for edge in ((edge[1], edge[2]), (edge[2], edge[1]))
        intersecting_triangles, collinear_segments, left, right = DT.locate_intersecting_triangles(tri, edge)
        @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), intersecting_triangles)
        @test all(!DT.is_none, [DT.triangle_line_segment_intersection(tri, T..., edge...) for T in intersecting_triangles])
        if typeof(true_triangles) <: AbstractVector || typeof(true_triangles) <: AbstractSet
            @test DT.compare_triangle_collections(true_triangles, intersecting_triangles)
        else
            @test any(V -> DT.compare_triangle_collections(intersecting_triangles, V), true_triangles)
        end
        @test allunique(intersecting_triangles)
        if typeof(true_collinear_segments) <: AbstractVector || typeof(true_collinear_segments) <: AbstractSet
            @test collinear_segments == true_collinear_segments
        else # Tuple of possibilities, in case the edge's endpoints have equal degree so that we could start at any point
            @test any(==(collinear_segments), true_collinear_segments)
        end
    end
    constrained_edges = get_interior_segments(tri)
    DT.add_edge!(constrained_edges, edge)
    _, collinear_segments = DT.locate_intersecting_triangles(tri, edge)
    DT.split_segment!(tri, edge, collinear_segments)
    !isnothing(current_constrained_edges) && @test compare_edge_vectors(constrained_edges, current_constrained_edges)
end

function get_random_vertices_and_constrained_edges(nverts1, nverts2, nedges, rng=Random.default_rng())
    ## To generate a random set of constrained edges, we get a random small triangulation, 
    ## and we just take the edges from that triangulation.
    points = [Tuple(rand(rng, 2)) for _ in 1:nverts1]
    tri = triangulate(points; rng)
    edges = Set{NTuple{2,Int}}()
    all_edges = collect(each_solid_edge(tri))
    iter = 0
    while length(edges) < nedges && iter < 10000
        S = DT.random_edge(rng, all_edges)
        push!(edges, S)
        iter += 1
    end
    ## Now get the rest of the points 
    append!(points, [Tuple(rand(rng, 2)) for _ in 1:(nverts2-nverts1)])
    return points, edges, vec(hcat(getindex.(edges, 1), getindex.(edges, 2))')
end

function second_shewchuk_example_constrained()
    p1 = [0.0, 0.0]
    p2 = [0.0, 20.0]
    p3 = [30.0, 0.0]
    p4 = [30.0, 20.0]
    p5 = [2.0, 16.0]
    p6 = [2.0, 6.0]
    p7 = [4.0, 12.0]
    p8 = [4.0, 2.0]
    p9 = [6.0, 14.0]
    p10 = [6.0, 18.0]
    p11 = [8.0, 12.0]
    p12 = [8.0, 4.0]
    p13 = [10.0, 18.0]
    p14 = [10.0, 2.0]
    p15 = [14.0, 16.0]
    p16 = [18.0, 16.0]
    p17 = [18.0, 4.0]
    p18 = [14.0, 4.0]
    p19 = [12.0, 16.0]
    p20 = [12.0, 2.0]
    p21 = [22.0, 18.0]
    p22 = [22.0, 12.0]
    p23 = [22.0, 10.0]
    p24 = [22.0, 4.0]
    p25 = [26.0, 0.0]
    p26 = [26.0, 16.0]
    p27 = [28.0, 18.0]
    p28 = [28.0, 4.0]
    C = [
        (1, 3),
        (3, 4),
        (4, 2),
        (2, 1),
        (5, 6),
        (7, 8),
        (10, 9),
        (11, 12),
        (13, 14),
        (15, 16),
        (17, 18),
        (19, 20),
        (21, 22),
        (23, 24),
        (25, 26),
        (27, 28)
    ]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
        p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24,
        p25, p26, p27, p28]
    return pts, Set(C)
end

function example_for_testing_add_point_on_constrained_triangulation()
    A = [0; 3]
    B = [4; 0]
    C = [7; 6]
    D = [-4; 8]
    E = [-2; 5]
    F = [-3; -2]
    G = [6; -4]
    H = [4; 4]
    I = [3; 7]
    J = [-3; 3]
    K = [1; 1]
    L = [5; 3]
    M = [2; 2]
    N = [1.6; 2]
    P = [A, B, C, D, E, F, G, H, I, J, K, L, M, N]
    return P, Set([(1, 2)])
end

function validate_statistics(tri::DT.Triangulation, stats=statistics(tri))
    ## Build up the array (see also test_iterators)
    I = DT.integer_type(tri)
    T = NTuple{3,I}
    E = NTuple{2,I}
    solid_triangles = T[]
    ghost_triangles = T[]
    all_triangles = T[]
    solid_vertices = Set{I}()
    ghost_vertices = Set{I}()
    all_vertices = Set{I}()
    solid_edges = E[]
    ghost_edges = E[]
    all_edges = E[]
    for T in each_triangle(tri)
        i, j, k = DT.triangle_vertices(T)
        push!(all_triangles, (i, j, k))
        if DT.is_ghost_triangle(i, j, k)
            push!(ghost_triangles, (i, j, k))
        else
            push!(solid_triangles, (i, j, k))
        end
        for (u, v) in DT.triangle_edges(i, j, k)
            push!(all_edges, (u, v))
            if DT.is_ghost_edge(u, v)
                push!(ghost_edges, (u, v))
            else
                push!(solid_edges, (u, v))
            end
        end
        for k in (i, j, k)
            push!(all_vertices, k)
            if DT.is_ghost_vertex(k)
                push!(ghost_vertices, k)
            else
                push!(solid_vertices, k)
            end
        end
    end
    push!(all_vertices, DT.all_ghost_vertices(tri)...)
    push!(ghost_vertices, DT.all_ghost_vertices(tri)...)
    for e in each_ghost_edge(tri)
        u, v = DT.edge_vertices(e)
        push!(ghost_edges, (min(u, v), max(u, v)))
        push!(all_edges, (min(u, v), max(u, v)))
    end
    for (i, e) in enumerate(all_edges)
        u, v = DT.edge_vertices(e)
        all_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(solid_edges)
        u, v = DT.edge_vertices(e)
        solid_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(ghost_edges)
        u, v = DT.edge_vertices(e)
        ghost_edges[i] = (min(u, v), max(u, v))
    end
    unique!(all_edges)
    unique!(solid_edges)
    unique!(ghost_edges)
    all_vertices = collect(all_vertices)
    solid_vertices = collect(solid_vertices)
    ghost_vertices = collect(ghost_vertices)
    sort!(all_vertices)
    sort!(solid_vertices)
    sort!(ghost_vertices)

    ## Build up the individual statistics 
    areas = Dict{T,Float64}()
    lengths = Dict{T,NTuple{3,Float64}}()
    circumcenters = Dict{T,NTuple{2,Float64}}()
    circumradii = Dict{T,Float64}()
    angles = Dict{T,NTuple{3,Float64}}()
    radius_edge_ratio = Dict{T,Float64}()
    edge_midpoints = Dict{T,NTuple{3,NTuple{2,Float64}}}()
    aspect_ratio = Dict{T,Float64}()
    inradius = Dict{T,Float64}()
    perimeter = Dict{T,Float64}()
    centroid = Dict{T,NTuple{2,Float64}}()
    offcenters = Dict{T,NTuple{2,Float64}}()
    sinks = Dict{T,NTuple{2,Float64}}()
    total_A = 0.0
    for T in each_solid_triangle(tri)
        u, v, w = DT.triangle_vertices(T)
        p, q, r = get_point(tri, u, v, w)
        p = [getx(p), gety(p)]
        q = [getx(q), gety(q)]
        r = [getx(r), gety(r)]
        areas[triangle_vertices(T)] = 0.5 * (p[1] * (q[2] - r[2]) + q[1] * (r[2] - p[2]) + r[1] * (p[2] - q[2]))
        total_A += areas[triangle_vertices(T)]
        ℓ1 = norm(q - p)
        ℓ2 = norm(r - q)
        ℓ3 = norm(r - p)
        ℓmin = min(ℓ1, ℓ2, ℓ3)
        ℓmax = max(ℓ1, ℓ2, ℓ3)
        ℓmed = ℓ1 + ℓ2 + ℓ3 - ℓmin - ℓmax
        ℓ1, ℓ2, ℓ3 = ℓmin, ℓmed, ℓmax
        lengths[triangle_vertices(T)] = (ℓ1, ℓ2, ℓ3)
        r′ = p - r
        s′ = q - r
        ox = r[1] + det([norm(r′)^2 r′[2]; norm(s′)^2 s′[2]]) / (4areas[triangle_vertices(T)])
        oy = r[2] + det([r′[1] norm(r′)^2; s′[1] norm(s′)^2]) / (4areas[triangle_vertices(T)])
        circumcenters[triangle_vertices(T)] = (ox, oy)
        circumradii[triangle_vertices(T)] = norm(r - collect(circumcenters[triangle_vertices(T)]))
        all_angles = [(norm(p - r)^2 + norm(q - r)^2 - norm(p - q)^2) / (2norm(p - r) * norm(q - r)) for (p, q, r) in ((p, q, r), (q, r, p), (r, p, q))]
        all_angles[all_angles.<-1.0] .= -1.0
        all_angles[all_angles.>1.0] .= 1.0
        all_angles = acos.(all_angles)
        sort!(all_angles)
        radius_edge_ratio[triangle_vertices(T)] = circumradii[triangle_vertices(T)] / ℓ1
        edge_midpoints[triangle_vertices(T)] = ((Tuple(0.5 * (p + q))), Tuple(0.5 * (q + r)), Tuple(0.5 * (r + p)))
        inradius[triangle_vertices(T)] = 2areas[triangle_vertices(T)] / (ℓ1 + ℓ2 + ℓ3)
        perimeter[triangle_vertices(T)] = ℓ1 + ℓ2 + ℓ3
        aspect_ratio[triangle_vertices(T)] = inradius[triangle_vertices(T)] / circumradii[triangle_vertices(T)]
        centroid[triangle_vertices(T)] = (1 / 3 * (p[1] + q[1] + r[1]), 1 / 3 * (p[2] + q[2] + r[2]))
        angles[triangle_vertices(T)] = Tuple(all_angles)
        offcenters[triangle_vertices(T)] = DT.triangle_offcenter(p, q, r)
        sinks[triangle_vertices(T)] = DT.triangle_sink(tri, T)
        @test radius_edge_ratio[triangle_vertices(T)] ≥ 1 / sqrt(3) - 0.1
        @test DT.get_radius_edge_ratio(stats, T) ≥ 1 / sqrt(3) - 0.1
        @test angles[triangle_vertices(T)][1] ≤ deg2rad(60) + 0.01
        @test DT.get_minimum_angle(stats, T) ≤ deg2rad(60) + 0.01
    end

    ## Now compare the statistics 
    for T in each_solid_triangle(tri)
        @test areas[triangle_vertices(T)] ≈ DT.get_area(stats, T) rtol = 1e-4 atol = 1e-4
        @test collect(lengths[triangle_vertices(T)]) ≈ collect(DT.get_lengths(stats, T)) rtol = 1e-4 atol = 1e-4
        @test collect(circumcenters[triangle_vertices(T)]) ≈ collect(DT.get_circumcenter(stats, T)) rtol = 1e-4 atol = 1e-4
        @test circumradii[triangle_vertices(T)] ≈ DT.get_circumradius(stats, T) rtol = 1e-4 atol = 1e-4
        @test radius_edge_ratio[triangle_vertices(T)] ≈ DT.get_radius_edge_ratio(stats, T) rtol = 1e-4 atol = 1e-4
        @test collect(collect.(edge_midpoints[triangle_vertices(T)])) ≈ collect(collect.(DT.get_edge_midpoints(stats, T))) rtol = 1e-4 atol = 1e-4
        @test aspect_ratio[triangle_vertices(T)] ≈ DT.get_aspect_ratio(stats, T) rtol = 1e-4 atol = 1e-4
        @test inradius[triangle_vertices(T)] ≈ DT.get_inradius(stats, T) rtol = 1e-4 atol = 1e-4
        @test perimeter[triangle_vertices(T)] ≈ DT.get_perimeter(stats, T) rtol = 1e-4 atol = 1e-4
        @test radius_edge_ratio[triangle_vertices(T)] ≈ 1 / (2sin(angles[triangle_vertices(T)][1])) rtol = 1e-4 atol = 1e-4
        @test (2sin(DT.get_minimum_angle(stats, T) / 2)^2 - 0.1 ≤ DT.get_aspect_ratio(stats, T) ≤ 2tan(DT.get_minimum_angle(stats, T) / 2) + 0.1)
        @test DT.get_radius_edge_ratio(stats, T) ≈ 1 / (2(sin(DT.get_minimum_angle(stats, T)))) rtol = 1e-4 atol = 1e-4
        @test areas[triangle_vertices(T)] ≈ inradius[triangle_vertices(T)] * 0.5perimeter[triangle_vertices(T)] rtol = 1e-4 atol = 1e-4
        @test DT.get_area(stats, T) ≈ DT.get_inradius(stats, T) * 0.5DT.get_perimeter(stats, T) rtol = 1e-4 atol = 1e-4
        @test collect(centroid[triangle_vertices(T)]) ≈ collect(DT.get_centroid(stats, T)) rtol = 1e-4 atol = 1e-4
        @test DT.get_angles(stats, T)[1] ≈ angles[triangle_vertices(T)][1] rtol = 1e-4 atol = 1e-4
        @test DT.get_angles(stats, T)[2] ≈ angles[triangle_vertices(T)][2] rtol = 1e-4 atol = 1e-4
        @test DT.get_angles(stats, T)[3] ≈ angles[triangle_vertices(T)][3] rtol = 1e-4 atol = 1e-4
        @test sum(DT.get_angles(stats, T)) ≈ π rtol = 1e-4 atol = 1e-4
        @test DT.get_minimum_angle(stats, T) ≈ angles[triangle_vertices(T)][1] rtol = 1e-4 atol = 1e-4
        @test DT.get_maximum_angle(stats, T) ≈ angles[triangle_vertices(T)][3] rtol = 1e-4 atol = 1e-4
        @test DT.get_minimum_angle(stats, T) ≈ DT.get_angles(stats, T)[1] rtol = 1e-4 atol = 1e-4
        @test DT.get_maximum_angle(stats, T) ≈ DT.get_angles(stats, T)[3] rtol = 1e-4 atol = 1e-4
        @test collect(offcenters[triangle_vertices(T)]) ≈ collect(DT.get_offcenter(stats, T)) rtol = 1e-4 atol = 1e-2
        @test collect(sinks[triangle_vertices(T)]) ≈ collect(DT.get_sink(stats, T)) rtol = 1e-4 atol = 1e-4
    end
    @test stats.individual_statistics == DT.get_individual_statistics(stats)
    @test stats.area ≈ DT.get_area(stats) rtol = 1e-4 atol = 1e-4
    @test stats.area ≈ total_A rtol = 1e-4 atol = 1e-4

    ## Test the number statistics 
    @test DT.num_vertices(stats) == length(all_vertices) == stats.num_vertices
    @test DT.num_solid_vertices(stats) == length(solid_vertices) == stats.num_solid_vertices
    @test DT.num_ghost_vertices(stats) == length(DT.all_ghost_vertices(tri)) == length(ghost_vertices) == stats.num_ghost_vertices
    @test DT.num_triangles(stats) == length(all_triangles) == stats.num_triangles
    @test DT.num_solid_triangles(stats) == length(solid_triangles) == stats.num_solid_triangles
    @test DT.num_ghost_triangles(stats) == length(ghost_triangles) == stats.num_ghost_triangles
    @test DT.num_edges(stats) == length(all_edges) == stats.num_edges
    @test DT.num_solid_edges(stats) == length(solid_edges) == stats.num_solid_edges
    @test DT.num_ghost_edges(stats) == length(ghost_edges) == stats.num_ghost_edges
    @test DT.num_boundary_segments(stats) == length(keys(get_boundary_edge_map(tri))) == stats.num_boundary_segments
    @test DT.num_interior_segments(stats) == num_edges(get_interior_segments(tri)) == stats.num_interior_segments
    @test DT.num_segments(stats) == num_edges(get_all_segments(tri)) == stats.num_segments
    @test DT.num_convex_hull_vertices(stats) == length(get_convex_hull_vertices(tri)) - 1 == stats.num_convex_hull_vertices

    ## Global statistics  
    smallest_angle = minimum([angles[triangle_vertices(T)][1] for T in each_solid_triangle(tri)])
    largest_angle = maximum([angles[triangle_vertices(T)][3] for T in each_solid_triangle(tri)])
    smallest_area = minimum([areas[triangle_vertices(T)] for T in each_solid_triangle(tri)])
    largest_area = maximum([areas[triangle_vertices(T)] for T in each_solid_triangle(tri)])
    smallest_radius_edge_ratio = minimum([radius_edge_ratio[triangle_vertices(T)] for T in each_solid_triangle(tri)])
    largest_radius_edge_ratio = maximum([radius_edge_ratio[triangle_vertices(T)] for T in each_solid_triangle(tri)])
    @test DT.get_smallest_angle(stats) ≈ smallest_angle rtol = 1e-2
    @test DT.get_largest_angle(stats) ≈ largest_angle rtol = 1e-2
    @test DT.get_smallest_area(stats) ≈ smallest_area rtol = 1e-2
    @test DT.get_largest_area(stats) ≈ largest_area rtol = 1e-2
    @test DT.get_smallest_radius_edge_ratio(stats) ≈ smallest_radius_edge_ratio rtol = 1e-2
    @test DT.get_largest_radius_edge_ratio(stats) ≈ largest_radius_edge_ratio rtol = 1e-2
    @test DT.get_smallest_radius_edge_ratio(stats) ≥ 1 / sqrt(3) - 0.1
    @test DT.get_smallest_angle(stats) ≤ deg2rad(60) + 0.01
end

function slow_encroachment_test(tri::DT.Triangulation)
    E = DT.edge_type(tri)
    I = DT.integer_type(tri)
    ch = Channel{Pair{E,Tuple{Bool,I}}}(Inf) # https://discourse.julialang.org/t/can-dicts-be-threadsafe/27172/17
    @sync for i in collect(each_solid_vertex(tri))
        Base.Threads.@spawn begin
            for j in each_solid_vertex(tri)
                if i < j
                    e = DT.construct_edge(E, i, j)
                    p, q = get_point(tri, i, j)
                    r2 = 0.25((getx(p) - getx(q))^2 + (gety(p) - gety(q))^2)
                    m = (0.5(getx(p) + getx(q)), 0.5(gety(p) + gety(q)))
                    flag = false
                    k_flag = 0
                    for k in each_solid_vertex(tri)
                        if i == j || i == k || j == k
                            continue
                        end
                        r = get_point(tri, k)
                        if (getx(r) - getx(m))^2 + (gety(r) - gety(m))^2 ≤ r2
                            flag = true
                            k_flag = k
                            break
                        end
                    end
                    put!(ch, e => (flag, k_flag))
                end
            end
        end
    end
    not_in_dt_encroached_edges = Dict{E,Tuple{Bool,I}}()
    in_dt_encroached_edges = Dict{E,Tuple{Bool,I}}()
    while !isempty(ch)
        e, (b, k) = take!(ch)
        if DT.edge_exists(tri, e) || DT.edge_exists(tri, DT.reverse_edge(e))
            in_dt_encroached_edges[e] = (b, k)
        else
            not_in_dt_encroached_edges[e] = (b, k)
        end
    end
    return in_dt_encroached_edges, not_in_dt_encroached_edges
end

function slow_encroachment_test_diametral_lens(tri::DT.Triangulation, lens_angle)
    E = DT.edge_type(tri)
    I = DT.integer_type(tri)
    ch = Channel{Pair{E,Tuple{Bool,I}}}(Inf) # https://discourse.julialang.org/t/can-dicts-be-threadsafe/27172/17
    @sync for i in collect(each_solid_vertex(tri))
        Base.Threads.@spawn begin
            for j in each_solid_vertex(tri)
                if i < j
                    e = DT.construct_edge(E, i, j)
                    p, q = get_point(tri, i, j)
                    _, _, lens = compute_diametral_lens(p, q, lens_angle)
                    unique!(lens)
                    push!(lens, lens[begin])
                    boundary_nodes = [(1:length(lens)); 1]
                    if DT.polygon_features(lens, boundary_nodes)[1] < 0.0
                        reverse!(lens)
                    end
                    flag = false
                    k_flag = 0
                    for k in each_solid_vertex(tri)
                        if i == j || i == k || j == k
                            continue
                        end
                        r = get_point(tri, k)
                        if DT.distance_to_polygon(r, lens, boundary_nodes) ≥ 0
                            flag = true
                            k_flag = k
                            break
                        end
                    end
                    put!(ch, e => (flag, k_flag))
                end
            end
        end
    end
    not_in_dt_encroached_edges = Dict{E,Tuple{Bool,I}}()
    in_dt_encroached_edges = Dict{E,Tuple{Bool,I}}()
    while !isempty(ch)
        e, (b, k) = take!(ch)
        if DT.edge_exists(tri, e) || DT.edge_exists(tri, DT.reverse_edge(e))
            in_dt_encroached_edges[e] = (b, k)
        else
            not_in_dt_encroached_edges[e] = (b, k)
        end
    end
    return in_dt_encroached_edges, not_in_dt_encroached_edges
end

function poor_triangulation_example()
    a = [0.0, 0.0] # 1
    b = [3.0, 0.0] # 2
    c = [6.0, 0.0] # 3
    d = [9.0, 0.0] # 4
    e = [12.0, 0.0] # 5
    f = [12.0, 4.0] # 6
    g = [10.0, 7.0] # 7
    h = [5.0, 6.0] # 8
    i = [2.0, 8.0] # 9
    j = [-1.62, 5.45] # 10
    k = [-3.14, 1.21] # 11
    ℓ = [-1.0, -3.0] # 12
    m = [5.0, -4.0] # 13
    n = [11.0, -3.0] # 14
    o = [2.0, 4.0] # 15
    p = [8.0, 4.0] # 16
    T1 = [10, 11, 1]
    T2 = [1, 11, 12]
    T3 = [1, 12, 2]
    T4 = [2, 12, 13]
    T5 = [2, 13, 3]
    T6 = [3, 13, 4]
    T7 = [4, 13, 14]
    T8 = [4, 14, 5]
    T9 = [6, 4, 5]
    T10 = [6, 3, 4]
    T11 = [6, 2, 3]
    T12 = [6, 16, 2]
    T13 = [16, 15, 2]
    T14 = [15, 1, 2]
    T15 = [15, 10, 1]
    T16 = [9, 10, 15]
    T17 = [9, 15, 8]
    T18 = [9, 8, 7]
    T19 = [7, 8, 16]
    T20 = [8, 15, 16]
    T21 = [7, 16, 6]
    pts = [a, b, c, d, e, f, g, h, i, j, k, ℓ, m, n, o, p]
    tri = DT.Triangulation(pts)
    for T in (T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21)
        add_triangle!(tri, T)
    end
    return tri
end

function Base.:(==)(stats1::DT.TriangulationStatistics, stats2::DT.TriangulationStatistics)
    for f in fieldnames(DT.TriangulationStatistics)
        if f ≠ :individual_statistics
            if getfield(stats1, f) ≠ getfield(stats2, f)
                return false
            end
        else
            indiv_dict1 = stats1.individual_statistics
            indiv_dict2 = stats2.individual_statistics
            if length(indiv_dict1) ≠ length(indiv_dict2)
                return false
            end
            for (T, v) in indiv_dict1
                V, flag = DT.contains_triangle(T, keys(indiv_dict1))
                if !flag
                    return false
                end
                if v ≠ indiv_dict1[V]
                    return false
                end
            end
        end
    end
    return true
end


## TODO: Implement a brute-force DT.VoronoiTessellation that we can compare with
function validate_tessellation(vorn::DT.VoronoiTessellation; check_convex=true, check_adjacent=true)
    tri = DT.get_triangulation(vorn)
    for (i, p) in DT.get_generators(vorn)
        flag = get_point(tri, i) == get_generator(vorn, i) == p
        if !flag
            println("Generator $i is not correct, mapped to $p.")
            return false
        end
    end
    flag = DT.get_triangulation(vorn) == vorn.triangulation
    if !flag
        println("DT.Triangulation is not correct.")
        return false
    end
    circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
    triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
    for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
        V = DT.sort_triangle(V)
        c = DT.get_triangle_to_circumcenter(vorn, V)
        c = get_polygon_point(vorn, c)
        i, j, k = triangle_vertices(V)
        p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
        cx, cy = DT.triangle_circumcenter(p, q, r)
        flag = cx == c[1] && cy == c[2]
        if !flag
            println("Circumcenter of $V is not correct, mapped to $c.")
            return false
        end
    end
    for (c, V) in circumcenter_to_triangle
        flag = DT.get_circumcenter_to_triangle(vorn, c) == V
        if !flag
            println("Circumcenter $c is not correct, mapped to $V.")
            return false
        end
        flag = DT.get_triangle_to_circumcenter(vorn, V) == c
        if !flag
            println("Triangle $V is not correct, mapped to $c.")
            return false
        end
    end
    for (V, c) in triangle_to_circumcenter
        flag = DT.get_circumcenter_to_triangle(vorn, c) == V
        if !flag # cocircular points disrupt these dicts somewhat
            _flag = c ∈ vorn.cocircular_circumcenters
            if !_flag
                println("Circumcenter $c is duplicated from another triangle, but not in cocircular_circumcenters.")
                return false
            end
            i, j, k = triangle_vertices(V)
            for i in (i, j, k)
                for e in each_edge(get_adjacent2vertex(DT.get_triangulation(vorn), i))
                    v, w = DT.edge_vertices(e)
                    V′ = DT.sort_triangle(DT.construct_triangle(DT.triangle_type(DT.get_triangulation(vorn)), i, v, w))
                    flag = flag || DT.get_triangle_to_circumcenter(vorn, V′) == c
                end
            end
        end
        if !flag
            println("Triangle $V is not correct, mapped to $c.")
            return false
        end
        flag = DT.get_triangle_to_circumcenter(vorn, V) == c
        if !flag
            println("Circumcenter $c is not correct, mapped to $V.")
            return false
        end
    end
    for i in each_polygon_index(vorn)
        A = get_area(vorn, i)
        flag = A ≥ 0
        if !flag
            println("Polygon $i has area $A.")
            return false
        end
        if !isfinite(A)
            flag = i ∈ DT.get_unbounded_polygons(vorn)
            if !flag
                println("Polygon $i has infinite area, but is not in the unbounded polygons.")
                return false
            end
        end
    end
    if check_adjacent
        for i in each_polygon_index(vorn)
            C = get_polygon(vorn, i)
            ne = DT.num_boundary_edges(C)
            for j in 1:ne
                u = get_boundary_nodes(C, j)
                v = get_boundary_nodes(C, j + 1)
                flag = get_adjacent(vorn, u, v) == get_adjacent(vorn, DT.construct_edge(DT.edge_type(DT.get_triangulation(vorn)), u, v)) == i
                if !flag
                    println("Polygon $i is not adjacent to points $u and $v.")
                    return false
                end
            end
        end
    end
    for i in each_polygon_index(vorn)
        if i ∉ DT.get_unbounded_polygons(vorn)
            verts = get_polygon(vorn, i)
            poly_points = get_polygon_point.(Ref(vorn), verts)
            flag = @views allunique(poly_points[begin:end-1])
            if !flag
                println("Polygon $i has repeated vertices.")
                return false
            end
            if check_convex
                _pts = unique(poly_points)
                ch = convex_hull(_pts)
                _poly_points = _pts[DT.get_vertices(ch)]
                flag = DT.circular_equality(collect.(poly_points), collect.(_poly_points), ≈)
                if !flag
                    println("Polygon $i is not convex.")
                    return false
                end
            end
        end
    end
    if isempty(DT.get_unbounded_polygons(vorn))
        A = 0.0
        for i in each_polygon_index(vorn)
            A += get_area(vorn, i)
        end
        @test isapprox(A, get_area(vorn.triangulation), rtol=1e-4)
    end
    return true
end

function _make_graph_from_adjacency(A, labels=Dict(axes(A, 1) .=> axes(A, 1)))
    g = DT.Graph{Int}()
    foreach(axes(A, 1)) do i
        push!(g.vertices, labels[i])
        g.neighbours[labels[i]] = Set{Int}()
    end
    for j in axes(A, 2)
        for i in axes(A, 1)
            flag = A[i, j]
            if flag == 1
                push!(g.vertices, labels[i])
                push!(g.vertices, labels[j])
                (labels[j], labels[i]) ∉ g.edges && push!(g.edges, (labels[i], labels[j]))
                if labels[i] ∉ keys(g.neighbours)
                    g.neighbours[labels[i]] = Set{Int}(labels[j])
                else
                    push!(g.neighbours[labels[i]], labels[j])
                end
                if labels[j] ∉ keys(g.neighbours)
                    g.neighbours[labels[j]] = Set{Int}(labels[i])
                else
                    push!(g.neighbours[labels[j]], labels[i])
                end
            end
        end
    end
    return g
end

function _slow_compute_concentric_shell_ternary_split_position(p, q)
    m = t -> p .+ t .* (q .- p)
    L = t -> norm(m(t) .- p)
    R = t -> norm(m(t) .- q)
    tₙ = n -> exp2(n) / norm(q .- p)
    n_upp = ceil(log2(norm(q .- p)))
    obj = n -> (L(tₙ(n)) - R(tₙ(n)))^2
    n_range = -64:n_upp
    obj_vals = obj.(n_range)
    t_vals = tₙ.(n_range)
    return t_vals[argmin(obj_vals)]
end

_orient(a, b, c) = det([a[1] a[2] 1; b[1] b[2] 1; c[1] c[2] 1])
function _validate_offcenter(p, q, r, β)
    c = DT.triangle_circumcenter(p, q, r)
    offcenter = DT.triangle_offcenter(p, q, r, c, β)
    p, q, r = DT.make_shortest_edge_first(p, q, r, DT.squared_triangle_lengths_and_smallest_index(p, q, r)[4])
    if DT.squared_triangle_lengths_and_smallest_index(p, q, r)[1] ≈ DT.squared_triangle_lengths_and_smallest_index(p, q, r)[2]
        p, q, r = DT.select_shortest_edge_for_offcenter(p, q, r, c, norm(p .- q)^2)
    end
    m = (p .+ q) ./ 2

    # Check orientations 
    @test _orient(m, c, offcenter) ≈ 0.0 atol = 1e-6
    @test _orient(p, q, offcenter) > 0

    # Check the radius-edge ratios 
    ρ = DT.triangle_radius_edge_ratio(p, q, offcenter)
    @test ρ ≈ β
end

function get_both_offcenter_triangle_candidates(p, q, r, β)
    h1 = (2β + sqrt(4β^2 - 1)) / 2
    h2 = 1 / (2 * sqrt(4β^2 - 1))
    L1, L2, L3, idx = DT.squared_triangle_lengths_and_smallest_index(p, q, r)
    L1 = sqrt(L1)
    p, q, r = DT.make_shortest_edge_first(p, q, r, idx)
    c = DT.triangle_circumcenter(p, q, r)
    dir = c .- (p .+ q) ./ 2
    dir = dir ./ norm(dir)
    offcenter1 = (p .+ q) ./ 2 .+ dir .* h1 .* L1
    offcenter2 = (p .+ q) ./ 2 .+ dir .* h2 .* L1
    return [p, q, offcenter1], [p, q, offcenter2]
end

function reflect_point_across_line(p, q, u)
    orig_p = p
    q = collect(q .- p)
    u = collect(u .- p)
    p = zeros(2)
    pq_θ = -atan(q[2], q[1])
    cwrot_θ = [cos(pq_θ) -sin(pq_θ); sin(pq_θ) cos(pq_θ)]
    rotu = cwrot_θ * u
    # pq is now a horizontal line on the x-axis. So, to reflect rotu about this line, just do 
    rotu = [rotu[1], -rotu[2]]
    # Now, rotate back
    rotu = cwrot_θ' * rotu
    u = rotu .+ orig_p
    return Tuple(u)
end
function compute_diametral_circle(p, q)
    m = (p .+ q) ./ 2
    θ = LinRange(0, 2π, 250)
    rad = norm(p .- q) ./ 2
    circle = tuple.(m[1] .+ rad * cos.(θ), m[2] .+ rad * sin.(θ))
    lower_circle = circle[1:125, :] |> vec
    upper_circle = circle[126:end, :] |> vec
    unique!(circle)
    push!(circle, circle[begin])
    return upper_circle, lower_circle, circle
end
function compute_diametral_lens(p, q, lens_angle)
    orig_p = p
    orig_q = q
    q = collect(q .- p)
    p = zeros(2)
    pq_θ = -atan(q[2], q[1])
    cwrot_θ = [cos(pq_θ) -sin(pq_θ); sin(pq_θ) cos(pq_θ)]
    p = cwrot_θ * p
    q = cwrot_θ * q
    m = (p .+ q) ./ 2
    pq_normal = (q[2] - p[2], q[1] - p[1])
    pq_normal = pq_normal ./ norm(pq_normal)
    a = norm(p .- q) / 2
    c = a / cosd(lens_angle)
    b = c * sind(lens_angle)
    u = m .+ b .* pq_normal
    circumcenter = DT.triangle_circumcenter(p, q, u)
    circumradius = DT.triangle_circumradius(p, q, u)
    p_angle = atan(p[2] - circumcenter[2], p[1] - circumcenter[1])
    q_angle = atan(q[2] - circumcenter[2], q[1] - circumcenter[1])
    θ = LinRange(p_angle, q_angle, 250)
    lens = hcat(circumcenter[1] .+ circumradius * cos.(θ), circumcenter[2] .+ circumradius * sin.(θ))
    reflected_lens = reflect_point_across_line.(Ref(p), Ref(q), eachrow(lens))
    reflected_lens = hcat(first.(reflected_lens), last.(reflected_lens))
    combined_lens = vcat(lens, reverse(reflected_lens))
    rerotated_and_shifted_lens = (cwrot_θ' * lens')'
    rerotated_and_shifted_reflected_lens = (cwrot_θ' * reflected_lens')'
    rerotated_and_shifted_combined_lens = (cwrot_θ' * combined_lens')'
    rerotated_and_shifted_lens = tuple.(rerotated_and_shifted_lens[:, 1] .+ orig_p[1], rerotated_and_shifted_lens[:, 2] .+ orig_p[2])
    rerotated_and_shifted_reflected_lens = tuple.(rerotated_and_shifted_reflected_lens[:, 1] .+ orig_p[1], rerotated_and_shifted_reflected_lens[:, 2] .+ orig_p[2])
    rerotated_and_shifted_combined_lens = vcat(rerotated_and_shifted_lens, reverse(rerotated_and_shifted_reflected_lens))
    unique!(rerotated_and_shifted_combined_lens)
    push!(rerotated_and_shifted_combined_lens, rerotated_and_shifted_combined_lens[begin])
    return rerotated_and_shifted_lens, rerotated_and_shifted_reflected_lens, rerotated_and_shifted_combined_lens
end
function get_points_in_diametral_circle(p, q)
    tri = triangulate(rand(2, 50))
    args = DT.RefinementArguments(tri; use_lens=false)
    θ = LinRange(0, 2π, 250)
    r = LinRange(0, norm(p .- q) / 2, 250)
    m = (p .+ q) ./ 2
    points = [(m[1] + r * cos(θ), m[2] + r * sin(θ)) for r in r, θ in θ] |> vec
    classifications = DT.encroaches_upon.(Ref(p), Ref(q), points, Ref(args))
    return points[classifications], points[.!classifications]
end
function get_points_in_diametral_lens(p, q, lens_angle)
    tri = triangulate(rand(2, 50))
    args = DT.RefinementArguments(tri; use_lens=true, min_angle=lens_angle)
    θ = LinRange(0, 2π, 250)
    r = LinRange(0, norm(p .- q) / 2, 250)
    m = (p .+ q) ./ 2
    points = [(m[1] + r * cos(θ), m[2] + r * sin(θ)) for r in r, θ in θ] |> vec
    classifications = DT.encroaches_upon.(Ref(p), Ref(q), points, Ref(args))
    return points[classifications], points[.!classifications]
end

function get_random_convex_polygon(points)
    tri = triangulate(points)
    S = get_convex_hull_vertices(tri)
    pop!(S) # Want S[begin] ≠ S[end]
    return S
end

function _lexicographically_sort_pair_vector(pairs)
    lt = (x, y) -> let (key_x, value_x) = x, (key_y, value_y) = y
        if value_x < value_y
            return true
        elseif value_x > value_y
            return false
        else
            return key_x < key_y
        end
    end
    return sort(pairs, lt=lt)
end
function _compare_pairs(pairs, dt_pairs)
    @test length(pairs) == length(dt_pairs)
    @test _lexicographically_sort_pair_vector(pairs) == _lexicographically_sort_pair_vector(dt_pairs)
end

function is_sink(tri::Triangulation, T)
    _, _, t = DT.triangle_angles(get_point(tri, T...)...)
    return t ≤ π / 2 || DT.is_boundary_triangle(tri, T)
end
function is_interior_sink(tri::Triangulation, T)
    return is_sink(tri, T) && !DT.is_boundary_triangle(tri, T)
end

function get_idx(tri, p, q, r)
    i = findfirst(==(p), get_points(tri))
    j = findfirst(==(q), get_points(tri))
    k = findfirst(==(r), get_points(tri))
    return i, j, k
end

⊢(x, y) = x == y && typeof(x) == typeof(y)

function manual_insertion_event_history(tri::Triangulation, orig_tri::Triangulation)
    # Initialise
    added_triangles = Set{DT.triangle_type(tri)}()
    deleted_triangles = Set{DT.triangle_type(tri)}()
    added_segments = Set{DT.edge_type(tri)}()
    deleted_segments = Set{DT.edge_type(tri)}()
    added_boundary_segments = Set{DT.edge_type(tri)}()
    deleted_boundary_segments = Set{DT.edge_type(tri)}()
    # Find added or deleted triangles
    for T in each_triangle(tri)
        _, flag = DT.contains_triangle(orig_tri, T)
        !flag && push!(added_triangles, T)
    end
    for T in each_triangle(orig_tri)
        _, flag = DT.contains_triangle(tri, T)
        !flag && push!(deleted_triangles, T)
    end
    # Find added or deleted segments 
    for e in DT.get_interior_segments(tri)
        flag = DT.contains_unoriented_edge(e, DT.get_interior_segments(orig_tri))
        !flag && push!(added_segments, e)
    end
    for e in DT.get_interior_segments(orig_tri)
        flag = DT.contains_unoriented_edge(e, DT.get_interior_segments(tri))
        !flag && push!(deleted_segments, e)
    end
    # Find added or delete boundary segments 
    for e in DT.each_boundary_edge(tri)
        flag = DT.contains_unoriented_edge(e, DT.each_boundary_edge(orig_tri))
        !flag && push!(added_boundary_segments, e)
    end
    for e in DT.each_boundary_edge(orig_tri)
        flag = DT.contains_unoriented_edge(e, DT.each_boundary_edge(tri))
        !flag && push!(deleted_boundary_segments, e)
    end
    # Done 
    return DT.InsertionEventHistory(added_triangles, deleted_triangles, added_segments, deleted_segments, added_boundary_segments, deleted_boundary_segments)
end

function validate_insertion_event_history(tri::Triangulation, orig_tri::Triangulation, history::DT.InsertionEventHistory)
    manual_history = manual_insertion_event_history(tri, orig_tri)
    @test length(manual_history.added_triangles) == length(history.added_triangles)
    @test length(manual_history.deleted_triangles) == length(history.deleted_triangles)
    @test length(manual_history.added_segments) == length(history.added_segments)
    @test length(manual_history.deleted_segments) == length(history.deleted_segments)
    @test length(manual_history.added_boundary_segments) == length(history.added_boundary_segments)
    @test length(manual_history.deleted_boundary_segments) == length(history.deleted_boundary_segments)
    @test DT.compare_triangle_collections(manual_history.added_triangles, history.added_triangles)
    @test DT.compare_unoriented_edge_collections(manual_history.added_segments, history.added_segments)
    @test DT.compare_unoriented_edge_collections(manual_history.added_boundary_segments, history.added_boundary_segments)
end

_approx_ispow2(x) = ispow2(x) || isapprox(log2(x), round(log2(x)), atol=1e-9, rtol=1e-9)

function compare_encroach_queues(args::DT.RefinementArguments, manual_enqueue)
    _manual_enqueue_pairs = collect(manual_enqueue)
    _args_queue_segments_pairs = collect(args.queue.segments)
    for pars in (_manual_enqueue_pairs, _args_queue_segments_pairs)
        for i in eachindex(pars)
            u, v = pars[i].first
            if u > v
                pars[i] = (v, u) => pars[i].second
            end
        end
    end
    _compare_pairs(_manual_enqueue_pairs, _args_queue_segments_pairs)
end

function compare_triangle_queues(args::DT.RefinementArguments, manual_enqueue)
    _manual_enqueue_pairs = collect(manual_enqueue)
    _args_queue_triangle_pairs = collect(args.queue.triangles)
    for pars in (_manual_enqueue_pairs, _args_queue_triangle_pairs)
        for i in eachindex(pars)
            T = pars[i].first
            pars[i] = DT.sort_triangle(T) => pars[i].second
        end
    end
    _compare_pairs(_manual_enqueue_pairs, _args_queue_triangle_pairs)
end

function is_conformal(tri::Triangulation)
    points = get_points(tri)
    segment_tree = DT.BoundaryRTree(points)
    for e in each_segment(tri)
        i, j = DT.edge_vertices(e)
        insert!(segment_tree, i, j)
    end
    for r in each_solid_vertex(tri)
        intersects = DT.get_intersections(segment_tree, r, cache_id=1) # can't use multithreading here
        c = get_point(tri, r)
        for box in intersects
            i, j = DT.get_edge(box)
            any(==(r), (i, j)) && continue
            a, b = get_point(tri, i, j)
            flag = DT.is_inside(DT.point_position_relative_to_diametral_circle(a, b, c))
            if flag
                #cert = test_visibility(tri, segment_tree, i, j, r)
                cert = DT.test_visibility(tri, i, j, r)
                flag = DT.is_visible(cert)
            end
            if flag
                println("Triangulation is not conformal, as vertex $r is in the diametral circle of the edge $((i, j))")
                return false
            end
        end
    end
    return true
end

function slow_triangle_assess(tri, args)
    good_T = NTuple{3,Int}[]
    bad_T = NTuple{3,Int}[]
    for T in each_solid_triangle(tri)
        u, v, w = T
        p, q, r = get_point(tri, u, v, w)
        t1, t2, t3 = DT.triangle_angles(p, q, r)
        t = min(t1, t2, t3) * 180 / π
        t′ = max(t1, t2, t3) * 180 / π
        A = DT.triangle_area(p, q, r)
        seditious = DT.is_triangle_seditious(tri, args, u, v, w, DT.squared_triangle_lengths_and_smallest_index(p, q, r)[4])
        nestled = DT.is_triangle_nestled(tri, T, DT.squared_triangle_lengths_and_smallest_index(p, q, r)[4])
        ρ, qual = DT.assess_triangle_quality(tri, args, T)
        @inferred DT.assess_triangle_quality(tri, args, T)
        condition = ((A > args.constraints.max_area) || (t < args.constraints.min_angle && ρ > args.constraints.max_radius_edge_ratio && !seditious && !nestled)) && A ≥ args.constraints.min_area && !args.constraints.custom_constraint(tri, T) && t′ ≤ args.constraints.max_angle
        @test ρ == DT.triangle_radius_edge_ratio(p, q, r)
        @test qual == condition
        if qual
            push!(bad_T, T)
        else
            push!(good_T, T)
        end
    end
    return good_T, bad_T
end

function slow_triangle_assess_queue(tri, args)
    good_T, bad_T = slow_triangle_assess(tri, args)
    queue = PriorityQueue{NTuple{3,Int},Float64}(Base.Order.Reverse)
    for T in bad_T
        queue[T] = DT.triangle_radius_edge_ratio(get_point(tri, T...)...)
    end
    return queue
end

function compute_φmin(tri)
    # φmin is the minimum angle between two adjoining segments 
    φmin = Inf
    locked = !DT.has_boundary_nodes(tri)
    if !DT.has_boundary_nodes(tri)
        lock_convex_hull!(tri)
    end
    for e in each_segment(tri)
        for η in each_segment(tri)
            if e ≠ η
                u, v = e
                x, y = η
                w = DT.get_shared_vertex(e, η)
                w == DT.∅ && continue
                if w == u && w == x
                    e1 = (u, v)
                    e2 = (x, y)
                    i, j, k = w, v, y
                elseif w == u && w == y
                    e1 = (u, v)
                    e2 = (y, x)
                    i, j, k = w, v, x
                elseif w == v && w == x
                    e1 = (v, u)
                    e2 = (x, y)
                    i, j, k = w, u, y
                else # if w == v && w == y
                    e1 = (v, u)
                    e2 = (y, x)
                    i, j, k = w, u, x
                end
                p, q, r = get_point(tri, i, j, k)
                s1 = q .- p
                s2 = r .- p
                rat = dot(s1, s2) / (norm(s1) * norm(s2))
                if rat > 1
                    rat = 1
                elseif rat < -1
                    rat = -1
                end
                θ = acosd(rat)
                φmin = min(φmin, θ)
            end
        end
    end
    locked && unlock_convex_hull!(tri)
    return φmin
end

function validate_refinement(tri, args; check_conformal=true, warn=true)
    ## Things to check:
    ##   1. All angle constraints are met, except for seditious and nestled triangles.
    ##   2. If !use_lens, the triangulation is conformal. 
    ##   3. All area constraints are met.
    ##   4. All custom constraints are met.
    ##   5. The triangulation is Delaunay.
    ##   6. All circumcenters are inside the domain, if !use_lens. 
    ##   7. There are no more points than there are maximum points.
    ##   8. There are no more encroached edges.
    ##   9. There are no more bad triangles.
    minθ = Inf
    maxθ = -Inf
    all_segments = Int[]
    for e in each_segment(tri)
        push!(all_segments, e...)
    end
    for T in each_solid_triangle(tri)
        is_in_corner = all(∈(all_segments), triangle_vertices(T)) # triangles in corners sometimes can't have their circumcenters inserted if doing use_lens, so their centroid gets inserted, but this breaks the constraint conditions sometimes for minimum angles.
        # Get some stats
        t1, t2, t3 = DT.triangle_angles(get_point(tri, T...)...)
        t1, t2, t3 = (t1, t2, t3) .* 180 ./ π
        A = DT.triangle_area(get_point(tri, T...)...)
        seditious = DT.is_triangle_seditious(tri, args, T..., DT.squared_triangle_lengths_and_smallest_index(get_point(tri, T...)...)[4])
        nestled = DT.is_triangle_nestled(tri, T, DT.squared_triangle_lengths_and_smallest_index(get_point(tri, T...)...)[4])
        _flag, steiner_point = DT.get_steiner_point(tri, args, T)
        show_print = true
        if !DT.is_none(_flag)
            warn && @warn "Steiner point found for triangle $T has precision issues. Proceeding with the tests, but will not return any results for this triangle."
            show_print = warn
        else
            # Only update extrema for good triangles 
            minθ = min(minθ, t1, t2, t3)
            maxθ = max(maxθ, t1, t2, t3)
        end
        V = find_triangle(tri, steiner_point)
        flag = DT.point_position_relative_to_triangle(tri, V, steiner_point)
        if DT.is_on(flag) && DT.is_ghost_triangle(V)
            V = DT.replace_ghost_triangle_with_boundary_triangle(tri, V)
        end
        if !(seditious || nestled || is_in_corner) && min(t1, t2, t3) < args.constraints.min_angle
            show_print && println("Triangle $T has minimum angle $(min(t1, t2, t3)), which is less than the minimum angle constraint $(args.constraints.min_angle).")
            DT.is_none(_flag) && return false
        end
        if max(t1, t2, t3) > args.constraints.max_angle
            show_print && println("Triangle $T has maximum angle $(max(t1, t2, t3)), which is greater than the maximum angle constraint $(args.constraints.max_angle).")
            return false
        end
        flag = args.constraints.min_area - 1e-16 ≤ A ≤ args.constraints.max_area
        if !flag
            show_print && println("Triangle $T has area $A, which is not in the range $(args.constraints.min_area) to $(args.constraints.max_area).")
            DT.is_none(_flag) && return false
        end
        flag = !args.constraints.custom_constraint(tri, T)
        if !flag
            show_print && println("Triangle $T violates the custom constraint.")
            DT.is_none(_flag) && return false
        end
        if !args.use_lens
            flag = DT.dist(tri, steiner_point) > -eps(Float64)
            if !flag
                show_print && println("The Steiner point associated with the triangle $T is a distance $(-DT.dist(tri,steiner_point)) away from the domain.")
                DT.is_none(_flag) && return false
            end
        end
        for (u, v) in DT.triangle_edges(T)
            if DT.contains_segment(tri, u, v)
                flag = DT.is_encroached(tri, args, (u, v))
                if flag
                    show_print && println("The edge ($u, $v) is encroached.")
                    DT.is_none(_flag) && return false
                end
            end
        end
        flag = DT.assess_triangle_quality(tri, args, T)[2]
        if flag && !seditious && !nestled && !is_in_corner
            show_print && println("Triangle $T is of poor quality.")
            DT.is_none(_flag) && return false
        end
        ρ = DT.triangle_radius_edge_ratio(get_point(tri, T...)...)
        flag_2 = flag == ((A > args.constraints.max_area) || (min(t1, t2, t3) < args.constraints.min_angle && ρ > args.constraints.max_radius_edge_ratio && !seditious && !nestled && !is_in_corner)) && A ≥ args.constraints.min_area && !args.constraints.custom_constraint(tri, T) && max(t1, t2, t3) ≤ args.constraints.max_angle
        if !flag_2 && !is_in_corner
            if flag
                show_print && println("Triangle $T was incorrectly marked as being of good quality.")
            else
                show_print && println("Triangle $T was incorrectly marked as being of poor quality.")
            end
            DT.is_none(_flag) && return false
        end
    end
    if !args.use_lens && check_conformal
        flag = is_conformal(tri)
        if !flag
            println("The triangulation is not conformal.")
            return false
        end
    end
    flag = validate_triangulation(tri)
    if !flag
        println("The triangulation is not Delaunay.")
        return false
    end
    flag = DT.num_solid_vertices(tri) > args.constraints.max_points
    if flag
        println("The number of points in the triangulation ($(DT.num_points(tri))) is greater than the maximum number of points allowed ($(args.constraints.max_points)).")
        return false
    end
    if !args.use_lens
        φmin = compute_φmin(tri)
        ρ = args.constraints.max_radius_edge_ratio
        τmin = sind(φmin) / sqrt(5 - 4cosd(φmin))
        τmax = 180 - 2asind(1 / (2ρ))
        flag = minθ ≥ τmin
        if !flag
            println("The minimum angle $(minθ) is not in the range $(τmin) to $(τmax).")
            return false
        end
        flag = maxθ ≤ τmax
        if !flag
            println("The maximum angle $(maxθ) is not in the range $(τmin) to $(τmax).")
            return false
        end
    end
    return true
end
validate_refinement(tri; check_conformal=true, warn=true, kwargs...) = validate_refinement(tri, DT.RefinementArguments(tri; kwargs...); warn, check_conformal)

function why_not_equal(tri1, tri2)
    !DT.has_ghost_triangles(tri1) && DT.has_ghost_triangles(tri2) && println("!has_ghost_triangles(tri1) && has_ghost_triangles(tri2)")
    !DT.has_ghost_triangles(tri2) && DT.has_ghost_triangles(tri1) && println("!has_ghost_triangles(tri2) && has_ghost_triangles(tri1)")
    DT.get_points(tri1) ≠ DT.get_points(tri2) && println("get_points(tri1) ≠ get_points(tri2)")
    !DT.compare_triangle_collections(DT.get_triangles(tri1), DT.get_triangles(tri2)) && println("!compare_triangle_collections(get_triangles(tri1), get_triangles(tri2))")
    !DT.compare_unoriented_edge_collections(DT.get_interior_segments(tri1), DT.get_interior_segments(tri2)) && println("!compare_unoriented_edge_collections(get_interior_segments(tri1), get_interior_segments(tri2))")
    !DT.compare_unoriented_edge_collections(DT.get_all_segments(tri1), DT.get_all_segments(tri2)) && println("!compare_unoriented_edge_collections(get_all_segments(tri1), get_all_segments(tri2))")
    DT.get_adjacent(tri1) ≠ get_adjacent(tri2) && println("get_adjacent(tri1) ≠ get_adjacent(tri2)")
    DT.get_adjacent2vertex(tri1) ≠ DT.get_adjacent2vertex(tri2) && println("get_adjacent2vertex(tri1) ≠ get_adjacent2vertex(tri2)")
    DT.get_graph(tri1) ≠ DT.get_graph(tri2) && println("get_graph(tri1) ≠ get_graph(tri2)")
    DT.get_boundary_edge_map(tri1) ≠ DT.get_boundary_edge_map(tri2) && println("get_boundary_edge_map(tri1) ≠ get_boundary_edge_map(tri2)")
    DT.get_ghost_vertex_map(tri1) ≠ DT.get_ghost_vertex_map(tri2) && println("get_ghost_vertex_map(tri1) ≠ get_ghost_vertex_map(tri2)")
    DT.get_ghost_vertex_ranges(tri1) ≠ DT.get_ghost_vertex_ranges(tri2) && println("get_ghost_vertex_ranges(tri1) ≠ get_ghost_vertex_ranges(tri2)")
    DT.get_convex_hull(tri1) ≠ DT.get_convex_hull(tri2) && println("get_convex_hull(tri1) ≠ get_convex_hull(tri2)")
    rep1 = DT.get_representative_point_list(tri1)
    rep2 = DT.get_representative_point_list(tri2)
    length(rep1) ≠ length(rep2) && println("length(rep1) ≠ length(rep2)")
    for i in 1:DT.num_curves(tri1)
        p1 = DT.get_representative_point_coordinates(tri1, i)
        p2 = DT.get_representative_point_coordinates(tri2, i)
        !([getx(p1), gety(p1)] ≈ [getx(p2), gety(p2)]) && println("!([getx(p1), gety(p1)] ≈ [getx(p2), gety(p2)]) for $i")
    end
    DT.get_polygon_hierarchy(tri1) ≠ DT.get_polygon_hierarchy(tri2) && println("get_polygon_hierarchy(tri1) ≠ get_polygon_hierarchy(tri2)")
    DT.get_boundary_nodes(tri1) ≠ DT.get_boundary_nodes(tri2) && println("get_boundary_nodes(tri1) ≠ get_boundary_nodes(tri2)")
end

#=
These weighted Delaunay examples come from the following MATLAB code, which exploits 
the relationship between the convex hull of the lifting map and the Delaunay triangulation.
1:78 are random triangulations, and 79:155 are triangulations of convex polygons.

For computing these triangulations, we first define the MATLAB functions:

    function [cells, lifted] = LiftingMap(points,weights)
    num = size(points, 1);
    lifted = zeros(num, 3);
    for i = 1:num
        p = points(i, :);
        lifted(i, :) = [p sum(p.^2) - weights(i)];
    end
    ch = convhull(lifted);
    sep_ch = reshape(lifted(ch, :), [size(ch) 3]);
    vecs = diff(sep_ch, 1, 2);
    n = cross(vecs(:, 1, :), vecs(:, 2, :), 3);
    downward_faces = n(:, :, 3) <= 0;
    cells = ch(downward_faces, :);
    end
    function save_to_file(points, weights, n)
    [lift_cells, ~] = LiftingMap(points,weights);
    res = [(1+size(points, 1)) * ones(1, 3); [points weights]; lift_cells];
    % so the points are in res(2:res(1), 1:2), the weights are in res(2:res(1), 3),
    % and the triangles are in res((res(1)+1):end, :).
    writematrix(res, fullfile(pwd, 'tempmats', strcat(string(n), '.txt')));
    end

Then, for the random triangulations, we use:

    digits(16)
    rng(123)
    ctr = 1;
    for pass = 1:3
        if pass == 1
            range = 4:20;
        elseif pass == 2
            range = 21:25:1500;
        else 
            range = 15000;
        end
        for n = range
            points = randn(n, 2);
            weights = randn(n, 1);
            s = 10*rand;
            u = rand;
            if u < 1/250
                weights = 0 * weights;
            elseif (1/250 <= u) && (u < 1/5)
                weights = weights.^2;
            elseif (1/5 <= u) && (u < 0.6)
                weights = s * weights;
            elseif (u <= 0.6) && (u < 0.65)
                sparse_prob = rand(n);
                weights(sparse_prob >= 0.95) = 0;
            end
            save_to_file(points, weights, ctr)
            ctr = ctr + 1;
        end
    end

The triangulations of convex polygons are computed with the code:

    digits(16)
    rng(123)
    ctr = 79;
    for pass = 1:3
        if pass == 1
            range = 6:20;
        elseif pass == 2
            range = 21:25:1500;
        else
            range = [100, 15000];
        end
        for n = range
            if pass < 3
                n0 = 0;
                while n0 <= 3
                    points = randn(n, 2);
                    convex_poly = points(convhull(points), :);
                    points = convex_poly(1:end-1, :);
                    n0 = size(points, 1);
                end
                weights = randn(n0, 1);
            else 
                theta = linspace(0, 2*pi, n);
                x = cos(theta);
                y = sin(theta);
                points = [x(:), y(:)];
                points = points(1:end-1, :);
                n0 = size(points, 1);
                weights = randn(n0, 1);
            end
            s = 10*rand;
            u = rand;
            if u < 1/250
                weights = 0 * weights;
            elseif (1/250 <= u) && (u < 1/5)
                weights = weights.^2;
            elseif (1/5 <= u) && (u < 0.6)
                weights = s * weights;
            elseif (u <= 0.6) && (u < 0.65)
                sparse_prob = rand(n0);
                weights(sparse_prob >= 0.95) = 0;
            end
            save_to_file(points, weights, ctr)
            ctr = ctr + 1;
        end
    end
=#
const WGT_DIR = joinpath(dirname(dirname(pathof(DelaunayTriangulation))), "test", "triangulation", "weighted_delaunay_mats")
const NUM_WEGT = length(readdir(WGT_DIR))
const NUM_CWEGT = NUM_WEGT - 78 # number of convex polygon weighted examples
function get_weighted_example(i)
    # get the points, weights, and triangles
    file = "$i.txt"
    mat = readdlm(joinpath(WGT_DIR, file), ',')
    num_pts = Int(mat[1, 1])
    points = mat[2:num_pts, 1:2]
    weights = mat[2:num_pts, 3]
    triangles = Int.(mat[num_pts+1:end, :])
    triangles = Set([Tuple(tri) for tri in eachrow(triangles)])
    _triangles = empty(triangles)
    for T in triangles
        i, j, k = T
        p, q, r = eachrow(points[[i, j, k], :])
        or = DT.triangle_orientation(p, q, r)
        if !DT.is_positively_oriented(or)
            push!(_triangles, (j, i, k))
        else
            push!(_triangles, (i, j, k))
        end
    end
    triangles = _triangles
    points = ElasticMatrix(points')

    # since some of the vertices might be submerged, we can't just get the 
    # convex hull of points to find the boundary (in cases of collinear points). instead, let's find 
    # all the edges which have only one adjoining triangle, and then sort them.
    d = Dict{NTuple{2,Int},Int}()
    for T in triangles
        i, j, k = T
        d[(i, j)] = k
        d[(j, k)] = i
        d[(k, i)] = j
    end
    boundary_nodes = Set{Int}()
    for (i, j) in keys(d)
        if !haskey(d, (j, i))
            push!(boundary_nodes, i, j)
        end
    end
    boundary_nodes = collect(boundary_nodes)
    DT.sort_convex_polygon!(boundary_nodes, points)
    push!(boundary_nodes, boundary_nodes[1])

    # now get all the submerged vertices 
    submerged_vertices = Set{Int}()
    nonsubmerged_vertices = Set{Int}()
    for T in triangles
        push!(nonsubmerged_vertices, T...)
    end
    for i in axes(points, 2)
        if i ∉ nonsubmerged_vertices
            push!(submerged_vertices, i)
        end
    end

    # return 
    tri = Triangulation(points, triangles, boundary_nodes; weights)
    unlock_convex_hull!(tri)
    return (tri=tri, submerged_vertices=submerged_vertices, nonsubmerged_vertices=nonsubmerged_vertices, weights=weights)
end
get_unstructured_weighted_example(i) = get_weighted_example(i)
function get_convex_polygon_weighted_example(i)
    tri, submerged_vertices, nonsubmerged_vertices, weights = get_weighted_example(78 + i)
    S = eachindex(weights)
    representative_point_list = DT.get_representative_point_list(tri)
    cx, cy = DT.mean_points(get_points(tri), S)
    representative_point_list[1] = DT.RepresentativeCoordinates(cx, cy, length(S))
    return (tri=tri, submerged_vertices=submerged_vertices, nonsubmerged_vertices=nonsubmerged_vertices, weights=weights, S=S)
end

function get_all_distances_to_witness_planes(tri, i)
    distances = Dict{NTuple{3,Int},Float64}()
    for T in each_solid_triangle(tri)
        T = DT.sort_triangle(T)
        δ = DT.get_distance_to_witness_plane(tri, i, T)
        distances[T] = δ
    end
    return distances
end

function get_nearest_power_point(tri, i)
    all_dists = Dict{Int,Float64}()
    for j in each_solid_vertex(tri)
        p = DT.get_point(tri, i)
        q = DT.get_point(tri, j)
        all_dists[j] = norm(p .- q)^2 - get_weight(tri, i) - get_weight(tri, j)
    end
    return argmin(all_dists)
end

⪧(a::Vector{Tuple}, b::Vector{Tuple}; kwargs...) = ⪧(collect.(a), collect.(b); kwargs...)
⪧(a::Vector{<:Union{Vector,Number}}, b::Vector; kwargs...) = isapprox(a, b; kwargs...)
⪧(a, b; kwargs...) = isapprox(collect(collect.(a)), collect(collect.(b)); kwargs...)

function closest_point_on_curve(c, p)
    t = LinRange(0, 1, 50000)
    t′ = argmin(t) do τ
        q = c(τ)
        norm(p .- q)
    end
    return t′, c(t′)
end

function slow_arc_length(c, t₁, t₂)
    t = LinRange(t₁, t₂, 15000)
    s = 0.0
    for i in 1:(length(t)-1)
        s += norm(c(t[i+1]) .- c(t[i]))
    end
    return s
end

function slow_total_absolute_curvature(c, t₁, t₂)
    t = LinRange(t₁, t₂, 1500)
    s = 0.0
    for i in 1:(length(t)-1)
        T₁ = DT.differentiate(c, t[i])
        T₂ = DT.differentiate(c, t[i+1])
        _rat = dot(T₁, T₂) / (norm(T₁) * norm(T₂))
        _dot = _rat > 1 ? 1 : _rat < -1 ? -1 : _rat
        θ = acos(_dot)
        s += abs(θ)
    end
    return s
end

function slow_get_segment(control_points, knots, alpha, tension, t)
    L = 0.0
    for i in 2:lastindex(control_points)
        L += norm(control_points[i] .- control_points[i-1])^alpha
    end
    is_closed = control_points[begin] == control_points[end]
    s = 0
    for outer s in eachindex(knots)
        knots[s] ≤ t ≤ knots[s+1] && break
    end
    points = NTuple{2,Float64}[]
    if s == 1
        if is_closed
            push!(points, control_points[end-1])
        else
            push!(points, DT.extend_left_control_point(control_points))
        end
    else
        push!(points, control_points[s-1])
    end
    push!(points, control_points[s], control_points[s+1])
    if s == length(knots) - 1
        if is_closed
            push!(points, control_points[begin+1])
        else
            push!(points, DT.extend_right_control_point(control_points))
        end
    else
        push!(points, control_points[s+2])
    end
    return DT.catmull_rom_spline_segment(points..., alpha, tension), s
end

function slow_eval_bspline(control_points, knots, t)
    if t == 0
        return control_points[begin] |> collect
    elseif t == 1
        return control_points[end] |> collect
    end
    val = Real[0, 0]
    order = length(knots) - length(control_points)
    a, b = knots[order], knots[length(control_points)+1]
    t = a + (b - a) * t
    for i in eachindex(control_points)
        val = val .+ control_points[i] .* slow_eval_bspline_basis(knots, i, order, t)
    end
    return val
end
function slow_eval_bspline_basis(knots, i, order, t)
    if order == 1
        return (knots[i] ≤ t < knots[i+1]) ? 1.0 : 0.0
    else
        coeff1 = (t - knots[i]) / (knots[i+order-1] - knots[i])
        coeff2 = (knots[i+order] - t) / (knots[i+order] - knots[i+1])
        coeff1 = isfinite(coeff1) ? coeff1 : 0.0
        coeff2 = isfinite(coeff2) ? coeff2 : 0.0
        return coeff1 * slow_eval_bspline_basis(knots, i, order - 1, t) + coeff2 * slow_eval_bspline_basis(knots, i + 1, order - 1, t)
    end
end

function slow_bezier_eval(points, t)
    points = collect.(points)
    n = length(points) - 1
    return collect(sum([binomial(n, i) .* (1 - t)^(n - i) .* t^i .* points[i+1] for i in 0:n]))
end

function flatten_boundary_nodes(points, boundary_nodes, segments=nothing)
    _points = NTuple{2,Float64}[]
    if DT.has_multiple_curves(boundary_nodes)
        nc = DT.num_curves(boundary_nodes)
        for i in 1:nc
            _boundary_nodes = get_boundary_nodes(boundary_nodes, i)
            ns = DT.num_sections(_boundary_nodes)
            for j in 1:ns
                __boundary_nodes = get_boundary_nodes(_boundary_nodes, j)
                if !(get_boundary_nodes(__boundary_nodes, 1) isa DT.AbstractParametricCurve)
                    n = DT.num_boundary_edges(__boundary_nodes)
                    for k in 1:(n+1)
                        push!(_points, get_point(points, get_boundary_nodes(__boundary_nodes, k)))
                    end
                else
                    append!(_points, get_boundary_nodes(__boundary_nodes, 1).(LinRange(0, 1, 250)))
                end
            end
            if !(get_boundary_nodes(get_boundary_nodes(_boundary_nodes, 1), 1) isa DT.AbstractParametricCurve)
                push!(_points, get_point(points, get_boundary_nodes(get_boundary_nodes(_boundary_nodes, 1), 1)))
            end
            push!(_points, (NaN, NaN))
        end
    elseif DT.has_multiple_sections(boundary_nodes)
        ns = DT.num_sections(boundary_nodes)
        for j in 1:ns
            __boundary_nodes = get_boundary_nodes(boundary_nodes, j)
            if !(get_boundary_nodes(__boundary_nodes, 1) isa DT.AbstractParametricCurve)
                n = DT.num_boundary_edges(__boundary_nodes)
                for k in 1:n
                    push!(_points, get_point(points, get_boundary_nodes(__boundary_nodes, k)))
                end
            else
                append!(_points, get_boundary_nodes(__boundary_nodes, 1).(LinRange(0, 1, 250)))
            end
        end
        if !(get_boundary_nodes(get_boundary_nodes(boundary_nodes, 1), 1) isa DT.AbstractParametricCurve)
            push!(_points, get_point(points, get_boundary_nodes(get_boundary_nodes(boundary_nodes, 1), 1)))
        else
            push!(_points, get_boundary_nodes(get_boundary_nodes(boundary_nodes, 1), 1).(0.0))
        end
    else
        if !(get_boundary_nodes(boundary_nodes, 1) isa DT.AbstractParametricCurve)
            n = DT.num_boundary_edges(boundary_nodes)
            for k in 1:n
                push!(_points, get_point(points, get_boundary_nodes(boundary_nodes, k)))
            end
            push!(_points, get_point(points, get_boundary_nodes(boundary_nodes, 1)))
        else
            append!(_points, get_boundary_nodes(boundary_nodes, 1).(LinRange(0, 1, 250)))
        end
    end
    if !isnothing(segments)
        push!(_points, (NaN, NaN))
        for e in segments
            i, j = e
            p, q = get_point(points, i, j)
            push!(_points, p, q, (NaN, NaN))
        end
    end
    push!(_points, (NaN, NaN))
    for p in DT.each_point(points)
        push!(_points, p)
        push!(_points, (NaN, NaN))
    end
    return _points
end

function Base.:(==)(rect1::SI.Rect, rect2::DT.BoundingBox)
    ab = rect1.low
    cd = rect1.high
    ab′ = DT.get_bl_corner(rect2)
    cd′ = DT.get_tr_corner(rect2)
    return (ab === ab′) && (cd === cd′)
end
Base.:(==)(rect2::DT.BoundingBox, rect1::SI.Rect) = rect1 == rect2
function Base.:(==)(rect1::SI.SpatialElem, rect2::DT.DiametralBoundingBox)
    id = rect1.val
    id′ = DT.get_edge(rect2)
    return (id == id′) && rect1.mbr == DT.get_bounding_box(rect2)
end
Base.:(==)(rect2::DT.DiametralBoundingBox, rect1::SI.SpatialElem) = rect1 == rect2
function Base.:(==)(tree1::SI.RTree, tree2::DT.RTree)
    si_rects, si_els = get_si_rectangles(tree1)
    dt_rects, dt_els = get_dt_rectangles(tree2)
    return (si_rects == dt_rects) && (si_els == dt_els)
end
Base.:(==)(tree2::DT.RTree, tree1::SI.RTree) = tree1 == tree2
Base.:(==)(int::DT.RTreeIntersectionIterator, int_si::SI.RTreeRegionQueryIterator) = (isempty(int) && isempty(int_si)) || collect(int) == collect(int_si)
Base.:(==)(int_si::SI.RTreeRegionQueryIterator, int::DT.RTreeIntersectionIterator) = int == int_si

function get_si_rectangles(tree::SI.RTree)
    rects = Any[]
    els = Any[]
    return get_si_rectangles!(rects, els, tree.root)
end
function get_dt_rectangles(tree::DT.RTree)
    rects = DT.BoundingBox[]
    els = DT.DiametralBoundingBox[]
    return get_dt_rectangles!(rects, els, DT.get_root(tree))
end
function get_si_rectangles!(rects, els, node::SI.Branch)
    rect = node.mbr
    push!(rects, rect)
    for child in node.children
        get_si_rectangles!(rects, els, child)
    end
    return identity.(rects), identity.(els)
end
function get_dt_rectangles!(rects, els, node::DT.Branch)
    rect = DT.get_bounding_box(node)
    push!(rects, rect)
    for child in DT.get_children(node)
        get_dt_rectangles!(rects, els, child)
    end
    return identity.(rects), identity.(els)
end
function get_si_rectangles!(rects, els, node::SI.Leaf)
    rect = node.mbr
    push!(rects, rect)
    for child in node.children
        push!(els, child)
    end
    return identity.(rects), identity.(els)
end
function get_dt_rectangles!(rects, els, node::DT.Leaf)
    rect = DT.get_bounding_box(node)
    push!(rects, rect)
    for child in DT.get_children(node)
        push!(els, child)
    end
    return identity.(rects), identity.(els)
end

to_lines(rect::DT.BoundingBox) = [(rect.x.a, rect.y.a), (rect.x.b, rect.y.a), (rect.x.b, rect.y.b), (rect.x.a, rect.y.b), (rect.x.a, rect.y.a)]
to_lines(rect::SI.Rect) =
    let a = rect.low[1], c = rect.low[2], b = rect.high[1], d = rect.high[2]
        [(a, c), (b, c), (b, d), (a, d), (a, c)]
    end
to_lines(rect::DT.DiametralBoundingBox) = to_lines(DT.get_bounding_box(rect))
to_lines(rect::SI.SpatialElem) = to_lines(rect.mbr)
to_edge(rect::DT.DiametralBoundingBox) = DT.get_edge(rect)
to_edge(rect::SI.SpatialElem) = rect.val
function get_plot_data(tree, points)
    if tree isa DT.RTree
        bounding_rectangles, id_rectangles = get_dt_rectangles(tree)
    else
        bounding_rectangles, id_rectangles = get_si_rectangles(tree)
    end
    _lines = NTuple{2,Float64}[]
    id_rectangle_lines = NTuple{2,Float64}[]
    bounding_rectangle_lines = NTuple{2,Float64}[]
    for rect in id_rectangles
        i, j = to_edge(rect)
        p, q = get_point(points, i, j)
        push!(_lines, p, q)
        append!(id_rectangle_lines, to_lines(rect))
        push!(id_rectangle_lines, (NaN, NaN))
    end
    for rect in bounding_rectangles
        append!(bounding_rectangle_lines, to_lines(rect))
        push!(bounding_rectangle_lines, (NaN, NaN))
    end
    return _lines, id_rectangle_lines, bounding_rectangle_lines
end
function plot_tree(tree, points)
    _lines, _id_rectangle_lines, _bounding_rectangle_lines = get_plot_data(tree, points)
    fig, ax, sc = lines(_id_rectangle_lines, color=:blue)
    linesegments!(ax, _lines, color=:black)
    lines!(ax, _bounding_rectangle_lines, color=:red)
    display(fig)
    return fig, ax, sc
end

function si_diametral_bounding_box(points, i, j)
    bbox = DT.bounding_box(points, i, j)
    bbox = SI.Rect(DT.get_bl_corner(bbox.bounding_box), DT.get_tr_corner(bbox.bounding_box))
    return bbox
end

function __inorder(tree)
    keys = Any[]
    root = tree.root
    __inorder!(keys, root)
    return identity.(keys)
end
function __inorder!(keys, node)
    left = node.leftChild
    __inorder!(keys, left)
    push!(keys, node.data)
    right = node.rightChild
    __inorder!(keys, right)
    return keys
end
function __inorder!(keys, node::Nothing)
    return keys
end

function get_child_from_tree(tree::DT.PolygonTree, i)
    return first(filter(subtree -> DT.get_index(subtree) == i, DT.get_children(tree)))
end

function traverse_tree(tree::DT.PolygonTree)
    height = DT.get_height(tree)
    index = DT.get_index(tree)
    parent_index = DT.has_parent(tree) ? DT.get_index(DT.get_parent(tree)) : 0
    tup = (height, index, parent_index)
    children = DT.get_children(tree)
    schildren = sort(OrderedSet(children), by=x -> DT.get_index(x))
    for child in schildren
        tup = (tup..., traverse_tree(child))
    end
    return tup
end
function traverse_tree(hierarchy::DT.PolygonHierarchy)
    dict = Dict()
    trees = DT.get_trees(hierarchy)
    strees = sort(OrderedDict(trees), by=x -> x[1])
    for (index, tree) in strees
        dict[index] = traverse_tree(tree)
    end
    return dict
end
function compare_trees(hierarchy1::DT.PolygonHierarchy, hierarchy2::DT.PolygonHierarchy)
    DT.get_bounding_boxes(hierarchy1) ≠ DT.get_bounding_boxes(hierarchy2) && return false
    DT.get_polygon_orientations(hierarchy1) ≠ DT.get_polygon_orientations(hierarchy2) && return false
    trav1 = traverse_tree(hierarchy1)
    trav2 = traverse_tree(hierarchy2)
    trav1 ≠ trav2 && return false
    return true
end

function plot_boundary_and_diametral_circles(points, boundary_nodes, segments=nothing)
    fpoints = flatten_boundary_nodes(points, boundary_nodes, segments)
    fig, ax, sc = lines(fpoints)
    scatter!(ax, points)
    if DT.has_multiple_curves(boundary_nodes)
        for i in 1:DT.num_curves(boundary_nodes)
            curve_nodes = get_boundary_nodes(boundary_nodes, i)
            for j in 1:DT.num_sections(curve_nodes)
                section_nodes = get_boundary_nodes(curve_nodes, j)
                for k in 1:DT.num_boundary_edges(section_nodes)
                    u, v = get_boundary_nodes(section_nodes, k), get_boundary_nodes(section_nodes, k + 1)
                    p, q = get_point(points, u, v)
                    c, r = DT.diametral_circle(p, q)
                    arc!(ax, c, r, 0, 2π, color=:red, linestyle=:dash)
                end
            end
        end
    elseif DT.has_multiple_sections(boundary_nodes)
        for i in 1:DT.num_sections(boundary_nodes)
            section_nodes = get_boundary_nodes(boundary_nodes, i)
            for j in 1:DT.num_boundary_edges(section_nodes)
                u, v = get_boundary_nodes(section_nodes, j), get_boundary_nodes(section_nodes, j + 1)
                p, q = get_point(points, u, v)
                c, r = DT.diametral_circle(p, q)
                arc!(ax, c, r, 0, 2π, color=:red, linestyle=:dash)
            end
        end
    else
        for i in 1:DT.num_boundary_edges(boundary_nodes)
            u, v = get_boundary_nodes(boundary_nodes, i), get_boundary_nodes(boundary_nodes, i + 1)
            p, q = get_point(points, u, v)
            c, r = DT.diametral_circle(p, q)
            arc!(ax, c, r, 0, 2π, color=:red, linestyle=:dash)
        end
    end
    if !isnothing(segments)
        for e in DT.each_edge(segments)
            i, j = DT.edge_vertices(e)
            p, q = get_point(points, i, j)
            c, r = DT.diametral_circle(p, q)
            arc!(ax, c, r, 0, 2π, color=:red, linestyle=:dash)
        end
    end
    display(fig)
    return fig, ax, sc
end
plot_boundary_and_diametral_circles(enricher::DT.BoundaryEnricher) = plot_boundary_and_diametral_circles(
    get_points(enricher),
    get_boundary_nodes(enricher),
    DT.get_segments(enricher)
)

function maximum_total_variation(points, boundary_nodes, boundary_curves)
    isempty(boundary_curves) && return 0.0
    max_Δθ = -Inf
    if DT.has_multiple_curves(boundary_nodes)
        ctr = 1
        for i in 1:DT.num_curves(boundary_nodes)
            curve_nodes = get_boundary_nodes(boundary_nodes, i)
            for j in 1:DT.num_sections(curve_nodes)
                section_nodes = get_boundary_nodes(curve_nodes, j)
                curve = boundary_curves[ctr]
                DT.is_piecewise_linear(curve) && continue
                for k in 1:DT.num_boundary_edges(section_nodes)
                    u, v = get_boundary_nodes(section_nodes, k), get_boundary_nodes(section_nodes, k + 1)
                    p, q = get_point(points, u, v)
                    t₁, t₂ = DT.get_inverse(curve, p), DT.get_inverse(curve, q)
                    Δθ = slow_total_absolute_curvature(curve, t₁, t₂)
                    if Δθ > max_Δθ
                        max_Δθ = Δθ
                    end
                end
                ctr += 1
            end
        end
    elseif DT.has_multiple_sections(boundary_nodes)
        for i in 1:DT.num_sections(boundary_nodes)
            section_nodes = get_boundary_nodes(boundary_nodes, i)
            curve = boundary_curves[i]
            DT.is_piecewise_linear(curve) && continue
            for j in 1:DT.num_boundary_edges(section_nodes)
                u, v = get_boundary_nodes(section_nodes, j), get_boundary_nodes(section_nodes, j + 1)
                p, q = get_point(points, u, v)
                t₁, t₂ = DT.get_inverse(curve, p), DT.get_inverse(curve, q)
                Δθ = slow_total_absolute_curvature(curve, t₁, t₂)
                if Δθ > max_Δθ
                    max_Δθ = Δθ
                end
            end
        end
    else
        curve = boundary_curves[1]
        DT.is_piecewise_linear(curve) && return 0.0
        for i in 1:DT.num_boundary_edges(boundary_nodes)
            u, v = get_boundary_nodes(boundary_nodes, i), get_boundary_nodes(boundary_nodes, i + 1)
            p, q = get_point(points, u, v)
            t₁, t₂ = DT.get_inverse(curve, p), DT.get_inverse(curve, q)
            if t₂ == 0.0
                t₂ = 1.0
            end
            Δθ = slow_total_absolute_curvature(curve, t₁, t₂)
            if Δθ > max_Δθ
                max_Δθ = Δθ
            end
        end
    end
    return max_Δθ
end

function ⪧(x::DT.CatmullRomSplineSegment, y::DT.CatmullRomSplineSegment; kwargs...)
    ⪧(x.a, y.a; kwargs...) || return false
    ⪧(x.b, y.b; kwargs...) || return false
    ⪧(x.c, y.c; kwargs...) || return false
    ⪧(x.d, y.d; kwargs...) || return false
    ⪧(x.p₁, y.p₁; kwargs...) || return false
    ⪧(x.p₂, y.p₂; kwargs...) || return false
    return true
end

function all_diametral_circles_are_empty(enricher::DT.BoundaryEnricher)
    points = get_points(enricher)
    tree = DT.BoundaryRTree(points)
    n = 0
    for e in DT.each_boundary_edge(enricher)
        n += 1
        i, j = DT.edge_vertices(e)
        insert!(tree, i, j)
    end
    if DT.has_segments(enricher)
        for e in DT.get_segments(enricher)
            n += 1
            i, j = DT.edge_vertices(e)
            insert!(tree, i, j)
        end
    end
    rat = n
    for i in DT.each_point_index(points)
        intersects = DT.get_intersections(tree, i)
        for box in intersects
            u, v = DT.get_edge(box)
            u, v = DT.reorient_edge(enricher, u, v)
            (u == i || v == i) && continue
            cert = DT.point_position_relative_to_diametral_circle(get_point(points, u, v, i)...)
            flag = DT.is_inside(cert) && DT.is_visible(DT.test_visibility(enricher, u, v, i))
            if flag
                p, q = get_point(points, u, v)
                t, Δθ, ct = DT.compute_split_position(enricher, u, v)
                flag = flag && !isnan(Δθ)
            end
            rat -= flag
        end
    end
    return rat / n # can't just look at rat/n == 1 since, due to precision issues, we might have some points just slightly in the circle
end

function all_points_are_inside(enricher::DT.BoundaryEnricher, orig_points, orig_boundary_nodes)
    orig_points, orig_boundary_nodes = deepcopy(orig_points), deepcopy(orig_boundary_nodes)
    boundary_curves, new_boundary_nodes = DT.convert_boundary_curves!(orig_points, orig_boundary_nodes, Int)
    new_points, new_boundary_nodes = DT.polygonise(orig_points, new_boundary_nodes, boundary_curves)
    points = get_points(enricher)
    hierarchy = DT.get_polygon_hierarchy(enricher)
    for p in DT.each_point(points)
        δ = DT.distance_to_polygon(p, new_points, new_boundary_nodes)
        δ < -1e-4 && return false
    end
    return true
end

function plot_small_angle_complexes(enricher)
    points, boundary_nodes, boundary_curves = get_points(enricher), get_boundary_nodes(enricher), DT.get_boundary_curves(enricher)
    complexes = DT.get_small_angle_complexes(points, boundary_nodes, boundary_curves)
    _points, _boundary_nodes = DT.polygonise(points, boundary_nodes, boundary_curves; n=2^10)
    fpoints = flatten_boundary_nodes(_points, _boundary_nodes)
    fig, ax, sc = lines(fpoints)
    colors = [:red, :green, :blue, :orange, :purple, :cyan, :magenta, :yellow]
    for (apex, _complexes) in complexes
        color_ctr = 1
        for complex in _complexes
            for member in DT.get_members(complex)
                next_edge = DT.get_next_edge(member)
                parent_curve = DT.get_parent_curve(member)
                curve = boundary_curves[parent_curve]
                p, q = get_point(points, apex, next_edge)
                if DT.is_piecewise_linear(curve)
                    lines!(ax, [p, q], linewidth=4, color=colors[color_ctr])
                else
                    t₁ = DT.get_inverse(curve, p)
                    t₂ = DT.get_inverse(curve, q)
                    t = LinRange(t₁, t₂, 1000)
                    lines!(ax, curve.(t), linewidth=4, color=colors[color_ctr])
                end
            end
            color_ctr = mod1(color_ctr + 1, length(colors))
        end
        p = get_point(points, apex)
        scatter!(ax, [p], color=colors[1], markersize=17)
    end
    display(fig)
    return fig
end

using DelaunayTriangulation:
    test_adjacent2vertex_map_matches_adjacent_map,
    test_state,
    test_adjacent_map_matches_adjacent2vertex_map,
    test_each_edge_has_two_incident_triangles,
    test_triangle_orientation,
    test_iterators

export validate_triangulation
export @_adj
export _make_graph_from_adjacency
export get_random_convex_polygon
export compare_trees
export second_shewchuk_example_constrained
export get_random_vertices_and_constrained_edges
export example_for_testing_add_point_on_constrained_triangulation
export sort_edge_vector
export complicated_geometry
export validate_refinement
export validate_statistics
export validate_tessellation
export compare_edge_vectors
export simple_geometry
export fixed_shewchuk_example_constrained
export ⊢
export is_sink
export is_conformal
export _validate_offcenter
export example_with_special_corners
export get_idx
export is_interior_sink
export _compare_pairs
export ⪧
export slow_total_absolute_curvature
export slow_arc_length
export slow_bezier_eval
export slow_eval_bspline
export slow_get_segment
export closest_point_on_curve
export get_dt_rectangles
export si_diametral_bounding_box
export __inorder
export get_child_from_tree
export traverse_tree
export poor_triangulation_example
export example_triangulation
export example_empty_triangulation
export shewchuk_example_constrained
export test_segment_triangle_intersections
export test_split_edges
export test_intersections
export validate_insertion_event_history
export _approx_ispow2
export _slow_compute_concentric_shell_ternary_split_position
export slow_triangle_assess_queue
export slow_triangle_assess
export compare_encroach_queues
export slow_encroachment_test
export compute_diametral_circle
export get_points_in_diametral_circle
export slow_encroachment_test_diametral_lens
export compare_triangle_queues
export flatten_boundary_nodes
export maximum_total_variation
export all_points_are_inside
export all_diametral_circles_are_empty
export compute_diametral_lens
export get_points_in_diametral_lens
export test_adjacent2vertex_map_matches_adjacent_map
export test_state
export test_adjacent_map_matches_adjacent2vertex_map
export test_each_edge_has_two_incident_triangles
export test_triangle_orientation
export test_iterators
end