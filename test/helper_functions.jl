using StaticArrays
using StableRNGs

const DT = DelaunayTriangulation

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
    T = Set{NTuple{3,Int64}}((Tuple(T) for T in eachrow(T)))
    outer = [[indexin([a, b, c, d, e, f, g, h, a], pts)...]]
    inner1 = [[indexin([ℓ, k, j, i, ℓ], pts)...]]
    inner2 = [[indexin([r, q, p], pts)...], [indexin([p, o, n, m, r], pts)...]]
    boundary_nodes = [outer, inner1, inner2]
    label_map = Dict(["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "ℓ",
        "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "z", "a1",
        "b1"] .=> pts)
    index_map = Dict(["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "ℓ",
        "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "z", "a1",
        "b1"] .=> each_point_index(pts))
    return Triangulation(pts, T, boundary_nodes), label_map, index_map
end

function validate_triangulation(tri::Triangulation;
    ignore_boundary_indices=false)
    # Tests aren't countered correctly when we use @async, 
    # so let's build thread-safe counters and then fake a 
    # testset later. 
    local passes = Threads.Atomic{Int64}(0)
    local fails = Threads.Atomic{Int64}(0)
    _tri = deepcopy(tri)
    DT.clear_empty_features!(_tri)

    ## Check that all triangles are positively oriented 
    @sync begin
        @async for T in each_triangle(_tri)
            cert = DT.triangle_orientation(_tri, T)
            DT.is_positively_oriented(cert) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
            @test DT.is_positively_oriented(cert)
        end

        ## Test that the triangulation is Delaunay 
        @async for T in each_triangle(_tri)
            for r in each_point_index(_tri)
                if r ∈ get_vertices(tri)
                    if !(ignore_boundary_indices && DT.is_ghost_triangle(T))
                        cert = DT.point_position_relative_to_circumcircle(_tri, T, r)
                        DT.is_inside(cert) && @show T, r
                        !DT.is_inside(cert) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                        @test !DT.is_inside(cert)
                    end
                end
            end
        end

        ## Test that every edge is incident to two triangles if it not a boundary edge, and one otherwise 
        @async for (ij, k) in get_adjacent(_tri)
            i, j = initial(ij), terminal(ij)
            vij = get_adjacent(tri, i, j)
            vji = get_adjacent(tri, j, i)
            if DT.is_boundary_edge(_tri, i, j)
                !DT.is_boundary_index(vij) && @show i, j, vij
                DT.is_boundary_index(vij) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test DT.is_boundary_index(vij)
            elseif DT.is_boundary_edge(_tri, j, i)
                !DT.is_boundary_index(vji) && @show j, i, vji
                DT.is_boundary_index(vji) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test DT.is_boundary_index(vji)
            else
                !(DT.edge_exists(vij) && DT.edge_exists(vji)) && @show i, j, vij, vji
                DT.edge_exists(vij) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test DT.edge_exists(vij) && DT.edge_exists(vji)
            end
        end

        ## Test the adjacent map 
        @async for T in each_triangle(_tri)
            u, v, w = indices(T)
            if !(ignore_boundary_indices && DT.is_ghost_triangle(T))
                all((get_adjacent(_tri, u, v) == w, get_adjacent(_tri, u, v) == w, get_adjacent(_tri, w, u) == v)) || @show u, v, w, T
                get_adjacent(_tri, u, v) == w ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test get_adjacent(_tri, u, v) == w
                get_adjacent(_tri, v, w) == u ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test get_adjacent(_tri, v, w) == u
                get_adjacent(_tri, w, u) == v ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test get_adjacent(_tri, w, u) == v
            end
        end

        ## Test the adjacent2vertex map 
        E = DT.edge_type(_tri)
        @async for T in each_triangle(_tri)
            u, v, w = indices(T)
            if !(ignore_boundary_indices && DT.is_ghost_triangle(T))
                DT.contains_edge(DT.construct_edge(E, v, w), get_adjacent2vertex(_tri, u)) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test DT.contains_edge(DT.construct_edge(E, v, w), get_adjacent2vertex(_tri, u))
                DT.contains_edge(DT.construct_edge(E, w, u), get_adjacent2vertex(_tri, v)) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test DT.contains_edge(DT.construct_edge(E, w, u), get_adjacent2vertex(_tri, v))
                DT.contains_edge(DT.construct_edge(E, u, v), get_adjacent2vertex(_tri, w)) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test DT.contains_edge(DT.construct_edge(E, u, v), get_adjacent2vertex(_tri, w))
            end
        end

        ## Check that the adjacent and adjacent2vertex maps are inverses 
        @async for (k, S) in get_adjacent2vertex(_tri)
            for ij in DT.each_edge(S)
                if !(ignore_boundary_indices && (DT.is_boundary_index(k) || DT.is_ghost_edge(ij)))
                    get_adjacent(_tri, ij) == k || @show i, j, k
                    get_adjacent(_tri, ij) == k ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                    @test get_adjacent(_tri, ij) == k
                end
            end
        end
        Base.Threads.@spawn for (ij, k) in get_adjacent(_tri)
            if !(ignore_boundary_indices && (DT.is_boundary_index(k) || DT.is_ghost_edge(ij)))
                DT.contains_edge(ij, get_adjacent2vertex(_tri, k)) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test DT.contains_edge(ij, get_adjacent2vertex(_tri, k))
            end
        end

        ## Check the graph 
        @async for (i, j) in get_edges(_tri)
            if DT.has_ghost_triangles(_tri)
                ij = DT.construct_edge(E, i, j)
                ji = DT.construct_edge(E, j, i)
                ij ∈ keys((get_adjacent ∘ get_adjacent)(_tri)) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test ij ∈ keys((get_adjacent ∘ get_adjacent)(_tri))
                ji ∈ keys((get_adjacent ∘ get_adjacent)(_tri)) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                @test ji ∈ keys((get_adjacent ∘ get_adjacent)(_tri))
            else
                if !DT.is_boundary_index(i) && !DT.is_boundary_index(j)
                    ij = DT.construct_edge(E, i, j)
                    ji = DT.construct_edge(E, j, i)
                    ij ∈ keys((get_adjacent ∘ get_adjacent)(_tri)) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                    @test ij ∈ keys((get_adjacent ∘ get_adjacent)(_tri))
                    ji ∈ keys((get_adjacent ∘ get_adjacent)(_tri)) ? Threads.atomic_add!(passes, 1) : Threads.atomic_add!(fails, 1)
                    @test ji ∈ keys((get_adjacent ∘ get_adjacent)(_tri))
                end
            end
        end
    end
    @testset "Triangulation validation" begin
        for _ in 1:passes[]
            @test true
        end
        for _ in 1:fails[]
            @test false
        end
    end
end

macro _adj(i, j, k)
    return :(($i, $j) => $k, ($j, $k) => $i, ($k, $i) => $j)
end

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
    T = Set{NTuple{3,Int64}}([
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
    DG = DT.Graph(DT.SimpleGraphs.relabel(DT.SimpleGraphs.UndirectedGraph(A), Dict(1:7 .=> [-1, (1:6)...])))
    adj = DT.Adjacent(DT.DataStructures.DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (1, 3) => 5, (3, 5) => 1,
            (4, 5) => DT.BoundaryIndex, (5, 2) => DT.BoundaryIndex,
            (2, 3) => DT.BoundaryIndex, (3, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex
        )
    ))
    adj2v = DT.Adjacent2Vertex(Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(5, 4), (3, 5), (6, 3), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3)]),
        3 => Set{NTuple{2,Int64}}([(1, 6), (5, 1), (2, 5)]),
        4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(1, 4), (3, 1)])
    ))
    tri = Triangulation(pts, T, adj, adj2v, DG,
        Int64[],
        DT.DataStructures.OrderedDict{Int64,Vector{Int64}}(),
        DT.DataStructures.OrderedDict{Int64,UnitRange{Int64}}(),
        Set{NTuple{2,Int64}}(),
        Set{NTuple{2,Int64}}(),
        ConvexHull(pts, [2, 6, 4, 5, 2]))
    return tri
end

function example_empty_triangulation()
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    pts = [p1, p2, p3]
    T = Set{NTuple{3,Int64}}([])
    A = zeros(Int64, 0, 0)
    DG = DT.Graph(DT.SimpleGraphs.UndirectedGraph(A))
    adj = DT.Adjacent(DT.DataStructures.DefaultDict(DT.DefaultAdjacentValue, Dict{NTuple{2,Int64},Int64}()))
    adj2v = DT.Adjacent2Vertex(Dict(DT.BoundaryIndex => Set{NTuple{2,Int64}}()))
    tri = Triangulation(pts, T, adj, adj2v, DG,
        Int64[],
        DT.DataStructures.OrderedDict{Int64,Vector{Int64}}(),
        DT.DataStructures.OrderedDict{Int64,UnitRange{Int64}}(),
        Set{NTuple{2,Int64}}(),
        Set{NTuple{2,Int64}}(),
        ConvexHull(pts, [2, 6, 4, 5, 2]))
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
        intersecting_triangles1, collinear_segments1, left1, right1 = DT.locate_intersecting_triangles(tri, e)
        intersecting_triangles2, collinear_segments2, left2, right2 = DT.locate_intersecting_triangles(
            e,
            get_points(tri),
            get_adjacent(tri),
            get_adjacent2vertex(tri),
            get_graph(tri),
            get_boundary_index_ranges(tri),
            get_boundary_map(tri),
            NTuple{3,Int64},
        )
        for (intersecting_triangles, collinear_segments) in zip(
            (intersecting_triangles1, intersecting_triangles2),
            (collinear_segments1, collinear_segments2)
        )
            @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), intersecting_triangles)
            @test all(!DT.is_none, [DT.triangle_line_segment_intersection(tri, T, e) for T in intersecting_triangles])
            @test DT.compare_triangle_collections(allT, intersecting_triangles)
            @test allunique(intersecting_triangles)
            if typeof(constrained_edges) <: AbstractVector
                @test collinear_segments == constrained_edges
            else # Tuple of possibilities, in case the edge's endpoints have equal degree so that we could start at any point
                @test any(==(collinear_segments), constrained_edges)
            end
        end
    end
end

function test_split_edges(tri, edge, current_constrained_edges)
    constrained_edges = get_constrained_edges(tri)
    DT.add_edge!(constrained_edges, edge)
    _, collinear_segments = DT.locate_intersecting_triangles(tri, edge)
    DT.split_constrained_edge!(tri, edge, collinear_segments)
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
    constrained_edges = get_constrained_edges(tri)
    for edge in ((edge[1], edge[2]), (edge[2], edge[1]))
        intersecting_triangles1, collinear_segments1, left1, right1 = DT.locate_intersecting_triangles(tri, edge)
        intersecting_triangles2, collinear_segments2, left2, right2 = DT.locate_intersecting_triangles(
            edge,
            get_points(tri),
            get_adjacent(tri),
            get_adjacent2vertex(tri),
            get_graph(tri),
            get_boundary_index_ranges(tri),
            get_boundary_map(tri),
            NTuple{3,Int64},
        )
        for (intersecting_triangles, collinear_segments) in zip(
            (intersecting_triangles1, intersecting_triangles2),
            (collinear_segments1, collinear_segments2)
        )
            @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), intersecting_triangles)
            @test all(!DT.is_none, [DT.triangle_line_segment_intersection(tri, T, edge) for T in intersecting_triangles])
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
    end
    constrained_edges = get_constrained_edges(tri)
    DT.add_edge!(constrained_edges, edge)
    _, collinear_segments = DT.locate_intersecting_triangles(tri, edge)
    DT.split_constrained_edge!(tri, edge, collinear_segments)
    @test compare_edge_vectors(constrained_edges, current_constrained_edges)
end
