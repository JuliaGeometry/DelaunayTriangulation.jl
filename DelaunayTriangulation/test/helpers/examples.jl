# Example geometries and triangulations used throughout the tests.

export @_adj, _make_graph_from_adjacency, get_random_convex_polygon, second_shewchuk_example_constrained
export get_random_vertices_and_constrained_edges, example_for_testing_add_point_on_constrained_triangulation
export complicated_geometry, simple_geometry, fixed_shewchuk_example_constrained, example_with_special_corners
export poor_triangulation_example, example_triangulation, example_empty_triangulation, shewchuk_example_constrained
export julia_logo_curves, tight_walls_boundary

function complicated_geometry()
    x1 = [
        collect(LinRange(0, 2, 4)),
        collect(LinRange(2, 2, 4)),
        collect(LinRange(2, 0, 4)),
        collect(LinRange(0, 0, 4)),
    ]
    y1 = [
        collect(LinRange(0, 0, 4)),
        collect(LinRange(0, 6, 4)),
        collect(LinRange(6, 6, 4)),
        collect(LinRange(6, 0, 4)),
    ]
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
    x4 = reverse(
        reverse.(
            [
            collect(LinRange(1, 1.5, 4)),
            collect(LinRange(1.5, 1.5, 4)),
            collect(LinRange(1.5, 1, 4)),
            collect(LinRange(1, 1, 4)),
        ],
        ),
    )
    y4 = reverse(
        reverse.(
            [
            collect(LinRange(2, 2, 4)),
            collect(LinRange(2, 5, 4)),
            collect(LinRange(5, 5, 4)),
            collect(LinRange(5, 2, 4)),
        ],
        ),
    )
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
    pts = [
        a, b, c, d, e, f, g, h, i, j, k, ℓ,
        m, n, o, p, q, r, s, t, u, v, w, z, a1, b1,
    ]
    T = [
        h a s
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
        r q d
    ]
    T = indexin(T, pts)
    T = Set{NTuple{3,Int}}((Tuple(T) for T in eachrow(T)))
    outer = [[indexin([a, b, c, d, e, f, g, h, a], pts)...]]
    inner1 = [[indexin([ℓ, k, j, i, ℓ], pts)...]]
    inner2 = [[indexin([r, q, p], pts)...], [indexin([p, o, n, m, r], pts)...]]
    boundary_nodes = [outer, inner1, inner2]
    label_map = Dict(
        [
            "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "ℓ",
            "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "z", "a1",
            "b1",
        ] .=> pts,
    )
    index_map = Dict(
        [
            "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "ℓ",
            "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "z", "a1",
            "b1",
        ] .=> DT.each_point_index(pts),
    )
    return DT.Triangulation(pts, T, boundary_nodes, delete_ghosts=true, predicates=DT.ExactKernel()), label_map, index_map
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
    T = Set{NTuple{3,Int}}(
        [
        (6, 3, 1),
        (3, 2, 5),
        (4, 1, 5),
        (4, 6, 1),
        (5, 1, 3),
    ],
    )
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
            (6, 4) => DT.𝒢,
        ),
    )
    adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.𝒢 => Set{NTuple{2,Int}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
            1 => Set{NTuple{2,Int}}([(5, 4), (3, 5), (6, 3), (4, 6)]),
            2 => Set{NTuple{2,Int}}([(5, 3)]),
            3 => Set{NTuple{2,Int}}([(1, 6), (5, 1), (2, 5)]),
            4 => Set{NTuple{2,Int}}([(1, 5), (6, 1)]),
            5 => Set{NTuple{2,Int}}([(4, 1), (1, 3), (3, 2)]),
            6 => Set{NTuple{2,Int}}([(1, 4), (3, 1)]),
        ),
    )
    rep = Dict{Int,DT.RepresentativeCoordinates{Int,Float64}}()
    tri = DT.Triangulation(
        pts, T, Int[], Set{NTuple{2,Int}}(), Set{NTuple{2,Int}}(),
        DT.ZeroWeight(), adj, adj2v, DG, (), DT.construct_boundary_edge_map(Int[]),
        Dict{Int,Vector{Int}}(), Dict{Int,UnitRange{Int}}(), DT.ConvexHull(pts, [2, 6, 4, 5, 2]), rep, DT.construct_polygon_hierarchy(pts), nothing,
        DT._build_cache(
            pts, Int, NTuple{2,Int}, NTuple{3,Int},
            Set{NTuple{2,Int}}, Set{NTuple{3,Int}}, DT.ZeroWeight(), Val(true),
        ),
        false,  # has_ghosts
        Dict{Int,Int}(),  # boundary_vertex_to_ghost
    )
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
    tri = DT.Triangulation(
        pts, T, Int[], Set{NTuple{2,Int}}(), Set{NTuple{2,Int}}(),
        DT.ZeroWeight(), adj, adj2v, DG, (), DT.construct_boundary_edge_map(Int[]),
        Dict{Int,Vector{Int}}(), Dict{Int,UnitRange{Int}}(), DT.ConvexHull(pts, [1, 2, 3, 1]), rep, DT.construct_polygon_hierarchy(pts), nothing,
        DT._build_cache(
            pts, Int, NTuple{2,Int}, NTuple{3,Int},
            Set{NTuple{2,Int}}, Set{NTuple{3,Int}}, DT.ZeroWeight(), Val(true),
        ),
        false,  # has_ghosts
        Dict{Int,Int}(),  # boundary_vertex_to_ghost
    )
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
    tri = triangulate(pts; rng, delete_ghosts=false, randomise=false, predicates=DT.ExactKernel())
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
    tri = triangulate(pts; rng, delete_ghosts=false, randomise=false, predicates=DT.ExactKernel())
    return tri
end

function get_random_vertices_and_constrained_edges(nverts1, nverts2, nedges, rng=Random.default_rng())
    ## To generate a random set of constrained edges, we get a random small triangulation, 
    ## and we just take the edges from that triangulation.
    points = [Tuple(rand(rng, 2)) for _ in 1:nverts1]
    tri = triangulate(points; predicates=DT.FastKernel(), rng)
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
        (27, 28),
    ]
    pts = [
        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
        p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24,
        p25, p26, p27, p28,
    ]
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

function get_random_convex_polygon(points)
    tri = triangulate(points; predicates=DT.ExactKernel(), rng=StableRNG(1))
    S = get_convex_hull_vertices(tri)
    pop!(S) # Want S[begin] ≠ S[end]
    return S
end

# A domain with tight walls and small angles, with a hole, as curves to pass to convert_boundary_points_to_indices.
function tight_walls_boundary()
    A = (0.0, 0.0)
    B = (0.0, 25.0)
    C = (5.0, 25.0)
    D = (5.0, 5.0)
    E = (10.0, 5.0)
    F = (10.0, 10.0)
    G = (25.0, 10.0)
    H = (25.0, 15.0)
    I = (10.0, 15.0)
    J = (10.0, 25.0)
    K = (45.0, 25.0)
    L = (45.0, 20.0)
    M = (40.0, 20.0)
    N = (40.0, 5.0)
    O = (45.0, 5.0)
    P = (45.0, 0.0)
    Q = (10.0, 0.0)
    R = (10.0, -5.0)
    S = (15.0, -5.0)
    T = (15.0, -10.0)
    U = (10.0, -10.0)
    V = (5.0, -10.0)
    W = (5.0, -5.0)
    Z = (5.0, 0.0)
    A1 = (5.0, 2.5)
    B1 = (10.0, 2.5)
    C1 = (38.0, 2.5)
    D1 = (38.0, 20.0)
    E1 = (27.0, 20.0)
    F1 = (27.0, 11.0)
    G1 = (27.0, 4.0)
    H1 = (2.0, 4.0)
    I1 = (2.0, 0.0)
    pts = [A, I1, H1, G1, F1, E1, D1, C1, B1, A1, Z, W, V, U, T, S, R, Q, P, O, N, M, L, K, J, I, H, G, F, E, D, C, B, A]
    J1 = (17.0603265896789, 7.623652007194)
    K1 = (14.8552854162067, 6.5423337394336)
    L1 = (16.6998871670921, 6.9875824379232)
    M1 = (16.0, 6.0)
    N1 = (16.9755173137761, 6.6483453343121)
    O1 = (17.0391242707032, 4.8885528593294)
    P1 = (17.4207660122657, 6.4575244635308)
    Q1 = (17.6327892020226, 4.9945644542079)
    R1 = (22.6789411182379, 6.1818943168468)
    S1 = (21.8096460402344, 6.4787267825065)
    T1 = (26.0, 8.0)
    U1 = (15.0673086059636, 9.086612016517)
    W1 = (15.0, 8.5)
    Z1 = (17.7913089332764, 8.3603005983396)
    inner_pts = [Z1, W1, U1, T1, S1, R1, Q1, P1, O1, N1, M1, L1, K1, J1, Z1]
    boundary_pts = [[pts], [inner_pts]]
    return boundary_pts
end

# The boundary of the Julia logo, as curves to pass to convert_boundary_points_to_indices.
function julia_logo_curves()
    C = (15.7109521325776, 33.244486807457)
    D = (14.2705719699703, 32.8530791545746)
    E = (14.3, 27.2)
    F = (14.1, 27.0)
    G = (13.7, 27.2)
    H = (13.4, 27.5)
    I = (13.1, 27.6)
    J = (12.7, 27.4)
    K = (12.5, 27.1)
    L = (12.7, 26.7)
    M = (13.1, 26.5)
    N = (13.6, 26.4)
    O = (14.0, 26.4)
    P = (14.6, 26.5)
    Q = (15.1983491346581, 26.8128534095401)
    R = (15.6, 27.6)
    S = (15.6952958264624, 28.2344688505621)
    T = (17.8088971520274, 33.1192363585346)
    U = (16.3058917649589, 33.0722674401887)
    V = (16.3215480710742, 29.7374742376305)
    W = (16.3841732955354, 29.393035503094)
    Z = (16.6190178872649, 28.9233463196351)
    A1 = (17.0417381523779, 28.5319386667527)
    B1 = (17.5114273358368, 28.3753756055997)
    C1 = (18.1376795804487, 28.3597192994844)
    D1 = (18.7169629067146, 28.5632512789833)
    E1 = (19.2805899268653, 28.8920337074045)
    F1 = (19.26493362075, 28.4536571361762)
    G1 = (20.6426885588962, 28.4223445239456)
    H1 = (20.689657477242, 33.1035800524193)
    I1 = (19.2805899268653, 33.0722674401887)
    J1 = (19.2962462329806, 29.7531305437458)
    K1 = (19.0614016412512, 29.393035503094)
    L1 = (18.7482755189452, 29.236472441941)
    M1 = (18.4508057027546, 29.1425346052493)
    N1 = (18.1689921926793, 29.3147539725175)
    O1 = (17.7932408459121, 29.6278800948235)
    P1 = (22.6466957416542, 35.4207133574833)
    Q1 = (21.2219718851621, 34.9979930923702)
    R1 = (21.2376281912774, 28.4693134422915)
    S1 = (22.6780083538847, 28.4380008300609)
    T1 = (24.5724213938357, 33.1975178891111)
    U1 = (23.3512295168425, 32.8530791545746)
    V1 = (23.3199169046119, 28.4380008300609)
    W1 = (24.6663592305274, 28.3753756055997)
    Z1 = (15.1942940307729, 35.4363696635986)
    A2 = (14.7246048473139, 35.3737444391374)
    B2 = (14.3645098066621, 35.1858687657538)
    C2 = (14.1766341332786, 34.8570863373326)
    D2 = (14.1140089088174, 34.3247719294125)
    E2 = (14.2705719699703, 33.8394264398383)
    F2 = (14.7246048473139, 33.6202381542241)
    G2 = (15.4604512347329, 33.6045818481088)
    H2 = (16.0, 34.0)
    I2 = (15.9771093365377, 34.6848669700643)
    J2 = (15.6170142958859, 35.2328376840997)
    K2 = (24.1653574348379, 35.4520259697138)
    L2 = (23.7739497819555, 35.4363696635986)
    M2 = (23.4608236596496, 35.2641502963303)
    N2 = (23.272947986266, 34.9040552556785)
    O2 = (23.1320412312284, 34.5909291333725)
    P2 = (23.1163849251131, 34.2151777866054)
    Q2 = (23.2886042923813, 33.8081138276077)
    R2 = (23.8209187003014, 33.6045818481088)
    S2 = (24.3062641898756, 33.5576129297629)
    T2 = (24.7602970672192, 33.8550827459536)
    U2 = (25.010797965064, 34.4656786844502)
    V2 = (24.8385785977957, 34.9666804801397)
    W2 = (24.5254524754898, 35.2641502963303)
    Z2 = (25.3708930057158, 37.4716894585871)
    A3 = (24.7916096794498, 37.3464390096648)
    B3 = (24.4471709449133, 36.9550313567823)
    C3 = (24.3062641898756, 36.5636237038999)
    D3 = (24.4941398632592, 35.9999966837492)
    E3 = (25.0264542711793, 35.5929327247515)
    F3 = (25.5587686790994, 35.5929327247515)
    F3 = (25.5587686790994, 35.5929327247515)
    G3 = (26.0, 36.0)
    H3 = (26.1380520053653, 36.5792800100152)
    I3 = (26.0, 37.0)
    J3 = (25.7466443524829, 37.2838137852036)
    K3 = (26.3885529032101, 35.4676822758291)
    L3 = (25.9814889442124, 35.3580881330221)
    M3 = (25.6840191280217, 35.1858687657538)
    N3 = (25.5274560668688, 34.9040552556785)
    O3 = (25.4961434546382, 34.5596165211419)
    P3 = (25.5274560668688, 34.246490398836)
    Q3 = (25.6683628219064, 33.8394264398383)
    R3 = (26.0284578625583, 33.6358944603394)
    S3 = (26.5451159643631, 33.6202381542241)
    T3 = (27.0, 34.0)
    U3 = (27.280962351782, 34.5596165211419)
    V3 = (27.0304614539373, 35.2171813779844)
    W3 = (26.1693646175959, 33.087923746304)
    Z3 = (26.0, 33.0)
    A4 = (25.5274560668688, 32.7278287056522)
    B4 = (25.2612988629087, 32.4147025833463)
    C4 = (25.1830173323322, 32.0702638488098)
    D4 = (25.2299862506781, 31.7727940326191)
    E4 = (25.6527065157911, 31.5222931347744)
    F4 = (26.2946150665183, 31.7258251142732)
    G4 = (26.5607722704784, 32.5086404200381)
    H4 = (27.1557119028596, 32.7434850117675)
    I4 = (27.6097447802033, 32.4929841139228)
    J4 = (27.6410573924338, 32.1015764610403)
    K4 = (27.7193389230103, 31.6005746653509)
    L4 = (27.437525412935, 31.4283552980826)
    M4 = (26.9834925355914, 31.2561359308143)
    N4 = (26.5764285765937, 31.0995728696614)
    O4 = (26.0441141686736, 30.7864467473554)
    P4 = (25.6527065157911, 30.5672584617413)
    Q4 = (25.3239240873699, 30.1915071149741)
    R4 = (25.1673610262169, 29.8783809926682)
    S4 = (25.1047358017558, 29.6122237887082)
    T4 = (25.0890794956405, 29.1895035235952)
    U4 = (25.2926114751393, 28.8294084829433)
    V4 = (25.6840191280217, 28.5632512789833)
    W4 = (26.1537083114806, 28.3753756055997)
    Z4 = (26.8269294744384, 28.391031911715)
    A5 = (27.4844943312809, 28.6102201973292)
    B5 = (27.7342002330051, 28.7239579596219)
    C5 = (27.7264126450755, 28.4202565942047)
    D5 = (29.1825559185446, 28.3922538389457)
    E5 = (29.1545531632856, 32.2146299318021)
    F5 = (29.000538009361, 32.5786657501693)
    G5 = (28.6785063238822, 32.9006974356481)
    H5 = (28.3144705055149, 33.0827153448317)
    I5 = (27.9084305542591, 33.2367304987563)
    J5 = (27.3343740714492, 33.3207387645334)
    K5 = (26.8303244767868, 33.2367304987563)
    L5 = (27.6564057569279, 30.786489413592)
    M5 = (27.6984098898165, 30.3944508399657)
    N5 = (27.6984098898165, 29.7363860913787)
    O5 = (27.5863988687804, 29.4143544059)
    P5 = (27.2643671833016, 29.2043337414573)
    Q5 = (26.9843396307114, 29.1763309861983)
    R5 = (26.6903107004917, 29.3163447624934)
    S5 = (26.5782996794556, 29.7503874690082)
    T5 = (26.7603175886393, 30.3384453294476)
    U5 = (27.3203726938197, 30.7024811478149)
    J_curve = [[C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, C]]
    U_curve = [[T, U, V, W, Z, A1, B1, C1, D1, E1, F1, G1, H1, I1, J1, K1, L1, M1, N1, O1, T]]
    L_curve = [[P1, Q1, R1, S1, P1]]
    I_curve = [[T1, U1, V1, W1, T1]]
    A_curve_outline = [
        [
            K5, W3, Z3, A4, B4, C4, D4, E4, F4, G4, H4, I4, J4, K4, L4, M4, N4,
            O4, P4, Q4, R4, S4, T4, U4, V4, W4, Z4, A5, B5, C5, D5, E5, F5, G5,
            H5, I5, J5, K5,
        ],
    ]
    A_curve_hole = [[L5, M5, N5, O5, P5, Q5, R5, S5, T5, U5, L5]]
    dot_1 = [[Z1, A2, B2, C2, D2, E2, F2, G2, H2, I2, J2, Z1]]
    dot_2 = [[Z2, A3, B3, C3, D3, E3, F3, G3, H3, I3, J3, Z2]]
    dot_3 = [[K2, L2, M2, N2, O2, P2, Q2, R2, S2, T2, U2, V2, W2, K2]]
    dot_4 = [[K3, L3, M3, N3, O3, P3, Q3, R3, S3, T3, U3, V3, K3]]
    curves = [J_curve, U_curve, L_curve, I_curve, A_curve_outline, A_curve_hole, dot_1, dot_2, dot_3, dot_4]
    return curves
end
