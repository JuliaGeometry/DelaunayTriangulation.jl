using StaticArrays

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

function validate_triangulation(tri::Triangulation)
    _tri = deepcopy(tri)
    DT.clear_empty_features!(_tri)

    ## Check that all triangles are positively oriented 
    @sync begin
        @async for T in each_triangle(_tri)
            cert = DT.triangle_orientation(_tri, T)
            @test DT.is_positively_oriented(cert)
        end

        ## Test that the triangulation is Delaunay 
        @async for T in each_triangle(_tri)
            for r in each_point_index(_tri)
                if r ∈ get_vertices(tri)
                    cert = DT.point_position_relative_to_circumcircle(_tri, T, r)
                    @test !DT.is_inside(cert)
                end
            end
        end

        ## Test that every edge is incident to two triangles if it not a boundary edge, and one otherwise 
        @async for (ij, k) in get_adjacent(_tri)
            i, j = initial(ij), terminal(ij)
            vij = get_adjacent(tri, i, j)
            vji = get_adjacent(tri, j, i)
            if DT.is_boundary_edge(_tri, i, j)
                @test DT.is_boundary_index(vij)
            elseif DT.is_boundary_edge(_tri, j, i)
                @test DT.is_boundary_index(vji)
            else
                @test DT.edge_exists(vij) && DT.edge_exists(vji)
            end
        end

        ## Test the adjacent map 
        @async for T in each_triangle(_tri)
            u, v, w = indices(T)
            @test get_adjacent(_tri, u, v) == w
            @test get_adjacent(_tri, v, w) == u
            @test get_adjacent(_tri, w, u) == v
        end

        ## Test the adjacent2vertex map 
        E = DT.edge_type(_tri)
        @async for T in each_triangle(_tri)
            u, v, w = indices(T)
            @test DT.contains_edge(DT.construct_edge(E, v, w), get_adjacent2vertex(_tri, u))
            @test DT.contains_edge(DT.construct_edge(E, w, u), get_adjacent2vertex(_tri, v))
            @test DT.contains_edge(DT.construct_edge(E, u, v), get_adjacent2vertex(_tri, w))
        end

        ## Check that the adjacent and adjacent2vertex maps are inverses 
        @async for (k, S) in get_adjacent2vertex(_tri)
            for ij in DT.each_edge(S)
                @test get_adjacent(_tri, ij) == k
            end
        end
        Base.Threads.@spawn for (ij, k) in get_adjacent(_tri)
            @test DT.contains_edge(ij, get_adjacent2vertex(_tri, k))
        end

        ## Check the graph 
        @async for (i, j) in get_edges(_tri)
            if DT.has_ghost_triangles(_tri)
                ij = DT.construct_edge(E, i, j)
                ji = DT.construct_edge(E, j, i)
                @test ij ∈ keys((get_adjacent ∘ get_adjacent)(_tri))
                @test ji ∈ keys((get_adjacent ∘ get_adjacent)(_tri))
            else
                if !DT.is_boundary_index(i) && !DT.is_boundary_index(j)
                    ij = DT.construct_edge(E, i, j)
                    ji = DT.construct_edge(E, j, i)
                    @test ij ∈ keys((get_adjacent ∘ get_adjacent)(_tri))
                    @test ji ∈ keys((get_adjacent ∘ get_adjacent)(_tri))
                end
            end
        end
    end
end

macro _adj(i, j, k)
    return :(($i, $j) => $k, ($j, $k) => $i, ($k, $i) => $j)
end
