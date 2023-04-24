using StaticArrays
using StableRNGs
using LinearAlgebra
using Random

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

function test_planarity(tri) # just check's Euler's formula. Doesn't guarantee planarity, though. Remember that Euler's formula includes the exterior face, which is covered by the boundary vertex.
    if DT.has_multiple_segments(tri)
        num_exterior_faces = DT.num_curves(tri)
        flag1 = (length ∘ collect ∘ each_solid_vertex)(tri) - (length ∘ collect ∘ each_solid_edge)(tri) + (length ∘ collect ∘ each_solid_triangle)(tri) + num_exterior_faces == 2
        if !flag1
            X = (length ∘ collect ∘ each_solid_vertex)(tri) - (length ∘ collect ∘ each_solid_edge)(tri) + (length ∘ collect ∘ each_solid_triangle)(tri) + num_exterior_faces
            println("Planarity test failed for the solid graph as the Euler characteristic was $X, not 2.")
            return false
        end
        nedges = length(keys(get_adjacent(get_adjacent(tri)))) ÷ 2 # when dealing with multiple boundary indices, num_edges(tri) is difficult as not every edge actually appears twice
        nverts = num_vertices(tri) - length(DT.all_boundary_indices(tri)) + DT.num_curves(tri) # add back in one boundary index for each exterior face 
        flag2 = nverts - nedges + num_triangles(tri) == 2
        if !flag2
            X = nverts - nedges + num_triangles(tri)
            println("Planarity test failed for the combined graph as the Euler characteristic was $X, not 2.")
            return false
        end
        flag3 = 2nedges == 3num_triangles(tri)
        if !flag3
            println("Planarity test failed as 2num_edges(tri), $nedges, did not equal 3num_triangles(tri), $(num_triangles(tri)).")
            return false
        end
    else
        flag = num_vertices(tri) - num_edges(tri) + num_triangles(tri) == 2
        if !flag
            X = num_vertices(tri) - num_edges(tri) + num_triangles(tri)
            println("Planarity test failed as the Euler characteristic was $X, not 2.")
            return false
        end
    end
    return true
end

function test_triangle_orientation(tri; check_ghost_triangle_orientation=true)
    for T in each_triangle(tri)
        cert = DT.triangle_orientation(tri, T)
        if DT.is_ghost_triangle(T) && !check_ghost_triangle_orientation
            continue
        end
        flag = DT.is_positively_oriented(cert)
        if !flag
            println("Orientation test failed for the triangle $T.")
            return false
        end
    end
    return true
end

function test_delaunay_criterion(tri; check_ghost_triangle_delaunay=true)
    for T in each_triangle(tri)
        if DT.is_ghost_triangle(T) && !check_ghost_triangle_delaunay
            continue
        end
        for r in each_solid_vertex(tri)
            cert = DT.point_position_relative_to_circumcircle(tri, T, r)
            if DT.is_inside(cert)
                ace = get_all_constrained_edges(tri)
                if DT.is_empty(ace)
                    flag = !DT.is_inside(cert)
                    if !flag
                        println("Delaunay criterion test failed for the triangle-vertex pair ($T, $r).")
                        return false
                    end
                else # This is extremely slow. Should probably get around to cleaning this up sometime.
                    i, j, k = DT.indices(T)
                    ## We need to see if a subsegment separates T from r. Let us walk from each vertex of T to r.
                    i_history = DT.PointLocationHistory{DT.triangle_type(tri),DT.edge_type(tri),DT.integer_type(tri)}()
                    j_history = DT.PointLocationHistory{DT.triangle_type(tri),DT.edge_type(tri),DT.integer_type(tri)}()
                    k_history = DT.PointLocationHistory{DT.triangle_type(tri),DT.edge_type(tri),DT.integer_type(tri)}()
                    DT.is_boundary_index(i) || jump_and_march(tri, get_point(tri, i); point_indices=nothing, m=nothing, try_points=nothing, k=r, store_history=Val(true), history=i_history)
                    DT.is_boundary_index(j) || jump_and_march(tri, get_point(tri, j); point_indices=nothing, m=nothing, try_points=nothing, k=r, store_history=Val(true), history=j_history)
                    DT.is_boundary_index(k) || jump_and_march(tri, get_point(tri, k); point_indices=nothing, m=nothing, try_points=nothing, k=r, store_history=Val(true), history=k_history)
                    i_walk = DT.is_boundary_index(i) ? Set{DT.triangle_type(tri)}() : i_history.triangles
                    j_walk = DT.is_boundary_index(j) ? Set{DT.triangle_type(tri)}() : j_history.triangles
                    k_walk = DT.is_boundary_index(k) ? Set{DT.triangle_type(tri)}() : k_history.triangles
                    all_edges = Set{NTuple{2,DT.integer_type(tri)}}()
                    E = DT.edge_type(tri)
                    for T in Iterators.flatten((each_triangle(i_walk), each_triangle(j_walk), each_triangle(k_walk)))
                        for e in DT.triangle_edges(T)
                            u, v = DT.edge_indices(e)
                            ee = DT.construct_edge(E, min(u, v), max(u, v))
                            push!(all_edges, (min(u, v), max(u, v)))
                        end
                    end
                    all_constrained_edges = Set{NTuple{2,DT.integer_type(tri)}}()
                    for e in each_edge(ace)
                        u, v = DT.edge_indices(e)
                        ee = DT.construct_edge(E, min(u, v), max(u, v))
                        push!(all_constrained_edges, (min(u, v), max(u, v)))
                    end
                    intersect!(all_edges, all_constrained_edges)
                    flags = [DT.line_segment_intersection_type(tri, initial(e), terminal(e), i, r) for e in each_edge(all_edges) for i in filter(!DT.is_boundary_index, (i, j, k))]
                    flag = !all(DT.is_none, flags) || isempty(flags)
                    if !flag
                        println("Delaunay criterion test failed for the triangle-vertex pair ($T, $r).")
                        return false
                    end
                end
            end
        end
    end
    return true
end

function test_each_edge_has_two_incident_triangles(tri)
    for e in each_edge(tri)
        i, j = DT.edge_indices(e)
        vᵢⱼ = get_adjacent(tri, i, j)
        vⱼᵢ = get_adjacent(tri, j, i)
        if DT.is_boundary_edge(tri, i, j)
            flag = DT.is_boundary_index(vᵢⱼ) && DT.edge_exists(vⱼᵢ)
            if !flag
                println("Test that each edge has two incident triangles failed for the edge $e.")
                return false
            end
        elseif DT.is_boundary_edge(tri, j, i)
            flag = DT.is_boundary_index(vⱼᵢ) && DT.edge_exists(vᵢⱼ)
            if !flag
                println("Test that each edge has two incident triangles failed for the edge $e.")
                return false
            end
        else
            flag = DT.edge_exists(vᵢⱼ) && DT.edge_exists(vⱼᵢ)
            if !flag
                println("Test that each edge has two incident triangles failed for the edge $e.")
                return false
            end
        end
    end
    DT.clear_empty_features!(tri)
    return true
end

function test_adjacent_map_matches_triangles(tri)
    for T in each_triangle(tri)
        u, v, w = DT.indices(T)
        flag1 = get_adjacent(tri, u, v) == w
        if !flag1
            println("Test that the adjacent map matches the triangle set failed for the triangle $T and edge ($u, $v).")
            return false
        end
        flag2 = get_adjacent(tri, v, w) == u
        if !flag2
            println("Test that the adjacent map matches the triangle set failed for the triangle $T and edge ($v, $w).")
            return false
        end
        flag3 = get_adjacent(tri, w, u) == v
        if !flag3
            println("Test that the adjacent map matches the triangle set failed for the triangle $T and edge ($w, $u).")
            return false
        end
    end
    DT.clear_empty_features!(tri)
    return true
end

function test_adjacent2vertex_map_matches_triangles(tri)
    E = DT.edge_type(tri)
    for T in each_triangle(tri)
        u, v, w = DT.indices(T)
        uv = DT.construct_edge(E, u, v)
        vw = DT.construct_edge(E, v, w)
        wu = DT.construct_edge(E, w, u)
        Su = get_adjacent2vertex(tri, u)
        Sv = get_adjacent2vertex(tri, v)
        Sw = get_adjacent2vertex(tri, w)
        flag1 = DT.contains_edge(vw, Su)
        if !flag1
            println("Test that the adjacent2vertex map matches the triangle set failed for the triangle $T and edge ($v, $w) for the edge set $Su from $u.")
            return false
        end
        flag2 = DT.contains_edge(wu, Sv)
        if !flag2
            println("Test that the adjacent2vertex map matches the triangle set failed for the triangle $T and edge ($w, $u) for the edge set $Sv from $v.")
            return false
        end
        flag3 = DT.contains_edge(uv, Sw)
        if !flag3
            println("Test that the adjacent2vertex map matches the triangle set failed for the triangle $T and edge ($u, $v) for the edge set $Sw from $w.")
            return false
        end
    end
    return true
end

function test_adjacent_map_matches_adjacent2vertex_map(tri)
    for (k, S) in get_adjacent2vertex(tri)
        for e in each_edge(S)
            flag1 = get_adjacent(tri, e) == k
            if !flag1
                println("Test that the adjacent map matches the adjacent2vertex map failed for the vertex $k and edge $e for the corresponding edge set $S, as get_adjacent(tri, $e) == $(get_adjacent(tri, e)) rather than $k.")
                return false
            end
            flag2 = DT.contains_edge(e, S)
            if !flag2
                println("Test that the adjacent map matches the adjacent2vertex map failed for the vertex $k and edge $e for the corresponding edge set $S, as the edge edge does not contain $e.")
                return false
            end
        end
    end
    DT.clear_empty_features!(tri)
    return true
end

function test_adjacent2vertex_map_matches_adjacent_map(tri)
    for (e, _) in get_adjacent(tri)
        k = get_adjacent(tri, e)
        S = get_adjacent2vertex(tri, k)
        flag = DT.contains_edge(e, S)
        if !flag
            println("Test that the adjacent2vertex map matches the adjacent map failed for the edge $e and adjacent vertex $k.")
            return false
        end
    end
    DT.clear_empty_features!(tri)
    return true
end

function test_graph_contains_all_vertices(tri)
    all_vertices = Set{Int64}()
    for T in each_triangle(tri)
        i, j, k = DT.indices(T)
        push!(all_vertices, i, j, k)
    end
    for i in all_vertices
        flag = i ∈ get_vertices(tri)
        if !flag
            println("Test that the graph contains all vertices failed as the vertex $i did not appear in the graph.")
            return false
        end
    end
    return true
end

function test_graph_matches_adjacent_map(tri)
    E = DT.edge_type(tri)
    adj_dict = get_adjacent(get_adjacent(tri))
    for e in get_edges(tri)
        i, j = DT.edge_indices(e)
        if DT.has_ghost_triangles(tri) || !DT.is_ghost_edge(i, j)
            eᵢⱼ = DT.construct_edge(E, i, j)
            eⱼᵢ = DT.construct_edge(E, j, i)
            flag = if !DT.has_multiple_segments(tri)
                eᵢⱼ ∈ keys(adj_dict) && eⱼᵢ ∈ keys(adj_dict)
            else
                DT.edge_exists(tri, i, j) && DT.edge_exists(tri, j, i)
            end
            if !flag
                println("Test that the graph matches the adjacent map failed for the edge $e.")
                return false
            end
        end
    end
    DT.clear_empty_features!(tri)
    return true
end

function test_adjacent_map_matches_graph(tri)
    for (e, k) in get_adjacent(tri)
        i, j = DT.edge_indices(e)
        flag1 = all(∈(get_neighbours(tri, i)), (j, k))
        if !flag1
            println("Test that the adjacent map matches the graph failed for the vertex $i, which defines the edge $e with adjacent vertex $k, as $j or $k were not in the neighbourhood of $i.")
        end
        flag2 = all(∈(get_neighbours(tri, j)), (k, i))
        if !flag2
            println("Test that the adjacent map matches the graph failed for the vertex $j, which defines the edge $e with adjacent vertex $k, as $k or $i were not in the neighbourhood of $j.")
        end
        flag3 = all(∈(get_neighbours(tri, k)), (i, j))
        if !flag3
            println("Test that the adjacent map matches the graph failed for the vertex $k, adjacent to the edge $e, as $i or $j were not in the neighbourhood of $k.")
        end
    end
    return true
end

function test_graph_matches_triangles(tri)
    for T in each_triangle(tri)
        i, j, k = DT.indices(T)
        flag1 = all(∈(get_neighbours(tri, i)), (j, k))
        if !flag1
            println("Test that the graph matches the triangle set failed for the triangle $T as $j or $k were not in the neighbourhood of $i.")
        end
        flag2 = all(∈(get_neighbours(tri, j)), (k, i))
        if !flag2
            println("Test that the graph matches the triangle set failed for the triangle $T as $k or $i were not in the neighbourhood of $j.")
        end
        flag3 = all(∈(get_neighbours(tri, k)), (i, j))
        if !flag3
            println("Test that the graph matches the triangle set failed for the triangle $T as $i or $j were not in the neighbourhood of $k.")
        end
    end
    return true
end

function test_constrained_edges(tri)
    for e in each_constrained_edge(tri)
        u, v = DT.edge_indices(e)
        flag = u ∈ get_neighbours(tri, v)
        if !flag
            println("Test that the constrained edge $e is in the triangulation failed.")
            return false
        end
        if !DT.contains_edge(e, get_constrained_edges(tri)) && !DT.contains_edge(DT.reverse_edge(e), get_constrained_edges(tri))
            flag = DT.contains_boundary_edge(tri, e) || DT.contains_boundary_edge(tri, DT.reverse_edge(e))
            if !flag
                println("Test that the constrained edge $e is in the triangulation passed, but it was not in the constrained_edge field and not in the boundary edges.")
                return false
            end
        end
    end
    return true
end

function test_boundary_edge_map_matches_boundary_nodes(tri::Triangulation)
    boundary_edge_map = DT.get_boundary_edge_map(tri)
    for (edge, pos) in boundary_edge_map
        flag = DT.get_boundary_nodes(DT.get_boundary_nodes(tri, pos[1]), pos[2]) == DT.initial(edge)
        if !flag
            println("Validation of the boundary edge map failed as the edge $edge, mapping to $pos, did not correspond to the correct edge in the boundary nodes.")
            return false
        end
    end
    return true
end

function test_boundary_nodes_matches_boundary_edge_map(tri::Triangulation)
    boundary_nodes = DT.get_boundary_nodes(tri)
    boundary_edge_map = DT.get_boundary_edge_map(tri)
    E = DT.edge_type(tri)
    if DT.has_multiple_curves(tri)
        nc = DT.num_curves(tri)
        for i in 1:nc
            curve_nodes = get_boundary_nodes(tri, i)
            ns = DT.num_segments(curve_nodes)
            for j in 1:ns
                segment_nodes = get_boundary_nodes(curve_nodes, j)
                ne = DT.num_boundary_edges(segment_nodes)
                for k in 1:ne
                    left_node = get_boundary_nodes(segment_nodes, k)
                    right_node = get_boundary_nodes(segment_nodes, k + 1)
                    edge = DT.construct_edge(E, left_node, right_node)
                    pos = DT.get_boundary_edge_map(tri, edge)
                    flag = pos == ((i, j), k)
                    if !flag
                        println("Validation of the boundary nodes failed as the edge $edge, located at $(((i, j), k)), was not correctly mapped to by the boundary edge map, instead giving $pos.")
                        return false
                    end
                end
            end
        end
    elseif DT.has_multiple_segments(tri)
        ns = DT.num_segments(tri)
        for i in 1:ns
            segment_nodes = get_boundary_nodes(tri, i)
            ne = DT.num_boundary_edges(segment_nodes)
            for j in 1:ne
                left_node = get_boundary_nodes(segment_nodes, j)
                right_node = get_boundary_nodes(segment_nodes, j + 1)
                edge = DT.construct_edge(E, left_node, right_node)
                pos = DT.get_boundary_edge_map(tri, edge)
                flag = pos == (i, j)
                if !flag
                    println("Validation of the boundary nodes failed as the edge $edge, located at $((i, j)), was not correctly mapped to by the boundary edge map, instead giving $pos.")
                    return false
                end
            end
        end
    else
        ne = DT.num_boundary_edges(boundary_nodes)
        for i in 1:ne
            left_node = get_boundary_nodes(boundary_nodes, i)
            right_node = get_boundary_nodes(boundary_nodes, i + 1)
            edge = DT.construct_edge(E, left_node, right_node)
            pos = DT.get_boundary_edge_map(tri, edge)
            flag = pos == (get_boundary_nodes(tri), i)
            if !flag
                println("Validation of the boundary nodes failed as the edge $edge, located at $((get_boundary_nodes(tri), i)), was not correctly mapped to by the boundary edge map, instead giving $pos.")
                return false
            end
        end
    end
    return true
end

function test_iterators(tri::Triangulation)
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
        i, j, k = DT.indices(T)
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
            if DT.is_boundary_index(k)
                push!(ghost_vertices, k)
            else
                push!(solid_vertices, k)
            end
        end
    end
    flag1 = allunique(all_triangles)
    flag2 = allunique(all_edges)
    flag3 = allunique(solid_triangles)
    flag4 = allunique(solid_edges)
    flag5 = allunique(ghost_triangles)
    flag6 = allunique(ghost_edges)
    if !flag1
        println("Validation of the triangle iterator failed as it returned duplicate triangles.")
        return false
    elseif !flag2
        println("Validation of the edge iterator failed as it returned duplicate edges.")
        return false
    elseif !flag3
        println("Validation of the solid triangle iterator failed as it returned duplicate triangles.")
        return false
    elseif !flag4
        println("Validation of the solid edge iterator failed as it returned duplicate edges.")
        return false
    elseif !flag5
        println("Validation of the ghost triangle iterator failed as it returned duplicate triangles.")
        return false
    elseif !flag6
        println("Validation of the ghost edge iterator failed as it returned duplicate edges.")
        return false
    end
    for (i, e) in enumerate(all_edges)
        u, v = DT.edge_indices(e)
        all_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(solid_edges)
        u, v = DT.edge_indices(e)
        solid_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(ghost_edges)
        u, v = DT.edge_indices(e)
        ghost_edges[i] = (min(u, v), max(u, v))
    end
    unique!(all_edges)
    unique!(solid_edges)
    unique!(ghost_edges)
    flag7 = length(all_triangles) == length(each_triangle(tri))
    flag8 = length(all_edges) == length(each_edge(tri))
    flag9 = length(solid_triangles) == length(each_solid_triangle(tri))
    flag10 = length(solid_edges) == length(each_solid_edge(tri))
    flag11 = length(ghost_triangles) == length(each_ghost_triangle(tri))
    flag12 = length(ghost_edges) == length(each_ghost_edge(tri))
    flag13 = length(all_vertices) == length(each_vertex(tri))
    flag14 = length(solid_vertices) == length(each_solid_vertex(tri))
    flag15 = length(ghost_vertices) == length(each_ghost_vertex(tri))
    if !flag7
        println("Validation of the triangle iterator failed as it returned a different number of triangles to the each_triangle iterator.")
        return false
    elseif !flag8
        println("Validation of the edge iterator failed as it returned a different number of edges to the each_edge iterator.")
        return false
    elseif !flag9
        println("Validation of the solid triangle iterator failed as it returned a different number of triangles to the each_solid_triangle iterator.")
        return false
    elseif !flag10
        println("Validation of the solid edge iterator failed as it returned a different number of edges to the each_solid_edge iterator.")
        return false
    elseif !flag11
        println("Validation of the ghost triangle iterator failed as it returned a different number of triangles to the each_ghost_triangle iterator.")
        return false
    elseif !flag12
        println("Validation of the ghost edge iterator failed as it returned a different number of edges to the each_ghost_edge iterator.")
        return false
    elseif !flag13
        println("Validation of the vertex iterator failed as it returned a different number of vertices to the each_vertex iterator.")
        return false
    elseif !flag14
        println("Validation of the solid vertex iterator failed as it returned a different number of vertices to the each_solid_vertex iterator.")
        return false
    elseif !flag15
        println("Validation of the ghost vertex iterator failed as it returned a different number of vertices to the each_ghost_vertex iterator.")
        return false
    end
    all_vertices = collect(all_vertices)
    solid_vertices = collect(solid_vertices)
    ghost_vertices = collect(ghost_vertices)
    sort!(all_vertices)
    sort!(solid_vertices)
    sort!(ghost_vertices)
    flag16 = DT.compare_triangle_collections(all_triangles, collect(each_triangle(tri)))
    flag17 = DT.compare_triangle_collections(solid_triangles, collect(each_solid_triangle(tri)))
    flag18 = DT.compare_triangle_collections(ghost_triangles, collect(each_ghost_triangle(tri)))
    flag19 = compare_edge_vectors(all_edges, each_edge(tri))
    flag20 = compare_edge_vectors(solid_edges, each_solid_edge(tri))
    flag21 = compare_edge_vectors(ghost_edges, each_ghost_edge(tri))
    flag22 = all_vertices == sort(collect(each_vertex(tri)))
    flag23 = solid_vertices == sort(collect(each_solid_vertex(tri)))
    flag24 = ghost_vertices == sort(collect(each_ghost_vertex(tri)))
    if !flag16
        println("Validation of the triangle iterator failed as it returned a different set of triangles to the each_triangle iterator.")
        return false
    elseif !flag17
        println("Validation of the solid triangle iterator failed as it returned a different set of triangles to the each_solid_triangle iterator.")
        return false
    elseif !flag18
        println("Validation of the ghost triangle iterator failed as it returned a different set of triangles to the each_ghost_triangle iterator.")
        return false
    elseif !flag19
        println("Validation of the edge iterator failed as it returned a different set of edges to the each_edge iterator.")
        return false
    elseif !flag20
        println("Validation of the solid edge iterator failed as it returned a different set of edges to the each_solid_edge iterator.")
        return false
    elseif !flag21
        println("Validation of the ghost edge iterator failed as it returned a different set of edges to the each_ghost_edge iterator.")
        return false
    elseif !flag22
        println("Validation of the vertex iterator failed as it returned a different set of vertices to the each_vertex iterator.")
        return false
    elseif !flag23
        println("Validation of the solid vertex iterator failed as it returned a different set of vertices to the each_solid_vertex iterator.")
        return false
    elseif !flag24
        println("Validation of the ghost vertex iterator failed as it returned a different set of vertices to the each_ghost_vertex iterator.")
        return false
    end
    return true
end

function validate_triangulation(_tri::Triangulation; check_planarity=true, check_ghost_triangle_orientation=true, check_ghost_triangle_delaunay=true) # doesn't work for non-convex. need to find a better way
    tri = deepcopy(_tri)
    DT.delete_ghost_triangles!(tri)
    DT.add_ghost_triangles!(tri)
    DT.clear_empty_features!(tri)
    return (!check_planarity || test_planarity(tri)) &&
           test_triangle_orientation(tri; check_ghost_triangle_orientation) &&
           test_delaunay_criterion(tri; check_ghost_triangle_delaunay) &&
           test_each_edge_has_two_incident_triangles(tri) &&
           test_adjacent2vertex_map_matches_triangles(tri) &&
           test_adjacent2vertex_map_matches_triangles(tri) &&
           test_adjacent_map_matches_adjacent2vertex_map(tri) &&
           test_adjacent2vertex_map_matches_adjacent_map(tri) &&
           test_graph_contains_all_vertices(tri) &&
           test_graph_matches_adjacent_map(tri) &&
           test_graph_matches_triangles(tri) &&
           test_adjacent_map_matches_graph(tri) &&
           test_constrained_edges(tri) &&
           test_boundary_edge_map_matches_boundary_nodes(tri) &&
           test_boundary_nodes_matches_boundary_edge_map(tri) &&
           test_iterators(tri)
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
    rep = DT.get_empty_representative_points()
    DT.compute_representative_points!(rep, pts, [2, 6, 4, 5, 2])
    tri = Triangulation(pts, T, adj, adj2v, DG,
        Int64[],
        DT.construct_boundary_edge_map(Int64[]),
        DT.DataStructures.OrderedDict{Int64,Vector{Int64}}(),
        DT.DataStructures.OrderedDict{Int64,UnitRange{Int64}}(),
        Set{NTuple{2,Int64}}(),
        Set{NTuple{2,Int64}}(),
        ConvexHull(pts, [2, 6, 4, 5, 2]),
        rep)
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
    rep = DT.get_empty_representative_points()
    DT.compute_representative_points!(rep, pts, [1, 2, 3, 1])
    tri = Triangulation(pts, T, adj, adj2v, DG,
        Int64[],
        DT.construct_boundary_edge_map(Int64[]),
        DT.DataStructures.OrderedDict{Int64,Vector{Int64}}(),
        DT.DataStructures.OrderedDict{Int64,UnitRange{Int64}}(),
        Set{NTuple{2,Int64}}(),
        Set{NTuple{2,Int64}}(),
        ConvexHull(pts, [1, 2, 3, 1]),
        rep)
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
            DT.get_representative_point_list(tri),
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
            DT.get_representative_point_list(tri),
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
    !isnothing(current_constrained_edges) && @test compare_edge_vectors(constrained_edges, current_constrained_edges)
end

function get_random_vertices_and_constrained_edges(nverts1, nverts2, nedges, rng=Random.default_rng())
    ## To generate a random set of constrained edges, we get a random small triangulation, 
    ## and we just take the edges from that triangulation.
    points = [Tuple(rand(rng, 2)) for _ in 1:nverts1]
    tri = triangulate(points; rng)
    edges = Set{NTuple{2,Int64}}()
    all_edges = collect(each_solid_edge(tri))
    iter = 0
    while length(edges) < nedges && iter < 10000
        S = DT.random_edge(all_edges, rng)
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

function validate_statistics(tri::Triangulation, stats=statistics(tri))
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
        i, j, k = DT.indices(T)
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
            if DT.is_boundary_index(k)
                push!(ghost_vertices, k)
            else
                push!(solid_vertices, k)
            end
        end
    end
    push!(all_vertices, DT.all_boundary_indices(tri)...)
    push!(ghost_vertices, DT.all_boundary_indices(tri)...)
    for e in each_ghost_edge(tri)
        u, v = DT.edge_indices(e)
        push!(ghost_edges, (min(u, v), max(u, v)))
        push!(all_edges, (min(u, v), max(u, v)))
    end
    for (i, e) in enumerate(all_edges)
        u, v = DT.edge_indices(e)
        all_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(solid_edges)
        u, v = DT.edge_indices(e)
        solid_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(ghost_edges)
        u, v = DT.edge_indices(e)
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
    total_A = 0.0
    for T in each_solid_triangle(tri)
        u, v, w = DT.indices(T)
        p, q, r = get_point(tri, u, v, w)
        p = [getx(p), gety(p)]
        q = [getx(q), gety(q)]
        r = [getx(r), gety(r)]
        areas[indices(T)] = 0.5 * (p[1] * (q[2] - r[2]) + q[1] * (r[2] - p[2]) + r[1] * (p[2] - q[2]))
        total_A += areas[indices(T)]
        ℓ1 = norm(q - p)
        ℓ2 = norm(r - q)
        ℓ3 = norm(r - p)
        ℓmin = min(ℓ1, ℓ2, ℓ3)
        ℓmax = max(ℓ1, ℓ2, ℓ3)
        ℓmed = ℓ1 + ℓ2 + ℓ3 - ℓmin - ℓmax
        ℓ1, ℓ2, ℓ3 = ℓmin, ℓmed, ℓmax
        lengths[indices(T)] = (ℓ1, ℓ2, ℓ3)
        r′ = p - r
        s′ = q - r
        ox = r[1] + det([norm(r′)^2 r′[2]; norm(s′)^2 s′[2]]) / (4areas[indices(T)])
        oy = r[2] + det([r′[1] norm(r′)^2; s′[1] norm(s′)^2]) / (4areas[indices(T)])
        circumcenters[indices(T)] = (ox, oy)
        circumradii[indices(T)] = norm(r - collect(circumcenters[indices(T)]))
        all_angles = [(norm(p - r)^2 + norm(q - r)^2 - norm(p - q)^2) / (2norm(p - r) * norm(q - r)) for (p, q, r) in ((p, q, r), (q, r, p), (r, p, q))]
        all_angles[all_angles.<-1.0] .= -1.0
        all_angles[all_angles.>1.0] .= 1.0
        all_angles = acos.(all_angles)
        sort!(all_angles)
        radius_edge_ratio[indices(T)] = circumradii[indices(T)] / ℓ1
        edge_midpoints[indices(T)] = ((Tuple(0.5 * (p + q))), Tuple(0.5 * (q + r)), Tuple(0.5 * (r + p)))
        inradius[indices(T)] = 2areas[indices(T)] / (ℓ1 + ℓ2 + ℓ3)
        perimeter[indices(T)] = ℓ1 + ℓ2 + ℓ3
        aspect_ratio[indices(T)] = inradius[indices(T)] / circumradii[indices(T)]
        centroid[indices(T)] = (1 / 3 * (p[1] + q[1] + r[1]), 1 / 3 * (p[2] + q[2] + r[2]))
        angles[indices(T)] = Tuple(all_angles)
        @test radius_edge_ratio[indices(T)] ≥ 1 / sqrt(3) - 0.1
        @test DT.get_radius_edge_ratio(stats, T) ≥ 1 / sqrt(3) - 0.1
        @test angles[indices(T)][1] ≤ deg2rad(60) + 0.01
        @test DT.get_minimum_angle(stats, T) ≤ deg2rad(60) + 0.01
    end

    ## Now compare the statistics 
    for T in each_solid_triangle(tri)
        @test areas[indices(T)] ≈ DT.get_area(stats, T) rtol = 1e-4 atol = 1e-4
        @test collect(lengths[indices(T)]) ≈ collect(DT.get_lengths(stats, T)) rtol = 1e-4 atol = 1e-4
        @test collect(circumcenters[indices(T)]) ≈ collect(DT.get_circumcenter(stats, T)) rtol = 1e-4 atol = 1e-4
        @test circumradii[indices(T)] ≈ DT.get_circumradius(stats, T) rtol = 1e-4 atol = 1e-4
        @test radius_edge_ratio[indices(T)] ≈ DT.get_radius_edge_ratio(stats, T) rtol = 1e-4 atol = 1e-4
        @test collect(collect.(edge_midpoints[indices(T)])) ≈ collect(collect.(DT.get_edge_midpoints(stats, T))) rtol = 1e-4 atol = 1e-4
        @test aspect_ratio[indices(T)] ≈ DT.get_aspect_ratio(stats, T) rtol = 1e-4 atol = 1e-4
        @test inradius[indices(T)] ≈ DT.get_inradius(stats, T) rtol = 1e-4 atol = 1e-4
        @test perimeter[indices(T)] ≈ DT.get_perimeter(stats, T) rtol = 1e-4 atol = 1e-4
        @test radius_edge_ratio[indices(T)] ≈ 1 / (2sin(angles[indices(T)][1])) rtol = 1e-4 atol = 1e-4
        @test (2sin(DT.get_minimum_angle(stats, T) / 2)^2 - 0.1 ≤ DT.get_aspect_ratio(stats, T) ≤ 2tan(DT.get_minimum_angle(stats, T) / 2) + 0.1)
        @test DT.get_radius_edge_ratio(stats, T) ≈ 1 / (2(sin(DT.get_minimum_angle(stats, T)))) rtol = 1e-4 atol = 1e-4
        @test areas[indices(T)] ≈ inradius[indices(T)] * 0.5perimeter[indices(T)] rtol = 1e-4 atol = 1e-4
        @test DT.get_area(stats, T) ≈ DT.get_inradius(stats, T) * 0.5DT.get_perimeter(stats, T) rtol = 1e-4 atol = 1e-4
        @test collect(centroid[indices(T)]) ≈ collect(DT.get_centroid(stats, T)) rtol = 1e-4 atol = 1e-4
        @test DT.get_angles(stats, T)[1] ≈ angles[indices(T)][1] rtol = 1e-4 atol = 1e-4
        @test DT.get_angles(stats, T)[2] ≈ angles[indices(T)][2] rtol = 1e-4 atol = 1e-4
        @test DT.get_angles(stats, T)[3] ≈ angles[indices(T)][3] rtol = 1e-4 atol = 1e-4
        @test sum(DT.get_angles(stats, T)) ≈ π rtol = 1e-4 atol = 1e-4
        @test DT.get_minimum_angle(stats, T) ≈ angles[indices(T)][1] rtol = 1e-4 atol = 1e-4
        @test DT.get_maximum_angle(stats, T) ≈ angles[indices(T)][3] rtol = 1e-4 atol = 1e-4
        @test DT.get_minimum_angle(stats, T) ≈ DT.get_angles(stats, T)[1] rtol = 1e-4 atol = 1e-4
        @test DT.get_maximum_angle(stats, T) ≈ DT.get_angles(stats, T)[3] rtol = 1e-4 atol = 1e-4
    end
    @test stats.individual_statistics == DT.get_individual_statistics(stats)
    @test stats.total_area ≈ DT.get_total_area(stats) rtol = 1e-4 atol = 1e-4
    @test stats.total_area ≈ total_A rtol = 1e-4 atol = 1e-4

    ## Test the number statistics 
    @test DT.num_vertices(stats) == length(all_vertices) == stats.num_vertices
    @test DT.num_solid_vertices(stats) == length(solid_vertices) == stats.num_solid_vertices
    @test DT.num_ghost_vertices(stats) == length(DT.all_boundary_indices(tri)) == length(ghost_vertices) == stats.num_ghost_vertices
    @test DT.num_triangles(stats) == length(all_triangles) == stats.num_triangles
    @test DT.num_solid_triangles(stats) == length(solid_triangles) == stats.num_solid_triangles
    @test DT.num_ghost_triangles(stats) == length(ghost_triangles) == stats.num_ghost_triangles
    @test DT.num_edges(stats) == length(all_edges) == stats.num_edges
    @test DT.num_solid_edges(stats) == length(solid_edges) == stats.num_solid_edges
    @test DT.num_ghost_edges(stats) == length(ghost_edges) == stats.num_ghost_edges
    @test DT.num_constrained_boundary_edges(stats) == length(keys(get_boundary_edge_map(tri))) == stats.num_constrained_boundary_edges
    @test DT.num_constrained_interior_edges(stats) == num_edges(get_constrained_edges(tri)) == stats.num_constrained_interior_edges
    @test DT.num_constrained_edges(stats) == num_edges(get_all_constrained_edges(tri)) == stats.num_constrained_edges
    @test DT.num_convex_hull_points(stats) == length(get_convex_hull_indices(tri)) - 1 == stats.num_convex_hull_points

    ## Global statistics  
    smallest_angle = minimum([angles[indices(T)][1] for T in each_solid_triangle(tri)])
    largest_angle = maximum([angles[indices(T)][3] for T in each_solid_triangle(tri)])
    smallest_area = minimum([areas[indices(T)] for T in each_solid_triangle(tri)])
    largest_area = maximum([areas[indices(T)] for T in each_solid_triangle(tri)])
    smallest_radius_edge_ratio = minimum([radius_edge_ratio[indices(T)] for T in each_solid_triangle(tri)])
    largest_radius_edge_ratio = maximum([radius_edge_ratio[indices(T)] for T in each_solid_triangle(tri)])
    @test DT.get_smallest_angle(stats) ≈ smallest_angle rtol = 1e-2
    @test DT.get_largest_angle(stats) ≈ largest_angle rtol = 1e-2
    @test DT.get_smallest_area(stats) ≈ smallest_area rtol = 1e-2
    @test DT.get_largest_area(stats) ≈ largest_area rtol = 1e-2
    @test DT.get_smallest_radius_edge_ratio(stats) ≈ smallest_radius_edge_ratio rtol = 1e-2
    @test DT.get_largest_radius_edge_ratio(stats) ≈ largest_radius_edge_ratio rtol = 1e-2
    @test DT.get_smallest_radius_edge_ratio(stats) ≥ 1 / sqrt(3) - 0.1
    @test DT.get_smallest_angle(stats) ≤ deg2rad(60) + 0.01
end

function slow_encroachment_test(tri::Triangulation)
    E = DT.edge_type(tri)
    I = DT.integer_type(tri)
    ch = Channel{Pair{E,Tuple{Bool,I}}}(Inf) # https://discourse.julialang.org/t/can-dicts-be-threadsafe/27172/17
    @sync for i in collect(each_solid_vertex(tri))
        Base.Threads.@spawn begin
            @show i
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
    tri = Triangulation(pts)
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


## TODO: Implement a brute-force VoronoiTessellation that we can compare with
function validate_tessellation(vorn::VoronoiTessellation; check_convex=true, check_adjacent=true)
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
        println("Triangulation is not correct.")
        return false
    end
    circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
    triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
    for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
        V = DT.rotate_triangle_to_standard_form(V)
        c = DT.get_triangle_to_circumcenter(vorn, V)
        c = get_polygon_point(vorn, c)
        i, j, k = indices(V)
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
            i, j, k = indices(V)
            for i in (i, j, k)
                for e in each_edge(get_adjacent2vertex(DT.get_triangulation(vorn), i))
                    v, w = DT.edge_indices(e)
                    V′ = DT.rotate_triangle_to_standard_form(DT.construct_triangle(DT.triangle_type(DT.get_triangulation(vorn)), i, v, w))
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
                _poly_points = _pts[get_indices(ch)]
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
        @test isapprox(A, get_total_area(vorn.triangulation), rtol=1e-4)
    end
    return true
end
