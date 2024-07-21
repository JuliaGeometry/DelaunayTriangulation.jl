#=
"""
    abstract type AbstractTriangulationState 

`AbstractTriangulationState`s are used for representing the 
state of a certain test of a triangulation. 

# Interface 
- Subtypes should have a `flag` field that is `true` if the state is positive, and `false` otherwise. 
    This allows `test_state` to return `true` and `false` according to the flag.
- Subtypes should implement `Base.summary(state)` which creates a string summarising the state. This is 
    used in `show` for printing the `state`.
"""
=#

abstract type AbstractTriangulationState end
test_state(state::AbstractTriangulationState) = state.flag
Base.show(io::IO, state::AbstractTriangulationState) = print(io, summary(state))

function compare_edge_vectors(E1, E2)
    E1s = sort_edge_vector(collect(E1))
    E2s = sort_edge_vector(collect(E2))
    return E1s == E2s
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

struct TriangleOrientationState <: AbstractTriangulationState
    flag::Bool
    bad_triangle::NTuple{3, Int}
    triangle_orientation::Symbol
end
TriangleOrientationState(tri) = test_triangle_orientation(tri)
function Base.summary(state::TriangleOrientationState)
    if test_state(state)
        return "All the triangles have positive orientation."
    else
        return "The triangle $(state.bad_triangle) is $(state.triangle_orientation) oriented."
    end
end
function test_triangle_orientation(tri)
    for T in each_solid_triangle(tri)
        cert = triangle_orientation(tri, T)
        flag = is_positively_oriented(cert)
        orientation = is_positively_oriented(cert) ? :positively : is_negatively_oriented(cert) ? :negatively : :degenerately
        !flag && return TriangleOrientationState(flag, Int.(triangle_vertices(T)), orientation)
    end
    return TriangleOrientationState(true, (∅, ∅, ∅), :positive)
end

struct DelaunayCriterionState <: AbstractTriangulationState
    flag::Bool
    bad_triangle::NTuple{3, Int}
    bad_vertex::Int
end
DelaunayCriterionState(tri) = test_delaunay_criterion(tri)
function Base.summary(state::DelaunayCriterionState)
    if test_state(state)
        return "All the triangles are Delaunay."
    else
        if state.bad_vertex !== ∅
            return "The Delaunay criterion does not hold for the triangle-vertex pair ($(state.bad_triangle), $(state.bad_vertex))."
        else
            return "The test of the Delaunay criterion failed as there was a BoundsError when testing the visibility."
        end
    end
end
function test_delaunay_criterion(tri)
    try
        points = get_points(tri)
        triangle_tree = BoundaryRTree(points)
        segment_tree = BoundaryRTree(points)
        failures = Tuple{triangle_type(tri), integer_type(tri)}[]
        for T in each_solid_triangle(tri)
            i, j, k = triangle_vertices(T)
            p, q, r = get_point(tri, i, j, k)
            c = triangle_circumcenter(p, q, r)
            cr = triangle_circumradius(p, q, r)
            xmin, xmax = getx(c) - cr, getx(c) + cr
            ymin, ymax = gety(c) - cr, gety(c) + cr
            any(!isfinite, (xmin, xmax, ymin, ymax)) && continue
            bbox = BoundingBox(xmin, xmax, ymin, ymax)
            dbx = DiametralBoundingBox(bbox, (i, j)) # just store (i, j) and get k via get_adjacent 
            insert!(triangle_tree.tree, dbx)
        end
        for e in each_segment(tri)
            i, j = edge_vertices(e)
            p, q = get_point(tri, i, j)
            px, py = getxy(p)
            qx, qy = getxy(q)
            xmin = min(px, qx)
            xmax = max(px, qx)
            ymin = min(py, qy)
            ymax = max(py, qy)
            bbox = BoundingBox(xmin, xmax, ymin, ymax)
            dbx = DiametralBoundingBox(bbox, (i, j))
            insert!(segment_tree.tree, dbx)
        end
        for r in shuffle(collect(each_solid_vertex(tri)))
            !isempty(failures) && break
            intersects = get_intersections(triangle_tree, r, cache_id = 1) # can't use multithreading here
            for box in intersects
                i, j = get_edge(box)
                k = get_adjacent(tri, i, j)
                any(==(r), (i, j, k)) && continue
                T = construct_triangle(triangle_type(tri), i, j, k)
                cert = point_position_relative_to_circumcircle(tri, T, r)
                c = triangle_centroid(get_point(tri, i, j, k)...)
                if is_inside(cert)
                    is_boundary_edge(tri, i, j) && is_right(point_position_relative_to_line(tri, i, j, r)) && continue # if it's outside of the domain relative to this edge, just continue
                    is_boundary_edge(tri, j, k) && is_right(point_position_relative_to_line(tri, j, k, r)) && continue
                    is_boundary_edge(tri, k, i) && is_right(point_position_relative_to_line(tri, k, i, r)) && continue
                    for (i, j) in triangle_edges(i, j, k)
                        contains_segment(tri, i, j) && continue # visibility is defined according to the relative interior of the simplex, which means that it's fine if a segment can see the vertex
                        if rand() < 1 / 2 # just testing both
                            cert = test_visibility(tri, i, j, r, shift = 0.01, attractor = c)
                        else
                            cert = test_visibility(tri, segment_tree, i, j, r, c)
                        end
                        flag = is_visible(cert)
                        flag && push!(failures, (T, r))
                        flag && break
                    end
                end
            end
        end
        if isempty(failures)
            return DelaunayCriterionState(true, (∅, ∅, ∅), ∅)
        else
            return DelaunayCriterionState(false, Int.(failures[1][1]), Int(failures[1][2]))
        end
    catch e
        e isa BoundsError && return DelaunayCriterionState(false, (∅, ∅, ∅), ∅)
        rethrow(e)
    end
end
function test_visibility(tri::Triangulation, segment_tree, i, j, k, centroid = nothing)
    if is_boundary_edge(tri, j, i)
        i, j = j, i
    end
    p, q, a = get_point(tri, i, j, k)
    if is_boundary_edge(tri, i, j)
        side_e = point_position_relative_to_line(p, q, a)
        is_right(side_e) && return Cert.Invisible
    end
    # Need to see if i or j is a boundary node without the other being a boundary node.
    # If this is the case, it's possible that one of them mistakenly sees k.
    for u in (i, j)
        flag, g = is_boundary_node(tri, u)
        if flag
            ℓ = get_left_boundary_node(tri, u, g)
            cert = point_position_relative_to_line(tri, ℓ, u, k)
            is_right(cert) && return Cert.Invisible
            ℓ = get_right_boundary_node(tri, u, g)
            cert = point_position_relative_to_line(tri, u, ℓ, k)
            is_right(cert) && return Cert.Invisible
        end
    end
    flags = falses(10)
    ts = LinRange(0.00001, 0.99999, 10)
    for idx in 1:10
        t = ts[idx]
        m = p .+ t .* (q .- p)
        # Need to push m towards the centroid so that we are not taking points right on the edge, 
        # since visibility is defined relative to the interior of the simplex. 
        if !isnothing(centroid)
            m = m .+ 0.01 .* (centroid .- m)
        end
        mx, my = getxy(m)
        ax, ay = getxy(a)
        xmin, ymin, xmax, ymax = min(mx, ax), min(my, ay), max(mx, ax), max(my, ay)
        bbox = BoundingBox(xmin, xmax, ymin, ymax)
        intersections = get_intersections(segment_tree, bbox; cache_id = 2)
        for box in intersections
            u, v = get_edge(box)
            p′, q′ = get_point(tri, u, v)
            cert = line_segment_intersection_type(m, a, p′, q′)
            flags[idx] = !has_no_intersections(cert) && !is_touching(cert)
            flags[idx] && break
        end
    end
    return all(flags) ? Cert.Invisible : Cert.Visible
end

struct EdgesHaveTwoIncidentTrianglesState <: AbstractTriangulationState
    flag::Bool
    bad_edge::NTuple{2, Int}
end
EdgesHaveTwoIncidentTrianglesState(tri) = test_each_edge_has_two_incident_triangles(tri)
function Base.summary(state::EdgesHaveTwoIncidentTrianglesState)
    if test_state(state)
        return "All the edges have two incident triangles."
    else
        return "The edge $(state.bad_edge) does not have two incident triangles."
    end
end
function test_each_edge_has_two_incident_triangles(tri)
    for e in each_edge(tri)
        i, j = edge_vertices(e)
        vᵢⱼ = get_adjacent(tri, i, j)
        vⱼᵢ = get_adjacent(tri, j, i)
        if is_boundary_edge(tri, j, i)
            flag = is_ghost_vertex(vᵢⱼ) && edge_exists(vⱼᵢ)
            !flag && return EdgesHaveTwoIncidentTrianglesState(flag, Int.(edge_vertices(e)))
        elseif is_boundary_edge(tri, i, j)
            flag = is_ghost_vertex(vⱼᵢ) && edge_exists(vᵢⱼ)
            !flag && return EdgesHaveTwoIncidentTrianglesState(flag, Int.(edge_vertices(e)))
        else
            flag = edge_exists(vᵢⱼ) && edge_exists(vⱼᵢ)
            !flag && return EdgesHaveTwoIncidentTrianglesState(flag, Int.(edge_vertices(e)))
        end
    end
    clear_empty_features!(tri)
    return EdgesHaveTwoIncidentTrianglesState(true, (∅, ∅))
end

struct AdjacentMapState <: AbstractTriangulationState
    flag::Bool
    bad_triangle::NTuple{3, Int}
    bad_edge::NTuple{2, Int}
    bad_vertex::Int
end
AdjacentMapState(tri) = test_adjacent_map_matches_triangles(tri)
function Base.summary(state::AdjacentMapState)
    if test_state(state)
        return "The adjacent map is correct."
    else
        return "The triangle $(state.bad_triangle) is incorrectly stored in the adjacent map, with the edge $(state.bad_edge) instead mapping to $(state.bad_vertex)."
    end
end
function test_adjacent_map_matches_triangles(tri)
    for T in each_triangle(tri)
        u, v, w = triangle_vertices(T)
        flag1 = get_adjacent(tri, u, v) == w
        !flag1 && return AdjacentMapState(flag1, Int.(triangle_vertices(T)), Int.((u, v)), get_adjacent(tri, u, v))
        flag2 = get_adjacent(tri, v, w) == u
        !flag2 && return AdjacentMapState(flag2, Int.(triangle_vertices(T)), Int.((v, w)), get_adjacent(tri, v, w))
        flag3 = get_adjacent(tri, w, u) == v
        !flag3 && return AdjacentMapState(flag3, Int.(triangle_vertices(T)), Int.((w, u)), get_adjacent(tri, w, u))
    end
    clear_empty_features!(tri)
    return AdjacentMapState(true, (∅, ∅, ∅), (∅, ∅), ∅)
end

struct Adjacent2VertexMapState <: AbstractTriangulationState
    flag::Bool
    bad_triangle::NTuple{3, Int}
    bad_edge::NTuple{2, Int}
    bad_vertex::Int
    bad_edge_set::Set{NTuple{2, Int}}
end
Adjacent2VertexMapState(tri) = test_adjacent2vertex_map_matches_triangles(tri)
function Base.summary(state::Adjacent2VertexMapState)
    if test_state(state)
        return "The adjacent2vertex map is correct."
    else
        return "The adjacent2vertex map is incorrect for the triangle $(state.bad_triangle), as the edge $(state.bad_edge) is not inside the vertex $(state.bad_vertex)'s associated edge set in the map, $(state.bad_edge_set)."
    end
end
function test_adjacent2vertex_map_matches_triangles(tri)
    E = edge_type(tri)
    for T in each_triangle(tri)
        u, v, w = triangle_vertices(T)
        for (u, v, w) in ((u, v, w), (v, w, u), (w, u, v))
            vw = construct_edge(E, v, w)
            Su = get_adjacent2vertex(tri, u)
            flag = contains_edge(vw, Su)
            if !flag
                _Su = Set{NTuple{2, Int}}(Int.(edge_vertices(e)) for e in each_edge(Su))
                return Adjacent2VertexMapState(flag, Int.(triangle_vertices(T)), Int.((v, w)), Int(u), _Su)
            end
        end
    end
    return Adjacent2VertexMapState(true, (∅, ∅, ∅), (∅, ∅), ∅, Set{NTuple{2, Int}}())
end

struct AdjacentMapAdjacent2VertexMapState <: AbstractTriangulationState
    flag::Bool
    bad_vertex::Int
    bad_vertex2::Int
    bad_edge_set::Set{NTuple{2, Int}}
    bad_edge::NTuple{2, Int}
end
AdjacentMapAdjacent2VertexMapState(tri) = test_adjacent_map_matches_adjacent2vertex_map(tri)
function Base.summary(state::AdjacentMapAdjacent2VertexMapState)
    if test_state(state)
        return "The adjacent2vertex map is consistent with the adjacent map."
    else
        return "The adjacent2vertex map is inconsistent with the adjacent map. The vertex $(state.bad_vertex) and edge $(state.bad_edge) are inconsistent with the associated edge set $(state.bad_edge_set) from the adjacent2vertex map, as $(state.bad_edge) maps to $(state.bad_vertex2) instead of $(state.bad_vertex)."
    end
end
function test_adjacent_map_matches_adjacent2vertex_map(tri)
    for (k, S) in get_adjacent2vertex(get_adjacent2vertex(tri))
        _S = Set{NTuple{2, Int}}(Int.(edge_vertices(e)) for e in each_edge(S))
        for e in each_edge(S)
            flag = get_adjacent(tri, e) == k
            !flag && return AdjacentMapAdjacent2VertexMapState(flag, Int(k), get_adjacent(tri, e), _S, Int.(edge_vertices(e)))
        end
    end
    clear_empty_features!(tri)
    return AdjacentMapAdjacent2VertexMapState(true, ∅, ∅, Set{NTuple{2, Int}}(), (∅, ∅))
end

struct Adjacent2VertexMapAdjacentMapState <: AbstractTriangulationState
    flag::Bool
    bad_vertex::Int
    bad_edge::NTuple{2, Int}
    in_adjacent2vertex::Bool
end
Adjacent2VertexMapAdjacentMapState(tri) = test_adjacent2vertex_map_matches_adjacent_map(tri)
function Base.summary(state::Adjacent2VertexMapAdjacentMapState)
    if test_state(state)
        return "The adjacent2vertex map is consistent with the adjacent map."
    else
        if !state.in_adjacent2vertex
            return "The adjacent2vertex map is inconsistent with the adjacent map. The edge $(state.bad_edge) is mapped to $(state.bad_vertex) but $(state.bad_vertex) is not a key in the adjacent2vertex map."
        else
            return "The adjacent2vertex map is inconsistent with the adjacent map. The edge $(state.bad_edge) is mapped to $(state.bad_vertex) but $(status.bad_edge) is not in the edge set that $(state.bad_vertex) maps to in the adjacent2vertex map."
        end
    end
end
function test_adjacent2vertex_map_matches_adjacent_map(tri)
    adj2v = get_adjacent2vertex(get_adjacent2vertex(tri))
    for (e, _) in get_adjacent(get_adjacent(tri))
        k = get_adjacent(tri, e)
        if haskey(adj2v, k)
            S = get_adjacent2vertex(tri, k)
            flag = contains_edge(e, S)
        else
            flag = false
        end
        !flag && return Adjacent2VertexMapAdjacentMapState(flag, Int(k), Int.(edge_vertices(e)), haskey(adj2v, k))
    end
    clear_empty_features!(tri)
    return Adjacent2VertexMapAdjacentMapState(true, ∅, (∅, ∅), true)
end

struct GraphState <: AbstractTriangulationState
    flag::Bool
    bad_vertex::Int
end
GraphState(tri) = test_graph_contains_all_vertices(tri)
function Base.summary(state::GraphState)
    if test_state(state)
        return "The graph correctly contains all vertices in the triangulation."
    else
        return "The graph does not include the vertex $(state.bad_vertex) despite it being in the triangulation."
    end
end
function test_graph_contains_all_vertices(tri)
    all_vertices = Set{Int}()
    for T in each_triangle(tri) # need a method that doesn't use the graph
        i, j, k = triangle_vertices(T)
        push!(all_vertices, i, j, k)
    end
    for i in all_vertices
        flag = has_vertex(tri, i)
        !flag && return GraphState(flag, i)
    end
    return GraphState(true, ∅)
end

struct GraphAdjacentMapState <: AbstractTriangulationState
    flag::Bool
    bad_edge::NTuple{2, Int}
end
GraphAdjacentMapState(tri) = test_graph_matches_adjacent_map(tri)
function Base.summary(state::GraphAdjacentMapState)
    if test_state(state)
        return "The graph's edges correctly matches the keys of the adjacent map."
    else
        return "The edge $(state.bad_edge) appears as an edge in the graph but it and its reverse are not both a key of the adjacent map."
    end
end
function test_graph_matches_adjacent_map(tri)
    E = edge_type(tri)
    adj_dict = get_adjacent(get_adjacent(tri))
    for e in each_edge(tri)
        i, j = edge_vertices(e)
        if has_ghost_triangles(tri) || !is_ghost_edge(i, j)
            eᵢⱼ = construct_edge(E, i, j)
            eⱼᵢ = construct_edge(E, j, i)
            flag = if !has_multiple_sections(tri)
                eᵢⱼ ∈ keys(adj_dict) && eⱼᵢ ∈ keys(adj_dict)
            else
                edge_exists(tri, i, j) && edge_exists(tri, j, i)
            end
            !flag && return GraphAdjacentMapState(flag, Int.(edge_vertices(e)))
        end
    end
    clear_empty_features!(tri)
    return GraphAdjacentMapState(true, (∅, ∅))
end

struct AdjacentMapGraphState <: AbstractTriangulationState
    flag::Bool
    bad_vertex::Int
    bad_edge::NTuple{2, Int}
    bad_edge_2::NTuple{2, Int}
    is_k::Bool
    adjacent_vertex::Int
end
AdjacentMapGraphState(tri) = test_adjacent_map_matches_graph(tri)
function Base.summary(state::AdjacentMapGraphState)
    if test_state(state)
        return "The adjacent map corrcetly matches the graph."
    else
        if state.is_k
            return "The vertex $(state.bad_vertex) is adjacent to the edge $(state.bad_edge) but either $(state.bad_edge_2[1]) or $(state.bad_edge_2[2]) fail to be in $(state.bad_vertex)'s neighbourhood."
        else
            return "The vertex $(state.bad_vertex), which defines the edge $(state.bad_edge) with adjacent vertex $(state.bad_vertex), is inconsistent with the graph which fails to include at least one of $(state.bad_edge_2[1]) and $(state.bad_edge_2[2]) in $(state.bad_vertex)'s neighbourhood."
        end
    end
end
function test_adjacent_map_matches_graph(tri)
    for (e, k) in get_adjacent(get_adjacent(tri))
        i, j = edge_vertices(e)
        flag1 = all(∈(get_neighbours(tri, i)), (j, k))
        !flag1 && return AdjacentMapGraphState(flag1, Int(i), Int.((i, j)), Int.((j, k)), false, Int(k))
        flag2 = all(∈(get_neighbours(tri, j)), (k, i))
        !flag2 && return AdjacentMapGraphState(flag2, Int(j), Int.((i, j)), Int.((k, i)), false, Int(k))
        flag3 = all(∈(get_neighbours(tri, k)), (i, j))
        !flag3 && return AdjacentMapGraphState(flag3, Int(k), Int.((i, j)), Int.((i, j)), true, Int(k))
    end
    return AdjacentMapGraphState(true, ∅, (∅, ∅), (∅, ∅), true, ∅)
end

struct GraphTrianglesState <: AbstractTriangulationState
    flag::Bool
    bad_triangle::NTuple{3, Int}
    bad_edge::NTuple{2, Int}
    bad_vertex::Int
end
GraphTrianglesState(tri) = test_graph_matches_triangles(tri)
function Base.summary(state::GraphTrianglesState)
    if test_state(state)
        return "The graph correctly matches the triangle set."
    else
        if state.bad_vertex == ∅
            return "The graph is inconsistent with the triangle set, as one of the vertices of $(state.bad_triangle) is not a vertex in the graph."
        else
            return "The graph is inconsistent with the triangle set. The triangle $(state.bad_triangle) is in the triangle set but either $(state.bad_edge[1]) or $(state.bad_edge[2]) are not in $(state.bad_vertex)'s neighbourhood."
        end
    end
end
function test_graph_matches_triangles(tri)
    for T in each_triangle(tri)
        try
            i, j, k = triangle_vertices(T)
            flag1 = all(∈(get_neighbours(tri, i)), (j, k))
            !flag1 && return GraphTrianglesState(flag1, Int.(triangle_vertices(T)), Int.((j, k)), Int(i))
            flag2 = all(∈(get_neighbours(tri, j)), (k, i))
            !flag2 && return GraphTrianglesState(flag2, Int.(triangle_vertices(T)), Int.((k, i)), Int(j))
            flag3 = all(∈(get_neighbours(tri, k)), (i, j))
            !flag3 && return GraphTrianglesState(flag3, Int.(triangle_vertices(T)), Int.((i, j)), Int(k))
        catch e
            e isa KeyError && return GraphTrianglesState(false, Int.(triangle_vertices(T)), (∅, ∅), ∅)
            rethrow(e)
        end
    end
    return GraphTrianglesState(true, (∅, ∅, ∅), (∅, ∅), ∅)
end

struct SegmentState <: AbstractTriangulationState
    flag::Bool
    segment::NTuple{2, Int}
    type::Int
end
SegmentState(tri) = test_segments(tri)
function Base.summary(state::SegmentState)
    if test_state(state)
        return "The set of segments is consistent with the triangulation."
    else
        if state.type == 0
            return "The segment $(state.segment) is not in the triangulation."
        elseif state.type == 1
            return "The segment $(state.segment) is in the triangulation, but it was not in the interior_segment field and not in the boundary edges."
        elseif state.type == 2
            return "The interior segments are not a subset of all segments."
        end
    end
end
function test_segments(tri)
    for e in each_segment(tri)
        flag = edge_exists(tri, e) || edge_exists(tri, reverse_edge(e))
        !flag && return SegmentState(flag, Int.(edge_vertices(e)), 0)
        flag = contains_edge(e, get_interior_segments(tri)) || contains_edge(reverse_edge(e), get_interior_segments(tri))
        if !flag
            flag = contains_boundary_edge(tri, e) || contains_boundary_edge(tri, reverse_edge(e))
            !flag && return SegmentState(flag, Int.(edge_vertices(e)), 1)
        end
    end
    for e in get_interior_segments(tri)
        flag = edge_exists(tri, e) || edge_exists(tri, reverse_edge(e))
        !flag && return SegmentState(flag, Int.(edge_vertices(e)), 0)
    end
    int_segs = Set{NTuple{2, Int}}()
    all_segs = Set{NTuple{2, Int}}()
    for e in get_interior_segments(tri)
        u, v = edge_vertices(e)
        push!(int_segs, minmax(u, v))
    end
    for e in each_segment(tri)
        u, v = edge_vertices(e)
        push!(all_segs, minmax(u, v))
    end
    flag = int_segs ⊆ all_segs
    !flag && return SegmentState(flag, (∅, ∅), 2)
    return SegmentState(true, (∅, ∅), 0)
end

struct BoundaryEdgeMapBoundaryNodesState <: AbstractTriangulationState
    flag::Bool
    bad_edge::NTuple{2, Int}
    bad_edge_2::NTuple{2, Int}
    bad_pos::Tuple
end
BoundaryEdgeMapBoundaryNodesState(tri) = test_boundary_edge_map_matches_boundary_nodes(tri)
function Base.summary(state::BoundaryEdgeMapBoundaryNodesState)
    if test_state(state)
        return "The boundary edge map is consistent with the boundary nodes."
    else
        return "The boundary edge map is inconsistent with the boundary nodes. The edge $(state.bad_edge) maps to $(state.bad_pos) in the boundary edge map but $(state.bad_pos) corresponds to the edge $(state.bad_edge_2)."
    end
end
function test_boundary_edge_map_matches_boundary_nodes(tri::Triangulation)
    boundary_edge_map = get_boundary_edge_map(tri)
    for (edge, pos) in boundary_edge_map
        u = get_boundary_nodes(get_boundary_nodes(tri, pos[1]), pos[2])
        v = get_boundary_nodes(get_boundary_nodes(tri, pos[1]), pos[2] + 1)
        flag = u == initial(edge)
        !flag && return BoundaryEdgeMapBoundaryNodesState(flag, Int.(edge_vertices(edge)), Int.((u, v)), pos)
    end
    return BoundaryEdgeMapBoundaryNodesState(true, (∅, ∅), (∅, ∅), ())
end

struct BoundaryNodesBoundaryEdgeMapState <: AbstractTriangulationState
    flag::Bool
    bad_edge::NTuple{2, Int}
    bad_pos::Tuple
    bad_pos_2::Tuple
end
BoundaryNodesBoundaryEdgeMapState(tri) = test_boundary_nodes_matches_boundary_edge_map(tri)
function Base.summary(state::BoundaryNodesBoundaryEdgeMapState)
    if test_state(state)
        return "The boundary nodes are consistent with the boundary edge map."
    else
        return "The boundary nodes are inconsistent with the boundary edge map. The edge $(state.bad_edge) maps to $(state.bad_pos_2) but is at $(state.bad_pos)."
    end
end
function test_boundary_nodes_matches_boundary_edge_map(tri::Triangulation)
    boundary_nodes = get_boundary_nodes(tri)
    E = edge_type(tri)
    if has_multiple_curves(tri)
        nc = num_curves(tri)
        for i in 1:nc
            curve_nodes = get_boundary_nodes(tri, i)
            ns = num_sections(curve_nodes)
            for j in 1:ns
                segment_nodes = get_boundary_nodes(curve_nodes, j)
                ne = num_boundary_edges(segment_nodes)
                for k in 1:ne
                    left_node = get_boundary_nodes(segment_nodes, k)
                    right_node = get_boundary_nodes(segment_nodes, k + 1)
                    edge = construct_edge(E, left_node, right_node)
                    pos = get_boundary_edge_map(tri, edge)
                    flag = pos == ((i, j), k)
                    !flag && return BoundaryNodesBoundaryEdgeMapState(flag, Int.(edge_vertices(edge)), pos, ((i, j), k))
                end
            end
        end
    elseif has_multiple_sections(tri)
        ns = num_sections(tri)
        for i in 1:ns
            segment_nodes = get_boundary_nodes(tri, i)
            ne = num_boundary_edges(segment_nodes)
            for j in 1:ne
                left_node = get_boundary_nodes(segment_nodes, j)
                right_node = get_boundary_nodes(segment_nodes, j + 1)
                edge = construct_edge(E, left_node, right_node)
                pos = get_boundary_edge_map(tri, edge)
                flag = pos == (i, j)
                !flag && return BoundaryNodesBoundaryEdgeMapState(flag, Int.(edge_vertices(edge)), pos, (i, j))
            end
        end
    else
        ne = num_boundary_edges(boundary_nodes)
        for i in 1:ne
            left_node = get_boundary_nodes(boundary_nodes, i)
            right_node = get_boundary_nodes(boundary_nodes, i + 1)
            edge = construct_edge(E, left_node, right_node)
            pos = get_boundary_edge_map(tri, edge)
            flag = pos == (get_boundary_nodes(tri), i)
            !flag && return BoundaryNodesBoundaryEdgeMapState(flag, Int.(edge_vertices(edge)), pos, (get_boundary_nodes(tri), i))
        end
    end
    return BoundaryNodesBoundaryEdgeMapState(true, (∅, ∅), (), ())
end

abstract type AbstractTriangulationIteratorState <: AbstractTriangulationState end
struct UniqueIteratorOutputState <: AbstractTriangulationIteratorState
    flag::Bool
    iterator_name::Symbol
    primitive_name::Symbol
end
function Base.summary(state::UniqueIteratorOutputState)
    if test_state(state)
        return "The $(state.iterator_name) iterator has no duplicate $(state.primitive_name)s."
    else
        return "The $(state.iterator_name) iterator has duplicate $(state.primitive_name)s."
    end
end
struct IteratorLengthState <: AbstractTriangulationIteratorState
    flag::Bool
    iterator_name::Symbol
    primitive_name::Symbol
    iterator_function_name::Symbol
    iterator_length::Int
    correct_length::Int
end
function Base.summary(state::IteratorLengthState)
    if test_state(state)
        return "The $(state.iterator_name) iterator correctly has the same length as the number of $(state.primitive_name)s in the triangulation."
    else
        return "The $(state.iterator_name) iterator has a length $(state.iterator_length), different to the number of $(state.primitive_name)s ($(state.correct_length)) returned from $(state.iterator_function_name)."
    end
end
struct IteratorOutputState <: AbstractTriangulationIteratorState
    flag::Bool
    iterator_name::Symbol
    primitive_name::Symbol
    iterator_function_name::Symbol
end
function Base.summary(state::IteratorOutputState)
    if test_state(state)
        return "The $(state.iterator_name) iterator returns the same set of $(state.primitive_name)s as returned from $(state.iterator_function_name)."
    else
        return "The $(state.iterator_name) iterator returns a different set of $(state.primitive_name)s as returned from $(state.iterator_function_name)."
    end
end
struct IteratorState <: AbstractTriangulationIteratorState
    flag::Bool
    unique_state::UniqueIteratorOutputState
    length_state::IteratorLengthState
    output_state::IteratorOutputState
    iterator_name::Symbol
end
function Base.summary(state::IteratorState)
    test_state(state) && return "The $(state.iterator_name) iterator is correctly constructed."
    !test_state(state.unique_state) && return summary(state.unique_state)
    !test_state(state.length_state) && return summary(state.length_state)
    return summary(state.output_state) # !test_state(state.output_state)
end
function IteratorState(
        unique_flag, length_flag, output_flag,
        iterator_length, correct_length,
        iterator_name, primitive_name, iterator_function_name,
    )
    flag = unique_flag && length_flag && output_flag
    return IteratorState(
        flag,
        UniqueIteratorOutputState(unique_flag, iterator_name, primitive_name),
        IteratorLengthState(length_flag, iterator_name, primitive_name, iterator_function_name, iterator_length, correct_length),
        IteratorOutputState(output_flag, iterator_name, primitive_name, iterator_function_name),
        iterator_name,
    )
end
function VertexIteratorState(length_flag, output_flag, iterator_length, correct_length)
    return IteratorState(true, length_flag, output_flag, iterator_length, correct_length, :vertex, :vertice, :each_vertex)
end
function SolidVertexIteratorState(length_flag, output_flag, iterator_length, correct_length)
    return IteratorState(true, length_flag, output_flag, iterator_length, correct_length, :solid_vertex, :vertice, :each_solid_vertex)
end
function GhostVertexIteratorState(length_flag, output_flag, iterator_length, correct_length)
    return IteratorState(true, length_flag, output_flag, iterator_length, correct_length, :ghost_vertex, :vertice, :each_ghost_vertex)
end
function TriangleIteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length)
    return IteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length, :triangle, :triangle, :each_triangle)
end
function SolidTriangleIteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length)
    return IteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length, :solid_triangle, :triangle, :each_solid_triangle)
end
function GhostTriangleIteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length)
    return IteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length, :ghost_triangle, :triangle, :each_ghost_triangle)
end
function EdgeIteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length)
    return IteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length, :edge, :edge, :each_edge)
end
function SolidEdgeIteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length)
    return IteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length, :solid_edge, :edge, :each_solid_edge)
end
function GhostEdgeIteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length)
    return IteratorState(unique_flag, length_flag, output_flag, iterator_length, correct_length, :ghost_edge, :edge, :each_ghost_edge)
end
function test_iterators(tri::Triangulation)
    I = integer_type(tri)
    T = NTuple{3, I}
    E = NTuple{2, I}
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
        i, j, k = triangle_vertices(T)
        push!(all_triangles, (i, j, k))
        if is_ghost_triangle(i, j, k)
            push!(ghost_triangles, (i, j, k))
        else
            push!(solid_triangles, (i, j, k))
        end
        for (u, v) in triangle_edges(i, j, k)
            push!(all_edges, (u, v))
            if is_ghost_edge(u, v)
                push!(ghost_edges, (u, v))
            else
                push!(solid_edges, (u, v))
            end
        end
        for k in (i, j, k)
            push!(all_vertices, k)
            if is_ghost_vertex(k)
                push!(ghost_vertices, k)
            else
                push!(solid_vertices, k)
            end
        end
    end
    unique_triangle_flag = allunique(all_triangles)
    unique_edge_flag = allunique(all_edges)
    unique_solid_triangle_flag = allunique(solid_triangles)
    unique_solid_edge_flag = allunique(solid_edges)
    unique_ghost_triangle_flag = allunique(ghost_triangles)
    unique_ghost_edge_flag = allunique(ghost_edges)
    for (i, e) in enumerate(all_edges)
        u, v = edge_vertices(e)
        all_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(solid_edges)
        u, v = edge_vertices(e)
        solid_edges[i] = (min(u, v), max(u, v))
    end
    for (i, e) in enumerate(ghost_edges)
        u, v = edge_vertices(e)
        ghost_edges[i] = (min(u, v), max(u, v))
    end
    unique!(all_edges)
    unique!(solid_edges)
    unique!(ghost_edges)
    triangle_correct_length, triangle_iterator_length = length(all_triangles), length(each_triangle(tri))
    edge_correct_length, edge_iterator_length = length(all_edges), length(each_edge(tri))
    solid_triangle_correct_length, solid_triangle_iterator_length = length(solid_triangles), length(each_solid_triangle(tri))
    solid_edge_correct_length, solid_edge_iterator_length = length(solid_edges), length(each_solid_edge(tri))
    ghost_triangle_correct_length, ghost_triangle_iterator_length = length(ghost_triangles), length(each_ghost_triangle(tri))
    ghost_edge_correct_length, ghost_edge_iterator_length = length(ghost_edges), length(each_ghost_edge(tri))
    vertices_correct_length, vertices_iterator_length = length(all_vertices), length(each_vertex(tri))
    solid_vertices_correct_length, solid_vertices_iterator_length = length(solid_vertices), length(each_solid_vertex(tri))
    ghost_vertices_correct_length, ghost_vertices_iterator_length = length(ghost_vertices), length(each_ghost_vertex(tri))
    triangle_length_flag = triangle_correct_length == triangle_iterator_length
    edge_length_flag = edge_correct_length == edge_iterator_length
    solid_triangle_length_flag = solid_triangle_correct_length == solid_triangle_iterator_length
    solid_edge_length_flag = solid_edge_correct_length == solid_edge_iterator_length
    ghost_triangle_length_flag = ghost_triangle_correct_length == ghost_triangle_iterator_length
    ghost_edge_length_flag = ghost_edge_correct_length == ghost_edge_iterator_length
    vertices_length_flag = vertices_correct_length == vertices_iterator_length
    solid_vertices_length_flag = solid_vertices_correct_length == solid_vertices_iterator_length
    ghost_vertices_length_flag = ghost_vertices_correct_length == ghost_vertices_iterator_length
    all_vertices = collect(all_vertices)
    solid_vertices = collect(solid_vertices)
    ghost_vertices = collect(ghost_vertices)
    sort!(all_vertices)
    sort!(solid_vertices)
    sort!(ghost_vertices)
    triangle_output_flag = compare_triangle_collections(all_triangles, each_triangle(tri))
    solid_triangle_output_flag = compare_triangle_collections(solid_triangles, each_solid_triangle(tri))
    ghost_triangle_output_flag = compare_triangle_collections(ghost_triangles, each_ghost_triangle(tri))
    edge_output_flag = compare_edge_vectors(all_edges, each_edge(tri))
    solid_edge_output_flag = compare_edge_vectors(solid_edges, each_solid_edge(tri))
    ghost_edge_output_flag = compare_edge_vectors(ghost_edges, each_ghost_edge(tri))
    vertex_output_flag = all_vertices == sort(collect(each_vertex(tri)))
    solid_vertex_output_flag = solid_vertices == sort(collect(each_solid_vertex(tri)))
    ghost_vertex_output_flag = ghost_vertices == sort(collect(each_ghost_vertex(tri)))
    return (
        VertexIteratorState(vertices_length_flag, vertex_output_flag, vertices_iterator_length, vertices_correct_length),
        SolidVertexIteratorState(solid_vertices_length_flag, solid_vertex_output_flag, solid_vertices_iterator_length, solid_vertices_correct_length),
        GhostVertexIteratorState(ghost_vertices_length_flag, ghost_vertex_output_flag, ghost_vertices_iterator_length, ghost_vertices_correct_length),
        TriangleIteratorState(unique_triangle_flag, triangle_length_flag, triangle_output_flag, triangle_iterator_length, triangle_correct_length),
        SolidTriangleIteratorState(unique_solid_triangle_flag, solid_triangle_length_flag, solid_triangle_output_flag, solid_triangle_iterator_length, solid_triangle_correct_length),
        GhostTriangleIteratorState(unique_ghost_triangle_flag, ghost_triangle_length_flag, ghost_triangle_output_flag, ghost_triangle_iterator_length, ghost_triangle_correct_length),
        EdgeIteratorState(unique_edge_flag, edge_length_flag, edge_output_flag, edge_iterator_length, edge_correct_length),
        SolidEdgeIteratorState(unique_solid_edge_flag, solid_edge_length_flag, solid_edge_output_flag, solid_edge_iterator_length, solid_edge_correct_length),
        GhostEdgeIteratorState(unique_ghost_edge_flag, ghost_edge_length_flag, ghost_edge_output_flag, ghost_edge_iterator_length, ghost_edge_correct_length),
    )
end

struct DuplicateSegmentsState <: AbstractTriangulationState
    flag::Bool
    interior::Bool
    bad_edge::NTuple{2, Int}
end
DuplicateSegmentsState(tri) = test_no_duplicate_segments(tri)
function Base.summary(state::DuplicateSegmentsState)
    if test_state(state)
        return "The segments contain no duplicates."
    else
        if state.interior
            return "The edge $(state.bad_edge) is duplicated in the interior segments field."
        else
            return "The edge $(state.bad_edge) is duplicated in the all segments field."
        end
    end
end
function test_no_duplicate_segments(tri)
    interior_segments = get_interior_segments(tri)
    all_segments = get_all_segments(tri)
    for e in each_edge(interior_segments)
        flag = !contains_edge(reverse_edge(e), interior_segments)
        !flag && return DuplicateSegmentsState(flag, true, Int.(edge_vertices(e)))
    end
    for e in each_edge(all_segments)
        flag = !contains_edge(reverse_edge(e), all_segments)
        !flag && return DuplicateSegmentsState(flag, false, Int.(edge_vertices(e)))
    end
    return DuplicateSegmentsState(true, true, (∅, ∅))
end

struct TriangulationState <: AbstractTriangulationState
    adjacent2vertex_adjacent_state::Adjacent2VertexMapAdjacentMapState
    adjacent2vertex_state::Adjacent2VertexMapState
    adjacent_adjacent2vertex_state::AdjacentMapAdjacent2VertexMapState
    adjacent_graph_state::AdjacentMapGraphState
    adjacent_state::AdjacentMapState
    boundary_edge_map_boundary_nodes_state::BoundaryEdgeMapBoundaryNodesState
    boundary_nodes_boundary_edge_map_state::BoundaryNodesBoundaryEdgeMapState
    delaunay_criterion_state::DelaunayCriterionState
    duplicate_segments_state::DuplicateSegmentsState
    edges_have_two_incident_triangles_state::EdgesHaveTwoIncidentTrianglesState
    graph_adjacent_state::GraphAdjacentMapState
    graph_state::GraphState
    graph_triangles_state::GraphTrianglesState
    segment_state::SegmentState
    triangle_orientation_state::TriangleOrientationState
    vertex_iterator_state::IteratorState
    solid_vertex_iterator_state::IteratorState
    ghost_vertex_iterator_state::IteratorState
    triangle_iterator_state::IteratorState
    solid_triangle_iterator_state::IteratorState
    ghost_triangle_iterator_state::IteratorState
    edge_iterator_state::IteratorState
    solid_edge_iterator_state::IteratorState
    ghost_edge_iterator_state::IteratorState
end

function TriangulationState(tri::Triangulation)
    has_ghosts = has_ghost_triangles(tri)
    delete_ghost_triangles!(tri)
    add_ghost_triangles!(tri)
    clear_empty_features!(tri)
    state = TriangulationState(
        Adjacent2VertexMapAdjacentMapState(tri),
        Adjacent2VertexMapState(tri),
        AdjacentMapAdjacent2VertexMapState(tri),
        AdjacentMapGraphState(tri),
        AdjacentMapState(tri),
        BoundaryEdgeMapBoundaryNodesState(tri),
        BoundaryNodesBoundaryEdgeMapState(tri),
        DelaunayCriterionState(tri),
        DuplicateSegmentsState(tri),
        EdgesHaveTwoIncidentTrianglesState(tri),
        GraphAdjacentMapState(tri),
        GraphState(tri),
        GraphTrianglesState(tri),
        SegmentState(tri),
        TriangleOrientationState(tri),
        test_iterators(tri)...,
    )
    !has_ghosts && delete_ghost_triangles!(tri)
    return state
end

function test_state(triangulation_state::TriangulationState)
    for f in fieldnames(typeof(triangulation_state))
        state = getfield(triangulation_state, f)
        !test_state(state) && return false
    end
    return true
end

function Base.show(io::IO, triangulation_state::TriangulationState)
    is_valid = test_state(triangulation_state)
    if is_valid
        # println(io, "The triangulation is correctly constructed.")
    else
        for f in fieldnames(typeof(triangulation_state))
            state = getfield(triangulation_state, f)
            !test_state(state) && println(io, state)
        end
    end
    return
end

"""
    validate_triangulation(tri::Triangulation; print_result=true) -> Bool 

Tests if `tri` is a valid `Triangulation`. Returns `true` if so, 
and `false` otherwise. If `print_result=true` and `tri` is not a 
valid triangulation, all the issues with `tri` will be printed.
"""
function validate_triangulation(tri::Triangulation; print_result = true)
    state = TriangulationState(tri)
    print_result && !test_state(state) && println(state)
    return test_state(state)
end
