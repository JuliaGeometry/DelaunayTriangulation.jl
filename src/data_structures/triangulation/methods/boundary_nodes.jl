"""
    num_curves(tri::Triangulation) -> Integer

Returns the number of curves in `tri`.
"""
num_curves(tri::Triangulation) = num_curves(get_boundary_nodes(tri))

"""
    num_sections(tri::Triangulation) -> Integer

Assuming `tri` only has one curve, returns the number of sections in `tri`.
"""
num_sections(tri::Triangulation) = num_sections(get_boundary_nodes(tri))

"""
    get_boundary_nodes(tri, mnℓ...) 

Given a triangulation `tri`, returns the specified component of the boundary nodes.
There are several forms for the methods:

1. `get_boundary_nodes(tri, m)`: If `tri` has multiple curves, this returns the `m`th curve. If `tri` has multiple sections, this returns the `m`th section. Otherwise, this returns the `m`th boundary node.
2. `get_boundary_nodes(tri, m, n)`: If `tri` has multiple curves, this returns the `n`th section of the `m`th curve. Otherwise, if `tri` has multiple sections, this returns the `n`th boundary node of the `m`th section.
3. `get_boundary_nodes(tri, (m, n))`: This is equivalent to `get_boundary_nodes(tri, m, n)`.
4. `get_boundary_nodes(tri::A, ::A)`: This just returns `boundary_nodes`.  
"""
get_boundary_nodes(tri::Triangulation, mnℓ...) = get_boundary_nodes(get_boundary_nodes(tri), mnℓ...)
get_boundary_nodes(tri::Triangulation, m::Integer) = get_boundary_nodes(get_boundary_nodes(tri), m) # method ambiguity
get_boundary_nodes(tri::Triangulation, (m, n)::NTuple{2,Integer}) = get_boundary_nodes(get_boundary_nodes(tri), (m, n)) # method ambiguity
get_boundary_nodes(tri::Triangulation, m::Integer, n::Integer) = get_boundary_nodes(get_boundary_nodes(tri), m, n) # method ambiguity

"""
    get_right_boundary_node(tri::Triangulation, k, ghost_vertex) -> Vertex

Returns the boundary node to the right of the vertex `k` in `tri`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `k`: The boundary vertex. 
- `ghost_vertex`: The ghost vertex associated with the boundary section that `k` is on. 

# Outputs 
- `r`: The vertex right of `k` on the boundary. 
"""
function get_right_boundary_node(tri::Triangulation, k, ghost_vertex)
    range = get_ghost_vertex_range(tri, ghost_vertex)
    for index in range
        i = get_adjacent(tri, k, index)
        i ≠ ∅ && return i
    end
    return get_adjacent(tri, k, ghost_vertex)
end

"""
    get_left_boundary_node(tri::Triangulation, k, ghost_vertex) -> Vertex

Returns the boundary node to the left of the vertex `k` in `tri`.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `k`: The boundary vertex.
- `ghost_vertex`: The ghost vertex associated with the boundary section that `k` is on.

# Outputs
- `ℓ`: The vertex left of `k` on the boundary.
"""
function get_left_boundary_node(tri::Triangulation, k, ghost_vertex)
    range = get_ghost_vertex_range(tri, ghost_vertex)
    for index in range
        i = get_adjacent(tri, index, k)
        i ≠ ∅ && return i
    end
    return get_adjacent(tri, ghost_vertex, k)
end

"""
    insert_boundary_node!(tri::Triangulation, pos, node)

Inserts a boundary node into `tri` at the specified position.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `pos`: The position to insert the node at, given as a `Tuple` so that `insert_boundary_node!(tri, pos, node)` is the same as `insert!(get_boundary_nodes(tri, pos[1]), pos[2], node)`.
- `node`: The node to insert.

# Outputs 
There are no outputs, but the boundary nodes of `tri` are updated in-place.
"""
insert_boundary_node!(tri::Triangulation, pos, node) = insert_boundary_node!(get_boundary_nodes(tri), pos, node)

"""
    delete_boundary_node!(tri::Triangulation, pos)

Deletes the boundary node at the specified position `pos` in `tri`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `pos`: The position to delete the node at, given as a `Tuple` so that `delete_boundary_node!(tri, pos)` is the same as `deleteat!(get_boundary_nodes(tri, pos[1]), pos[2])`.

# Outputs
There are no outputs, but the boundary nodes of `tri` are updated in-place.
"""
delete_boundary_node!(tri::Triangulation, pos) = delete_boundary_node!(get_boundary_nodes(tri), pos)

"""
    split_boundary_edge!(tri::Triangulation, ij, node)
    split_boundary_edge!(tri::Triangulation, i, j, node)

Splits the boundary edge `edge` in `tri` at the edge `(i, j)`.

See also [`merge_boundary_edge!`](@ref).
"""
function split_boundary_edge!(tri::Triangulation, edge, node)
    i, j = edge_vertices(edge)
    return split_boundary_edge!(tri, i, j, node)
end
function split_boundary_edge!(tri::Triangulation, i, j, node)
    pos = get_boundary_edge_map(tri, i, j)
    new_pos = (pos[1], pos[2] + 1)
    bnn = get_boundary_edge_map(tri)
    insert_boundary_node!(tri, new_pos, node)
    E = edge_type(tri)
    boundary_nodes = get_boundary_nodes(tri)
    split_boundary_edge_map!(bnn, boundary_nodes, pos, i, j)
    segments = get_all_segments(tri)
    delete_unoriented_edge!(segments, construct_edge(E, i, j))
    !contains_edge(construct_edge(E, node, i), segments) && add_edge!(segments, construct_edge(E, i, node))
    !contains_edge(construct_edge(E, j, node), segments) && add_edge!(segments, construct_edge(E, node, j))
    # Sometimes, the user might accidentally also put an interior segment in that forms part of the boundary. Let's remove it. 
    interior_segments = get_interior_segments(tri)
    if contains_unoriented_edge(construct_edge(E, i, j), interior_segments)
        delete_unoriented_edge!(interior_segments, construct_edge(E, i, j))
    end
    return tri
end

"""
    merge_boundary_edge!(tri::Triangulation, ij, node)
    merge_boundary_edge!(tri::Triangulation, i, j, node)

Merges the edges `(i, node)` and `(node, j)` into a single edge `(i, j)`, i.e. does the inverse operation to [`split_boundary_edge!`](@ref).
"""
function merge_boundary_edge!(tri::Triangulation, edge, node)
    i, j = edge_vertices(edge)
    return merge_boundary_edge!(tri, i, j, node)
end
function merge_boundary_edge!(tri::Triangulation, i, j, node)
    pos = get_boundary_edge_map(tri, i, node)
    node_pos = (pos[1], pos[2] + 1)
    bnn = get_boundary_edge_map(tri)
    delete_boundary_node!(tri, node_pos)
    E = edge_type(tri)
    delete!(bnn, construct_edge(E, i, node))
    delete!(bnn, construct_edge(E, node, j))
    e = construct_edge(E, i, j)
    bnn[e] = pos
    nodes = get_boundary_nodes(tri, pos[1])
    ne = num_boundary_edges(nodes)
    for k in pos[2]:ne
        u = get_boundary_nodes(nodes, k)
        v = get_boundary_nodes(nodes, k + 1)
        e = construct_edge(E, u, v)
        bnn[e] = (pos[1], k)
    end
    segments = get_all_segments(tri)
    delete_unoriented_edge!(segments, construct_edge(E, i, node))
    delete_unoriented_edge!(segments, construct_edge(E, node, j))
    !contains_edge(construct_edge(E, j, i), segments) && add_edge!(segments, construct_edge(E, i, j))
    return tri
end

"""
    split_boundary_edge_at_collinear_segments!(tri::Triangulation, collinear_segments)

Splits a boundary edge into pieces defined by `collinear_segments`. In particular, if `r = collinear_segments` and

    u = initial(r[1])
    v = terminal(r[end]),

then the boundary edge is `(u, v)` and the edges are split so that all segments in `collinear_segments` appear instead.
"""
function split_boundary_edge_at_collinear_segments!(tri::Triangulation, collinear_segments)
    # To understand this function, consider for example a segment (u, v) such that the points r₁, …, r₅ 
    # are collinear with u and v. We want to split the boundary edge (u, v) into 5 edges, each of which
    # is collinear with the points r₁, …, r₅. We assume that u < r₁ < ⋯ < r₅ < v. With this setup, we can first split 
    # (u, v) at r₁, then split (r₁, v) at r₂, and so on: 
    #   u -------------------------------- v
    #   u -- r₁ -------------------------- v -> split(u, v, r₁) 
    #   u -- r₁ -- r₂ -------------------- v -> split(r₁, v, r₂)
    #   u -- r₁ -- r₂ -- r₃ -------------- v -> split(r₂, v, r₃)
    #   u -- r₁ -- r₂ -- r₃ -- r₄ -------- v -> split(r₃, v, r₄)
    #   u -- r₁ -- r₂ -- r₃ -- r₄ -- r₅ -- v -> split(r₄, v, r₅)
    v = terminal(last(collinear_segments))
    for k in (firstindex(collinear_segments)):(lastindex(collinear_segments)-1)
        segment = collinear_segments[k]
        rₖ₋₁, rₖ = edge_vertices(segment)
        split_boundary_edge!(tri, rₖ₋₁, v, rₖ)
        if is_curve_bounded(tri)
            enricher = get_boundary_enricher(tri)
            split_boundary_edge!(enricher, rₖ₋₁, v, rₖ, Val(false))
        end
    end
    return tri
end

"""
    get_all_boundary_nodes(tri::Triangulation) -> Set{Vertex}

Returns the set of all boundary vertices in `tri`, in no specific order.
"""
function get_all_boundary_nodes(tri::Triangulation)
    boundary_sections = (values ∘ get_ghost_vertex_map)(tri)
    I = integer_type(tri)
    all_nodes = Set{I}()
    if has_boundary_nodes(tri)
        for section_index in boundary_sections
            bn_nodes = get_boundary_nodes(tri, section_index)
            nedges = num_boundary_edges(bn_nodes)
            for node_idx in 1:(nedges+1)
                vᵢ = get_boundary_nodes(bn_nodes, node_idx)
                push!(all_nodes, vᵢ)
            end
        end
    end
    return all_nodes
end

"""
    is_positively_oriented(tri::Triangulation, curve_index) -> Bool

Tests if the `curve_index`th curve in `tri` is positively oriented, returning `true` if so and `false` otherwise.
"""
function is_positively_oriented(tri::Triangulation, curve_index::Integer)
    points = get_points(tri)
    if has_boundary_nodes(tri)
        boundary_nodes = get_boundary_nodes(tri)
        if has_multiple_curves(boundary_nodes)
            curve_boundary_nodes = get_boundary_nodes(boundary_nodes, curve_index)
        else
            curve_boundary_nodes = boundary_nodes
        end
    else
        curve_boundary_nodes = get_convex_hull_vertices(tri)
    end
    area = polygon_features(points, curve_boundary_nodes)[1]
    return area > 0.0
end

"""
    contains_boundary_edge(tri::Triangulation, ij) -> Bool 
    contains_boundary_edge(tri::Triangulation, i, j) -> Bool

Returns `true` if the boundary edge `(i, j)` is in `tri`, and `false` otherwise. Orientation matters here.
"""
contains_boundary_edge(tri::Triangulation, e) = e ∈ keys(get_boundary_edge_map(tri))
function contains_boundary_edge(tri::Triangulation, i, j)
    E = edge_type(tri)
    e = construct_edge(E, i, j)
    return contains_boundary_edge(tri, e)
end