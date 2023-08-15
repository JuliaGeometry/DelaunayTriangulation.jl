"""
    has_multiple_curves(tri::Triangulation)

Returns `has_multiple_curves(get_boundary_nodes(tri))`, testing if `tri`'s boundary is 
comprised of multiple curves.
"""
@inline has_multiple_curves(tri::Triangulation) = has_multiple_curves(get_boundary_nodes(tri))

"""
    has_multiple_segments(tri::Triangulation)

Returns `has_multiple_segments(get_boundary_nodes(tri))`, testing if `tri`'s boundary is 
comprised of multiple segments.
"""
@inline has_multiple_segments(tri::Triangulation) = has_multiple_segments(get_boundary_nodes(tri))

"""
    num_curves(tri::Triangulation)

Returns `num_curves(get_boundary_nodes(tri))`, the number of curves in `tri`'s boundary.
"""
@inline num_curves(tri::Triangulation) = num_curves(get_boundary_nodes(tri))

"""
    num_segments(tri::Triangulation)

Returns `num_segments(get_boundary_nodes(tri))`, the number of segments in `tri`'s boundary.
"""
@inline num_segments(tri::Triangulation) = num_segments(get_boundary_nodes(tri))

"""
    get_boundary_nodes(tri::Triangulation, mnℓ...)

Returns `get_boundary_nodes(get_boundary_nodes(tri), mnℓ...)`.
"""
@inline get_boundary_nodes(tri::Triangulation, mnℓ...) = get_boundary_nodes(get_boundary_nodes(tri), mnℓ...)

"""
    map_boundary_index(tri::Triangulation, i)

Returns `map_boundary_index(get_boundary_map(tri), i)`, mapping the boundary index `i` to
the corresponding index in the boundary nodes.
"""
@inline map_boundary_index(tri::Triangulation, i) = map_boundary_index(get_boundary_map(tri), i)

"""
    get_curve_index(tri::Triangulation, i)

Returns `get_curve_index(get_boundary_map(tri), i)`, the curve corresponding to the boundary
index `i`.
"""
@inline get_curve_index(tri::Triangulation, i) = get_curve_index(get_boundary_map(tri), i)

"""
    get_segment_index(tri::Triangulation, i)

Returns `get_segment_index(get_boundary_map(tri), i)`, the segment corresponding to the boundary 
index `i`.
"""
@inline get_segment_index(tri::Triangulation, i) = get_segment_index(get_boundary_map(tri), i)

"""
    num_outer_boundary_segments(tri::Triangulation)

Returns `num_outer_boundary_segments(get_boundary_nodes(tri))`, returns the number of boundary segments on the 
outermost boundary.
"""
@inline num_outer_boundary_segments(tri::Triangulation) = num_outer_boundary_segments(get_boundary_nodes(tri))

"""
    get_right_boundary_node(tri::Triangulation, k, boundary_index)

Returns `get_right_boundary_node(get_adjacent(tri), k, boundary_index, get_boundary_index_ranges(tri), Val(has_multiple_segments(tri)))`, 
the node to the right of the boundary node `k`. It is assumed that `k` is on the part of the boundary of `tri` with index `boundary_index`.
"""
@inline get_right_boundary_node(tri::Triangulation, k, boundary_index) = get_right_boundary_node(get_adjacent(tri), k, boundary_index, get_boundary_index_ranges(tri), Val(has_multiple_segments(tri)))

"""
    get_left_boundary_node(tri::Triangulation, k, boundary_index)

Returns `get_left_boundary_node(get_adjacent(tri), k, boundary_index, get_boundary_index_ranges(tri), Val(has_multiple_segments(tri)))`,
the node to the left of the boundary node `k`. It is assumed that `k` is on the part of the boundary of `tri` with index `boundary_index`.
"""
@inline get_left_boundary_node(tri::Triangulation, k, boundary_index) = get_left_boundary_node(get_adjacent(tri), k, boundary_index, get_boundary_index_ranges(tri), Val(has_multiple_segments(tri)))

"""
    get_boundary_index_range(tri::Triangulation, i)

Returns `map_boundary_index(get_boundary_index_ranges(tri), i)`, the list of boundary indices belonging to the curve corresponding 
to the boundary index `i`.
"""
@inline get_boundary_index_range(tri::Triangulation, i) = map_boundary_index(get_boundary_index_ranges(tri), i)

"""
    get_boundary_edge_map(tri::Triangulation, ij)
    get_boundary_edge_map(tri::Triangulation, i, j)

Returns `get_boundary_edge_map(get_boundary_nodes(tri), ij)` or `get_boundary_edge_map(get_boundary_nodes(tri), construct_edge(edge_type(tri), i, j))`,
returning a `Tuple` `pos` such that `get_boundary_nodes(get_boundary_nodes(tri, pos[1]), pos[2])` is the edge `ij = (i, j)`.
"""
@inline get_boundary_edge_map(tri::Triangulation, ij) = map_boundary_index(get_boundary_edge_map(tri), ij)
@inline get_boundary_edge_map(tri::Triangulation, i, j) = get_boundary_edge_map(tri, construct_edge(edge_type(tri), i, j))

"""
    insert_boundary_node!(tri::Triangulation, pos, node)

Calls `insert_boundary_node!(get_boundary_nodes(tri), pos, node)`, splitting the boundary node at `pos` at `node`.
"""
@inline insert_boundary_node!(tri::Triangulation, pos, node) = insert_boundary_node!(get_boundary_nodes(tri), pos, node)

"""
    delete_boundary_node!(tri::Triangulation, pos)

Calls `delete_boundary_node!(get_boundary_nodes(tri), pos)`, deleting the boundary node at `pos`.
"""
@inline delete_boundary_node!(tri::Triangulation, pos) = delete_boundary_node!(get_boundary_nodes(tri), pos)

"""
    split_boundary_edge!(tri::Triangulation, edge, node)
    split_boundary_edge!(tri::Triangulation, i, j, node)

Given a triangulation `tri`, an `edge = (i, j)`, and a node on `edge`, this function splits the boundary edge into two edges  
`edge = (i, node)` and `(node, j)`, updating `tri.boundary_nodes` and `tri.boundary_edge_map` accordingly.

See also [`split_boundary_edge_at_collinear_segments!`](@ref).
"""
@inline function split_boundary_edge!(tri::Triangulation, edge, node)
    i, j = edge_indices(edge)
    split_boundary_edge!(tri, i, j, node)
    return nothing
end
@inline function split_boundary_edge!(tri::Triangulation, i, j, node)
    pos = get_boundary_edge_map(tri, i, j)
    new_pos = (pos[1], pos[2] + 1)
    bnn = get_boundary_edge_map(tri)
    insert_boundary_node!(tri, new_pos, node)
    E = edge_type(tri)
    delete!(bnn, construct_edge(E, i, j))
    # When we add in a new node, we still have to modify the whole set of nodes to the right, shifting their index by 1.
    # I wonder if we should use a better data structure this, but alas.
    nodes = get_boundary_nodes(tri, pos[1])
    ne = num_boundary_edges(nodes)
    for k in pos[2]:ne
        u = get_boundary_nodes(nodes, k)
        v = get_boundary_nodes(nodes, k + 1)
        e = construct_edge(E, u, v)
        bnn[e] = (pos[1], k)
    end
    return nothing
end

"""
    merge_boundary_edge!(tri::Triangulation, edge)
    merge_boundary_edge!(tri::Triangulation, i, j)

This is the inverse of [`split_boundary_edge!`](@ref). Given a triangulation `tri` and an edge `edge = (i, j)`,
assumed to be a boundary node separated by a single `node` between them, deletes `node` and merges the edge,
updating `tri.boundary_nodes` and `tri.boundary_edge_map` accordingly.
"""
function merge_boundary_edge!(tri::Triangulation, edge, node)
    i, j = edge_indices(edge)
    merge_boundary_edge!(tri, i, j, node)
    return nothing
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
end

"""
    split_boundary_edge_at_collinear_segments!(tri::Triangulation, collinear_segments)

Given a triangulation `tri` and a list of collinear segments `collinear_segments`, assumed to represent a 
boundary edge `(initial(first(collinear_segments)), terminal(last(collinear_segments)))`, splits 
the edge at the `collinear_segments` and updates `tri.boundary_nodes` and `tri.boundary_edge_map` accordingly.

See also [`split_boundary_edge!`](@ref).
"""
@inline function split_boundary_edge_at_collinear_segments!(tri::Triangulation, collinear_segments)
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
        rₖ₋₁ = initial(segment)
        rₖ = terminal(segment)
        split_boundary_edge!(tri, rₖ₋₁, v, rₖ)
    end
    return nothing
end

"""
    contains_boundary_edge(tri::Triangulation, e)
    contains_boundary_edge(tri::Triangulation, i, j)

Tests if the triangulation `tri` has the constrained boundary edge `e = (i, j)`,
returning `e ∈ keys(get_boundary_edge_map(tri))`.
"""
@inline contains_boundary_edge(tri::Triangulation, e) = e ∈ keys(get_boundary_edge_map(tri))
@inline function contains_boundary_edge(tri::Triangulation, i, j)
    E = edge_type(tri)
    e = construct_edge(E, i, j)
    return contains_boundary_edge(tri, e)
end

"""
    merge_constrained_edges(bn_map, boundary_nodes, constrained_edges::Es)
    merge_constrained_edges(tri::Triangulation, bn_map=get_boundary_map(tri))

Merges the boundary edges, defined via `boundary_nodes` or `get_boundary_nodes(tri)`,
and the `constrained_edges` into a single collection. `bn_map` is used to iterate 
over all the boundary edges.
"""
function merge_constrained_edges(bn_map, boundary_nodes, constrained_edges::Es) where {Es}
    all_constrained = initialise_edges(Es)
    E = edge_type(Es)
    for segment_index in values(bn_map)
        bn_nodes = get_boundary_nodes(boundary_nodes, segment_index)
        nedges = num_boundary_edges(bn_nodes)
        for edge_idx in 1:nedges
            vᵢ = get_boundary_nodes(bn_nodes, edge_idx)
            vᵢ₊₁ = get_boundary_nodes(bn_nodes, edge_idx + 1)
            e = construct_edge(E, vᵢ, vᵢ₊₁)
            add_edge!(all_constrained, e)
        end
    end
    for e in each_edge(constrained_edges)
        add_edge!(all_constrained, e)
    end
    return all_constrained
end
merge_constrained_edges(tri::Triangulation, bn_map=get_boundary_map(tri)) = merge_constrained_edges(bn_map, get_boundary_nodes(tri), get_constrained_edges(tri))

"""
    get_all_boundary_nodes(tri::Triangulation)

Returns a `Set` of all the boundary nodes in the triangulation `tri`.
"""
function get_all_boundary_nodes(tri::Triangulation)
    bn_map = get_boundary_map(tri)
    all_nodes = Set{integer_type(tri)}()
    if has_boundary_nodes(tri)
        for segment_index in values(bn_map)
            bn_nodes = get_boundary_nodes(tri, segment_index)
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
    all_boundary_indices(tri::Triangulation)

Returns all the boundary indices in the triangulation `tri`.
"""
all_boundary_indices(tri::Triangulation) = keys(get_boundary_index_ranges(tri))

"""
    is_positively_oriented(tri::Triangulation, curve_index)

Tests if the curve with index `curve_index` in the triangulation `tri` is positively oriented.
"""
function is_positively_oriented(tri::Triangulation, curve_index)
    points = get_points(tri)
    if has_boundary_nodes(tri)
        boundary_nodes = get_boundary_nodes(tri)
        if has_multiple_curves(boundary_nodes)
            curve_boundary_nodes = get_boundary_nodes(boundary_nodes, curve_index)
        else
            curve_boundary_nodes = boundary_nodes
        end
    else
        curve_boundary_nodes = get_convex_hull_indices(tri)
    end
    area = polygon_features(points, curve_boundary_nodes)[1]
    return area > 0.0
end