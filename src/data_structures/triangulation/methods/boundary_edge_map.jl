"""
    get_boundary_edge_map(tri::Triangulation, ij) 
    get_boundary_edge_map(tri::Triangulation, i, j)

Returns the value from the key `(i, j)` in the boundary edge map of `tri`. The returned value is a `Tuple` 
`(position, index)` so that `boundary_nodes = get_boundary_nodes(tri, position)` are the boundary nodes associated 
with the section that `(i, j)` resides on, and `i = get_boundary_nodes(boundary_nodes, index)` and 
`j = get_boundary_nodes(boundary_nodes, index + 1)`.
"""
get_boundary_edge_map(tri::Triangulation, ij) = get_boundary_edge_map(tri)[ij]
get_boundary_edge_map(tri::Triangulation, i, j) = get_boundary_edge_map(tri, construct_edge(edge_type(tri), i, j))


"""
    each_boundary_edge(tri::Triangulation) -> KeySet

Returns an iterator over the boundary edges of `tri`, in no specific order.
"""
each_boundary_edge(tri::Triangulation) = keys(get_boundary_edge_map(tri))


"""
    split_boundary_edge_map!(boundary_edge_map, boundary_nodes, pos) 

After splitting an edge starting at `pos` on the boundary, updates the `boundary_edge_map` to reflect the new
boundary edges. See [`split_boundary_edge!`](@ref).
"""
function split_boundary_edge_map!(boundary_edge_map::Dict{E, T}, boundary_nodes, pos, i, j) where {E, T}
    e = construct_edge(E, i, j)
    delete!(boundary_edge_map, e)
    nodes = get_boundary_nodes(boundary_nodes, pos[1])
    ne = num_boundary_edges(nodes)
    for k in pos[2]:ne
        u = get_boundary_nodes(nodes, k)
        v = get_boundary_nodes(nodes, k + 1)
        e = construct_edge(E, u, v)
        boundary_edge_map[e] = (pos[1], k)
    end
    return boundary_edge_map
end