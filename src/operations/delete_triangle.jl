"""
    delete_triangle!(tri::Triangulation, T)
    delete_triangle!(tri::Ts, u::Integer, v::Integer, w::Integer) where {Ts<:Triangulation}

Given a triangle `T`, or a triple of integers `u, v, w`, deletes the triangle from 
the [`Triangulation`](@ref) type `tri`, updating the [`Adjacent`](@ref), 
[`Adjacent2Vertex`](@ref), [`Graph`](@ref), and triangles fields.

No values are returned.
"""
function delete_triangle!(tri::Ts, u::Integer, v::Integer,
                          w::Integer) where {Ts<:Triangulation}
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    triangles = get_triangles(tri)
    delete_triangle!(adj, u, v, w)
    delete_triangle!(adj2v, u, v, w)
    delete_triangle!(graph, u, v, w)
    delete_triangle!(triangles, u, v, w)
    vu_bnd = is_boundary_edge(v, u, adj)
    uw_bnd = is_boundary_edge(u, w, adj)
    wv_bnd = is_boundary_edge(w, v, adj)
    vu_exists = edge_exists(v, u, adj)
    uw_exists = edge_exists(u, w, adj)
    wv_exists = edge_exists(w, v, adj)
    protect_vu_neighbour = vu_exists && !vu_bnd
    protect_uw_neighbour = uw_exists && !uw_bnd
    protect_wv_neighbour = wv_exists && !wv_bnd
    !protect_vu_neighbour && delete_neighbour!(graph, u, v)
    !protect_uw_neighbour && delete_neighbour!(graph, w, u)
    !protect_wv_neighbour && delete_neighbour!(graph, v, w)
    return nothing
end
function delete_triangle!(tri::Triangulation, T)
    u, v, w = indices(T)
    delete_triangle!(tri, u, v, w)
    return nothing
end
