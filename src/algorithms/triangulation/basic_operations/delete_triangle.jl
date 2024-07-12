"""
    delete_triangle!(tri::Triangulation, u, v, w; protect_boundary=false, update_ghost_edges=false)
    delete_triangle!(tri::Triangulation, T; protect_boundary=false, update_ghost_edges=false)

Deletes the triangle `T = (u, v, w)` from `tri`. This won't delete the points from `tri`, it will just update the fields 
so that its non-existence in the triangulation is known.

# Arguments 
- `tri::Triangulation`: The triangulation to delete the triangle from.
- `u, v, w`: The vertices of the triangle to delete.

# Keyword Arguments 
- `protect_boundary=false`: If `true`, then the boundary edges will not be updated. Otherwise, `delete_boundary_edges_single!`, `delete_boundary_edges_double!`, or `delete_boundary_edges_triple!` will be called depending on the number of boundary edges in the triangle.
- `update_ghost_edges=false`: If `true`, then the ghost edges will be updated. Otherwise, the ghost edges will not be updated. Will only be used if `protect_boundary=false`.

# Outputs 
There are no outputs as `tri` is updated in-place.
"""
function delete_triangle!(
        tri::Ts, u::Integer, v::Integer, w::Integer;
        protect_boundary = false, update_ghost_edges = false,
    ) where {Ts <: Triangulation}
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    triangles = get_triangles(tri)
    delete_triangle!(adj, u, v, w)
    delete_triangle!(adj2v, u, v, w)
    # delete_triangle!(graph, u, v, w)
    delete_triangle!(triangles, u, v, w)
    vu_bnd = is_boundary_edge(tri, u, v)
    uw_bnd = is_boundary_edge(tri, w, u)
    wv_bnd = is_boundary_edge(tri, v, w)
    num_bnd_edges = !protect_boundary ? count((vu_bnd, uw_bnd, wv_bnd)) : 0 # If we are protecting the boundary, setting this count to zero will do that 
    vu_exists = edge_exists(tri, v, u) # Could also write this as vu_is_bnd || edge_exists(v, u, adj), but this way makes the code simpler.
    uw_exists = edge_exists(tri, u, w)
    wv_exists = edge_exists(tri, w, v)
    protect_vu_neighbour = vu_exists && !vu_bnd
    protect_uw_neighbour = uw_exists && !uw_bnd
    protect_wv_neighbour = wv_exists && !wv_bnd
    !protect_vu_neighbour && delete_neighbour!(graph, u, v)
    !protect_uw_neighbour && delete_neighbour!(graph, w, u)
    !protect_wv_neighbour && delete_neighbour!(graph, v, w)
    if num_bnd_edges == 1
        delete_boundary_edges_single!(u, v, w, vu_bnd, uw_bnd, wv_bnd, triangles, adj, adj2v, graph, update_ghost_edges)
    elseif num_bnd_edges == 2
        delete_boundary_edges_double!(u, v, w, vu_bnd, uw_bnd, wv_bnd, triangles, adj, adj2v, graph, update_ghost_edges)
    elseif (num_bnd_edges == 3 || num_triangles(tri) == 0) && !protect_boundary # If this is true, then we just have a single triangle that we have deleted, and thus num_bnd_edges should be 3
        delete_boundary_edges_triple!(u, v, w, triangles, adj, adj2v, graph, update_ghost_edges)
    end
    return nothing
end
function delete_triangle!(tri::Triangulation, T; protect_boundary = false, update_ghost_edges = false)
    u, v, w = triangle_vertices(T)
    delete_triangle!(tri, u, v, w; protect_boundary, update_ghost_edges)
    return nothing
end


function delete_boundary_edges_single!(u, v, w, vu_bnd, uw_bnd, wv_bnd, triangles, adj::Adjacent{I, E}, adj2v, graph, update_ghost_edges) where {I, E}
    g = I(𝒢)
    u, v, w = choose_uvw(vu_bnd, wv_bnd, uw_bnd, u, v, w)
    delete_adjacent!(adj, v, u)
    delete_adjacent2vertex!(adj2v, g, v, u)
    add_adjacent!(adj, v, w, g)
    add_adjacent!(adj, w, u, g)
    add_adjacent2vertex!(adj2v, g, v, w)
    add_adjacent2vertex!(adj2v, g, w, u)
    add_neighbour!(graph, g, w)
    if update_ghost_edges
        ## Delete vu∂, add vw∂ and wu∂. See some of the comments in add_boundary_edges_double! for more info.
        # Ghost edges
        add_adjacent!(adj, w, g, v)
        add_adjacent!(adj, g, v, w)
        add_adjacent!(adj, u, g, w)
        add_adjacent!(adj, g, w, u)
        delete_adjacent2vertex!(adj2v, v, u, g)
        delete_adjacent2vertex!(adj2v, u, g, v)
        add_adjacent2vertex!(adj2v, v, w, g)
        add_adjacent2vertex!(adj2v, w, g, v)
        add_adjacent2vertex!(adj2v, w, u, g)
        add_adjacent2vertex!(adj2v, u, g, w)
        # Ghost triangles
        add_triangle!(triangles, v, w, g)
        add_triangle!(triangles, w, u, g)
        delete_triangle!(triangles, v, u, g)
    end
    return nothing
end


function delete_boundary_edges_double!(u, v, w, vu_bnd, uw_bnd, wv_bnd, triangles, adj::Adjacent{I, E}, adj2v, graph, update_ghost_edges) where {I, E}
    g = I(𝒢)
    u, v, w = choose_uvw(!vu_bnd, !wv_bnd, !uw_bnd, u, v, w)
    delete_adjacent!(adj, u, w)
    delete_adjacent!(adj, w, v)
    delete_adjacent2vertex!(adj2v, g, u, w)
    delete_adjacent2vertex!(adj2v, g, w, v)
    add_adjacent!(adj, u, v, g)
    add_adjacent2vertex!(adj2v, g, u, v)
    delete_neighbour!(graph, g, w)
    if update_ghost_edges
        ## Delete uw∂ and wv∂ and add uv∂. See some of the comments in add_boundary_edges_single! for more info.
        # Ghost edges
        add_adjacent!(adj, v, g, u)
        add_adjacent!(adj, g, u, v)
        delete_adjacent!(adj, w, g)
        delete_adjacent!(adj, g, w)
        delete_adjacent2vertex!(adj2v, u, w, g)
        delete_adjacent2vertex!(adj2v, w, g, u)
        delete_adjacent2vertex!(adj2v, w, v, g)
        delete_adjacent2vertex!(adj2v, v, g, w)
        add_adjacent2vertex!(adj2v, u, v, g)
        add_adjacent2vertex!(adj2v, v, g, u)
        # Ghost triangles
        add_triangle!(triangles, u, v, g)
        delete_triangle!(triangles, u, w, g)
        delete_triangle!(triangles, w, v, g)
    end
    return nothing
end


function delete_boundary_edges_triple!(u, v, w, triangles, adj::Adjacent{I, E}, adj2v, graph, update_ghost_edges) where {I, E}
    g = I(𝒢)
    delete_adjacent!(adj, w, v)
    delete_adjacent!(adj, v, u)
    delete_adjacent!(adj, u, w)
    delete_adjacent2vertex!(adj2v, g, w, v)
    delete_adjacent2vertex!(adj2v, g, v, u)
    delete_adjacent2vertex!(adj2v, g, u, w)
    delete_neighbour!(graph, g, u, v, w)
    if update_ghost_edges
        ## Delete vu∂, wv∂, and uw∂. See some of the comments in add_boundary_edges_triple! for more info.
        # Ghost edges
        delete_adjacent!(adj, u, g)
        delete_adjacent!(adj, g, v)
        delete_adjacent!(adj, v, g)
        delete_adjacent!(adj, g, w)
        delete_adjacent!(adj, w, g)
        delete_adjacent!(adj, g, u)
        delete_adjacent2vertex!(adj2v, v, u, g)
        delete_adjacent2vertex!(adj2v, u, g, v)
        delete_adjacent2vertex!(adj2v, w, v, g)
        delete_adjacent2vertex!(adj2v, v, g, w)
        delete_adjacent2vertex!(adj2v, u, w, g)
        delete_adjacent2vertex!(adj2v, w, g, u)
        # Ghost triangles 
        delete_triangle!(triangles, v, u, g)
        delete_triangle!(triangles, w, v, g)
        delete_triangle!(triangles, u, w, g)
    end
    return nothing
end