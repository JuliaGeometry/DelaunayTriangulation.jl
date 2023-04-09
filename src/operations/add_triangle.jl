"""
    add_triangle!(tri::Triangulation, u, v, w; protect_boundary=false, update_ghost_edges=false)
    add_triangle!(tri::Triangulation, T; protect_boundary=false, update_ghost_edges=false)

Given a triangle `T = (u, v, w)`, adds the triangle into the triangulation `tri`. To update the ghost triangles directly, 
set `update_ghost_edges=true`. The other parts of the boundary information will be handled, though, unless you set `protect_boundary=true`.
"""
function add_triangle!(tri::Ts, u::Integer, v::Integer, w::Integer;
    protect_boundary=false, update_ghost_edges=false) where {Ts<:Triangulation}
    ## Add the necessary triangles
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    triangles = get_triangles(tri)
    # Need to check for boundaries before we replace the adjacencies with the new triangles
    uv_bnd = is_boundary_edge(u, v, adj)
    vw_bnd = is_boundary_edge(v, w, adj)
    wu_bnd = is_boundary_edge(w, u, adj)
    # Now we can add
    add_triangle!(adj, u, v, w)
    add_triangle!(adj2v, u, v, w)
    add_triangle!(graph, u, v, w)
    add_triangle!(triangles, u, v, w)
    ## Now consider the boundary 
    num_bnd_edges = protect_boundary ? 0 : count((uv_bnd, vw_bnd, wu_bnd))
    if num_bnd_edges == 1
        add_boundary_edges_single!(u, v, w, uv_bnd, vw_bnd, wu_bnd, triangles, adj, adj2v, graph, update_ghost_edges)
    elseif num_bnd_edges == 2
        add_boundary_edges_double!(u, v, w, uv_bnd, vw_bnd, wu_bnd, triangles, adj, adj2v, graph, update_ghost_edges)
    elseif (num_bnd_edges == 3 || num_triangles(tri) == 1) && !protect_boundary # If this is true, then we just have a single triangle, and thus num_bnd_edges should be 3
        add_boundary_edges_triple!(u, v, w, triangles, adj, adj2v, graph, update_ghost_edges)
    end
    return nothing
end
function add_triangle!(tri::Triangulation, T; protect_boundary=false, update_ghost_edges=false)
    u, v, w = indices(T)
    add_triangle!(tri, u, v, w; protect_boundary, update_ghost_edges)
    return nothing
end

function add_boundary_edges_single!(u, v, w, uv_bnd, vw_bnd, wu_bnd, triangles, adj::Adjacent{I,E}, adj2v, graph, update_ghost_edges) where {I,E}
    g = I(BoundaryIndex)
    # Here, we are adding two ghost triangles uwg and wvg, where g is the ghost vertex, coming from 
    # the two new boundary edges uw and wv. 
    u, v, w = choose_uvw(uv_bnd, vw_bnd, wu_bnd, u, v, w)
    add_adjacent!(adj, u, w, g)
    add_adjacent!(adj, w, v, g)
    add_adjacent2vertex!(adj2v, g, u, w)
    add_adjacent2vertex!(adj2v, g, w, v)
    delete_adjacent2vertex!(adj2v, g, u, v)
    add_neighbour!(graph, g, w) # Note that u and v are already in graph[g]
    if update_ghost_edges
        # Need to make sure the edges vg and gu do not accidentally get deleted. 
        add_adjacent!(adj, w, g, u)
        add_adjacent!(adj, g, u, w)
        add_adjacent!(adj, v, g, w)
        add_adjacent!(adj, g, w, v)
        add_adjacent2vertex!(adj2v, u, w, g)
        add_adjacent2vertex!(adj2v, w, g, u)
        add_adjacent2vertex!(adj2v, w, v, g)
        add_adjacent2vertex!(adj2v, v, g, w)
        delete_adjacent2vertex!(adj2v, u, v, g)
        delete_adjacent2vertex!(adj2v, v, g, u)
        add_triangle!(triangles, u, w, g)
        add_triangle!(triangles, w, v, g)
        delete_triangle!(triangles, u, v, g)
    end
end

function add_boundary_edges_double!(u, v, w, uv_bnd, vw_bnd, wu_bnd, triangles, adj::Adjacent{I,E}, adj2v, graph, update_ghost_edges) where {I,E}
    g = I(BoundaryIndex)
    # Here, we are only adding a single ghost triangle vug, where g is the ghost vertex, 
    # coming from the new boundary edge vu. 
    u, v, w = choose_uvw(!uv_bnd, !vw_bnd, !wu_bnd, u, v, w)
    add_adjacent!(adj, v, u, g)
    add_adjacent2vertex!(adj2v, g, v, u)
    delete_adjacent2vertex!(adj2v, g, v, w)
    delete_adjacent2vertex!(adj2v, g, w, u)
    delete_neighbour!(graph, g, w) # u and v are still in graph[g] 
    if update_ghost_edges
        # Make sure we delete the two original ghost triangles vwg and wug, and also make sure 
        # the edges gv and ug do not get accidentally deleted. 
        add_adjacent!(adj, u, g, v)
        add_adjacent!(adj, g, v, u)
        delete_adjacent!(adj, w, g)
        delete_adjacent!(adj, g, w)
        add_adjacent2vertex!(adj2v, v, u, g)
        add_adjacent2vertex!(adj2v, u, g, v)
        delete_adjacent2vertex!(adj2v, w, u, g)
        delete_adjacent2vertex!(adj2v, u, g, w)
        delete_adjacent2vertex!(adj2v, v, w, g)
        delete_adjacent2vertex!(adj2v, w, g, v)
        add_triangle!(triangles, v, u, g)
        delete_triangle!(triangles, v, w, g)
        delete_triangle!(triangles, w, u, g)
    end
end

function add_boundary_edges_triple!(u, v, w, triangles, adj::Adjacent{I,E}, adj2v, graph, update_ghost_edges) where {I,E}
    g = I(BoundaryIndex)
    # Here, we are adding three ghost triangles uwg, wvg, and vug, where g is the ghost vertex.
    add_adjacent!(adj, v, u, g)
    add_adjacent!(adj, u, w, g)
    add_adjacent!(adj, w, v, g)
    add_adjacent2vertex!(adj2v, g, v, u)
    add_adjacent2vertex!(adj2v, g, u, w)
    add_adjacent2vertex!(adj2v, g, w, v)
    add_neighbour!(graph, g, u, v, w)
    if update_ghost_edges
        # Add in the triangles themselves 
        add_adjacent!(adj, u, g, v)
        add_adjacent!(adj, g, v, u)
        add_adjacent!(adj, v, g, w)
        add_adjacent!(adj, g, w, v)
        add_adjacent!(adj, w, g, u)
        add_adjacent!(adj, g, u, w)
        add_adjacent2vertex!(adj2v, v, u, g)
        add_adjacent2vertex!(adj2v, u, g, v)
        add_adjacent2vertex!(adj2v, w, v, g)
        add_adjacent2vertex!(adj2v, v, g, w)
        add_adjacent2vertex!(adj2v, u, w, g)
        add_adjacent2vertex!(adj2v, w, g, u)
        add_triangle!(triangles, v, u, g)
        add_triangle!(triangles, u, w, g)
        add_triangle!(triangles, w, v, g)
    end
end