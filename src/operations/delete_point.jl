"""
    delete_point!(tri::Triangulation, point; rng::AbstractRNG = Random.default_rng())
       
Deletes the point `point` into the triangulation.

This function will not update the convex hull - if you need it to 
be corrected, you could use e.g. [`convex_hull!`](@ref). You should 
also be careful with using this if you have deleted ghost triangles.

Points are deleted by using Chew's algorithm for triangulating the cavity 
left over when a point is deleted. No specialisation is made for low-degree 
vertices, but maybe one day we could implement e.g.
https://doi.org/10.1016/j.comgeo.2010.10.001.

# Arguments 
- `tri::Triangulation`: The triangulation. 
- `point`: The index of the point in `get_points(tri)` to delete.

# Keyword Argumens 
- `rng::AbstractRNG = Random.default_rng()`: RNG to use.

# Outputs 
There are no outputs, but `tri` is updated in-place.
"""
function delete_point!(tri::Triangulation, point; rng::AbstractRNG=Random.default_rng())
    ## The first step is to get the polygon that surrounds the point. 
    S = get_surrounding_polygon(tri, point)
    is_bnd, boundary_index = is_boundary_node(tri, point)
    check_orientation = Val(is_bnd) # Val(any(is_boundary_index, S)) - see the comment in add_point_convex_triangulation!
    ## Now delete all the neighbouring triangles, and remove other information 
    neighbouring_edges = get_adjacent2vertex(tri, point)
    for uv in each_edge(neighbouring_edges)
        u = initial(uv)
        v = terminal(uv)
        delete_triangle!(tri, point, u, v; protect_boundary=true)
    end
    ## Triangulate the cavity
    triangulate_convex!(tri, S; rng, update_ghost_edges=false, check_orientation)
    delete_adjacent2vertex!(tri, point)
    delete_vertex!(tri, point)
    ## 
    # There's an issue with non-convex cavities. I'm not 100% sure what the fix is,
    # but it should just naturally work inside the convex triangulation code itself.
    # There's an associated issue for this, but basically the problem is that triangles that 
    # do itself exist away from the cavity, made by concave indents, are temporarily replaced.
    # These "temporary" replacements are not fixed later. So, here we will just force this fix 
    # by looping over edges and looking at the triangle pointing away from the cavity.
    # I imagine this issue is related to the problem in a constrained Delaunay triangulation, 
    # wherein a similar treatment does actually require a separate data structure:
    #=
    "Be forewarned that CavityCDT, unlike ConvexDT, cannot use the same triangulation data
    structure as the triangulation in which the segment s is being inserted, because CavityCDT sometimes temporarily creates triangles that conflict with those outside the cavity.
    CavityCDT requires the use of a separate, empty triangulation data structure, and the final
    triangles must subsequently be added to the main triangulation."
    =#
    # (p. 76, Delaunay Mesh Triangulation, Cheng, Dey, and Shewchuk 2013),
    # but I should be able to leverage something here to avoid this entirely.
    # This is probably nice anyway and less bug-prone since we can also fix boundary edges here.
    fix_edges_after_deletion!(tri, S, boundary_index, check_orientation)
    return nothing
end

# check_ghost_edges comes from check_orientation above, which is true only if we are deleting a boundary node, 
# which is the only occasion where we should check for ghost edges at all. 
# We need the boundary index argument to know what boundary index to use for the ghost triangles we add.
function fix_edges_after_deletion!(tri, S, boundary_index, check_ghost_edges::Val{T}) where {T}
    push!(S, S[begin])
    for i in firstindex(S):(lastindex(S)-1)
        @show get_adjacent(tri, 12, 5)
        u = S[i]
        v = S[i+1]
        w = get_adjacent(tri, v, u)
        if edge_exists(w)
            add_adjacent!(tri, u, w, v)
            add_adjacent!(tri, w, v, u)
            add_neighbour!(tri, u, w)
            # Here, we might need to correct for any missing ghost edges. We detect this by seeing if an edge has only one 
            # incident triangle, meaning one of the adjacent vertices is DefaultAdjacentValue
            if T && !edge_exists(tri, u, v)
                @show u, v, boundary_index
                add_triangle!(tri, u, v, boundary_index; protect_boundary=true)
                # Now that we've added this triangle, it is possible that the next point in the 
                # sequence is also missing a ghost triangle. This next point will be some 
                # point in S, but it's not necessarily the *next* point in S. The best 
                # way I know to handle this currently is to take the intersection of S 
                # with the neighbourhood of the point v (not u since we're taking the end 
                # of the edge) with S, and then to check if any of these edges 
                # have no adjacent vertex.
                W = intersect(get_neighbours(tri, v), S)
                for w in W
                    if !edge_exists(tri, v, w)
                        add_triangle!(tri, v, w, boundary_index; protect_boundary=true)
                    end
                end
            end
        end
    end
    @show get_adjacent(tri, 12, 5)
end