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
    nn = num_neighbours(tri, point)
    @show nn
    is_bnd, boundary_index = is_boundary_node(tri, point)
    check_orientation = Val(is_bnd) # Val(any(is_boundary_index, S)) - see the comment in add_point_convex_triangulation!
    if is_bnd
        left_bnd = get_left_boundary_node(tri, point, boundary_index)
        right_bnd = get_right_boundary_node(tri, point, boundary_index)
    else
        I = integer_type(tri)
        left_bnd = I(DefaultAdjacentValue)
        right_bnd = I(DefaultAdjacentValue)
    end
    ## As a zeroth step, we take a special case where the point being deleted is a boundary pint with only two neighbours (i.e. nn == 3, of them being BoundaryIndex) 
    if is_bnd && nn == 3
        _delete_point_2!(tri, point, boundary_index, left_bnd, right_bnd)
        return nothing
    end
    ## The first step is to get the polygon that surrounds the point. 
    S = get_surrounding_polygon(tri, point)
    ## Before we do any deleting, we need to handle the case where the point 
    ## we are deleting is a boundary point that is collinear with its neighbouring 
    ## boundary nodes. This is for reasoning similar to what we do in the Bowyer-Watson 
    ## algorithm, where this point, when deleted, does not get the boundary information 
    ## correctly updated.
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
    ## The first step is to repair the relationships for those triangles on the other side of the cavity 
    fix_edges_after_deletion!(tri, S)
    ## The next step, if we have deleted a boundary node (as determined via `check_ghost_edges::Val{T}`), is to 
    ## ensure we have correctly repaired the ghost information for the updated boundary. The need for this step 
    ## comes from noting that the "cavity" that we delete includes the ghost vertex, which complicates things.
    fix_ghost_edges_after_deletion!(tri, S, boundary_index, check_orientation, left_bnd, right_bnd)
    return nothing
end

function fix_edges_after_deletion!(tri, S)
    push!(S, S[begin])
    for i in firstindex(S):(lastindex(S)-1)
        u = S[i]
        v = S[i+1]
        w = get_adjacent(tri, v, u)
        if edge_exists(w)
            add_adjacent!(tri, u, w, v)
            add_adjacent!(tri, w, v, u)
            add_neighbour!(tri, u, w)
            #= 
            This can't be done here because we first need to fix what we can of the ghost triagnles in fix_ghost_edges_after_deletion! (called after 
            this function). So, we do it at the end of fix_ghost_edges_after_deletion! instead.
            if !edge_exists(tri, u, v)
                add_triangle!(tri, u, v, boundary_index; protect_boundary=true)
            end
            =#
        end
    end
end

# check_ghost_edges comes from check_orientation above, which is true only if we are deleting a boundary node, 
# which is the only occasion where we should check for ghost edges at all. 
# We need the boundary index argument to know what boundary index to use for the ghost triangles we add.
function fix_ghost_edges_after_deletion!(tri, S, boundary_index, check_ghost_edges::Val{T}, left_bnd, right_bnd) where {T}
    if T
        ## The idea is as follows: When we delete the boundary point, the two neighbouring points on the boundary 
        ## should still remain on the boundary. There are two cases to consider:
        ##  1. These two boundary nodes are now connected.
        ##  2. Extra points have been added onto the boundary.
        ## We can identify this case by seeing if any edges no longer exist for nodes connected 
        ## to these neighbouring boundary nodes. This is done the two following lines.
        bad_right_edges = filter(u -> !edge_exists(tri, right_bnd, u), get_neighbours(tri, right_bnd))
        bad_left_edges = filter(u -> !edge_exists(tri, u, left_bnd), get_neighbours(tri, left_bnd))
        ## Now, we can check if we have case (1) by seeing if any of the sets above are non-empty. The 
        ## edges that fail to give an adjacency must be new boundary edges. When checking this, we need 
        ## to directly delete the non-existent triangle, since sometimes it can exist.
        ## There is one case, though, where this triangle is allowed to exist - when the edge 
        ## itself is not a boundary edge.
        if any(≠(0), (length(bad_left_edges), length(bad_right_edges)))
            if !is_boundary_edge(tri, right_bnd, left_bnd)
                ## If this is not a boundary edge, we might still have added a ghost triangle (right_bnd, left_bnd, boundary_index), so this still needs to be 
                ## removed. When do we do this removal, though, we will accidentally remove a correct adjacency. So, let's store it and correct it after deletion.
                w = get_adjacent(tri, right_bnd, left_bnd)
                delete_triangle!(tri, right_bnd, left_bnd, boundary_index; protect_boundary=true)
                add_adjacent!(tri, right_bnd, left_bnd, w)
            else
                delete_triangle!(tri, right_bnd, left_bnd, boundary_index; protect_boundary=true)
            end
        end
        ## Let us now add in the new ghost triangles. 
        if length(bad_right_edges) ≠ 0
            bad_right_edge = first(bad_right_edges)
            add_triangle!(tri, right_bnd, bad_right_edge, boundary_index; protect_boundary=true)
        end
        if length(bad_left_edges) ≠ 0
            bad_left_edge = first(bad_left_edges)
            add_triangle!(tri, bad_left_edge, left_bnd, boundary_index; protect_boundary=true)
        end

        ## Now, it is still possible that we have some remaining edges that have not been corrected. We do this here.
        for i in firstindex(S):(lastindex(S)-1)
            u = S[i]
            v = S[i+1]
            @show get_adjacent(tri, u, v), edge_exists(tri, v, u) && !edge_exists(tri, u, v), u, v
            if edge_exists(tri, v, u) && !edge_exists(tri, u, v)
                add_triangle!(tri, u, v, boundary_index; protect_boundary=true)
            end
        end
    end
end

# special case for deleting a boundary point with two neighbours
function _delete_point_2!(tri::Triangulation, point, boundary_index, left_bnd, right_bnd)
    delete_triangle!(tri, point, right_bnd, left_bnd; protect_boundary=false, update_ghost_edges=true)
    #add_triangle!(tri, boundary_index, left_bnd, right_bnd; protect_boundary=true)
    delete_adjacent2vertex!(tri, point)
    delete_vertex!(tri, point)
    return nothing
end