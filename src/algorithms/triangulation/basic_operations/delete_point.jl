"""
    get_surrounding_polygon(tri::Triangulation, u; skip_ghost_vertices=false) -> Vector 

Returns the counter-clockwise sequence of neighbours of `u` in `tri`.

# Arguments 
- `tri::Triangulation`: [`Triangulation`](@ref).
- `u`: The vertex.

# Keyword Arguments
- `skip_ghost_vertices=false`: Whether to skip ghost vertices in the returned polygon. 

# Outputs
- `S`: The surrounding polygon. This will not be circular, meaning `S[begin] ≠ S[end]`.
   In case `u` is an exterior ghost vertex, the returned polygon is a clockwise list of vertices for 
   the associated boundary curve. If you do not have ghost triangles and you try to get the surrounding polygon
   of a boundary vertex, then this function may return an invalid polygon.
"""
function get_surrounding_polygon(tri::Triangulation, u; skip_ghost_vertices=false)
    return copy(get_surrounding_polygon!(tri, u; skip_ghost_vertices))
end

function get_surrounding_polygon!(tri::Triangulation, u; skip_ghost_vertices=false)
    cache = get_cache(tri)
    S = get_surrounding_polygon(cache)
    neighbouring_vertices = get_neighbours(tri, u)
    v = first(neighbouring_vertices)
    if !is_ghost_vertex(u)
        k = num_neighbours(tri, u)
        nbnd = count(is_ghost_vertex, neighbouring_vertices)
        if nbnd > 0
            # If we are on a boundary that has multiple ghost vertices, 
            # k will not be counted correctly.
            k = k - nbnd + 1 # plus 1 because we do still want one ghost vertex 
        end
    else
        # Need to be careful when there are multiple ghost vertices - get_neighbours would only 
        # return the section for a given ghost vertices, so we need to check all. The nodes do not uniquely 
        # appear in a single segment, so we cannot just keep adding to k above. To keep things 
        # a bit simple, we just take a union. This could be improved.
        ghost_vertex_range = get_ghost_vertex_range(tri, u)
        _neighbouring_vertices = copy(neighbouring_vertices)
        for vertex in ghost_vertex_range
            union!(_neighbouring_vertices, get_neighbours(tri, vertex))
        end
        k = length(_neighbouring_vertices)
    end
    resize!(S, k)
    S[1] = v
    for i in 2:k
        v = get_adjacent(tri, u, v)
        S[i] = v
    end
    skip_ghost_vertices && filter!(!is_ghost_vertex, S)
    return S
end

struct InvalidVertexDeletionError{I} <: Exception
    v::I
    is_segment::Bool
    is_ghost::Bool
end
function Base.showerror(io::IO, e::InvalidVertexDeletionError)
    err = "Tried to delete the vertex $(e.v) which "
    if e.is_segment
        err *= "adjoins a segment of the triangulation."
    elseif e.is_ghost
        err *= "is a ghost vertex."
    else
        err *= "forms part of the boundary of the triangulation."
    end
    err *= " This is not allowed - only interior vertices not adjoining a segment may be deleted."
    return print(io, err)
end

"""
    check_delete_point_args(tri::Triangulation, vertex, S) -> Bool

Checks that the vertex `vertex` can be deleted from the triangulation `tri`. Returns `true` if so, 
and throws an `InvalidVertexDeletionError` otherwise. This will occur if:

1. `vertex` is a boundary node of `tri`.
2. `vertex` is a ghost vertex of `tri`.
3. `vertex` adjoins a segment of `tri`.
"""
function check_delete_point_args(tri::Triangulation, vertex, S)
    is_boundary_node(tri, vertex)[1] && throw(InvalidVertexDeletionError(vertex, false, false))
    is_ghost_vertex(vertex) && throw(InvalidVertexDeletionError(vertex, false, true))
    any(v -> contains_segment(tri, vertex, v), S) && throw(InvalidVertexDeletionError(vertex, true, false))
    return true
end

"""
    fix_edges_after_deletion!(tri::Triangulation, S) 

Ensures that the edges in `S` surrounding a deleted vertex of `tri` are correctly updated.
"""
function fix_edges_after_deletion!(tri, S)
    for i in firstindex(S):(lastindex(S)-1)
        u = S[i]
        v = S[i+1]
        w = get_adjacent(tri, v, u)
        if edge_exists(w)
            add_adjacent!(tri, u, w, v)
            add_adjacent!(tri, w, v, u)
            add_neighbour!(tri, u, w)
        end
    end
    return tri
end

"""
    delete_point!(tri::Triangulation, vertex; kwargs...)

Deletes the `vertex` of `tri`, retriangulating the cavity formed by the surrounding polygon of `vertex` using 
[`triangulate_convex`](@ref).

It is not possible to delete vertices that are on the boundary, are ghost vertices, or adjoin a segment of `tri`. 
See also [`check_delete_point_args`](@ref).

!!! warn "Point deletion"

    This function will not actually delete the corresponding coordinates from `get_points(tri)`, nor will it remove 
    the associated weight from `get_weights(tri)`.

# Arguments 
- `tri::Triangulation`: [`Triangulation`](@ref).
- `vertex`: The vertex to delete.

# Keyword Arguments
- `store_event_history=Val(false)`: Whether to store the event history of the triangulation from deleting the point. 
- `event_history=nothing`: The event history of the triangulation from deleting the point. Only updated if `store_event_history` is true, in which case it needs to be an [`InsertionEventHistory`](@ref) object.
- `rng::AbstractRNG=Random.default_rng()`: The random number generator to use for the triangulation.
"""
function delete_point!(tri::Triangulation, vertex;
    store_event_history=Val(false),
    event_history=nothing,
    rng::AbstractRNG=Random.default_rng())
    cache = get_cache(tri)
    empty!(cache)
    convex_tri = get_triangulation(cache)
    S = get_surrounding_polygon!(tri, vertex)
    check_delete_point_args(tri, vertex, S)
    neighbouring_edges = get_adjacent2vertex(tri, vertex)
    trit = triangle_type(tri)
    for uv in each_edge(neighbouring_edges) # note that we are mutating this iterator during iteration
        u, v = edge_vertices(uv)
        delete_triangle!(tri, vertex, u, v; protect_boundary=true)
        is_true(store_event_history) && delete_triangle!(event_history, construct_triangle(trit, vertex, u, v))
    end
    triangulate_convex!(convex_tri, S; rng) # fill in the cavity
    for T in each_solid_triangle(convex_tri)
        u, v, w = triangle_vertices(T)
        add_triangle!(tri, u, v, w; protect_boundary=true, update_ghost_edges=false)
        is_true(store_event_history) && add_triangle!(event_history, T)
    end
    delete_adjacent2vertex!(tri, vertex)
    delete_vertex!(tri, vertex)
    push!(S, S[begin])
    fix_edges_after_deletion!(tri, S)
    return tri
end
