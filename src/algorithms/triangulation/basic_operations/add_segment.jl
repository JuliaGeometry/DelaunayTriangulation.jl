"""
    optimise_edge_order(tri::Triangulation, segment) -> Edge 

Optimises the orientation of `segment` for inserting it into the triangulation. 

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `segment`: The segment to arrange. 

# Outputs 
- `e`: If `segment` is a boundary edge, then `e = segment`, Otherwise, `e = sort_edge_by_degree(tri, segment)` so that `initial(e)` has the smaller degree of the two vertices.
"""
optimise_edge_order(tri::Triangulation, segment) = contains_boundary_edge(tri, segment) ? segment : sort_edge_by_degree(tri, segment) # We want to preserve the orientation of the boundary edges as it makes it easier to split them later if we run into collinearities

"""
    fix_edge_order_after_rotation!(tri::Triangulation, segment, e) 

Fixes the edge order in `get_interior_segments(tri)` after `segment` was rotated by [`optimise_edge_order`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `segment`: The segment that was arranged.
- `e`: The arranged segment from [`optimise_edge_order`](@ref).

# Outputs 
There is no output, but `tri` will be updated so that `e` is in `get_interior_segments(tri)` instead of `segment`.
"""
function fix_edge_order_after_rotation!(tri::Triangulation, segment, e)
    interior_segments = get_interior_segments(tri)
    if initial(e) ≠ initial(segment) && contains_edge(segment, interior_segments) # If we switch the edge, let's also change it in constrained_edges
        delete_edge!(interior_segments, segment)
        add_edge!(interior_segments, e)
    end
    return tri
end

"""
    add_segment_to_list!(tri::Triangulation, e)

Adds `e` to `get_interior_segments(tri)` and `get_all_segments(tri)` if it, or `reverse_edge(e)`, is not already in the sets.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `e`: The edge to add.

# Outputs
There is no output, but `tri` will be updated so that `e` is in `get_interior_segments(tri)` and `get_all_segments(tri)`.
"""
function add_segment_to_list!(tri::Triangulation, e)
    interior_segments = get_interior_segments(tri)
    all_segments = get_all_segments(tri)
    if !contains_unoriented_edge(e, interior_segments) && !contains_boundary_edge(tri, e)
        add_edge!(interior_segments, e)
    end
    if !contains_unoriented_edge(e, all_segments)
        add_edge!(all_segments, e)
    end
    return tri
end

"""
    add_segment!(tri::Triangulation, segment; predicates::AbstractPredicateKernel=Adaptive(), rng::Random.AbstractRNG=Random.default_rng())
    add_segment!(tri::Triangulation, i, j; predicates::AbstractPredicateKernel=Adaptive(), rng::Random.AbstractRNG=Random.default_rng())

Adds `segment = (i, j)` to `tri`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `segment`: The segment to add. The second method uses `(i, j)` to represent the segment instead.

# Keyword Arguments
- `predicates::AbstractPredicateKernel=Adaptive()`: Method to use for computing predicates. Can be one of [`Fast`](@ref), [`Exact`](@ref), and [`Adaptive`](@ref). See the documentation for a further discussion of these methods.
- `rng::AbstractRNG=Random.default_rng()`: The RNG object.
# Outputs 
There is no output, but `tri` will be updated so that it now contains `segment`.
"""
function add_segment!(tri::Triangulation, segment; predicates::AbstractPredicateKernel=Adaptive(), rng::Random.AbstractRNG=Random.default_rng())
    e = optimise_edge_order(tri, segment)
    fix_edge_order_after_rotation!(tri, segment, e)
    add_segment_to_list!(tri, e)
    unoriented_edge_exists(tri, e) && return tri # Don't need to add edges that already appear
    intersecting_triangles, collinear_segments, left_cavity, right_cavity = locate_intersecting_triangles(tri, e, !contains_boundary_edge(tri, segment), rng, predicates)
    flag = process_collinear_segments!(tri, e, collinear_segments; predicates, rng)
    flag && return tri
    delete_intersected_triangles!(tri, intersecting_triangles)
    cache = get_cache(tri)
    for cavity in (left_cavity, right_cavity)
        empty!(cache)
        tri_cavity = get_triangulation(cache)
        tri_fan = get_triangulation_2(cache)
        marked_vertices = get_marked_vertices(cache)
        fan_triangles = get_fan_triangles(cache)
        triangulate_cavity_cdt!(tri_cavity, cavity, tri_fan, marked_vertices, fan_triangles; rng, predicates)
        add_new_triangles!(tri, tri_cavity)
    end
    return tri
end
add_segment!(tri::Triangulation, i, j; predicates::AbstractPredicateKernel=Adaptive(), rng=Random.default_rng()) = add_segment!(tri, construct_edge(edge_type(tri), i, j); predicates, rng)