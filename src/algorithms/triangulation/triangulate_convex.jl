"""
    triangulate_convex(points, S; delete_ghosts=false, delete_empty_features=true, rng=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel(), kwargs...) -> Triangulation

Triangulates the convex polygon `S`.

# Arguments 
- `points`: The point set corresponding to the vertices in `S`.
- `S`: A convex polygon represented as a vector of vertices. The vertices should be given in counter-clockwise order, and must not be circular so that `S[begin] ≠ S[end]`.

# Keyword Arguments
- `delete_ghosts=false`: If `true`, the ghost triangles are deleted after triangulation. 
- `delete_empty_features=true`: If `true`, the empty features are deleted after triangulation.
- `rng=Random.default_rng()`: The random number generator used to shuffle the vertices of `S` before triangulation.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `kwargs...`: Additional keyword arguments passed to `Triangulation`.

# Output 
- `tri::Triangulation`: The triangulated polygon. 
"""
function triangulate_convex(
        points, S;
        rng::Random.AbstractRNG = Random.default_rng(),
        delete_ghosts = false,
        delete_empty_features = true,
        predicates::AbstractPredicateKernel = AdaptiveKernel(),
        kwargs...,
    )
    tri = Triangulation(points; kwargs...)
    triangulate_convex!(tri, S; predicates, rng)
    postprocess_triangulate_convex!(tri, S; delete_ghosts, delete_empty_features)
    return tri
end

"""
    triangulate_convex!(tri::Triangulation, S; rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())

Triangulates the convex polygon `S` in-place into `tri`.

# Arguments
- `tri::Triangulation`: The triangulation to be modified.
- `S`: A convex polygon represented as a vector of vertices. The vertices should be given in counter-clockwise order, and must not be circular so that `S[begin] ≠ S[end]`.

# Keyword Arguments
- `rng=Random.default_rng()`: The random number generator used to shuffle the vertices of `S` before triangulation.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs 
There is no output, as `tri` is updated in-place. This function does not do any post-processing, e.g. deleting any ghost triangles. This is done by 
[`triangulate_convex`](@ref) or [`postprocess_triangulate_convex!`](@ref).
"""
function triangulate_convex!(tri::Triangulation, S; rng::Random.AbstractRNG = Random.default_rng(), predicates::AbstractPredicateKernel = AdaptiveKernel())
    list = ShuffledPolygonLinkedList(S; rng)
    delete_vertices_in_random_order!(list, tri, rng, predicates)
    u, v, w = get_triplet(list, 1)
    add_triangle!(tri, u, v, w; protect_boundary = true, update_ghost_edges = false)
    for i in 4:list.k
        u, v, w = get_triplet(list, i)
        add_point_convex_triangulation!(tri, u, v, w, S, predicates)
    end
    return tri
end

"""
    delete_vertices_in_random_order!(list::Triangulation, tri::ShuffledPolygonLinkedList, rng, predicates::AbstractPredicateKernel)

Deletes the vertices of the polygon represented by `list` in random order, done by switching 
the pointers of the linked list. Only three vertices will survive. If these these three vertices are 
collinear, then the deletion is attempted again after reshuffling the vertices.

# Arguments 
- `list::ShuffledPolygonLinkedList`: The linked list representing the polygon to be deleted.
- `tri::Triangulation`: The [`Triangulation`](@ref). 
- `rng::Random.AbstractRNG`: The random number generator used to shuffle the vertices of `S`. 
- `predicates::AbstractPredicateKernel`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
There is no output, as `list` is modified in-place.
"""
function delete_vertices_in_random_order!(list::ShuffledPolygonLinkedList, tri::Triangulation, rng::Random.AbstractRNG, predicates::AbstractPredicateKernel)
    for i in list.k:-1:4
        delete_vertex!(list, i)
    end
    # Check if the three surviving vertices that survived are collinear. Note that we do not 
    #    check for negative orientation here, since the polygon is provided in a counter-clockwise 
    #    order already.
    u, v, w = get_triplet(list, 1)
    a, b, c = get_point(tri, u, v, w)
    degenerate_cert = triangle_orientation(predicates, a, b, c)
    if !is_positively_oriented(degenerate_cert)
        reset!(list; rng)
        delete_vertices_in_random_order!(list, tri, rng, predicates)
    end
    return tri
end

"""
    add_point_convex_triangulation!(tri::Triangulation, u, v, w, S, predicates::AbstractPredicateKernel=AdaptiveKernel())

Adds the point `u` into the triangulation `tri`.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `u`: The vertex to add.
- `v`: The vertex next to `u`.
- `w`: The vertex previous to `u`.
- `S`: The set of vertices of the polygon.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
There is no output, as `tri` is modified in-place.

# Extended help 
This function forms part of Chew's algorithm for triangulating a convex polygon. There are some important points 
to make.

1. Firstly, checking that `x = get_adjacent(tri, w, v)` is needed to prevent the algorithm from exiting the polygon.
  This is important in case this algorithm is used as part of [`delete_point!`](@ref). When you are just triangulating
  a convex polygon by itself, this checked is the same as checking `edge_exists(tri, w, v)`.
2. For this method to be efficient, the set `x ∈ S` must be `O(1)` time. This is why we use a `Set` type for `S`.
3. The algorithm is recursive, recursively digging further through the polygon to find non-Delaunay edges to adjoins with `u`.
"""
function add_point_convex_triangulation!(tri::Triangulation, u, v, w, S, predicates::AbstractPredicateKernel = AdaptiveKernel())
    x = get_adjacent(tri, w, v)
    if edge_exists(x) && is_inside(point_position_relative_to_circumcircle(predicates, tri, u, v, w, x; cache = get_incircle_cache(tri)))
        # uvw and wvx are not Delaunay 
        delete_triangle!(tri, w, v, x; protect_boundary = true, update_ghost_edges = false)
        add_point_convex_triangulation!(tri, u, v, x, S, predicates)
        add_point_convex_triangulation!(tri, u, x, w, S, predicates)
    else
        # vw is a Delaunay edge 
        add_triangle!(tri, u, v, w; protect_boundary = true, update_ghost_edges = false)
    end
    return tri
end

"""
    postprocess_triangulate_convex!(tri::Triangulation, S; delete_ghosts, delete_empty_features)

Postprocesses the completed triangulation `tri` of the convex polygon `S`.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `S`: The vertices of the convex polygon, as in [`triangulate_convex`](@ref).

# Keyword Arguments
- `delete_ghosts=false`: If `true`, the ghost triangles are deleted after triangulation.
- `delete_empty_features=true`: If `true`, the empty features are deleted after triangulation.

# Output 
There are no output, as `tri` is modified in-place. The postprocessing that is done is:

1. The convex hull of `tri` is set to `S`.
2. The ghost triangles are deleted if `delete_ghosts=true`.
3. The empty features are deleted if `delete_empty_features=true`.
4. The representative points are set to the centroid of `S`.
"""
function postprocess_triangulate_convex!(tri::Triangulation, S; delete_ghosts, delete_empty_features)
    hull = get_convex_hull_vertices(tri)
    append!(hull, S)
    push!(hull, S[begin])
    I = integer_type(tri)
    for i in firstindex(S):(lastindex(S) - 1)
        u = S[i]
        v = S[i + 1]
        if !delete_ghosts
            add_triangle!(tri, v, u, I(𝒢); protect_boundary = true, update_ghost_edges = false)
        else
            # Still want the ghost vertex in the graph, adjacent2vertex map, and the adjacent map.
            add_neighbour!(tri, I(𝒢), u)
            add_adjacent2vertex!(tri, I(𝒢), v, u)
            add_adjacent!(tri, v, u, I(𝒢))
        end
    end
    u = S[end]
    v = S[begin]
    if !delete_ghosts
        add_triangle!(tri, v, u, I(𝒢); protect_boundary = true, update_ghost_edges = false)
    else
        add_neighbour!(tri, I(𝒢), u)
        add_adjacent2vertex!(tri, I(𝒢), v, u)
        add_adjacent!(tri, v, u, I(𝒢))
    end
    delete_empty_features && clear_empty_features!(tri)
    empty_representative_points!(tri)
    cx, cy = mean_points(get_points(tri), S)
    representative_point_list = get_representative_point_list(tri)
    representative_point_list[1] = RepresentativeCoordinates(cx, cy, length(S))
    return tri
end
