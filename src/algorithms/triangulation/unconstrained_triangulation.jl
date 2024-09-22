"""
    get_insertion_order(points, randomise, skip_points, ::Type{I}, rng) where {I} -> Vector{I}
    get_insertion_order(tri::Triangulation, randomise, skip_points, rng) -> Vector{I}

Gets the insertion order for points into a triangulation. 

# Arguments 
- `points`: The points to insert.
- `randomise`: If `true`, then the insertion order is randomised. Otherwise, the insertion order is the same as the order of the points.
- `skip_points`: The points to skip.
- `I::Type{I}`: The type of the vertices.
- `rng::Random.AbstractRNG`: The random number generator to use.

# Output 
- `order`: The order to insert the points in.

!!! warning "Mutation of `order`"

    This `order` might be mutated (by `circshift!`) in [`get_initial_triangle`](@ref).
"""
function get_insertion_order(points, randomise, skip_points, ::Type{I}, rng) where {I}
    point_indices = each_point_index(points)
    collected_point_indices = convert(Vector{I}, point_indices)
    randomise && Random.shuffle!(rng, collected_point_indices)
    setdiff!(collected_point_indices, skip_points)
    return collected_point_indices
end
function get_insertion_order(tri::Triangulation, randomise, skip_points, rng)
    return get_insertion_order(get_points(tri), randomise, skip_points, integer_type(tri), rng)
end

"""
    get_initial_triangle(tri::Triangulation, insertion_order, itr=0) -> Triangle 

Gets the initial triangle for the Bowyer-Watson algorithm. 

# Arguments 
- `tri`: The triangulation.
- `insertion_order`: The insertion order of the points. See [`get_insertion_order`](@ref).
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `itr=0`: To avoid issues with degenerate triangles and infinite loops, this counts the number of times `insertion_order` had to be shifted using `circshift!` to find an initial non-degenerate triangle.

# Output
- `initial_triangle`: The initial triangle.
"""
function get_initial_triangle(tri::Triangulation, insertion_order, predicates::AbstractPredicateKernel = AdaptiveKernel(), itr = 0)
    i, j, k = @view insertion_order[1:3] # insertion_order got converted into a Vector, so indexing is safe 
    initial_triangle = construct_positively_oriented_triangle(tri, i, j, k, predicates)
    i, j, k = triangle_vertices(initial_triangle)
    degenerate_cert = triangle_orientation(predicates, tri, i, j, k)
    if length(insertion_order) > 3 && (is_degenerate(degenerate_cert) || check_precision(triangle_area(tri, (i, j, k)))) && itr ≤ length(insertion_order) # Do not get stuck in an infinite loop if there are just three points, the three of them being collinear. The itr ≤ length(insertion_order) is needed because, if all the points are collinear, the loop could go on forever 
        @static if VERSION ≥ v"1.8.1"
            circshift!(insertion_order, -1)
        else
            _insertion_order = circshift(insertion_order, -1)
            copyto!(insertion_order, _insertion_order)
        end
        return get_initial_triangle(tri, insertion_order, predicates, itr + 1)
    end
    return initial_triangle
end

"""
    initialise_bowyer_watson!(tri::Triangulation, insertion_order, predicates::AbstractPredicateKernel=AdaptiveKernel()) -> Triangulation

Initialises the Bowyer-Watson algorithm.

# Arguments
- `tri`: The triangulation.
- `insertion_order`: The insertion order of the points. See [`get_insertion_order`](@ref).
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Output
`tri` is updated in place to contain the initial triangle from which the Bowyer-Watson algorithm starts.
"""
function initialise_bowyer_watson!(tri::Triangulation, insertion_order, predicates::AbstractPredicateKernel = AdaptiveKernel())
    I = integer_type(tri)
    initial_triangle = get_initial_triangle(tri, insertion_order, predicates)
    u, v, w = triangle_vertices(initial_triangle)
    g = I(𝒢)
    add_triangle!(tri, u, v, w; protect_boundary = true, update_ghost_edges = false)
    add_triangle!(tri, v, u, g; protect_boundary = true, update_ghost_edges = false)
    add_triangle!(tri, w, v, g; protect_boundary = true, update_ghost_edges = false)
    add_triangle!(tri, u, w, g; protect_boundary = true, update_ghost_edges = false)
    new_representative_point!(tri, I(1))
    for i in triangle_vertices(initial_triangle)
        p = get_point(tri, i)
        update_centroid_after_addition!(tri, I(1), p)
    end
    return tri
end

"""
    get_initial_search_point(tri::Triangulation, num_points, new_point, insertion_order, num_sample_rule::F, rng, try_last_inserted_point) where {F} -> Vertex 

For a given iteration of the Bowyer-Watson algorithm, finds the point to start the point location with [`find_triangle`](@ref) at.

# Arguments
- `tri`: The triangulation.
- `num_points`: The number of points currently in the triangulation.
- `new_point`: The point to insert.
- `insertion_order`: The insertion order of the points. See [`get_insertion_order`](@ref).
- `num_sample_rule::F`: The rule to use to determine the number of points to sample. See [`default_num_samples`](@ref) for the default. 
- `rng::Random.AbstractRNG`: The random number generator to use.
- `try_last_inserted_point`: If `true`, then the last inserted point is also considered as the start point. 

# Output
- `initial_search_point`: The vertex to start the point location with [`find_triangle`](@ref) at.
"""
function get_initial_search_point(tri::Triangulation, num_points, new_point, insertion_order, num_sample_rule::F, rng, try_last_inserted_point) where {F}
    num_currently_inserted = num_points + 3 - 1     # + 3 for the points already inserted
    last_inserted_point_index = insertion_order[num_currently_inserted]
    currently_inserted_points = each_solid_vertex(tri) # We can't just do something like insertion_order[1:num_currently_inserted] because, for weighted triangulations, not all previous points may still be in the triangulation if they are submerged
    m = num_sample_rule(num_currently_inserted)
    try_points = try_last_inserted_point ? (last_inserted_point_index,) : (oftype(last_inserted_point_index, ∅),)
    initial_search_point = select_initial_point(tri, new_point; m, point_indices = currently_inserted_points, rng, try_points)
    return initial_search_point
end

"""
    add_point_bowyer_watson!(tri::Triangulation, new_point, initial_search_point::I, rng::Random.AbstractRNG=Random.default_rng(), update_representative_point=true, store_event_history=Val(false), event_history=nothing, peek::P=Val(false), predicates::AbstractPredicateKernel=AdaptiveKernel()) -> Triangle 

Adds `new_point` into `tri`. 

# Arguments 
- `tri`: The triangulation.
- `new_point`: The point to insert.
- `initial_search_point::I`: The vertex to start the point location with [`find_triangle`](@ref) at. See [`get_initial_search_point`](@ref).
- `rng::Random.AbstractRNG`: The random number generator to use.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `update_representative_point=true`: If `true`, then the representative point is updated. See [`update_centroid_after_addition!`](@ref).
- `store_event_history=Val(false)`: If `true`, then the event history from the insertion is stored. 
- `event_history=nothing`: The event history to store the event history in. Should be an [`InsertionEventHistory`](@ref) if `store_event_history` is `true`, and `false` otherwise.
- `peek=Val(false)`: Whether to actually add `new_point` into `tri`, or just record into `event_history` all the changes that would occur from its insertion.

# Output 
- `V`: The triangle in `tri` containing `new_point`.

# Extended help 
This function works as follows:

1. First, the triangle containing the new point, `V`, is found using [`find_triangle`](@ref).
2. Once the triangle is found, we call into `add_point_bowyer_watson_and_process_after_found_triangle` to properly insert the point. 
3. Inside `add_point_bowyer_watson_and_process_after_found_triangle`, we first call into `add_point_bowyer_watson_after_found_triangle` to add the point into the cavity. We then 
   call into `add_point_bowyer_watson_onto_segment` to make any changes necessary incase the triangulation is constrained and `new_point` lies on a segment, since the depth-first search of the triangles containing `new_point` in its circumcenter
   must be performed on each side of the segment that `new_point` lies on.

The function [`add_point_bowyer_watson_dig_cavities!`](@ref) is the main workhorse of this function from `add_point_bowyer_watson_after_found_triangle`. See its docstring for the details. 
"""
function add_point_bowyer_watson!(tri::Triangulation, new_point, initial_search_point::I, rng::Random.AbstractRNG = Random.default_rng(), predicates::AbstractPredicateKernel = AdaptiveKernel(), update_representative_point = true, store_event_history = Val(false), event_history = nothing, peek::P = Val(false)) where {I, P}
    _new_point = is_true(peek) ? new_point : I(new_point)
    q = get_point(tri, _new_point)
    V = find_triangle(tri, q; predicates, m = nothing, point_indices = nothing, try_points = nothing, k = initial_search_point, rng)
    if is_weighted(tri)
        cert = point_position_relative_to_circumcircle(predicates, tri, V, _new_point; cache = get_orient3_cache(tri)) # redirects to point_position_relative_to_witness_plane
        is_outside(cert) && return V # If the point is submerged, then we don't add it
    end
    flag = point_position_relative_to_triangle(predicates, tri, V, q)
    add_point_bowyer_watson_and_process_after_found_triangle!(tri, _new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek, predicates)
    return V
end

function add_point_bowyer_watson_and_process_after_found_triangle!(tri::Triangulation, new_point, V, q, flag, update_representative_point = true, store_event_history = Val(false), event_history = nothing, peek::P = Val(false), predicates::AbstractPredicateKernel = AdaptiveKernel()) where {P}
    I = integer_type(tri)
    _new_point = is_true(peek) ? new_point : I(new_point)
    add_point_bowyer_watson_after_found_triangle!(tri, _new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek, predicates)
    add_point_bowyer_watson_onto_segment!(tri, _new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek, predicates)
    return V
end

function add_point_bowyer_watson_after_found_triangle!(tri::Triangulation, new_point, V, q, flag, update_representative_point = true, store_event_history = Val(false), event_history = nothing, peek::P = Val(false), predicates::AbstractPredicateKernel = AdaptiveKernel()) where {P}
    if is_ghost_triangle(V) && is_constrained(tri)
        # When we have a constrained boundary edge, we don't want to walk into its interior. So let's just check this case now.
        V = sort_triangle(V)
        u, v, w = triangle_vertices(V) # w is the ghost vertex
        if !contains_boundary_edge(tri, v, u)
            return add_point_bowyer_watson_dig_cavities!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek, predicates)
        else
            return tri
        end
    else
        return add_point_bowyer_watson_dig_cavities!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek, predicates)
    end
end

"""
    add_point_bowyer_watson_dig_cavities!(tri::Triangulation, new_point::N, V, q, flag, update_representative_point=true, store_event_history=Val(false), event_history=nothing, peek::F=Val(false), predicates::AbstractPredicateKernel=AdaptiveKernel()) where {N,F} 

Deletes all the triangles in `tri` whose circumcircle contains `new_point`. This leaves behind a polygonal cavity, whose boundary edges are then connected to `new_point`, restoring the Delaunay property from `new_point`'s insertion.

# Arguments 
- `tri`: The [`Triangulation`](@ref).
- `new_point::N`: The point to insert.
- `V`: The triangle in `tri` containing `new_point`.
- `q`: The point to insert.
- `flag`: The position of `q` relative to `V`. See [`point_position_relative_to_triangle`](@ref).
- `update_representative_point=true`: If `true`, then the representative point is updated. See [`update_centroid_after_addition!`](@ref).
- `store_event_history=Val(false)`: If `true`, then the event history from the insertion is stored.
- `event_history=nothing`: The event history to store the event history in. Should be an [`InsertionEventHistory`](@ref) if `store_event_history` is `true`, and `false` otherwise.
- `peek=Val(false)`: Whether to actually add `new_point` into `tri`, or just record into `event_history` all the changes that would occur from its insertion.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Output
There are no changes, but `tri` is updated in-place.

# Extended help 
This function works as follows: 

1. To dig the cavity, we call [`dig_cavity!`](@ref) on each edge of `V`, stepping towards the adjacent triangles to excavate the cavity recursively.
2. Once the cavity has been excavated, extra care is needed in case `is_on(flag)`, meaning `new_point` is on one of the edges of `V`. In particular, extra care is needed if:
   `is_on(flag) && (is_boundary_triangle(tri, V) || is_ghost_triangle(V) && !is_boundary_node(tri, new_point)[1])`.
   The need for this check is in case `new_point` is on a boundary edge already exists, since we need to fix the associated ghost edges. For example, a boundary edge `(i, k)` might have been 
   split into the edges `(i, j)` and `(j, k)`, which requires that the ghost triangle `(k, i, g)` be split into `(j, i, g)` and `(k, j, g)`, where `g` is the ghost vertex. This part of the function 
   will fix this case. The need for `is_ghost_triangle(V) && !is_boundary_node(tri, new_point)[1]` is in case the ghost edges were already correctly added. Nothing happens if the edge of `V` that `new_point` 
   is on is not the boundary edge.
"""
function add_point_bowyer_watson_dig_cavities!(tri::Triangulation, new_point::N, V, q, flag, update_representative_point = true, store_event_history = Val(false), event_history = nothing, peek::F = Val(false), predicates::AbstractPredicateKernel = AdaptiveKernel()) where {N, F}
    _new_point = is_true(peek) ? num_points(tri) + 1 : new_point # If we are peeking, then we need to use the number of points in the triangulation as the index for the new point since we don't actually insert the point
    i, j, k = triangle_vertices(V)
    ℓ₁ = get_adjacent(tri, j, i)
    ℓ₂ = get_adjacent(tri, k, j)
    ℓ₃ = get_adjacent(tri, i, k)
    !is_true(peek) && delete_triangle!(tri, V; protect_boundary = true, update_ghost_edges = false)
    is_true(store_event_history) && delete_triangle!(event_history, V)
    dig_cavity!(tri, new_point, i, j, ℓ₁, flag, V, store_event_history, event_history, peek, predicates)
    dig_cavity!(tri, new_point, j, k, ℓ₂, flag, V, store_event_history, event_history, peek, predicates)
    dig_cavity!(tri, new_point, k, i, ℓ₃, flag, V, store_event_history, event_history, peek, predicates)
    if is_on(flag) && (is_boundary_triangle(tri, V) || is_ghost_triangle(V) && !is_boundary_node(tri, new_point)[1])
        # ^ Need to fix the ghost edges if the point is added onto an existing boundary edge. Note that the last 
        #   condition is in case the ghost edges were already correctly added.
        # This isn't done using split_edge! since there are some special cases to consider with constraints, and also because we need to use peek here.
        e = find_edge(predicates, tri, V, new_point)
        u, v = edge_vertices(e)
        is_bnd = is_boundary_edge(tri, u, v) || is_boundary_edge(tri, v, u)
        if is_bnd # If the edge is not itself a boundary edge, no need to worry.
            if !is_boundary_edge(tri, v, u)
                u, v = v, u
            end
            g = get_adjacent(tri, u, v)
            if !is_true(peek)
                delete_triangle!(tri, v, u, new_point; protect_boundary = true, update_ghost_edges = false)
                delete_triangle!(tri, u, v, g; protect_boundary = true, update_ghost_edges = false)
                add_triangle!(tri, new_point, v, g; protect_boundary = true, update_ghost_edges = false)
                add_triangle!(tri, u, new_point, g; protect_boundary = true, update_ghost_edges = false)
            end
            if is_true(store_event_history)
                trit = triangle_type(tri)
                delete_triangle!(event_history, construct_triangle(trit, u, v, g))
                add_triangle!(event_history, construct_triangle(trit, _new_point, v, g))
                add_triangle!(event_history, construct_triangle(trit, u, _new_point, g))
            end
            if is_constrained(tri) && is_bnd # If we don't do this here now, then when we try and do it later we will get a KeyError since we've already modified the boundary edge but we wouldn't have updated the segment fields
                # We also only do this if contains_boundary_edge since we want to assume that (u, v) does not appear in interior_segments
                if contains_boundary_edge(tri, u, v)
                    if !is_true(peek)
                        split_boundary_edge!(tri, u, v, new_point)
                        if is_curve_bounded(tri)
                            enricher = get_boundary_enricher(tri)
                            split_boundary_edge!(enricher, u, v, new_point, Val(false))
                        end
                    end
                    is_true(store_event_history) && split_boundary_edge!(event_history, u, v, _new_point)
                elseif contains_boundary_edge(tri, v, u)
                    if !is_true(peek)
                        split_boundary_edge!(tri, v, u, new_point)
                        if is_curve_bounded(tri)
                            enricher = get_boundary_enricher(tri)
                            split_boundary_edge!(enricher, v, u, new_point, Val(false))
                        end
                    end
                    is_true(store_event_history) && split_boundary_edge!(event_history, v, u, _new_point)
                end
            end
        end
    end
    I = integer_type(tri)
    update_representative_point && !is_true(peek) && update_centroid_after_addition!(tri, I(1), q) # How do we efficiently determine which curve to update for a given q when adding into an existing triangulation?
    return tri
end

"""
    enter_cavity(tri::Triangulation, r, i, j, ℓ, predicates::AbstractPredicateKernel=AdaptiveKernel()) -> Bool

Determines whether to enter the cavity in `tri` through the edge `(i, j)` when inserting `r` into the triangulation.

# Arguments 
- `tri`: The [`Triangulation`](@ref).
- `r`: The new point being inserted.
- `i`: The first vertex of the edge `(i, j)`.
- `j`: The second vertex of the edge `(i, j)`.
- `ℓ`: The vertex adjacent to `(j, i)`, so that the triangle being stepped into is `(j, i, ℓ)`.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Output
- `true` if the cavity should be entered, and `false` otherwise. See also [`dig_cavity!`](@ref) and [`point_position_relative_to_circumcircle`](@ref).
"""
function enter_cavity(tri::Triangulation, r, i, j, ℓ, predicates::AbstractPredicateKernel = AdaptiveKernel())
    contains_segment(tri, i, j)  && return false
    if is_ghost_vertex(ℓ)
        cert = point_position_relative_to_circumcircle(predicates, tri, j, i, ℓ, r; cache = get_incircle_cache(tri))
    else
        cert = point_position_relative_to_circumcircle(predicates, tri, r, i, j, ℓ; cache =  get_incircle_cache(tri))
    end
    return is_weighted(tri) ? !is_outside(cert) : is_inside(cert)
end

"""
    dig_cavity!(tri::Triangulation, r, i, j, ℓ, flag, V, store_event_history=Val(false), event_history=nothing, peek::F=Val(false), predicates::AbstractPredicateKernel=AdaptiveKernel()) where {F}

Excavates the cavity in `tri` through the edge `(i, j)`, stepping towards the adjacent triangles to excavate the cavity recursively, eliminating all 
triangles containing `r` in their circumcircle. 

# Arguments 
- `tri`: The [`Triangulation`](@ref).
- `r`: The new point being inserted. 
- `i`: The first vertex of the edge `(i, j)`.
- `j`: The second vertex of the edge `(i, j)`.
- `ℓ`: The vertex adjacent to `(j, i)`, so that the triangle being stepped into is `(j, i, ℓ)`.
- `flag`: The position of `r` relative to `V`. See [`point_position_relative_to_triangle`](@ref).
- `V`: The triangle in `tri` containing `r`.
- `store_event_history=Val(false)`: If `true`, then the event history from the insertion is stored.
- `event_history=nothing`: The event history to store the event history in. Should be an [`InsertionEventHistory`](@ref) if `store_event_history` is `true`, and `false` otherwise.
- `peek=Val(false)`: Whether to actually add `new_point` into `tri`, or just record into `event_history` all the changes that would occur from its insertion.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Output
There are no changes, but `tri` is updated in-place.

# Extended help
This function works as follows:

1. First, we check if `ℓ` is $∅, which would imply that the triangle `(j, i, ℓ)` doesn't exist as it has already been deleted. If this is the check, we return and stop digging. 
2. If the edge `(i, j)` is not a segment, `r` is inside of the circumcenter of `(j, i, ℓ)`, and `ℓ` is not a ghost vertex, then we can step forward into the next triangles. In particular, 
   we delete the triangle `(j, i, ℓ)` from `tri` and then call `dig_cavity!` again on two edges of `(j, i, ℓ)` other than `(j, i)`.
3. If we did not do step 2, then this means that we are on the edge of the polygonal cavity; this also covers the case that `ℓ` is a ghost vertex, i.e. we are at the boundary of the triangulation.
   There are two cases to consider here. Firstly, if `(i, j)` is not a segment, we just need to add the triangle `(r, i, j)` into `tri` to connect `r` to this edge of the polygonal cavity. If 
   `(i, j)` is a segment, then the situation is more complicated. In particular, `r` being on an edge of `V` might imply that we are going to add a degenerate triangle `(r, i, j)` into `tri`, 
   and so this needs to be avoided. So, we check if `is_on(flag) && contains_segment(tri, i, j)` and, if the edge that `r` is on is `(i, j)`, we add the triangle `(r, i, j)`. Otherwise, we do nothing.
"""
function dig_cavity!(tri::Triangulation, r, i, j, ℓ, flag, V, store_event_history = Val(false), event_history = nothing, peek::F = Val(false), predicates::AbstractPredicateKernel = AdaptiveKernel()) where {F}
    if !edge_exists(ℓ)
        # The triangle has already been deleted in this case.
        return tri
    end
    _r = is_true(peek) ? num_points(tri) + 1 : r # If we are peeking, then we need to use the number of points in the triangulation as the index for the new point since we don't actually insert the point
    if enter_cavity(tri, r, i, j, ℓ, predicates)
        ℓ₁ = get_adjacent(tri, ℓ, i)
        ℓ₂ = get_adjacent(tri, j, ℓ)
        !is_true(peek) && delete_triangle!(tri, j, i, ℓ; protect_boundary = true, update_ghost_edges = false)
        dig_cavity!(tri, r, i, ℓ, ℓ₁, flag, V, store_event_history, event_history, peek, predicates)
        dig_cavity!(tri, r, ℓ, j, ℓ₂, flag, V, store_event_history, event_history, peek, predicates)
        if is_true(store_event_history)
            trit = triangle_type(tri)
            delete_triangle!(event_history, construct_triangle(trit, j, i, ℓ))
        end
    else
        # If we are here, then this means that we are on an edge of the polygonal cavity. 
        # Note that this also covers the is_ghost_vertex(ℓ) case.
        if is_on(flag) && contains_segment(tri, i, j)
            # When we have a segment (i, j) and we add a point r that lies on the edge, 
            # we can run into an issue where we add a degenerate triangle (i, j, r). We need to 
            # check this. There is probably a much smarter way to do this, but this should only 
            # be done very rarely anyway, so I'm not too concerned about the performance hit here. 
            e = find_edge(predicates, tri, V, r)
            u, v = edge_vertices(e)
            if u == i && v == j
                return tri
            else
                !is_true(peek) && add_triangle!(tri, r, i, j; protect_boundary = true, update_ghost_edges = false)
                if is_true(store_event_history)
                    trit = triangle_type(tri)
                    add_triangle!(event_history, construct_triangle(trit, _r, i, j))
                end
            end
        else
            !is_true(peek) && add_triangle!(tri, r, i, j; protect_boundary = true, update_ghost_edges = false)
            if is_true(store_event_history)
                trit = triangle_type(tri)
                add_triangle!(event_history, construct_triangle(trit, _r, i, j))
            end
        end
    end
    return tri
end

function add_point_bowyer_watson_onto_segment!(tri::Triangulation, new_point, V, q, flag, update_representative_point = true, store_event_history = Val(false), event_history = nothing, peek::P = Val(false), predicates::AbstractPredicateKernel = AdaptiveKernel()) where {P}
    _new_point = is_true(peek) ? num_points(tri) + 1 : new_point # If we are peeking, then we need to use the number of points in the triangulation as the index for the new point since we don't actually insert the point
    if is_on(flag) && is_constrained(tri)
        # If the point we are adding appears on a segment, then we perform the depth-first search 
        # on each side the segment. We also need to update the segment lists.
        e = find_edge(predicates, tri, V, new_point)
        u, v = edge_vertices(e)
        if contains_segment(tri, u, v)
            is_bnd = (contains_boundary_edge(tri, u, v) || contains_boundary_edge(tri, v, u))
            # First, let us search on the other side of the segment, provided we are not on a boundary edge. 
            w = get_adjacent(tri, v, u)
            T = triangle_type(tri)
            V = construct_triangle(T, v, u, w)
            add_point_bowyer_watson_dig_cavities!(tri, new_point, V, q, flag, update_representative_point, store_event_history, event_history, peek, predicates)
            # Now, we need to replace the segment by this new segment.
            E = edge_type(tri)
            interior_segments = get_interior_segments(tri)
            all_segments = get_all_segments(tri)
            if !is_true(peek)
                for edges in (interior_segments, all_segments)
                    delete_unoriented_edge!(edges, construct_edge(E, u, v))
                    !contains_edge(construct_edge(E, new_point, u), edges) && !is_bnd && add_edge!(edges, construct_edge(E, u, new_point))
                    !contains_edge(construct_edge(E, v, new_point), edges) && !is_bnd && add_edge!(edges, construct_edge(E, new_point, v))
                end
            end
            if is_true(store_event_history) && !is_bnd
                !contains_edge(construct_edge(E, v, u), event_history.deleted_segments) && delete_edge!(event_history, construct_edge(E, u, v))
                !contains_edge(construct_edge(E, _new_point, u), event_history.added_segments) && add_edge!(event_history, construct_edge(E, u, _new_point))
                !contains_edge(construct_edge(E, v, _new_point), event_history.added_segments) && add_edge!(event_history, construct_edge(E, _new_point, v))
            end
        end
        if contains_boundary_edge(tri, u, v)
            !is_true(peek) && split_boundary_edge!(tri, u, v, new_point)
            is_true(store_event_history) && split_boundary_edge!(event_history, u, v, _new_point)
        elseif contains_boundary_edge(tri, v, u)
            !is_true(peek) && split_boundary_edge!(tri, v, u, new_point)
            is_true(store_event_history) && split_boundary_edge!(event_history, v, u, _new_point)
        end
    end
    return tri
end

"""
    unconstrained_triangulation!(tri::Triangulation; kwargs...)

Computes the unconstrained Delaunay triangulation of the points in `tri`.

# Arguments
- `tri`: The triangulation.

# Keyword Arguments 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `randomise=true`: If `true`, then the insertion order is randomised. Otherwise, the insertion order is the same as the order of the points.
- `skip_points=()`: The vertices to skip. 
- `num_sample_rule::M=default_num_samples`: The rule to use to determine the number of points to sample. See [`default_num_samples`](@ref) for the default.
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator to use.
- `insertion_order=get_insertion_order(tri, randomise, skip_points, rng)`: The insertion order of the points. See [`get_insertion_order`](@ref).

# Outputs 
There is no output, but `tri` is updated in-place.
"""
function unconstrained_triangulation!(
        tri::Triangulation;
        predicates::AbstractPredicateKernel = AdaptiveKernel(),
        randomise = true,
        try_last_inserted_point = true,
        skip_points = (),
        num_sample_rule::M = default_num_samples,
        rng::Random.AbstractRNG = Random.default_rng(),
        insertion_order = get_insertion_order(tri, randomise, skip_points, rng),
    ) where {M}
    initialise_bowyer_watson!(tri, insertion_order, predicates)
    remaining_points = @view insertion_order[(begin + 3):end]
    for (num_points, new_point) in enumerate(remaining_points)
        initial_search_point = get_initial_search_point(tri, num_points, new_point, insertion_order, num_sample_rule, rng, try_last_inserted_point)
        add_point_bowyer_watson!(tri, new_point, initial_search_point, rng, predicates)
    end
    convex_hull!(tri; predicates, reconstruct = false)
    return tri
end
