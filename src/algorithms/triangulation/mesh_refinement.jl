"""
    refine!(tri::Triangulation; kwargs...) 

Refines the given [`Triangulation`](@ref) `tri` to meet the given quality constraints.

See the documentation for more information about mesh refinement, e.g. convergence issues and issues with small input-angles.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref) to refine.

# Keyword Arguments
- `min_angle=30.0`: The minimum angle constraint, in degrees.
- `max_angle=180.0`: The maximum angle constraint, in degrees. 

!!! danger "Maximum angle constraints"

    Maximum angle constraints are not currently implemented.

- `min_area=get_area(tri) / 1e9`: The minimum area constraint.
- `max_area=typemax(number_type(tri))`: The maximum area constraint.
- `max_points=max(1_000, num_solid_vertices(tri))^2`: The maximum number of vertices allowed in the triangulation. Note that this refers to [`num_solid_vertices`](@ref), not the amount returned by [`num_points`](@ref).
- `seditious_angle=20.0`: The angle at which a triangle is considered seditious, in degrees. See [`is_triangle_seditious`](@ref).
- `custom_constraint=(tri, T) -> false`: A custom constraint function that takes a [`Triangulation`](@ref) and a triangle, and returns `true` if the triangle should be refined and `false` otherwise.
- `use_circumcenter=true`: Whether to insert circumcenters for refining a triangle or generalised Steiner points.

!!! danger "Generalised Steiner points"

    Generalised Steiner points are not yet implemented. Thus, this argument must be `true` (and the `steiner_scale` keyword below is ignored).

- `use_lens=true`: Whether to use the diametral lens or the diametral circle for checking encroachment.
- `steiner_scale=0.999`: The perturbation factor to use for generalised Steiner points if `use_circumcenter=false`. (Not currently used - see above.)
- `rng=Random.default_rng()`: The random number generator to use in case it is needed during point location.
- `concavity_protection=false`: Whether to use concavity protection or not for [`find_triangle`](@ref). Most likely not needed, but may help in pathological cases.

# Output 
The triangulation is refined in-place.

!!! warning "Duplicate points and unused points" 

    During refinement, points are often deleted, which may often lead to points in `get_points(tri)` that do not 
    appear anywhere in the triangulation. (This is why we recommend e.g. [`each_solid_vertex`](@ref) over [`each_point`](@ref).)
    Similarly, since points are deleted, when two triangles have a common circumcenter it might happen (if they are near an input segment)
    that a point is duplicated inside `get_points(tri)`, in case one circumcenter was deleted previously.
"""
function refine!(tri::Triangulation; kwargs...) # WARNING: Minimum error might be dropped below slightly from a triangle being split that was only just about the minimum area.
    is_weighted(tri) && throw(ArgumentError("You can not use mesh refinement on a weighted triangulation."))
    args = RefinementArguments(tri; kwargs...)
    refine!(tri, args)
    finalise!(tri, args)
    return tri
end

function refine!(tri::Triangulation, args::RefinementArguments)
    enqueue_all_encroached_segments!(args, tri)
    split_all_encroached_segments!(tri, args)
    enqueue_all_bad_triangles!(args, tri)
    while keep_iterating(tri, args)
        refine_itr!(tri, args)
    end
    # finalise!(tri, args) # don't finalise here, can be dangerous when reusing this method on an unconstrained triangulation
    return tri
end

"""
    refine_itr!(tri::Triangulation, args::RefinementArguments)

Performs a single iteration of the refinement algorithm.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref) to refine.
- `args::RefinementArguments`: The [`RefinementArguments`](@ref) for the refinement.

# Output
The triangulation is refined in-place.
"""
function refine_itr!(tri::Triangulation, args::RefinementArguments)
    T, ρ = popfirst_triangle!(args.queue)
    u, v, w = triangle_vertices(T)
    if !is_ghost_triangle(T) && get_adjacent(tri, u, v) == w # Need the last part in case T was already split and thus no longer exists 
        success = split_triangle!(tri, args, T)
        if is_encroachment_failure(success)
            !haskey(args.queue, T) && (args.queue[T] = ρ)
            split_all_encroached_segments!(tri, args)
        elseif is_successful_insertion(success)
            assess_added_triangles!(args, tri)
        end # the other case is is_precision_failure(success), in which case nothing is done.
    end
    return tri
end

"""
    enqueue_all_encroached_segments!(args::RefinementArguments, tri::Triangulation)

Enqueues all encroached segments in the triangulation into `args.queue`.
"""
function enqueue_all_encroached_segments!(args::RefinementArguments, tri::Triangulation)
    for e in each_segment(tri)
        if !haskey(args.queue, e) && is_encroached(tri, args, e)
            ℓ² = edge_length_sqr(tri, e)
            args.queue[e] = ℓ²
        end
    end
    return args
end

"""
    enqueue_all_bad_triangles!(args::RefinementArguments, tri::Triangulation)

Enqueues all bad triangles in the triangulation into `args.queue`.
"""
function enqueue_all_bad_triangles!(args::RefinementArguments, tri::Triangulation)
    for T in each_solid_triangle(tri)
        ρ, flag = assess_triangle_quality(tri, args, T)
        flag && (args.queue[T] = ρ)
    end
    return args
end

"""
    finalise!(tri::Triangulation, args::RefinementArguments)

Finalises the triangulation after refinement, e.g. by deleting ghost triangles and unlocking the convex hull if needed.
"""
function finalise!(tri::Triangulation, args::RefinementArguments)
    !args.had_ghosts && delete_ghost_triangles!(tri)
    args.locked_convex_hull && unlock_convex_hull!(tri; reconstruct = true)
    return tri
end

"""
    split_triangle!(tri::Triangulation, args::RefinementArguments, T) -> Certificate 

Splits a bad triangle `T` of `tri` to improve its quality.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref) to split a triangle of.
- `args::RefinementArguments`: The [`RefinementArguments`](@ref) for the refinement.
- `T`: The triangle to split.

# Output
- `cert`: A [`Certificate`](@ref) indicating whether the split was successful or not. In particular, returns one of:
    - `Cert.SuccessfulInsertion`: The triangle was split successfully.
    - `Cert.EncroachmentFailure`: The triangle was not split successfully as the newly inserted point encroached upon a segment.
    - `Cert.PrecisionFailure`: The triangle was not split successfully due to precision issues.
"""
function split_triangle!(tri::Triangulation, args::RefinementArguments, T)
    empty!(args.events)
    precision_cert, c = get_steiner_point(tri, args, T)
    is_precision_failure(precision_cert) && return precision_cert
    V, flag = locate_steiner_point(tri, args, T, c)
    c′, V′ = check_for_invisible_steiner_point(tri, V, T, flag, c)
    check_steiner_point_precision(tri, V′, c) && return Cert.PrecisionFailure
    push_point!(tri, c′)
    new_point = num_points(tri)
    segment_flag = check_for_steiner_point_on_segment(tri, V, V′, new_point, flag)
    segment_flag && push!(args.offcenter_split_list, new_point)
    add_point_bowyer_watson_and_process_after_found_triangle!(tri, new_point, V′, c′, flag, false, Val(true), args.events, Val(false))
    any_encroached = enqueue_newly_encroached_segments!(args, tri)
    if any_encroached
        segment_flag && pop!(args.offcenter_split_list, new_point)
        undo_insertion!(tri, args.events)
        return Cert.EncroachmentFailure
    else
        return Cert.SuccessfulInsertion
    end
end

"""
    get_steiner_point(tri::Triangulation, args::RefinementArguments, T) -> Certificate, Point

Computes the Steiner point for a triangle `T` of `tri` to improve its quality in [`split_triangle!`](@ref).

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref) to split a triangle of.
- `args::RefinementArguments`: The [`RefinementArguments`](@ref) for the refinement.
- `T`: The triangle to split.

# Output
- `precision_flag`: A [`Certificate`](@ref) which is `Cert.PrecisionFailure` if the Steiner point could not be computed due to precision issues, and `Cert.None` otherwise.
- `c`: The Steiner point. If `is_precision_failure(precision_flag)`, then this is just an arbitrary point of `T` to ensure type stability.
"""
function get_steiner_point(tri::Triangulation, args::RefinementArguments, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    A² = squared_triangle_area(p, q, r)
    A = max(zero(A²), sqrt(A²))
    check_precision(A) && return Cert.PrecisionFailure, q # the point q is just for type stability with the return 
    c = triangle_circumcenter(p, q, r, A)
    if !args.use_circumcenter
        # TODO: Implement generalised Steiner points so that c above is translated appropriately.
    end
    check_steiner_point_precision(tri, T, c) && return Cert.PrecisionFailure, c
    return Cert.None, c # nothing went wrong yet
end

"""
    check_steiner_point_precision(tri::Triangulation, T, c) -> Bool

Checks if the Steiner point `c` of a triangle `T` of `tri` can be computed without precision issues, returning `true` if there are precision issues and `false` otherwise.
"""
function check_steiner_point_precision(tri::Triangulation, T, c)
    is_ghost_triangle(T) && return true # we aren't supposed to get ghost triangles, unless we are _just_ off the boundary due to precision errors
    i, j, k = triangle_vertices(T)
    (px, py), (qx, qy), (rx, ry) = get_point(tri, i, j, k)
    cx, cy = getxy(c)
    rel_err_flag_1 = check_relative_precision(px, cx) && check_relative_precision(py, cy)
    rel_err_flag_2 = check_relative_precision(qx, cx) && check_relative_precision(qy, cy)
    rel_err_flag_3 = check_relative_precision(rx, cx) && check_relative_precision(ry, cy)
    rel_err_flag = rel_err_flag_1 || rel_err_flag_2 || rel_err_flag_3
    abs_err_flag_1 = check_absolute_precision(px, cx) && check_absolute_precision(py, cy)
    abs_err_flag_2 = check_absolute_precision(qx, cx) && check_absolute_precision(qy, cy)
    abs_err_flag_3 = check_absolute_precision(rx, cx) && check_absolute_precision(ry, cy)
    abs_err_flag = abs_err_flag_1 || abs_err_flag_2 || abs_err_flag_3
    return rel_err_flag || abs_err_flag
end

"""
    get_init_for_steiner_point(tri::Triangulation, T) -> Vertex 

Gets the initial vertex to start the search for the Steiner point of a triangle `T` of `tri` in [`get_steiner_point`](@ref). The 
initial vertex is chosen so that it is opposite the longest edge.
"""
function get_init_for_steiner_point(tri::Triangulation, T)
    # We want to start at the vertex associated with the largest angle, noting that the circumcenter is always opposite this vertex.
    # This vertex can be found by considering the vertex not belonging to the longest edge, since the largest angle is opposite the longest edge.
    # This will help make it more likely that we find the circumcenters that are outside of the triangle, ensuring that we also stop at barriers appropriately.
    i, j, k = triangle_vertices(T)
    ℓij² = edge_length_sqr(tri, i, j)
    ℓjk² = edge_length_sqr(tri, j, k)
    ℓki² = edge_length_sqr(tri, k, i)
    max_edge = max(ℓij², ℓjk², ℓki²)
    init = max_edge == ℓij² ? k : max_edge == ℓjk² ? i : j
    return init
end

"""
    locate_steiner_point(tri::Triangulation, args::RefinementArguments, T, c) -> Triangle, Cert 

Locates the Steiner point `c` of a triangle `T` of `tri` in [`get_steiner_point`](@ref). The Steiner point is located by walking from the initial vertex `init` to `c` using [`find_triangle`](@ref).

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref) to split a triangle of.
- `args::RefinementArguments`: The [`RefinementArguments`](@ref) for the refinement.
- `T`: The triangle that the Steiner point is from.
- `c`: The Steiner point.

# Output
- `V`: The triangle that the Steiner point is in.
- `flag`: A [`Certificate`](@ref) which is `Cert.On` if the Steiner point is on the boundary of `V`, `Cert.Outside` if the Steiner point is outside of `V`, and `Cert.Inside` if the Steiner point is inside of `V`.
"""
function locate_steiner_point(tri::Triangulation, args::RefinementArguments, T, c)
    flag = point_position_relative_to_triangle(tri, T, c)
    !is_outside(flag) && return T, flag # T is never a ghost triangle, so don't worry about checking is_on(flag) here
    init = get_init_for_steiner_point(tri, T)
    V, _ = find_triangle(tri, c; m = nothing, point_indices = nothing, try_points = nothing, k = init, args.rng, args.concavity_protection, use_barriers = Val(true))
    flag = point_position_relative_to_triangle(tri, V, c)
    if is_ghost_triangle(V) && is_on(flag)
        V = replace_ghost_triangle_with_boundary_triangle(tri, V)
    end
    return V, flag
end

# This function is needed for two reasons: Firstly, when using diametral lens, boundary triangles might have their circumcenters on the boundary, and so we cannot split the triangle to improve its quality.
# Secondly, if we have to walk past a segment to insert the circumcenter, then we will not improve the triangle's quality which is the whole point of split_triangle!. 
# So, to get around this, whenever this occurs we will use the triangle's centroid rather than the proposed Steiner point. The "invisible" part refers to the "visible" definition of a CDT.
# It's not sufficient to just try and classify and segments adjoining the triangle as encroached, as it could be that the triangle 
# itself has no encroached edges (but it will be near a segment, but searching for that directly is not worth it).
"""
    check_for_invisible_steiner_point(tri::Triangulation, V, T, flag, c) -> Point, Triangle

Determines if the Steiner point `c`'s insertion will not affect the quality of `T`, and if so instead changes `c` to be `T`'s centroid.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref) to split a triangle of.
- `V`: The triangle that the Steiner point is in.
- `T`: The triangle that the Steiner point is from.
- `flag`: A [`Certificate`](@ref) which is `Cert.On` if the Steiner point is on the boundary of `V`, `Cert.Outside` if the Steiner point is outside of `V`, and `Cert.Inside` if the Steiner point is inside of `V`.
- `c`: The Steiner point.

# Output
- `c′`: The Steiner point to use instead of `c`, which is `T`'s centroid if `c` is not suitable.
- `V′`: The triangle that the Steiner point is in, which is `T` if `c` is not suitable.
"""
function check_for_invisible_steiner_point(tri::Triangulation, V, T, flag, c)
    !is_outside(flag) && !is_ghost_triangle(V) && return c, V # don't need to check if the point is on the solid edge of a ghost triangle, since we already do that check previously via is_on(flag) in locate_steiner_point
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    c′ = triangle_centroid(p, q, r)
    return c′, T
end

"""
    check_for_steiner_point_on_segment(tri::Triangulation, V, V′, new_point, flag) -> Bool

Checks if the Steiner point with vertex `new_point` is on a segment. If so, then its vertex is pushed into the offcenter-split list from `args`,
indicating that it should no longer be regarded as a free vertex (see [`is_free`](@ref)).

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `V`: The triangle that the Steiner point was originally in prior to [`check_for_invisible_steiner_point`](@ref).
- `V′`: The triangle that the Steiner point is in.
- `new_point`: The vertex associated with the Steiner point.
- `flag`: A [`Certificate`](@ref) which is `Cert.On` if the Steiner point is on the boundary of `V`, `Cert.Outside` if the Steiner point is outside of `V`, and `Cert.Inside` if the Steiner point is inside of `V`.

# Output
- `onflag`: Whether the Steiner point is on a segment or not.
"""
function check_for_steiner_point_on_segment(tri::Triangulation, V, V′, new_point, flag)
    !compare_triangles(V, V′) && return false # the point is a centroid, so it won't be on a segment
    if is_on(flag)
        e = find_edge(tri, V, new_point)
        u, v = edge_vertices(e)
        contains_segment(tri, u, v) && return true
    end
    return false
end

"""
    enqueue_newly_encroached_segments!(args::RefinementArguments, tri::Triangulation) -> Bool

Enqueues all segments that are newly encroached upon after a point insertion into the triangulation into `args.queue`.

# Arguments 
- `args::RefinementArguments`: The [`RefinementArguments`](@ref) for the refinement.
- `tri::Triangulation`: The [`Triangulation`](@ref) to enqueue newly encroached segments of.

# Output
- `any_encroached`: Whether any segments were newly encroached upon.
"""
function enqueue_newly_encroached_segments!(args::RefinementArguments, tri::Triangulation)
    any_encroached = false
    steiner_vertex = num_points(tri)
    neighbouring_edges = get_adjacent2vertex(tri, steiner_vertex)
    for edge in each_edge(neighbouring_edges)
        if contains_segment(tri, edge) && is_encroached(tri, args, edge)
            any_encroached = true
            ℓ² = edge_length_sqr(tri, edge)
            args.queue[edge] = ℓ²
        end
    end
    return any_encroached
end

const MIDPOINT_TOLERANCE = 1.0e-6

"""
    compute_split_position(tri::Triangulation, args::RefinementArguments, e) -> NTuple{2, Float}

Computes the position to split a segment `e` of `tri` at in [`split_subsegment!`](@ref). 

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref) to split a segment of.
- `args::RefinementArguments`: The [`RefinementArguments`](@ref) for the refinement.
- `e`: The segment to split.

# Output
- `mx, my`: The position to split the segment at.

This point is computed according to a set of rules:

1. If `e` is not a subsegment, meaning it is an input segment, then its midpoint is returned.
2. If `e` is a subsegment and the segment adjoins two other distinct segments (one for each vertex) at an acute angle, as determined by 
   [`segment_vertices_adjoin_other_segments_at_acute_angle`](@ref), then the point is returned so that `e` can be split such that one of the new subsegments has a power-of-two 
   length between 1/4 and 1/2 of the length of `e`, computed using [`compute_concentric_shell_quarternary_split_position`](@ref).
3. If `e` is a subsegment and the segment adjoins one other segment at an acute angle, as determined by 
   [`segment_vertices_adjoin_other_segments_at_acute_angle`](@ref), then the point is returned so that `e` can be split such that one of the new subsegments has a power-of-two
    length between 1/3 and 2/3 of the length of `e`, computed using [`compute_concentric_shell_ternary_split_position`](@ref).
4. Otherwise, the midpoint is returned.
"""
function compute_split_position(tri::Triangulation, args::RefinementArguments, e)
    # We don't consider curve bounded domains here, since they might possibly require a 
    # set of split points. Thus, for type stability reasons, we just use the piecewise linear version and 
    # do something different for curve bounded domains inside split_subsegment!.
    return _compute_split_position_piecewise_linear(tri, args, e)
end
function _compute_split_position_piecewise_linear(tri::Triangulation, args::RefinementArguments, e)
    p, q = get_point(tri, initial(e), terminal(e))
    px, py = getxy(p)
    qx, qy = getxy(q)
    mx, my = midpoint(px, qx), midpoint(py, qy)
    F = number_type(tri)
    if !is_subsegment(args, e)
        push!(args.midpoint_split_list, num_points(tri) + 1)
        # Split at the midpoint when splitting for the first time
        return mx, my
    else
        num_adjoin, adjoin_vert = segment_vertices_adjoin_other_segments_at_acute_angle(tri, e)
        if num_adjoin == 2
            # In this case, t is always ≤ 1/2. We can choose the vertex that this split is relative 
            # to arbitrarily, but just leaving this t as it is means that we get different results for e and reverse_edge(e).
            # So, we choose the smallest vertex.
            if initial(e) < terminal(e)
                t = compute_concentric_shell_quarternary_split_position(p, q)
            else
                t = compute_concentric_shell_quarternary_split_position(q, p)
                t = one(t) - t
            end
        elseif num_adjoin == 1
            t = compute_concentric_shell_ternary_split_position(p, q)
            if adjoin_vert ≠ initial(e)
                t = one(t) - t
            end
        else # No need to worry about ping-pong encroachment
            t = one(F) / 2
        end
        if abs(t - 1 / 2) < MIDPOINT_TOLERANCE # close enough to the midpoint, so just bisect - every third split on a segment should be a bisection anyway. 
            push!(args.midpoint_split_list, num_points(tri) + 1)
            return mx, my
        else
            push!(args.offcenter_split_list, num_points(tri) + 1)
            ax, ay = px + t * (qx - px), py + t * (qy - py)
            return ax, ay
        end
    end
end

"""
    compute_concentric_shell_ternary_split_position(p, q) -> Float64

Returns the value of `t ∈ [0, 1]` that gives the most balanced ternary split of the segment `pq`, so that one of the 
segments has a power-of-two length and both segments have lengths between 1/3 and 2/3 of the length of `pq`.
"""
function compute_concentric_shell_ternary_split_position(p, q)
    ℓ = dist(p, q)
    s = balanced_power_of_two_ternary_split(ℓ)
    t = s / ℓ
    return t
end

"""
    compute_concentric_shell_quarternary_split_position(p, q) -> Float64

Returns the value of `t ∈ [0, 1]` that gives the most balanced quarternary split of the segment `pq`, so that one of the
segments has a power-of-two length between 1/4 and 1/2 of the length of `pq`.
"""
function compute_concentric_shell_quarternary_split_position(p, q)
    ℓ = dist(p, q)
    s = balanced_power_of_two_quarternary_split(ℓ)
    t = s / ℓ
    return t
end

"""
    balanced_power_of_two_ternary_split(ℓ) -> Float64

Returns the value of `s ∈ [0, ℓ]` that gives the most balanced ternary split of the segment `pq`, so `s` is a power-of-two 
and `s ∈ [ℓ / 3, 2ℓ / 3]`.
"""
function balanced_power_of_two_ternary_split(ℓ)
    # This is the same as doing rₛₕₑₗₗ = 2^k, k = ⌊log₂(rₘₐₓ)⌋, rₘₐₓ = 2ℓ/3.
    # This version is actually faster, though, by about 4x.
    balanced_split = one(ℓ)
    while ℓ / 3 > balanced_split
        balanced_split = 2balanced_split
    end
    while 2ℓ / 3 < balanced_split
        balanced_split = balanced_split / 2
    end
    return balanced_split
end

"""
    balanced_power_of_two_quarternary_split(ℓ) -> Float
    
Returns the value of `s ∈ [0, ℓ]` that gives the most balanced quarternary split of the segment `pq`, so `s` is a power-of-two
and `s ∈ [ℓ / 4, ℓ / 2]`.
"""
function balanced_power_of_two_quarternary_split(ℓ)
    balanced_split = one(ℓ)
    while ℓ / 4 > balanced_split
        balanced_split = 2balanced_split
    end
    while ℓ / 2 < balanced_split
        balanced_split = balanced_split / 2
    end
    return balanced_split
end

"""
    segment_vertices_adjoin_other_segments_at_acute_angle(tri::Triangulation, e) -> Int, Vertex

Determines if the vertices of a segment `e` of `tri` adjoin other segments at an acute angle.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `e`: The segment.

# Output
- `num_adjoin`: The number of vertices of `e` that adjoin other segments at an acute angle.
- `adjoin_vert`: The vertex of `e` that adjoins another segment at an acute angle if `num_adjoin == 1`, and `∅` otherwise.
"""
function segment_vertices_adjoin_other_segments_at_acute_angle(tri::Triangulation, e)
    u, v = edge_vertices(e)
    p, q = get_point(tri, u, v)
    I = integer_type(tri)
    other_segments = (I(∅), I(∅))
    # Check u first 
    for w in get_neighbours(tri, u)
        if w ≠ v && contains_segment(tri, u, w)
            r = get_point(tri, w)
            angle = opposite_angle(r, q, p)
            if is_acute(angle)
                other_segments = (w, I(∅))
                break
            end
        end
    end
    # Now check v  
    for w in get_neighbours(tri, v)
        if w ≠ u && w ≠ other_segments[1] && contains_segment(tri, v, w) # should be adjoining OTHER segments, not both the same segment. This case is just a triangle anyway. For the case of a triangle nestled in the corner of a small input angle, see is_triangle_nestled
            r = get_point(tri, w)
            angle = opposite_angle(r, p, q)
            if is_acute(angle)
                other_segments = (other_segments[1], w)
                break
            end
        end
    end
    # Now classify 
    a, b = other_segments
    num_adjoin, adjoin_vert = if a ≠ I(∅) && b ≠ I(∅)
        2, I(∅)
    elseif a ≠ I(∅) || b ≠ I(∅)
        1, a ≠ I(∅) ? u : v
    else
        0, I(∅)
    end
    return num_adjoin, adjoin_vert
end

"""
    is_encroached(tri::Triangulation, args::RefinementArguments, edge) -> Bool

Determines if a segment `edge` of `tri` is encroached upon.

See also [`encroaches_upon`](@ref).
"""
function is_encroached(tri::Triangulation, args::RefinementArguments, edge)
    !unoriented_edge_exists(tri, edge) && return true
    is_ghost_edge(edge) && return false
    u, v = edge_vertices(edge)
    p, q = get_point(tri, u, v)
    w = get_adjacent(tri, u, v)
    r = get_point(tri, w)
    !is_ghost_vertex(w) && encroaches_upon(p, q, r, args) && return true
    w′ = get_adjacent(tri, v, u)
    r′ = get_point(tri, w′)
    !is_ghost_vertex(w′) && encroaches_upon(p, q, r′, args) && return true
    return false
end

"""
    encroaches_upon(p, q, r, args::RefinementArguments) -> Bool

Determines if a point `r` encroaches upon a segment `pq` according to the [`RefinementArguments`](@ref).

See also [`point_position_relative_to_diametral_circle`](@ref) and [`point_position_relative_to_diametral_lens`](@ref).
"""
function encroaches_upon(p, q, r, args::RefinementArguments)
    if !args.use_lens
        d = point_position_relative_to_diametral_circle(p, q, r)
    else
        d = point_position_relative_to_diametral_lens(p, q, r, args.constraints.min_angle)
    end
    return !is_outside(d)
end

"""
    check_split_subsegment_precision(mx, my, p, q) -> Bool

Checks if there are precision issues related to the computed split position `(mx, my)` of a segment `(p, q)`, 
returning `true` if so and `false` otherwise.
"""
function check_split_subsegment_precision(mx, my, p, q)
    px, py = getxy(p)
    qx, qy = getxy(q)
    abs_err_flag_1 = check_absolute_precision(mx, px) && check_absolute_precision(my, py)
    abs_err_flag_2 = check_absolute_precision(mx, qx) && check_absolute_precision(my, qy)
    abs_err_flag = abs_err_flag_1 || abs_err_flag_2
    return abs_err_flag
end

"""
    split_subsegment!(tri::Triangulation, args::RefinementArguments, e) 

Splits a subsegment `e` of `tri` at a position determined by [`compute_split_position`](@ref); for curve-bounded domains, 
the position is determined by [`split_subcurve!`](@ref). After the split, [`assess_triangle_quality`](@ref) is used to 
find any new bad quality triangles. Before splitting, all free vertices in the segment's diametral circle 
are deleted using [`delete_free_vertices_around_subsegment!`](@ref).
"""
function split_subsegment!(tri::Triangulation, args::RefinementArguments, e)
    if is_curve_bounded(tri)
        return _split_subsegment_curve_bounded!(tri, args, e)
    else
        return _split_subsegment_piecewise_linear!(tri, args, e)
    end
end

"""
    _split_subsegment_piecewise_linear!(tri::Triangulation, args::RefinementArguments, e)

Splits a subsegment `e` of `tri` at a position determined by [`compute_split_position`](@ref) for piecewise linear domains. See
[`split_subsegment!`](@ref).
"""
function _split_subsegment_piecewise_linear!(tri::Triangulation, args::RefinementArguments, e)
    u, v = edge_vertices(e)
    p, q = get_point(tri, u, v)
    mx, my = compute_split_position(tri, args, e)
    # We can run into precision issues (e.g. I found this was needed when reproducing Figure 6.14 in the Delaunay 
    # Mesh Generation book). Let's ensure that the computed splitting point isn't just one of p or q.
    check_split_subsegment_precision(mx, my, p, q) && (pop!(args.midpoint_split_list); pop!(args.offcenter_split_list); return tri)
    empty!(args.events)
    delete_free_vertices_around_subsegment!(tri, args, e)
    push_point!(tri, mx, my)
    I = integer_type(tri)
    r = I(num_points(tri))
    complete_split_edge_and_legalise!(tri, u, v, r, Val(true), args.events)
    assess_added_triangles!(args, tri)
    return tri
end

"""
    _split_subsegment_curve_bounded!(tri::Triangulation, args::RefinementArguments, e)

Splits a subsegment `e` of `tri` at a position determined by [`split_subcurve!`](@ref) for curve-bounded domains. See 
[`split_subsegment!`](@ref). See also `_split_subsegment_curve_bounded_standard!` and `_split_subsegment_curve_bounded_small_angle!`,
as well as the original functions `_split_subsegment_curve_standard!` and `_split_subcurve_complex!`, respectively,
used during boundary enrichment.
"""
function _split_subsegment_curve_bounded!(tri::Triangulation, args::RefinementArguments, e)
    u, v = edge_vertices(e)
    enricher = get_boundary_enricher(tri)
    flag, apex, complex_id, _ = is_small_angle_complex_member(enricher, u, v)
    if !flag
        return _split_subsegment_curve_bounded_standard!(tri, args, e)
    else
        return _split_subsegment_curve_bounded_small_angle!(tri, args, e, apex, complex_id)
    end
end
function _split_subsegment_curve_bounded_standard!(tri::Triangulation, args::RefinementArguments, e)
    enricher = get_boundary_enricher(tri)
    i, j = edge_vertices(e)
    t, Δθ, ct = compute_split_position(enricher, i, j)
    cx, cy = getxy(ct)
    p, q = get_point(tri, i, j)
    precision_flag = check_split_subsegment_precision(cx, cy, p, q) || isnan(Δθ)
    if precision_flag
        return tri
    else
        empty!(args.events)
        delete_free_vertices_around_subsegment!(tri, args, e)
        push_point!(tri, cx, cy)
        I = integer_type(tri)
        r = I(num_points(tri))
        is_interior = is_segment(enricher, i, j)
        complete_split_edge_and_legalise!(tri, i, j, r, Val(true), args.events)
        split_edge!(enricher, i, j, r, Val(false), Val(false), is_interior)
        assess_added_triangles!(args, tri)
        push!(args.midpoint_split_list, r)
        return tri
    end
end
function _split_subsegment_curve_bounded_small_angle!(tri::Triangulation, args::RefinementArguments, e, apex, complex_id)
    # All members of the complex get split at their intersection with the same circular shell centered at apex 
    enricher = get_boundary_enricher(tri)
    complexes = get_small_angle_complexes(enricher, apex)
    complex = complexes[complex_id]
    points = get_points(tri)
    emin = get_minimum_edge_length(complex, points)
    circle_radius = balanced_power_of_two_ternary_split(emin)
    members = get_members(complex)
    p = get_point(points, apex)
    E = edge_type(tri)
    for (member_id, member) in enumerate(members)
        ct = _compute_split_position_complex(enricher, apex, member, circle_radius)
        next_edge = get_next_edge(member)
        q = get_point(points, next_edge)
        cx, cy = getxy(ct)
        precision_flag = check_split_subsegment_precision(cx, cy, p, q)
        precision_flag && continue
        empty!(args.events)
        push_point!(tri, cx, cy)
        I = integer_type(tri)
        r = I(num_points(tri))
        i, j = reorient_edge(enricher, apex, next_edge)
        e′ = construct_edge(E, i, j)
        delete_free_vertices_around_subsegment!(tri, args, e′)
        is_interior = is_segment(enricher, i, j)
        complete_split_edge_and_legalise!(tri, i, j, r, Val(true), args.events)
        split_edge!(enricher, i, j, r, Val(false), Val(false), is_interior)
        replace_next_edge!(enricher, apex, complex_id, member_id, r)
        assess_added_triangles!(args, tri)
        push!(args.midpoint_split_list, r)
    end
    return tri
end

"""
    delete_free_vertices_around_subsegment!(tri::Triangulation, args::RefinementArguments, e)

Deletes all free vertices (i.e., Steiner points) contained in the diametral circle of `e`.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `args::RefinementArguments`: The [`RefinementArguments`](@ref).
- `e`: The segment.

# Output
The free vertices are deleted from `tri` in-place.
"""
function delete_free_vertices_around_subsegment!(tri::Triangulation, args::RefinementArguments, e)
    for e′ in (e, reverse_edge(e))
        u, v = edge_vertices(e′)
        p, q = get_point(tri, u, v)
        w = get_adjacent(tri, e′)
        r = get_point(tri, w)
        while is_free(args, w) && !is_outside(point_position_relative_to_diametral_circle(p, q, r))
            delete_point!(tri, w; store_event_history = Val(true), event_history = args.events, args.rng)
            w = get_adjacent(tri, e′)
            r = get_point(tri, w)
        end
    end
    return tri
end

"""
    split_all_encroached_segments!(tri::Triangulation, args::RefinementArguments)

Splits all encroached segments of `tri` according to [`split_subsegment!`](@ref) until no more encroached segments exist in `args.queue`.
"""
function split_all_encroached_segments!(tri::Triangulation, args::RefinementArguments)
    while keep_splitting(tri, args)
        e, _ = popfirst_segment!(args.queue)
        !unoriented_edge_exists(tri, e) && continue
        split_subsegment!(tri, args, e)
    end
    return tri
end

"""
    is_triangle_seditious(tri::Triangulation, args, u, v, w, smallest_idx) -> Bool

Determines if a triangle `uvw` of `tri` is seditious according to the [`RefinementArguments`](@ref).

See also [`is_triangle_nestled`](@ref).

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `args`: The [`RefinementArguments`](@ref).
- `u, v, w`: The vertices of the triangle.
- `smallest_idx`: The index of the smallest edge of the triangle, so that `1` means `uv` is the smallest edge, `2` means `vw` is the smallest edge, and `3` means `wu` is the smallest edge.

# Output
- `flag`: Whether the triangle is seditious.

A triangle is seditious if it is nestled in the corner of a small input angle, or if it is nestled in the corner of a small input angle and the shortest edge is seditious. Here, 'small' is defined by 
`args.constraints.seditious_angle`.
"""
function is_triangle_seditious(tri::Triangulation, args, u, v, w, smallest_idx)
    i, j, k = make_shortest_edge_first(u, v, w, smallest_idx)
    # !is_segment_vertex(args, k) && return false
    cosine_scale² = cosd(args.constraints.seditious_angle)^2
    for (a, b) in ((i, j), (j, i))
        c = get_adjacent(tri, a, b)
        !is_segment_vertex(args, c) && continue # only apply seditious check to adjacent vertices which are input vertices. in this case, this means that c is a segment vertex (this check is thus the same as doing !is_free(args, c))
        (!contains_segment(tri, a, c) || !contains_segment(tri, b, c)) && continue
        ℓac² = edge_length_sqr(tri, a, c)
        ℓbc² = edge_length_sqr(tri, b, c)
        #if (a, b) == (j, i) # in this case, we don't know that (i, j) is the shortest edge of the triangle (j, i, c), since we only tested for (i, j, c) == (i, j, k)
        #    ℓab = dist(p, q)
        #    if (ℓac < ℓab || ℓbc < ℓab) 
        #        return false
        #    end
        #end
        # Actually, don't need the above check. It's just about asking if the shortest edge is seditious - an edge can be seditious without being the shortest, it just isn't normally checked is all.
        !check_seditious_precision(ℓac², ℓbc²) && continue
        if args.constraints.seditious_angle > 0.0
            p, q, r = get_point(tri, a, b, c)
            px, py = getxy(p)
            qx, qy = getxy(q)
            rx, ry = getxy(r)
            ddot = (px - rx) * (qx - rx) + (py - ry) * (qy - ry)
            ddot < 0 && continue
            ℓacℓbc² = ℓac² * ℓbc²
            ddot^2 > cosine_scale² * ℓacℓbc² && return true
        else
            return true
        end
    end
    # If we did want to actually check that they're on the same concentric shell, 
    # we would compute the length of each segment and see if log2 of those lengths are 
    # both about the same and both approximately integers.
    return false
end

"""
    check_seditious_precision(ℓrp, ℓrq) -> Bool

Checks if there are precision issues related to the seditiousness of a triangle, returning `true` if so and `false` otherwise.
"""
function check_seditious_precision(ℓrp, ℓrq)
    ratio_flag = check_ratio_precision(ℓrp, ℓrq)
    rel_err_flag = check_relative_precision(ℓrp, ℓrq)
    abs_err_flag = check_absolute_precision(ℓrp, ℓrq)
    zero_flag = iszero(ℓrp) || iszero(ℓrq)
    return ratio_flag || rel_err_flag || abs_err_flag || zero_flag
end

"""
    is_triangle_nestled(tri::Triangulation, T, idx) -> Bool

Determines if a triangle `T` of `tri` is nestled in the corner of a small input angle.

See also [`is_triangle_seditious`](@ref).

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `T`: The triangle.
- `idx`: The index of the smallest edge of the triangle, so that `1` means `uv` is the smallest edge, `2` means `vw` is the smallest edge, and `3` means `wu` is the smallest edge.

# Output
- `flag`: Whether the triangle is nestled in the corner of a small input angle.

A triangle is nestled in the corner of a small input angle if it is nestled in the corner of a small input angle and the shortest edge is seditious. The size of the angle is not checked by this function,
and is instead determined by [`assess_triangle_quality`](@ref).
"""
function is_triangle_nestled(tri::Triangulation, T, idx) # nestled: a triangle is nestled if it is nestled in the corner of a small input angle. This function does NOT check the angle.
    u, v, w = triangle_vertices(T)
    i, j, k = make_shortest_edge_first(u, v, w, idx) # small angles are opposite the shortest edge
    e_ki_seg = contains_segment(tri, k, i)
    e_kj_seg = contains_segment(tri, k, j)
    return e_ki_seg && e_kj_seg
end

"""
    assess_triangle_quality(tri::Triangulation, args::RefinementArguments, T) -> Float64, Bool

Assesses the quality of a triangle `T` of `tri` according to the [`RefinementArguments`](@ref).

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `args::RefinementArguments`: The [`RefinementArguments`](@ref).
- `T`: The triangle.

# Output
- `ρ`: The radius-edge ratio of the triangle.
- `flag`: Whether the triangle is bad quality.

A triangle is bad quality if it does not meet the area constraints, violates the custom constraint, or if it is skinny but neither seditious or nestled.
"""
function assess_triangle_quality(tri::Triangulation, args::RefinementArguments, T)
    haskey(args.queue, T) && return args.queue[T], true
    if is_ghost_triangle(T)
        T = replace_ghost_triangle_with_boundary_triangle(tri, T)
    end
    # First, get the radius-edge ratio.
    u, v, w = triangle_vertices(T)
    p, q, r = get_point(tri, u, v, w)
    ℓmin², ℓmed², ℓmax², idx = squared_triangle_lengths_and_smallest_index(p, q, r)
    # A² = squared_triangle_area(ℓmin², ℓmed², ℓmax²)
    A² = squared_triangle_area(p, q, r) # Don't want to use the above approach since we want to allow for squared_triangle_area to try and consider more accurate methods if necessary
    A = sqrt(A²)
    cr = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    ρ = triangle_radius_edge_ratio(cr, sqrt(ℓmin²))
    # Next, check the area constraints.
    A < args.constraints.min_area && return ρ, false
    A > args.constraints.max_area && return ρ, true
    # Now check the custom constraints.
    violates_custom_constraint(args.constraints, tri, T) && return ρ, true
    # Now check if the triangle is skinny and can be refined.
    is_skinny = ρ > args.constraints.max_radius_edge_ratio
    !is_skinny && return ρ, false
    is_seditious = is_triangle_seditious(tri, args, u, v, w, idx)
    is_nestled = is_triangle_nestled(tri, T, idx)
    # Now we know that the triangle is skinny, so we can return depending on whether it is safe to split it 
    return ρ, !(is_seditious || is_nestled)
end

"""
    assess_added_triangles!(args::RefinementArguments, tri::Triangulation)

Assesses the quality of all triangles in `args.events.added_triangles` according to [`assess_triangle_quality`](@ref), and enqueues any bad quality triangles into `args.queue`.
"""
function assess_added_triangles!(args::RefinementArguments, tri::Triangulation)
    E = edge_type(tri)
    for T in each_added_triangle(args.events)
        if !is_ghost_triangle(T)
            u, v, w = triangle_vertices(T)
            if get_adjacent(tri, u, v) == w # It's not guaranteed that each triangle in events.added_triangles _actually_ still exists, edge flipping could have deleted it
                ρ, flag = assess_triangle_quality(tri, args, T)
                flag && !haskey(args.queue, T) && (args.queue[T] = ρ)
                for e in triangle_edges(T)
                    ee = construct_edge(E, initial(e), terminal(e)) # Need to convert e to the correct edge type since triangle_edges returns Tuples
                    if contains_segment(tri, ee) && !haskey(args.queue, ee) && is_encroached(tri, args, ee)
                        ℓ² = edge_length_sqr(tri, ee)
                        args.queue[ee] = ℓ²
                    end
                end
            end
        end
    end
    return args
end
