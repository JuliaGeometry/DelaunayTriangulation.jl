"""
    refine!(tri::Triangulation;
        min_area=1e-9get_total_area(tri),
        max_area=typemax(number_type(tri)),
        max_radius_edge_ratio=nothing,
        max_points=typemax(Int64),
        min_angle=nothing,
        rng::AbstractRNG=Random.default_rng(),
        maxiters=100,
        lock_convex_hull=!has_boundary_nodes(tri),
        exterior_curve_index=1
        )

Refine the triangulation `tri` according to the given refinement targets.

# Arguments 
- `tri::Triangulation`: The triangulation to refine.

# Keyword Arguments
- `min_area=1e-9get_total_area(tri)`

The minimum permissible area of a triangle. Any circumcenter insertion attempts into a triangle with area below this minimum will be prevented, and the triangle will never be queued for refinement. This is useful for preventing the mesh from becoming too fine in certain regions and preventing convergence.
- `max_area=typemax(number_type(tri))`

The maximum permissible area of a triangle.
- `max_radius_edge_ratio=nothing`

The maximum allowed radius-edge ratio (the ratio of the circumradius to the shortest edge length) of a triangle. 
This ratio is related to a `min_angle` constraint by the formula `min_angle = asin(0.5 / max_radius_edge_ratio)`.
If `max_radius_edge_ratio=nothing`, then it is computed from the provided `min_angle` if it is not nothing, 
or we set `max_radius_edge_ratio = 1`. This can also be a function of the form `f(T, p, q, r, ρ)`,
where `T` is the triangle's indices, `p, q, r` its coordinates, and `ρ` is the triangle's radius-edge 
ratio. This function should return `true` if it should be refined, and `false` otherwise.

The maximum radius-edge ratio cannot be less than `1/√3`. If it is, it is replaced with `1/√3`.
- `min_angle=nothing`

The minimum allowed angle of a triangle, given in degrees. If `max_radius_edge_ratio` is not `nothing`, then this is ignored. Otheriwse, 
we set `max_radius_edge_ratio = 1 / (2sind(min_angle))`.

If `min_angle > 33.9°, a warning is given as convergence of the algorithm may struggle to converge in this case. Additionally, if `min_angle < 0.0` or `min_angle > 60.00005`, then `min_angle` is replaced with `30°` (`max_radius_edge_ratio = 1`).
You should also be cautious that, for some domains (e.g. ones with many small angles), an even smaller angle might be needed to accommodate floating point arithmetic.

!!! note 

    A minimum angle constraint of `θₘᵢₙ°` corresponds to a maximum angle constraint of `180° - 2θₘᵢₙ`, e.g. `θₘᵢₙ° = 30°` gives `θₘₐₓ° = 120°`.

- `max_points=typemax(Int64)`

How many points are allowed to be in the triangulation before terminating.
- `rng::AbstractRNG=Random.default_rng()`

The random number generator to use.
- `maxiters=5000`

The maximum number of iterations to perform.
- `lock_convex_hull=!has_boundary_nodes(tri)`

Whether to lock the convex hull of the triangulation. If `true`, then the convex hull edges are put into the constrained edges of the triangulation. This 
is only relevant if the triangulation has no boundary nodes.

!!! danger 

    Not enabling this option can lead to the triangulation being refined in a way that 
    continually adds points that go out to infinity. Do so at your own peril.

- `exterior_curve_index=1`

The curve (or curves) corresponding to the outermost boundary. 

# Outputs 

The `TriangulationStatistics` of `tri` are returned. The triangulation is refined in-place.

See also [`RefinementTargets`](@ref).

!!! note 

    It is not strictly guaranteed that the refinement targets will be met. For example, precision 
    issues may cause small angles from near-degenerate triangles. Despite this, the algorithm will
    still typically give no angle larger than `180° - 2θ°`, where `θ° = asind[1 / (2max_radius_edge_ratio)]`
    is the minimum angle constraint associated with `max_radius_edge_ratio`. Moreover, domains with 
    small angles (between constrained segments) will cause issues, forcing small angles to be 
    introduced into the triangulation to accommodate them, but only between the constrained 
    segments with small angles. If you are finding issues with small angles, try increasing the
    `max_radius_edge_ratio` (or reduce `min_angle`). 

    Basically, be aware that arithmetic is not perfected and floating point arithmetic 
    only has so much precision to offer.
"""
function refine!(tri::Triangulation;
    min_area=1e-9get_total_area(tri),
    max_area=typemax(number_type(tri)),
    max_radius_edge_ratio=nothing,
    max_points=typemax(Int64),
    min_angle=nothing,
    rng::AbstractRNG=Random.default_rng(),
    maxiters=100_000,
    lock_convex_hull=!has_boundary_nodes(tri),
    exterior_curve_index=1
)
    tri, queue, events, targets, has_ghosts, segment_list = initialise_refine(tri; min_area, max_area, max_radius_edge_ratio, max_points, min_angle, lock_convex_hull)
    _refine_all!(tri, queue, events, targets, segment_list, exterior_curve_index, maxiters, rng)
    finalise_refine(tri, has_ghosts, lock_convex_hull)
    stats = statistics(tri)
    return stats
end

function initialise_refine(tri::Triangulation;
    min_area=zero(number_type(tri)),
    max_area=typemax(number_type(tri)),
    max_radius_edge_ratio=nothing,
    max_points=typemax(Int64),
    min_angle=nothing,
    lock_convex_hull=!has_boundary_nodes(tri))
    has_ghosts = has_ghost_triangles(tri)
    !has_ghosts && add_ghost_triangles!(tri)
    lock_convex_hull && lock_convex_hull!(tri)
    targets = RefinementTargets(; min_area, max_area, max_radius_edge_ratio, max_points, min_angle)
    queue = initialise_refinement_queue(tri, targets)
    events = initialise_event_history(tri)
    E = edge_type(tri)
    segment_list = Set{E}()
    for e in each_constrained_edge(tri)
        add_edge!(segment_list, e)
    end
    return tri, queue, events, targets, has_ghosts, segment_list
end

function finalise_refine(tri::Triangulation, has_ghosts, lock_convex_hull)
    !has_ghosts && delete_ghost_triangles!(tri)
    if lock_convex_hull
        idx = get_convex_hull_indices(tri) # these indices could have changed after refinement, so let's update the convex hull, which is now get_boundary_nodes(tri)
        bn = get_boundary_nodes(tri)
        resize!(idx, num_boundary_edges(bn) + 1)
        copyto!(idx, bn)
        unlock_convex_hull!(tri)
    end
    return nothing
end

function _refine_all!(tri::Triangulation, queue::RefinementQueue, events::InsertionEventHistory, targets::RefinementTargets, segment_list, exterior_curve_index, maxiters, rng::AbstractRNG=Random.default_rng())
    iters = 1
    split_all_encroached_segments!(tri, queue, events, targets, segment_list)
    while !isempty(queue) && !compare_points(targets, num_points(tri)) && iters ≤ maxiters
        ρ = peek_triangle_ρ(queue)
        T = triangle_dequeue!(queue)
        iters = _refine_itr!(tri, queue, events, targets, T, ρ, segment_list, exterior_curve_index, iters, rng)
    end
    return nothing
end

function _refine_itr!(tri::Triangulation, queue::RefinementQueue, events::InsertionEventHistory, targets::RefinementTargets, T, ρ, segment_list, exterior_curve_index, iters, rng::AbstractRNG=Random.default_rng())
    u, v, w = indices(T)
    if !is_ghost_triangle(T) && get_adjacent(tri, u, v) == w
        iters += 1
        success = split_triangle!(tri, queue, events, T, exterior_curve_index, rng)
        if is_encroachment_failure(success)
            split_all_encroached_segments!(tri, queue, events, targets, segment_list)
            triangle_enqueue!(queue, T, ρ) # Re-enqueue triangle that the circumcenter came from 
        elseif is_successful_insertion(success)
            assess_added_triangles!(tri, queue, events, targets)
        end
    end
    return iters
end