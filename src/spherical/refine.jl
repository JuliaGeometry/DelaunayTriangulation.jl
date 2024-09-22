function _quality_statistics(tri::SphericalTriangulation, args::RefinementArguments, T)
    u, v, w = triangle_vertices(T)
    A = triangle_area(tri, T)
    cr = triangle_circumradius(tri, T, A)
    points = get_points(tri)
    p, q, r = _getindex(points, u), _getindex(points, v), _getindex(points, w) # can't do get_point since DelaunayTriangulation projects it
    d1, d2, d3 = spherical_distance(p, q), spherical_distance(q, r), spherical_distance(r, p)
    mind = min(d1, d2, d3)
    idx = mind == d1 ? 1 : mind == d2 ? 2 : 3
    ρ = triangle_radius_edge_ratio(cr, mind)
    return ρ, A, idx
end

function triangle_sink(tri::SphericalTriangulation, args...) # just avoid computing this in statistics()
    V = number_type(tri)
    return (V(NaN), V(NaN))
end 

function RefinementArguments(
    tri::SphericalTriangulation;
    min_angle = 0.0, 
    max_angle = 180.0,
    min_area = number_type(tri)(4π / 1.0e9),
    max_area = typemax(number_type(tri)),
    max_points = max(1_000, num_solid_vertices(tri))^2,
    seditious_angle = -Inf,
    custom_constraint = (_tri, T) -> false,
    use_circumcenter = true,
    use_lens = true,
    steiner_scale = 0.999,
    rng = Random.default_rng(),
    concavity_protection = false,
    predicates::AbstractPredicateKernel = AdaptiveKernel()
)
    # Ignoring angle constraints. Spherical refinement uses centroids instead of Ruppert's algorithm.
    has_ghosts = has_ghost_triangles(tri)
    if !has_ghosts 
        throw(ArgumentError("Spherical triangulations must have ghost triangles for refinement to work. Please avoid using solidify! for this."))
    end
    constraints = RefinementConstraints(;
        min_angle = 0.0,
        max_angle = 180.0,
        min_area,
        max_area,
        max_points,
        seditious_angle=-Inf,
        custom_constraint)
    queue, events, segment_list, segment_vertices, midpoint_split_list, offcenter_split_list, min_steiner_vertex = _build_queues(tri)
    return RefinementArguments(
        queue,
        constraints,
        events,
        min_steiner_vertex,
        segment_list,
        segment_vertices,
        midpoint_split_list,
        offcenter_split_list,
        use_circumcenter,
        use_lens,
        steiner_scale,
        false,
        true,
        rng,
        concavity_protection,
        predicates
    )
end

function get_steiner_point(tri::SphericalTriangulation, args::RefinementArguments, T)
    c = triangle_centroid(tri, T)
    check_steiner_point_precision(tri, T, c) && return Certificate.PrecisionFailure, c
    return Certificate.None, c
end

function enqueue_triangle!(tri::SphericalTriangulation, queue::RefinementQueue, _, T)
    A = triangle_area(tri, T) # Instead of using radius-edge ratios, spherical triangulations use areas since we are doing centroid refinement
    return setindex!(queue, A, T)
end

function finalise!(tri::Triangulation, args::RefinementArguments)
    #unlock_convex_hull!(tri; reconstruct = true)
    #lock_convex_hull!(tri; rng = args.rng, predicates = args.predicates)
    convex_hull!(tri; reconstruct = false, predicates = args.predicates)
    return tri
end

function _each_solid_triangle(tri::SphericalTriangulation)
    return each_triangle(tri)
end

function _is_ghost_triangle(_::SphericalTriangulation, T)
    return false
end

function locate_steiner_point(tri::SphericalTriangulation, args::RefinementArguments, T, c)
    return T, Certificate.Inside # centroids are always inside the triangle
end