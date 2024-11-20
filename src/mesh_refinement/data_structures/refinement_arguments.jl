struct RefinementArguments{Q, C, H, I, E, R, T, P <: AbstractPredicateKernel}
    queue::Q
    constraints::C
    events::H
    min_steiner_vertex::I
    segment_list::Set{E}
    segment_vertices::Set{I}
    midpoint_split_list::Set{I}
    offcenter_split_list::Set{I}
    use_circumcenter::Bool
    use_lens::Bool
    steiner_scale::T
    locked_convex_hull::Bool
    had_ghosts::Bool
    rng::R
    concavity_protection::Bool
    predicates::P
end

function RefinementArguments(
        tri::Triangulation;
        min_angle = 30.0,
        max_angle = 180.0,
        min_area = number_type(tri)(get_area(tri) / 1.0e9),
        max_area = typemax(number_type(tri)),
        max_points = max(1_000, num_solid_vertices(tri))^2,
        seditious_angle = 20.0,
        custom_constraint = (_tri, T) -> false,
        use_circumcenter = true, # TODO: When we implement generalised Steiner points, change this default to FALSE.
        use_lens = true,
        steiner_scale = 0.999,
        rng = Random.default_rng(),
        concavity_protection = false,
        predicates::AbstractPredicateKernel = AdaptiveKernel(),
    )
    if !use_circumcenter
        throw(ArgumentError("Generalised Steiner points are not yet implemented."))
    end
    has_ghosts = has_ghost_triangles(tri)
    !has_ghosts && add_ghost_triangles!(tri)
    lock_convex_hull = !has_boundary_nodes(tri)
    lock_convex_hull && lock_convex_hull!(tri; rng, predicates)
    constraints = RefinementConstraints(;
        min_angle,
        max_angle,
        min_area,
        max_area,
        max_points,
        seditious_angle,
        custom_constraint,
    )
    queue = RefinementQueue(tri)
    events = InsertionEventHistory(tri)
    E = edge_type(tri)
    I = integer_type(tri)
    segment_list = Set{E}()
    segment_vertices = Set{I}()
    sizehint!(segment_list, num_edges(get_all_segments(tri)))
    sizehint!(segment_vertices, num_edges(get_all_segments(tri)) ÷ 2)
    for e in each_segment(tri)
        push!(segment_list, e)
        push!(segment_vertices, initial(e), terminal(e))
    end
    midpoint_split_list = Set{I}()
    offcenter_split_list = Set{I}()
    sizehint!(midpoint_split_list, num_edges(get_all_segments(tri)))
    sizehint!(offcenter_split_list, 2num_edges(get_all_segments(tri)))
    min_steiner_vertex = I(num_points(tri) + 1)
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
        lock_convex_hull,
        has_ghosts,
        rng,
        concavity_protection,
        predicates,
    )
end

@inline is_free(args::RefinementArguments, u) = !(u < args.min_steiner_vertex) && !(is_midpoint_split(args, u) || is_offcenter_split(args, u))
@inline keep_iterating(tri::Triangulation, args::RefinementArguments) = !isempty(args.queue) && num_solid_vertices(tri) < args.constraints.max_points
@inline keep_splitting(tri::Triangulation, args::RefinementArguments) = has_segments(args.queue) && num_solid_vertices(tri) < args.constraints.max_points
@inline is_midpoint_split(args::RefinementArguments, u) = u ∈ args.midpoint_split_list
@inline is_offcenter_split(args::RefinementArguments, u) = u ∈ args.offcenter_split_list
@inline is_subsegment(args::RefinementArguments, e) = !contains_unoriented_edge(e, args.segment_list)
@inline is_segment_vertex(args::RefinementArguments, u) = u ∈ args.segment_vertices
