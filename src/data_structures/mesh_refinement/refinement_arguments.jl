"""
    RefinementArguments{Q,C,H,I,E,T,R,P<:AbstractPredicateKernel}

A struct for storing arguments for mesh refinement.

# Fields 
- `queue::Q`: The [`RefinementQueue`](@ref).
- `constraints::C`: The [`RefinementConstraints`](@ref).
- `events::H`: The [`InsertionEventHistory`](@ref).
- `min_steiner_vertex::I`: The minimum vertex of a Steiner point. All vertices greater than or equal to this can be considered as Steiner vertices. 
- `segment_list::Set{E}`: The set of segments in the triangulation before refinement.
- `segment_vertices::Set{I}`: Set of vertices that are vertices of segments in `segment_list`.
- `midpoint_split_list::Set{I}`: Set of vertices that are centre-splits of encroached edges.
- `offcenter_split_list::Set{I}`: Set of vertices that are off-centre splits of encroached edges.
- `use_circumcenter::Bool`: Whether to use circumcenters for Steiner points, or the more general approach of [Erten and Üngör (2009)](https://doi.org/10.1109/ISVD.2009.32).
- `use_lens::Bool`: Whether to use diametral lens (`true`) or diametral circles (`false`) for the defining encroachment.
- `steiner_scale::T`: The factor by which to scale the Steiner points closer to the triangle's shortest edge.
- `locked_convex_hull::Bool`: Whether the convex hull of the triangulation had to be locked for refinement. 
- `had_ghosts::Bool`: Whether the triangulation initially had ghost triangles or not.
- `rng::R`: The random number generator.
- `concavity_protection::Bool`: Whether to use concavity protection or not for [`find_triangle`](@ref). Most likely not needed, but may help in pathological cases.
- `predicates::P<:AbstractPredicateKernel`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Constructors 
In addition to the default constructor, we provide 

    RefinementArguments(tri::Triangulation; kwargs...)

for constructing this struct. This constructor will lock the convex hull and add ghost triangles to `tri` if
needed ([`refine!`](@ref) will undo these changes once the refinement is finished))
"""
struct RefinementArguments{Q,C,H,I,E,R,T,P<:AbstractPredicateKernel}
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

"""
    RefinementArguments(tri::Triangulation; kwargs...) -> RefinementArguments

Initialises the [`RefinementArguments`](@ref) for the given [`Triangulation`](@ref), `tri`. The `kwargs...`
match those from [`refine!`](@ref).

!!! warning "Mutation"

    If `tri` has no ghost triangles, it will be mutated so that it has them. Similarly, if the triangulation 
    has no constrained boundary, then the convex hull will be locked so that it is treated as a constrained 
    boundary. These changes will be undone in [`refine!`](@ref) once the refinement is finished.
"""
function RefinementArguments(tri::Triangulation;
    min_angle=30.0,
    max_angle=180.0,
    min_area=number_type(tri)(get_area(tri) / 1e9),
    max_area=typemax(number_type(tri)),
    max_points=max(1_000, num_solid_vertices(tri))^2,
    seditious_angle=20.0,
    custom_constraint=(_tri, T) -> false,
    use_circumcenter=true, # TODO: When we implement generalised Steiner points, change this default to FALSE.
    use_lens=true,
    steiner_scale=0.999,
    rng=Random.default_rng(),
    concavity_protection=false,
    predicates::AbstractPredicateKernel=AdaptiveKernel())
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
        custom_constraint
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
        predicates
    )
end

"""
    is_free(args::RefinementArguments, u) -> Bool

Returns `true` if `u` is a free vertex, and `false` otherwise. A free vertex is a Steiner vertex (meaning a non-input vertex) that is not part of a segment or subsegment.
"""
is_free(args::RefinementArguments, u) = !(u < args.min_steiner_vertex) && !(is_midpoint_split(args, u) || is_offcenter_split(args, u))

"""
    keep_iterating(tri::Triangulation, args::RefinementArguments) -> Bool

Returns `true` if the refinement should continue, and `false` otherwise. The check is based on 
whether the [`RefinementQueue`](@ref) is empty or not, and whether the number of points in the triangulation
is less than or equal to the maximum number of points allowed by the [`RefinementConstraints`](@ref).
"""
keep_iterating(tri::Triangulation, args::RefinementArguments) = !isempty(args.queue) && num_solid_vertices(tri) < args.constraints.max_points

"""
    keep_splitting(tri::Triangulation, args::RefinementArguments) -> Bool

Returns `true` if more encroached segments need to be split, and `false` otherwise. The check is based on
whether the segment queue in the [`RefinementQueue`](@ref) of `args` is empty or not, and whether the number of points in the triangulation
is less than or equal to the maximum number of points allowed by the [`RefinementConstraints`](@ref) in `args`.
"""
keep_splitting(tri::Triangulation, args::RefinementArguments) = has_segments(args.queue) && num_solid_vertices(tri) < args.constraints.max_points

"""
    is_midpoint_split(args::RefinementArguments, u) -> Bool

Returns `true` if `u` is a midpoint of an encroached segment, and `false` otherwise.
"""
is_midpoint_split(args::RefinementArguments, u) = u ∈ args.midpoint_split_list

"""
    is_offcenter_split(args::RefinementArguments, u) -> Bool

Returns `true` if `u` is an off-centre split of an encroached segment, and `false` otherwise.
"""
is_offcenter_split(args::RefinementArguments, u) = u ∈ args.offcenter_split_list

"""
    is_subsegment(args::RefinementArguments, u, v) -> Bool

Returns `true` if the edge `(u, v)` is a subsegment, and `false` otherwise.
"""
is_subsegment(args::RefinementArguments, e) = !contains_unoriented_edge(e, args.segment_list)

"""
    is_segment_vertex(args::RefinementArguments, u) -> Bool

Returns `true` if `u` is a vertex of a segment, and `false` otherwise. Note that this excludes vertices of subsegments.
"""
is_segment_vertex(args::RefinementArguments, u) = u ∈ args.segment_vertices
