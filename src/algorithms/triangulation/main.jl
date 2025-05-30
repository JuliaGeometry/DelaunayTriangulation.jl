"""
    postprocess_triangulate!(tri; delete_ghosts=false, delete_empty_features=true, recompute_representative_points=true)

Postprocesses the triangulation `tri` after it has been constructed using [`triangulate`](@ref). This includes:

- Deleting ghost triangles using [`delete_ghost_triangles!`](@ref) if `delete_ghosts` is `true`.
- Clearing empty features using [`clear_empty_features!`](@ref) if `delete_empty_features` is `true`.
- Recomputing the representative points using [`compute_representative_points!`](@ref) if `recompute_representative_points` is `true`.
"""
function postprocess_triangulate!(tri; delete_ghosts = false, delete_empty_features = true, recompute_representative_points = true)
    delete_ghosts && delete_ghost_triangles!(tri)
    delete_empty_features && clear_empty_features!(tri)
    recompute_representative_points && compute_representative_points!(tri)
    return tri
end

"""
    triangulate(points; segments=nothing, boundary_nodes=nothing, kwargs...) -> Triangulation

Computes the Delaunay triangulation of `points`, and then the constrained Delaunay triangulation if any of `segments` and `boundary_nodes` are not `nothing`.

# Arguments 
- `points`: The points to triangulate. (This might get mutated for curve-bounded domains.)

# Keyword Arguments 
For the keyword arguments below, you may like to review the extended help as some of the arguments carry certain warnings.

- `segments=nothing`: The segments to include in the triangulation. If `nothing`, then no segments are included.
- `boundary_nodes=nothing`: The boundary nodes to include in the triangulation. If `nothing`, then no boundary nodes are included, and the convex hull of `points` remains as the triangulation. These boundary nodes 
   should match the specification given in [`check_args`](@ref) if a boundary is provided as a set of vertices, meaning the boundary is a piecewise linear curve. To specify a curve-bounded domain, you should 
   follow the same specification, but use [`AbstractParametricCurve`](@ref)s to fill out the vector, and any piecewise linear section should still be provided as a sequence of vertices. 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `weights=ZeroWeight{number_type(points)}()`: The weights to use for the triangulation. By default, the triangulation is unweighted. The weights can also be provided as a vector, with the `i`th weight referring to the `i`th vertex, or more generally any object that defines [`get_weight`](@ref). 
- `IntegerType=Int`: The integer type to use for the triangulation. This is used for representing vertices.
- `EdgeType=isnothing(segments) ? NTuple{2,IntegerType} : (edge_type ∘ typeof)(segments)`: The edge type to use for the triangulation. 
- `TriangleType=NTuple{3,IntegerType}`: The triangle type to use for the triangulation.
- `EdgesType=isnothing(segments) ? Set{EdgeType} : typeof(segments)`: The type to use for storing the edges of the triangulation.
- `TrianglesType=Set{TriangleType}`: The type to use for storing the triangles of the triangulation.
- `randomise=true`: Whether to randomise the order in which the points are inserted into the triangulation. This is done using [`get_insertion_order`](@ref).
- `delete_ghosts=false`: Whether to delete the ghost triangles after the triangulation is computed. This is done using [`delete_ghost_triangles!`](@ref).
- `delete_empty_features=true`: Whether to delete empty features after the triangulation is computed. This is done using [`clear_empty_features!`](@ref).
- `try_last_inserted_point=true`: Whether to try the last inserted point first when inserting points into the triangulation. 
- `skip_points=()`: The points to skip when inserting points into the triangulation. 
   Note that, for curve-bounded domains, `skip_points` is ignored when using [`enrich_boundary!`](@ref).
- `num_sample_rule=default_num_samples`: A function mapping a number of points `n` to a number of samples `m` to use for sampling the initial points during the point location step of the algorithm within [`find_triangle`](@ref).
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator.
- `insertion_order::Vector=get_insertion_order(points, randomise, skip_points, IntegerType, rng)`: The insertion order to use for inserting points into the triangulation. This is ignored if you are defining a curve-bounded domain.
- `recompute_representative_points=true`: Whether to recompute the representative points after the triangulation is computed. This is done using [`compute_representative_points!`](@ref). 
- `delete_holes=true`: Whether to delete holes after the triangulation is computed. This is done using [`delete_holes!`](@ref).
- `check_arguments=true`: Whether to check the arguments `points` and `boundary_nodes` are valid. This is done using [`check_args`](@ref).
- `polygonise_n=4096`: Number of points to use for polygonising the boundary when considering the polygon hierarchy for a curve-bounded domain using [`polygonise`](@ref). See [`triangulate_curve_bounded`](@ref).
- `coarse_n=0`: Number of points to use for initialising a curve-bounded domain. See [`triangulate_curve_bounded`](@ref). (A value of `0` means the number of points is chosen automatically until the diametral circles of all edges are empty.)

# Outputs
- `tri::Triangulation`: The triangulation.

!!! note "Type stability"

    The output from this function is currently not type stable. In particular, the inferred type is only `Triangulation` without any other information.
    If you are depending on the output from `triangulate` inside some other function, you should consider putting the output behind a function barrier;
    information about using function barriers is given [here](https://docs.julialang.org/en/v1/manual/performance-tips/#kernel-functions).

# Extended help

Here are some warnings to consider for some of the arguments.

- `points`

!!! warning "Mutation"

    For curve-bounded domains, `points` may get mutated to include the endpoints of the provided curves, and when inserting 
    Steiner points to split segments or refine boundaries.

!!! danger "Floating point precision"

    If your points are defined using non-`Float64` coordinates, you may run into precision issues that 
    lead to problems with robustness. The consequences of this could be potentially catastrophic, leading to 
    infinite loops for example. If you do encounter such issues, consider converting your coordinates to `Float64`.

- `segments`

!!! warning "Segments outside of the domain"

    When segments are outside of the domain, are if they are not entirely contained with the domain, you may run into issues - especially for curve-bounded domains.
    It is your responsibility to ensure that the segments are contained within the domain.

!!! warning "Mutation"

    The `segments` may get mutated in two ways: (1) Segments may get rotated so that `(i, j)` becomes `(j, i)`. (2) If there are
    segments that are collinear with other segments, then they may get split into chain of non-overlapping connecting segments (also see below). For 
    curve-bounded domains, segments are also split so that no subsegment's diametral circle contains any other point.

!!! warning "Intersecting segments"

    Currently, segments that intersect in their interiors (this excludes segments that only intersect by sharing a vertex) cause problems for triangulating. While there 
    is some support for collinear segments that lie on top of each other (they get split automatically), this is not the case for segments that intersect in their interiors.
    Moreover, this automatic splitting should not be heavily relied upon, and for curve-bounded domains you should not rely on it at all as it causes problems during the enrichment phase
    from [`enrich_boundary!`](@ref).

- `boundary_nodes`

!!! warning "Points outside of boundary curves"

    While for standard domains with piecewise linear boundaries (or no boundaries) it is fine for points to be 
    outside of the domain (they just get automatically deleted if needed), they may cause problems for curve-bounded domains.
    Please ensure that all your points are inside the curve-bounded domain if you are providing curves in `boundary_nodes`.

!!! warning "Aliasing"

    For curve-bounded domains, the `boundary_nodes` in the resulting [`Triangulation`](@ref) will not be aliased with the input boundary nodes.

!!! note "Refinement"

    For curve-bounded domains, note that the triangulation produced from this function is really just an initial coarse discretisation of the true curved boundaries. You will need to 
    refine further, via [`refine!`](@ref), to improve the discretisation, or increase `coarse_n` below. See also [`polygonise`](@ref) for a more direct approach to discretising a boundary (which 
    might not give as high-quality meshes as you can obtain from [`refine!`](@ref) though, note).
"""
function triangulate(
        points::P;
        segments = nothing,
        boundary_nodes = nothing,
        predicates::AbstractPredicateKernel = AdaptiveKernel(),
        weights = ZeroWeight{number_type(P)}(),
        IntegerType::Type{I} = Int,
        EdgeType::Type{E} = isnothing(segments) ? NTuple{2, IntegerType} : edge_type(typeof(segments)),
        TriangleType::Type{V} = NTuple{3, IntegerType},
        EdgesType::Type{Es} = isnothing(segments) ? Set{EdgeType} : typeof(segments),
        TrianglesType::Type{Ts} = Set{TriangleType},
        randomise = true,
        delete_ghosts = false,
        delete_empty_features = true,
        try_last_inserted_point = true,
        skip_points = (),
        num_sample_rule::M = default_num_samples,
        rng::Random.AbstractRNG = Random.default_rng(),
        insertion_order::Vector = get_insertion_order(points, randomise, skip_points, IntegerType, rng),
        recompute_representative_points = true,
        delete_holes = true,
        check_arguments = true,
        full_polygon_hierarchy::H = nothing,
        boundary_enricher = nothing,
        boundary_curves = (),
        polygonise_n = 4096,
        coarse_n = 0,
        enrich = false,
    ) where {P, I, E, V, Es, Ts, M, H}
    check_config(points, weights, segments, boundary_nodes, predicates)
    if enrich || (isempty(boundary_curves) && is_curve_bounded(boundary_nodes)) # If boundary_curves is not empty, then we are coming from triangulate_curve_bounded
        return triangulate_curve_bounded(points; segments, boundary_nodes, predicates, weights, IntegerType, EdgeType, TriangleType, EdgesType, TrianglesType, randomise, delete_ghosts, delete_empty_features, try_last_inserted_point, skip_points, num_sample_rule, rng, insertion_order, recompute_representative_points, delete_holes, check_arguments, polygonise_n, coarse_n)
    end
    if isnothing(full_polygon_hierarchy)
        full_polygon_hierarchy = construct_polygon_hierarchy(points, boundary_nodes; IntegerType)
    end
    skip_points_set = Set{IntegerType}(skip_points)
    n = length(skip_points_set)
    check_arguments && check_args(points, boundary_nodes, full_polygon_hierarchy, boundary_curves; skip_points = skip_points_set)
    if length(skip_points_set) > n 
        setdiff!(insertion_order, skip_points_set) 
    end
    tri = Triangulation(points; IntegerType, EdgeType, TriangleType, EdgesType, TrianglesType, weights, boundary_curves, boundary_enricher, build_cache = Val(true))
    return _triangulate!(
        tri, segments, boundary_nodes, predicates, randomise, try_last_inserted_point, skip_points_set, num_sample_rule, rng, insertion_order,
        recompute_representative_points, delete_holes, full_polygon_hierarchy, delete_ghosts, delete_empty_features,
    )
end

function check_config(points, weights, segments, boundary_nodes, kernel)
    (number_type(points) == Float32 && kernel != AdaptiveKernel()) && @warn "Using non-Float64 coordinates may cause issues. If you run into problems, consider using Float64 coordinates." maxlog = 1
    is_constrained = !(isnothing(segments) || isempty(segments)) || !(isnothing(boundary_nodes) || !has_boundary_nodes(boundary_nodes))
    # is_weighted(weights) && is_constrained && throw(ArgumentError("You cannot compute a constrained triangulation with weighted points."))
    # The above is still not implemented, but it errors when someone uses lock_convex_hull! (like in centroidal_smooth). So, just ignore it for now I guess.
    return nothing
end

function _triangulate!(
        tri::Triangulation, segments, boundary_nodes, predicates, randomise, try_last_inserted_point, skip_points, num_sample_rule, rng, insertion_order,
        recompute_representative_points, delete_holes, full_polygon_hierarchy, delete_ghosts, delete_empty_features,
    )
    unconstrained_triangulation!(tri; predicates, randomise, try_last_inserted_point, skip_points, num_sample_rule, rng, insertion_order)
    _tri = if !(isnothing(segments) || isempty(segments)) || !(isnothing(boundary_nodes) || !has_boundary_nodes(boundary_nodes))
        constrained_triangulation!(tri, segments, boundary_nodes, predicates, full_polygon_hierarchy; rng, delete_holes)
    else
        tri
    end
    postprocess_triangulate!(_tri; delete_ghosts, delete_empty_features, recompute_representative_points)
    return _tri
end

"""
    retriangulate(tri::Triangulation, points=get_points(tri); kwargs...)

Retriangulates the triangulation `tri` using the points `points`, returning a new [`Triangulation`](@ref). 

# Arguments
- `tri::Triangulation`: The triangulation to retriangulate.
- `points=get_points(tri)`: The points to use for retriangulating the triangulation. By default, this is simply `get_points(tri)`.

# Keyword Arguments
- `skip_points=Set(filter(i -> !has_vertex(tri, i), each_point_index(tri)))`: The points to skip when inserting points into the triangulation.
- `kwargs...`: Extra keyword arguments passed to `triangulate`. Other keyword arguments, like `segments` and `boundary_nodes`, 
   are automatically passed from the fields of `tri`, but may be overridden by passing the corresponding keyword arguments.
"""
function retriangulate(
        tri::T, points = get_points(tri);
        segments = get_interior_segments(tri),
        boundary_nodes = get_boundary_nodes(tri),
        skip_points = Set(filter(i -> !has_vertex(tri, i), each_point_index(tri))),
        kwargs...,
    ) where {T <: Triangulation}
    return triangulate(
        points;
        segments,
        boundary_nodes,
        weights = get_weights(tri),
        IntegerType = integer_type(tri),
        EdgeType = edge_type(tri),
        TriangleType = triangle_type(tri),
        EdgesType = edges_type(tri),
        TrianglesType = triangles_type(tri),
        check_arguments = false,
        skip_points,
        kwargs...,
    )::T
end
