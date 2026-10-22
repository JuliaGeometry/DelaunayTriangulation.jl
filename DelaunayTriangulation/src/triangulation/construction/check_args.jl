"""
    check_args(points, boundary_nodes, hierarchy::PolygonHierarchy, boundary_curves = (); skip_points = Set{Int}(), segments = nothing, weights = nothing) -> Bool

Check that the arguments `points`, `boundary_nodes`, `segments` and `weights` to [`triangulate`](@ref), and a constructed
[`PolygonHierarchy`](@ref) given by `hierarchy`, are valid. In particular, the function checks:

- The dimension of the points. If the points are not 2D, a warning is thrown.
- The coordinates of the points. If they are integers, an `IntegerCoordinatesError` is thrown, and if any coordinate of a point that is not in `skip_points` is not finite (e.g. `NaN` or `Inf`), a `NonFiniteCoordinatesError` is thrown.
- The points are all unique. If they are not, a warning is thrown and the indices of the duplicates are merged into `skip_points`. `-0.0` and `0.0` are treated as the same coordinate.
- There are at least three points. If there are not, an `InsufficientPointsError` is thrown.
- If `weights` is a vector, there is a weight for each point. If there are fewer weights than points, an `InconsistentWeightsError` is thrown. (Extra weights are allowed, e.g. for points that will be added later.)
- If `segments` are provided, each segment joins two different vertices that are indices of `points`. If not, an `InvalidSegmentError` is thrown.
- If the points are weighted, there are no `segments` or `boundary_nodes`, since only unconstrained weighted triangulations are supported. If there are, a `WeightedConstrainedTriangulationError` is thrown.

If any duplicate points are found, the indices of the duplicates are merged into `skip_points` in-place. If any of these
duplicates are also in `boundary_nodes`, a `DuplicateBoundaryNodeError` is thrown since they would not be in the triangulation.

If `boundary_nodes` are provided, meaning [`has_boundary_nodes`](@ref), then the function also checks:

- If the boundary nodes are all indices of `points`. If not, an `InvalidBoundaryNodeError` is thrown.
- If the boundary curves all connect consistently. Here, this means that each section of a boundary curve ends at the start of the next boundary section;
   for contiguous boundary curves, this means that the start and end boundary nodes are the same.
- If the orientation of the boundary curves are all consistent. This means that the curves are all positively oriented relative to the domain,
   so that e.g. the exterior boundary curves are all counter-clockwise (relative to just themselves), the next exterior-most curves inside those 
   exteriors are all clockwise (again, relative to just themselves), and so on.

!!! danger "Intersecting boundaries"

    Another requirement for [`triangulate`](@ref) is that none of the boundaries intersect in their interior, which also prohibits 
    interior self-intersections. This is NOT checked. Similarly, segments should not intersect in their interior, which is not checked.
"""
function check_args(points, boundary_nodes, hierarchy, boundary_curves = (); skip_points = Set{Int}(), segments = nothing, weights = nothing)
    check_dimension(points)
    has_valid_coordinates(points, skip_points)
    has_unique_points!(skip_points, points, boundary_nodes)
    has_enough_points(points)
    has_consistent_weights(points, weights)
    has_valid_segments(points, segments)
    is_unconstrained_if_weighted(weights, segments, boundary_nodes)
    has_bnd = has_boundary_nodes(boundary_nodes)
    if has_bnd
        has_valid_boundary_nodes(points, boundary_nodes)
        has_consistent_connections(boundary_nodes)
        has_consistent_orientations(hierarchy, boundary_nodes, is_curve_bounded(boundary_curves))
    end
    return true
end

struct InsufficientPointsError{P} <: Exception
    points::P
end
struct DuplicateBoundaryNodeError{I} <: Exception
    vertex::I
    duplicate_of::I
end
struct InconsistentConnectionError{I, J} <: Exception
    curve_index::I
    section_index₁::I
    section_index₂::I
    vertex₁::J
    vertex₂::J
end
struct InconsistentOrientationError{I} <: Exception
    index::I
    should_be_positive::Bool
    is_sectioned::Bool
    is_curve_bounded::Bool
end
struct IntegerCoordinatesError{T} <: Exception end
struct NonFiniteCoordinatesError{I, P} <: Exception
    index::I
    point::P
end
struct InconsistentWeightsError <: Exception
    num_weights::Int
    num_points::Int
end
struct InvalidSegmentError{E} <: Exception
    segment::E
    num_points::Int
end
struct InvalidBoundaryNodeError{I} <: Exception
    vertex::I
    num_points::Int
end
struct WeightedConstrainedTriangulationError <: Exception
    has_segments::Bool
    has_boundary_nodes::Bool
end
function Base.showerror(io::IO, err::InsufficientPointsError)
    points = err.points
    print(io, "InsufficientPointsError: The provided point set has ", num_points(points), " points, but triangulations require at least three points.")
    return io
end
function Base.showerror(io::IO, err::DuplicateBoundaryNodeError)
    print(io, "DuplicateBoundaryNodeError: The boundary node ", err.vertex, " is a duplicate of vertex ", err.duplicate_of, ", so it would be skipped. Boundary nodes must refer to unique points, so either remove the repeated point from the boundary or use vertex ", err.duplicate_of, " instead.")
    return io
end
function Base.showerror(io::IO, err::InconsistentConnectionError)
    print(io, "InconsistentConnectionError: ")
    if !iszero(err.section_index₁) && !(isone(err.section_index₁) && isone(err.section_index₂))
        print(io, "Section ", err.section_index₁)
        !iszero(err.curve_index) && print(io, " of curve ", err.curve_index)
        print(io, " ends at vertex ", err.vertex₁, " but the next section, section ", err.section_index₂, ", starts at vertex ", err.vertex₂, ".")
    else
        is_contiguous = iszero(err.curve_index)
        print(io, "The ")
        !is_contiguous ? print(io, "boundary curve with index ", err.curve_index) : print(io, "boundary")
        v₁, v₂ = isone(err.section_index₁) ? (err.vertex₁, err.vertex₂) : (err.vertex₂, err.vertex₁)
        print(io, " ends in vertex ", v₁, " but starts at vertex ", v₂, ".")
    end
    return io
end
function Base.showerror(io::IO, err::InconsistentOrientationError)
    print(io, "InconsistentOrientationError: ")
    suggestion = err.is_sectioned ? "reverse(reverse.(curve))" : "reverse(curve)"
    str = " You may be able to fix this by passing the curve as $suggestion."
    if err.is_curve_bounded
        # Only show this longer message if this part of the boundary could be defined by an AbstractParametricCurve. 
        # It's hard to detect if the curve is indeed defined by an AbstractParametricCurve since the curve could be defined 
        # by a combination of multiple AbstractParametricCurves and possibly a PiecewiseLinear part. Thus, the above advice
        # might not be wrong.
        str2 = "\nIf this curve is defined by an AbstractParametricCurve, you may instead need to reverse the order of the control points defining" *
            " the sections of the curve; the `positive` keyword may also be of interest for CircularArcs and EllipticalArcs. Alternatively, for individual" *
            " AbstractParametricCurves, note that `reverse` can be used to reverse the orientation of the curve directly instead of the control points."
        str *= str2
    end
    sign = err.should_be_positive ? "positive" : "negative"
    sign2 = err.should_be_positive ? "negative" : "positive"
    print(io, "The orientation of the boundary curve with index ", err.index, " should be $sign, but it is $sign2.", str)
    return io
end
function Base.showerror(io::IO, ::IntegerCoordinatesError{T}) where {T}
    print(io, "IntegerCoordinatesError: The provided points have coordinates of type ", T, ", but integer coordinates are not supported since, for example, centroids and circumcenters are computed using the same number type as the points. Convert the coordinates to floating point numbers first, e.g. to Float64.")
    return io
end
function Base.showerror(io::IO, err::NonFiniteCoordinatesError)
    print(io, "NonFiniteCoordinatesError: The point with index ", err.index, " has coordinates ", err.point, ", but all coordinates must be finite. If this point should not be triangulated, add its index to `skip_points`.")
    return io
end
function Base.showerror(io::IO, err::InconsistentWeightsError)
    print(io, "InconsistentWeightsError: There are ", err.num_weights, " weights but ", err.num_points, " points. There must be a weight for each point, with the `i`th weight being the weight of the `i`th point.")
    return io
end
function Base.showerror(io::IO, err::InvalidSegmentError)
    u, v = edge_vertices(err.segment)
    print(io, "InvalidSegmentError: ")
    if u == v
        print(io, "The segment ", (u, v), " joins the vertex ", u, " to itself.")
    else
        print(io, "The segment ", (u, v), " has a vertex that is not the index of a point. There are ", err.num_points, " points, so the vertices of a segment must be between 1 and ", err.num_points, ".")
    end
    return io
end
function Base.showerror(io::IO, err::WeightedConstrainedTriangulationError)
    constraints = err.has_segments && err.has_boundary_nodes ? "segments and boundary nodes" : err.has_segments ? "segments" : "boundary nodes"
    print(io, "WeightedConstrainedTriangulationError: The points are weighted and ", constraints, " were provided, but only unconstrained weighted triangulations are supported.")
    return io
end
function Base.showerror(io::IO, err::InvalidBoundaryNodeError)
    print(io, "InvalidBoundaryNodeError: The boundary node ", err.vertex, " is not the index of a point. There are ", err.num_points, " points, so the boundary nodes must be between 1 and ", err.num_points, ".")
    return io
end

function check_dimension(points)
    valid = is_planar(points) 
    if !valid 
        @warn "The provided points are not in the plane. All but the first two coordinates of each point will be ignored." maxlog=1
    end
    return valid
end

function has_valid_coordinates(points, skip_points = ())
    T = number_type(points)
    T <: Integer && throw(IntegerCoordinatesError{T}())
    for i in each_point_index(points)
        i ∈ skip_points && continue
        p = get_point(points, i)
        x, y = getxy(p)
        !(isfinite(x) && isfinite(y)) && throw(NonFiniteCoordinatesError(i, p))
    end
    return true
end

function has_consistent_weights(points, weights)
    # Weights can be any object that defines get_weight, and only vectors have a length to check. Extra weights 
    # are fine, e.g. for points that are pushed into points later.
    !(weights isa AbstractVector) && return true
    nw, np = length(weights), num_points(points)
    nw < np && throw(InconsistentWeightsError(nw, np))
    return true
end

function is_unconstrained_if_weighted(weights, segments, boundary_nodes)
    (isnothing(weights) || !is_weighted(weights)) && return true
    has_segments = !(isnothing(segments) || isempty(segments))
    has_bnd = has_boundary_nodes(boundary_nodes)
    (has_segments || has_bnd) && throw(WeightedConstrainedTriangulationError(has_segments, has_bnd))
    return true
end

function has_valid_segments(points, segments)
    isnothing(segments) && return true
    n = num_points(points)
    for e in each_edge(segments)
        u, v = edge_vertices(e)
        (u == v || !(1 ≤ u ≤ n) || !(1 ≤ v ≤ n)) && throw(InvalidSegmentError(e, n))
    end
    return true
end

function has_valid_boundary_nodes(points, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        for k in 1:num_curves(boundary_nodes)
            has_valid_boundary_nodes(points, get_boundary_nodes(boundary_nodes, k))
        end
    elseif has_multiple_sections(boundary_nodes)
        for k in 1:num_sections(boundary_nodes)
            has_valid_boundary_nodes(points, get_boundary_nodes(boundary_nodes, k))
        end
    else
        n = num_points(points)
        for i in 1:(num_boundary_edges(boundary_nodes) + 1)
            v = get_boundary_nodes(boundary_nodes, i)
            !(1 ≤ v ≤ n) && throw(InvalidBoundaryNodeError(v, n))
        end
    end
    return true
end

function has_unique_points!(skip_points, points, boundary_nodes = nothing)
    all_unique = points_are_unique(points)
    if !all_unique 
        dup_seen = find_duplicate_points(points)
        if WARN_ON_DUPES[]
            io = IOBuffer()
            println(io, "There were duplicate points. Only one of each duplicate will be used, and all other duplicates will be skipped. The indices of the duplicates are:")
        end
        for (p, ivec) in dup_seen 
            for j in 2:lastindex(ivec)
                push!(skip_points, ivec[j])
            end
            if WARN_ON_DUPES[]
                println(io, "  ", p, " at indices ", ivec)
            end
        end
        if WARN_ON_DUPES[]
            println(io, "To suppress this warning, call `DelaunayTriangulation.toggle_warn_on_dupes!()`.")
        end
        if WARN_ON_DUPES[]
            @warn String(take!(io))
        end
        has_boundary_nodes(boundary_nodes) && has_no_duplicate_boundary_nodes(boundary_nodes, dup_seen)
    end
    return true
end

function has_no_duplicate_boundary_nodes(boundary_nodes, dup_seen)
    # Only the first of each duplicate point is kept, so the others can't be boundary nodes. See Issue #220.
    duplicate_of = Dict(ivec[j] => ivec[1] for (_, ivec) in dup_seen for j in 2:lastindex(ivec))
    for section_index in values(construct_ghost_vertex_map(boundary_nodes))
        for v in each_boundary_node(get_boundary_nodes(boundary_nodes, section_index))
            u = get(duplicate_of, v, nothing)
            !isnothing(u) && throw(DuplicateBoundaryNodeError(v, oftype(v, u)))
        end
    end
    return true
end

function has_enough_points(points)
    has_enough = num_points(points) ≥ 3
    !has_enough && throw(InsufficientPointsError(points))
    return true
end

function has_consistent_orientations(hierarchy::PolygonHierarchy, boundary_nodes, is_curve_bounded)
    # Since trees start at height zero, the heights 0, 2, 4, ... must be positive, and 1, 3, 5, ... must be negative.
    for (_, tree) in get_trees(hierarchy)
        has_consistent_orientations(tree, hierarchy, boundary_nodes, is_curve_bounded)
    end
    return true
end
function has_consistent_orientations(tree::PolygonTree, hierarchy::PolygonHierarchy, boundary_nodes, is_curve_bounded)
    height = get_height(tree)
    index = get_index(tree)
    pos_orientation = get_polygon_orientation(hierarchy, index)
    is_sectioned = has_multiple_curves(boundary_nodes) || has_multiple_sections(boundary_nodes)
    if iseven(height)
        !pos_orientation && throw(InconsistentOrientationError(index, true, is_sectioned, is_curve_bounded))
    else
        pos_orientation && throw(InconsistentOrientationError(index, false, is_sectioned, is_curve_bounded))
    end
    for child in get_children(tree)
        has_consistent_orientations(child, hierarchy, boundary_nodes, is_curve_bounded)
    end
    return true
end

function has_consistent_connections(boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return has_consistent_connections_multiple_curves(boundary_nodes)
    elseif has_multiple_sections(boundary_nodes)
        return has_consistent_connections_multiple_sections(boundary_nodes)
    else
        return has_consistent_connections_contiguous(boundary_nodes)
    end
end
function has_consistent_connections_multiple_curves(boundary_nodes)
    nc = num_curves(boundary_nodes)
    for k in 1:nc
        curve_nodes = get_boundary_nodes(boundary_nodes, k)
        has_consistent_connections_multiple_sections(curve_nodes, nc > 1 ? k : 0)
    end
    return true
end
function has_consistent_connections_multiple_sections(boundary_nodes, curve_index = 0)
    ns = num_sections(boundary_nodes)
    sectionⱼ₋₁ = get_boundary_nodes(boundary_nodes, 1)
    nn = num_boundary_edges(sectionⱼ₋₁) + 1
    for j in 2:ns
        sectionⱼ = get_boundary_nodes(boundary_nodes, j)
        vⱼ₋₁ⁿ = get_boundary_nodes(sectionⱼ₋₁, nn)
        vⱼ¹ = get_boundary_nodes(sectionⱼ, 1)
        vⱼ₋₁ⁿ ≠ vⱼ¹ && throw(InconsistentConnectionError(curve_index, j - 1, j, vⱼ₋₁ⁿ, vⱼ¹))
        sectionⱼ₋₁ = sectionⱼ
        nn = num_boundary_edges(sectionⱼ) + 1
    end
    sectionⱼ = get_boundary_nodes(boundary_nodes, 1)
    vⱼ₋₁ⁿ = get_boundary_nodes(sectionⱼ₋₁, nn)
    vⱼ¹ = get_boundary_nodes(sectionⱼ, 1)
    vⱼ₋₁ⁿ ≠ vⱼ¹ && throw(InconsistentConnectionError(curve_index, ns, 1, vⱼ₋₁ⁿ, vⱼ¹))
    return true
end
function has_consistent_connections_contiguous(boundary_nodes)
    nn = num_boundary_edges(boundary_nodes) + 1
    v₁ = get_boundary_nodes(boundary_nodes, 1)
    vₙ = get_boundary_nodes(boundary_nodes, nn)
    v₁ ≠ vₙ && throw(InconsistentConnectionError(0, 0, 0, v₁, vₙ))
    return true
end
