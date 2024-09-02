"""
    check_args(points, boundary_nodes, hierarchy::PolygonHierarchy, boundary_curves = ()) -> Bool 

Check that the arguments `points` and `boundary_nodes` to [`triangulate`](@ref), and a constructed 
[`PolygonHierarchy`](@ref) given by `hierarchy`, are valid. In particular, the function checks:

- The dimension of the points. If the points are not 2D, a warning is thrown.
- The points are all unique. If they are not, a `DuplicatePointsError` is thrown.
- There are at least three points. If there are not, an `InsufficientPointsError` is thrown.

If `boundary_nodes` are provided, meaning [`has_boundary_nodes`](@ref), then the function also checks:

- If the boundary curves all connect consistently. Here, this means that each section of a boundary curve ends at the start of the next boundary section;
   for contiguous boundary curves, this means that the start and end boundary nodes are the same.
- If the orientation of the boundary curves are all consistent. This means that the curves are all positively oriented relative to the domain,
   so that e.g. the exterior boundary curves are all counter-clockwise (relative to just themselves), the next exterior-most curves inside those 
   exteriors are all clockwise (again, relative to just themselves), and so on.

!!! danger "Intersecting boundaries"

    Another requirement for [`triangulate`](@ref) is that none of the boundaries intersect in their interior, which also prohibits 
    interior self-intersections. This is NOT checked. Similarly, segments should not intersect in their interior, which is not checked.
"""
function check_args(points, boundary_nodes, hierarchy, boundary_curves = ())
    check_dimension(points)
    has_unique_points(points)
    has_enough_points(points)
    has_bnd = has_boundary_nodes(boundary_nodes)
    if has_bnd
        has_consistent_connections(boundary_nodes)
        has_consistent_orientations(hierarchy, boundary_nodes, is_curve_bounded(boundary_curves))
    end
    return true
end

struct DuplicatePointsError{P} <: Exception
    points::P
end
struct InsufficientPointsError{P} <: Exception
    points::P
end
struct InconsistentConnectionError{I, J} <: Exception
    curve_index::I
    segment_index₁::I
    segment_index₂::I
    vertex₁::J
    vertex₂::J
end
struct InconsistentOrientationError{I} <: Exception
    index::I
    should_be_positive::Bool
    is_sectioned::Bool
    is_curve_bounded::Bool
end
function Base.showerror(io::IO, err::DuplicatePointsError)
    points = err.points
    dup_seen = find_duplicate_points(points)
    println(io, "DuplicatePointsError: There were duplicate points. The following points are duplicated:")
    n = length(dup_seen)
    ctr = 1
    for (p, ivec) in dup_seen
        if ctr < n
            println(io, "    ", p, " at indices ", ivec)
        else
            print(io, "    ", p, " at indices ", ivec, ".")
        end
        ctr += 1
    end
    return io
end
function Base.showerror(io::IO, err::InsufficientPointsError)
    points = err.points
    print(io, "InsufficientPointsError: The provided point set has ", num_points(points), " points, but triangulations require at least three points.")
    return io
end
function Base.showerror(io::IO, err::InconsistentConnectionError)
    print(io, "InconsistentConnectionError: ")
    if !iszero(err.segment_index₁) && !(isone(err.segment_index₁) && isone(err.segment_index₂))
        print(io, "Segment ", err.segment_index₁)
        !iszero(err.curve_index) && print(io, " of curve ", err.curve_index)
        print(io, " ends at vertex ", err.vertex₁, " but the next segment, segment ", err.segment_index₂, ", starts at vertex ", err.vertex₂, ".")
    else
        is_contiguous = iszero(err.curve_index)
        print(io, "The ")
        !is_contiguous ? print(io, "boundary curve with index ", err.curve_index) : print(io, "boundary")
        v₁, v₂ = isone(err.segment_index₁) ? (err.vertex₁, err.vertex₂) : (err.vertex₂, err.vertex₁)
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
        # might nto be wrong.
        str2 = "\nIf this curve is defined by an AbstractParametricCurve, you may instead need to reverse the order of the control points defining" *
            " the sections of the curve; the `positive` keyword may also be of interest for CircularArcs and EllipticalArcs."
        str *= str2
    end
    sign = err.should_be_positive ? "positive" : "negative"
    sign2 = err.should_be_positive ? "negative" : "positive"
    print(io, "The orientation of the boundary curve with index ", err.index, " should be $sign, but it is $sign2.", str)
    return io
end

function check_dimension(points)
    valid = is_planar(points) 
    if !valid 
        @warn "The provided points are not in the plane. All but the first two coordinates will be ignored." maxlog=1
    end
    return valid
end

function has_unique_points(points)
    all_unique = points_are_unique(points)
    !all_unique && throw(DuplicatePointsError(points))
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
    segmentⱼ₋₁ = get_boundary_nodes(boundary_nodes, 1)
    nn = num_boundary_edges(segmentⱼ₋₁) + 1
    for j in 2:ns
        segmentⱼ = get_boundary_nodes(boundary_nodes, j)
        vⱼ₋₁ⁿ = get_boundary_nodes(segmentⱼ₋₁, nn)
        vⱼ¹ = get_boundary_nodes(segmentⱼ, 1)
        vⱼ₋₁ⁿ ≠ vⱼ¹ && throw(InconsistentConnectionError(curve_index, j - 1, j, vⱼ₋₁ⁿ, vⱼ¹))
        segmentⱼ₋₁ = segmentⱼ
        nn = num_boundary_edges(segmentⱼ) + 1
    end
    segmentⱼ = get_boundary_nodes(boundary_nodes, 1)
    vⱼ₋₁ⁿ = get_boundary_nodes(segmentⱼ₋₁, nn)
    vⱼ¹ = get_boundary_nodes(segmentⱼ, 1)
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
