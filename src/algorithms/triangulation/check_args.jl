"""
    check_args(points, boundary_nodes, segments, hierarchy::PolygonHierarchy, boundary_curves = (); skip_points = Set{Int}()) -> Bool 

Check that the arguments `points` and `boundary_nodes` to [`triangulate`](@ref), and a constructed 
[`PolygonHierarchy`](@ref) given by `hierarchy`, are valid. In particular, the function checks:

- The dimension of the points. If the points are not 2D, a warning is thrown.
- The points are all unique. If they are not, a warning is thrown and the indices of the duplicates are merged into `skip_points`.
- There are at least three points. If there are not, an `InsufficientPointsError` is thrown.

If any duplicate points are found, the indices of the duplicates are merged into `skip_points` in-place.

If `boundary_nodes` are provided, meaning [`has_boundary_nodes`](@ref), then the function also checks:

- If the boundary curves all connect consistently. Here, this means that each section of a boundary curve ends at the start of the next boundary section;
   for contiguous boundary curves, this means that the start and end boundary nodes are the same.
- If the orientation of the boundary curves are all consistent. This means that the curves are all positively oriented relative to the domain,
   so that e.g. the exterior boundary curves are all counter-clockwise (relative to just themselves), the next exterior-most curves inside those 
   exteriors are all clockwise (again, relative to just themselves), and so on.

The arguments `boundary_nodes` and `segments` are also used when checking for duplicate points. Any duplicate points that are also referenced in `boundary_nodes`
and `segments` are updated so the vertex refers to the first instance of the duplicate point. This is done in-place, so that the original `boundary_nodes` and `segments` are modified.

!!! danger "Mutation"

    If indeed duplicate points are found, the function modifies the `boundary_nodes` and `segments` in-place. This means that the original
    `boundary_nodes` and `segments` are modified, and the original points are not modified. The indices of the duplicates are merged into `skip_points` in-place.

!!! danger "Intersecting boundaries"

    Another requirement for [`triangulate`](@ref) is that none of the boundaries intersect in their interior, which also prohibits 
    interior self-intersections. This is NOT checked. Similarly, segments should not intersect in their interior, which is not checked.
"""
function check_args(points, boundary_nodes, segments, hierarchy, boundary_curves=(); skip_points=Set{Int}())
    check_dimension(points)
    has_unique_points!(skip_points, points, boundary_nodes, segments)
    has_enough_points(points)
    if has_boundary_nodes(boundary_nodes)
        has_consistent_connections(boundary_nodes)
        has_consistent_orientations(hierarchy, boundary_nodes, is_curve_bounded(boundary_curves))
    end
    return true
end

struct InsufficientPointsError{P} <: Exception
    points::P
end
struct InconsistentConnectionError{I,J} <: Exception
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

function check_dimension(points)
    valid = is_planar(points)
    if !valid
        @warn "The provided points are not in the plane. All but the first two coordinates of each point will be ignored." maxlog = 1
    end
    return valid
end

function has_unique_points!(skip_points, points, boundary_nodes, segments)
    all_unique = points_are_unique(points)
    if !all_unique
        dup_seen = find_duplicate_points(points)
        if WARN_ON_DUPES[]
            io = IOBuffer()
            println(io, "There were duplicate points. Only one of each duplicate will be used (the first vertex encountered in the order that follows), and all other duplicates will be skipped. The indices of the duplicates are:")
        end
        for (p, ivec) in dup_seen
            # Skip all but the first duplicate point for each point.
            for j in (firstindex(ivec)+1):lastindex(ivec)
                push!(skip_points, ivec[j])
            end
            if WARN_ON_DUPES[]
                println(io, "  ", p, " at indices ", ivec)
            end
            # We need to be careful about the duplicate points in the boundary nodes and segments. 
            # https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/220
            ref = first(ivec)
            if !(isnothing(segments) || num_edges(segments) == 0)
                # We can't delete the segments while iterating, so we need to build a list 
                E = edge_type(segments)
                segment_map = Dict{E,E}()
                for e in each_edge(segments)
                    u, v = edge_vertices(e)
                    if u ∈ ivec
                        segment_map[e] = construct_edge(E, ref, v)
                    elseif v ∈ ivec
                        segment_map[e] = construct_edge(E, u, ref)
                    end
                end
                for (e, new_e) in segment_map
                    delete_edge!(segments, e)
                    add_edge!(segments, new_e)
                end
            end 
            if has_boundary_nodes(boundary_nodes)
                if has_multiple_curves(boundary_nodes)
                    for k in 1:num_curves(boundary_nodes)
                        curve = get_boundary_nodes(boundary_nodes, k)
                        for j in 1:num_sections(curve)
                            segment = get_boundary_nodes(curve, j)
                            for i in 1:(num_boundary_edges(segment)+1)
                                v = get_boundary_nodes(segment, i)
                                if v ∈ ivec
                                    set_boundary_node!(boundary_nodes, ((k, j), i), ref)
                                end
                            end
                        end
                    end
                elseif has_multiple_sections(boundary_nodes)
                    for j in 1:num_sections(boundary_nodes)
                        segment = get_boundary_nodes(boundary_nodes, j)
                        for i in 1:(num_boundary_edges(segment)+1)
                            v = get_boundary_nodes(segment, i)
                            if v ∈ ivec
                                set_boundary_node!(boundary_nodes, (j, i), ref)
                            end
                        end
                    end
                else
                    for i in 1:(num_boundary_edges(boundary_nodes)+1)
                        v = get_boundary_nodes(boundary_nodes, i)
                        if v ∈ ivec
                            set_boundary_node!(boundary_nodes, (boundary_nodes, i), ref)
                        end
                    end
                end
            end
        end
        if WARN_ON_DUPES[]
            println(io, "To suppress this warning, call `DelaunayTriangulation.toggle_warn_on_dupes!()`.")
        end
        if WARN_ON_DUPES[]
            @warn String(take!(io))
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
function has_consistent_connections_multiple_sections(boundary_nodes, curve_index=0)
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
