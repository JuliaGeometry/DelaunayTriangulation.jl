@doc """
    is_curve_bounded(tri::Triangulation) -> Bool 
    is_curve_bounded(boundary_nodes) -> Bool

Returns `true` if `tri` is curve bounded, and `false` otherwise; similarly for the 
`boundary_nodes` method.
"""
is_curve_bounded
@inline is_curve_bounded(::Tuple{}) = false
@inline is_curve_bounded(::Any) = false
@inline is_curve_bounded(::Number) = false
@inline is_curve_bounded(boundary_curves::Tuple) = any(is_curve_bounded, boundary_curves)
@inline is_curve_bounded(tri::Triangulation) = is_curve_bounded(get_boundary_curves(tri))
@inline function is_curve_bounded(boundary_nodes::AbstractVector)::Bool
    if has_multiple_curves(boundary_nodes)
        nc = num_curves(boundary_nodes)
        for curve_index in 1:nc
            is_curve_bounded(get_boundary_nodes(boundary_nodes, curve_index))::Bool && return true
        end
        return false
    elseif has_multiple_sections(boundary_nodes)
        ns = num_sections(boundary_nodes)
        for section_index in 1:ns
            is_curve_bounded(get_boundary_nodes(boundary_nodes, section_index))::Bool && return true
        end
        return false
    else
        return has_boundary_nodes(boundary_nodes) && is_curve_bounded(get_boundary_nodes(boundary_nodes, 1))::Bool
    end
end

"""
    to_boundary_curves(points, boundary_nodes) -> NTuple{N, AbstractParametricCurve} where N 

Returns the set of boundary curves associated with `boundary_nodes` and `points`. 
"""
@inline function to_boundary_curves(points, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return _to_boundary_curves_multiple_curves(points, boundary_nodes)
    elseif has_multiple_sections(boundary_nodes)
        return _to_boundary_curves_multiple_sections(points, boundary_nodes)
    else
        return _to_boundary_curves_single_curve(points, boundary_nodes)
    end
end
@inline function _to_boundary_curves_single_curve(points, boundary_nodes)
    boundary_node = get_boundary_nodes(boundary_nodes, 1)
    if is_curve_bounded(boundary_node)
        return (boundary_node,)
    else
        return (PiecewiseLinear(points, boundary_nodes),)
    end
end
@inline function _to_boundary_curves_multiple_sections(points, boundary_nodes, section = 1, boundary_curves = ())
    if section > num_sections(boundary_nodes)
        return boundary_curves
    else
        section_nodes = get_boundary_nodes(boundary_nodes, section)
        new_boundary_curves = _to_boundary_curves_single_curve(points, section_nodes)
        return _to_boundary_curves_multiple_sections(points, boundary_nodes, section + 1, (boundary_curves..., new_boundary_curves...))
    end
end
@inline function _to_boundary_curves_multiple_curves(points, boundary_nodes, curve = 1, boundary_curves = ())
    if curve > num_curves(boundary_nodes)
        return boundary_curves
    else
        curve_nodes = get_boundary_nodes(boundary_nodes, curve)
        new_boundary_curves = _to_boundary_curves_multiple_sections(points, curve_nodes)
        return _to_boundary_curves_multiple_curves(points, boundary_nodes, curve + 1, (boundary_curves..., new_boundary_curves...))
    end
end

"""
    convert_boundary_curves!(points, boundary_nodes, IntegerType) -> (NTuple{N, AbstractParametricCurve} where N, Vector)

Converts the provided `points` and `boundary_nodes` into a set of boundary curves and modified boundary nodes suitable for 
triangulation. In particular:

1. The function gets `boundary_curves` from [`to_boundary_curves`](@ref).
2. `boundary_nodes` is replaced with a set of initial boundary nodes (from [`get_skeleton`](@ref)). These nodes come from evaluating each boundary curve at `t = 0` and `t = 1`. In the case of a piecewise linear boundary,
   the vertices are copied directly. Note that not all control points of a [`CatmullRomSpline`](@ref) (which [`is_interpolating`](@ref)) will be added - only those at `t = 0` and `t = 1`.
3. The `points` are modified to include the new boundary nodes. If a point is already in `points`, it is not added again.  
    
# Arguments 
- `points`: The point set. This is modified in place with the new boundary points.
- `boundary_nodes`: The boundary nodes to be converted. This is not modified in place.
- `IntegerType`: The type of integer to use for the boundary nodes.

# Output 
- `boundary_curves`: The boundary curves associated with `boundary_nodes`.
- `boundary_nodes`: The modified boundary nodes.
"""
@inline function convert_boundary_curves!(points, boundary_nodes, ::Type{I}) where {I <: Integer}
    boundary_curves = to_boundary_curves(points, boundary_nodes)
    !is_curve_bounded(boundary_curves) && return boundary_curves, boundary_nodes
    new_boundary_nodes = get_skeleton(boundary_nodes, I)
    if has_multiple_curves(boundary_nodes)
        _convert_boundary_curves_multiple_curves!(points, boundary_nodes, boundary_curves, new_boundary_nodes)
    elseif has_multiple_sections(boundary_nodes)
        _convert_boundary_curves_multiple_sections!(points, boundary_nodes, boundary_curves, new_boundary_nodes)
    else
        _convert_boundary_curves_contiguous!(points, boundary_nodes, boundary_curves, 1, new_boundary_nodes)
    end
    return boundary_curves, new_boundary_nodes
end
@inline function _convert_boundary_curves_multiple_curves!(points, boundary_nodes, boundary_curves, new_boundary_nodes)
    ctr = 1
    for curve_index in 1:num_curves(boundary_nodes)
        curve_nodes = get_boundary_nodes(boundary_nodes, curve_index)
        new_curve_nodes = get_boundary_nodes(new_boundary_nodes, curve_index)
        for section_index in 1:num_sections(curve_nodes)
            section_nodes = get_boundary_nodes(curve_nodes, section_index)
            new_section_nodes = get_boundary_nodes(new_curve_nodes, section_index)
            _convert_boundary_curves_contiguous!(points, section_nodes, boundary_curves, ctr, new_section_nodes)
            ctr += 1
        end
    end
    return nothing
end
@inline function _convert_boundary_curves_multiple_sections!(points, boundary_nodes, boundary_curves, new_boundary_nodes)
    for section_index in 1:num_sections(boundary_nodes)
        section_nodes = get_boundary_nodes(boundary_nodes, section_index)
        new_section_nodes = get_boundary_nodes(new_boundary_nodes, section_index)
        _convert_boundary_curves_contiguous!(points, section_nodes, boundary_curves, section_index, new_section_nodes)
    end
    return nothing
end
@inline function _convert_boundary_curves_contiguous!(points, boundary_nodes, boundary_curves, curve_index, new_boundary_nodes)
    if is_piecewise_linear(boundary_curves, curve_index)
        n = num_boundary_edges(boundary_nodes)
        for i in 1:(n + 1)
            v = get_boundary_nodes(boundary_nodes, i)
            insert_boundary_node!(new_boundary_nodes, (new_boundary_nodes, i), v)
        end
    else
        _points = (eval_boundary_curve(boundary_curves, curve_index, 0.0), eval_boundary_curve(boundary_curves, curve_index, 1.0))
        ctr = 1
        for p in _points
            i = find_point_index(points, p)
            if i == ∅
                push_point!(points, p)
                i = num_points(points)
            end
            insert_boundary_node!(new_boundary_nodes, (new_boundary_nodes, ctr), i)
            ctr += 1
        end
    end
    return nothing
end
