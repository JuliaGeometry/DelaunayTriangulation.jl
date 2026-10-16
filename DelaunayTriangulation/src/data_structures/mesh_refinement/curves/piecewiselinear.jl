"""
    PiecewiseLinear <: AbstractParametricCurve

Struct for representing a piecewise linear curve. This curve should not be 
interacted with or constructed directly. It only exists so that it can be 
an [`AbstractParametricCurve`](@ref). Instead, triangulations use this curve to 
know that its `boundary_nodes` field should be used instead.

!!! danger "Existing methods"

    This struct does have fields, namely `points` and `boundary_nodes` (and boundary_nodes should be a contiguous section). These are only used so that 
    we can use this struct in [`angle_between`](@ref) easily. In particular, we need to allow 
    for evaluating this curve at `t=0` and at `t=1`, and similarly for differentiating the curve at `t=0` 
    and at `t=1`. For this, we have defined, letting `L` be a `PiecewiseLinear` curve, `L(0)` to return the first point 
    on the curve, and the last point otherwise (meaning `L(h)` is constant for `h > 0`), and similarly for differentiation.
    Do NOT rely on the implementation of these methods.
"""
struct PiecewiseLinear{P,V} <: AbstractParametricCurve
    points::P
    boundary_nodes::V
end
Base.show(io::IO, ::PiecewiseLinear) = print(io, "PiecewiseLinear()")
Base.show(io::IO, ::MIME"text/plain", L::PiecewiseLinear) = Base.show(io, L)
Base.:(==)(L1::PiecewiseLinear, L2::PiecewiseLinear) = get_points(L1) == get_points(L2) && get_boundary_nodes(L1) == get_boundary_nodes(L2)

Base.copy(L::PiecewiseLinear) = _plcopy(L) # so we can check the aliasing when we copy from a triangulation
_plcopy(L::PiecewiseLinear; points=copy(get_points(L)), boundary_nodes=copy(get_boundary_nodes(L))) = PiecewiseLinear(points, boundary_nodes)


is_piecewise_linear(::PiecewiseLinear) = true
get_points(pl::PiecewiseLinear) = pl.points
get_boundary_nodes(pl::PiecewiseLinear) = pl.boundary_nodes
function (L::PiecewiseLinear)(t) # ONLY FOR EVALUATING AT THE ENDPOINTS.
    points = get_points(L)
    boundary_nodes = get_boundary_nodes(L)
    if iszero(t)
        u = get_boundary_nodes(boundary_nodes, 1)
    else
        n = num_boundary_edges(boundary_nodes)
        u = get_boundary_nodes(boundary_nodes, n + 1)
    end
    p = get_point(points, u)
    return p
end
function differentiate(L::PiecewiseLinear, t)
    points = get_points(L)
    boundary_nodes = get_boundary_nodes(L)
    if iszero(t)
        u = get_boundary_nodes(boundary_nodes, 1)
        v = get_boundary_nodes(boundary_nodes, 2)
    else
        n = num_boundary_edges(boundary_nodes)
        u = get_boundary_nodes(boundary_nodes, n)
        v = get_boundary_nodes(boundary_nodes, n + 1)
    end
    p, q = get_point(points, u), get_point(points, v)
    px, py = getxy(p)
    qx, qy = getxy(q)
    return (qx - px, qy - py)
end

function get_circle_intersection(L::PiecewiseLinear, t₁, t₂, r)
    points = get_points(L)
    boundary_nodes = get_boundary_nodes(L)
    if iszero(t₁)
        u, v = get_boundary_nodes(boundary_nodes, 1), get_boundary_nodes(boundary_nodes, 2)
    else
        n = num_boundary_edges(boundary_nodes)
        u, v = get_boundary_nodes(boundary_nodes, n + 1), get_boundary_nodes(boundary_nodes, n)
    end
    p, q = get_point(points, u, v)
    Ls = LineSegment(p, q)
    return get_circle_intersection(Ls, 0.0, 1.0, r)
end