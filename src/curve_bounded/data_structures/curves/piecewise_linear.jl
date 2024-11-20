struct PiecewiseLinear{P,V} <: AbstractParametricCurve
    points::P
    boundary_nodes::V
end
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