function polygon_features(points, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return polygon_features_multiple_curves(points, boundary_nodes)
    elseif has_multiple_sections(boundary_nodes)
        return polygon_features_multiple_segments(points, boundary_nodes)
    else
        return polygon_features_single_segment(points, boundary_nodes)
    end
end
function polygon_features_single_segment(points, boundary_nodes; scale = Val(true))
    F = number_type(points)
    cx = zero(F)
    cy = zero(F)
    a = zero(F)
    n_edge = num_boundary_edges(boundary_nodes)
    vᵢ = get_boundary_nodes(boundary_nodes, 1)
    pᵢ = get_point(points, vᵢ)
    xᵢ, yᵢ = getxy(pᵢ)
    for j in 2:(n_edge + 1)
        vᵢ₊₁ = get_boundary_nodes(boundary_nodes, j)
        pᵢ₊₁ = get_point(points, vᵢ₊₁)
        xᵢ₊₁, yᵢ₊₁ = getxy(pᵢ₊₁)
        area_contribution = xᵢ * yᵢ₊₁ - xᵢ₊₁ * yᵢ
        cx += (xᵢ + xᵢ₊₁) * area_contribution
        cy += (yᵢ + yᵢ₊₁) * area_contribution
        a += area_contribution
        vᵢ, pᵢ, xᵢ, yᵢ = vᵢ₊₁, pᵢ₊₁, xᵢ₊₁, yᵢ₊₁
    end
    if is_true(scale)
        return a / 2, (cx / (3a), cy / (3a))
    else
        return a, (cx, cy)
    end
end
function polygon_features_multiple_segments(points, boundary_nodes)
    F = number_type(points)
    cx = zero(F)
    cy = zero(F)
    a = zero(F)
    ns = num_sections(boundary_nodes)
    for i in 1:ns
        bn = get_boundary_nodes(boundary_nodes, i)
        sub_a, (sub_cx, sub_cy) = polygon_features_single_segment(
            points, bn;
            scale = Val(false),
        )
        a += sub_a
        cx += sub_cx
        cy += sub_cy
    end
    return a / 2, (cx / (3a), cy / (3a))
end
function polygon_features_multiple_curves(points, boundary_nodes)
    F = number_type(points)
    cx = zero(F)
    cy = zero(F)
    a = zero(F)
    nc = num_curves(boundary_nodes)
    for i in 1:nc
        bn = get_boundary_nodes(boundary_nodes, i)
        sub_a, (sub_cx, sub_cy) = polygon_features_multiple_segments(points, bn)
        cx += sub_a * sub_cx
        cy += sub_a * sub_cy
        a += sub_a
    end
    return a, (cx / a, cy / a)
end

function squared_distance_to_segment(x₁, y₁, x₂, y₂, x, y)
    qp₁_x = x - x₁
    qp₁_y = y - y₁
    p₁p₂_x = x₂ - x₁
    p₁p₂_y = y₂ - y₁
    denom = dist_sqr((x₁, y₁), (x₂, y₂))
    t = (qp₁_x * p₁p₂_x + qp₁_y * p₁p₂_y) / denom
    ts = min(max(t, zero(t)), one(t)) # https://math.stackexchange.com/a/330329/861404
    intersect_x = x₁ + ts * p₁p₂_x
    intersect_y = y₁ + ts * p₁p₂_y
    δ² = dist_sqr((x, y), (intersect_x, intersect_y))
    return δ²
end

function distance_to_polygon(q, points, boundary_nodes)
    if has_multiple_curves(boundary_nodes)
        return distance_to_polygon_multiple_curves(q, points, boundary_nodes)
    elseif has_multiple_sections(boundary_nodes)
        return distance_to_polygon_multiple_segments(q, points, boundary_nodes)
    else
        return distance_to_polygon_single_segment(q, points, boundary_nodes)
    end
end
function distance_to_polygon_single_segment(q, points, boundary_nodes, is_in_outer = false; return_sqrt = Val(true))
    x, y = getxy(q)
    F = number_type(points)
    dist = typemax(F)
    n_edge = num_boundary_edges(boundary_nodes)
    vᵢ = get_boundary_nodes(boundary_nodes, 1)
    pᵢ = get_point(points, vᵢ)
    xᵢ, yᵢ = getxy(pᵢ)
    for j in 2:(n_edge + 1)
        vᵢ₊₁ = get_boundary_nodes(boundary_nodes, j)
        pᵢ₊₁ = get_point(points, vᵢ₊₁)
        xᵢ₊₁, yᵢ₊₁ = getxy(pᵢ₊₁)
        if (yᵢ₊₁ > y) ≠ (yᵢ > y)
            intersect_x = (xᵢ - xᵢ₊₁) * (y - yᵢ₊₁) / (yᵢ - yᵢ₊₁) + xᵢ₊₁
            if x < intersect_x
                is_in_outer = !is_in_outer
            end
        end
        new_dist = squared_distance_to_segment(xᵢ, yᵢ, xᵢ₊₁, yᵢ₊₁, x, y)
        dist = new_dist < dist ? new_dist : dist
        vᵢ, pᵢ, xᵢ, yᵢ = vᵢ₊₁, pᵢ₊₁, xᵢ₊₁, yᵢ₊₁
    end
    dist = is_true(return_sqrt) ? sqrt(dist) : dist
    return is_in_outer ? dist : -dist
end
function distance_to_polygon_multiple_segments(q, points, boundary_nodes, is_in_outer = -one(number_type(points)); return_sqrt = Val(true))
    F = number_type(points)
    dist = typemax(F)
    ns = num_sections(boundary_nodes)
    for i in 1:ns
        bn = get_boundary_nodes(boundary_nodes, i)
        new_dist = distance_to_polygon_single_segment(q, points, bn, is_in_outer == one(F); return_sqrt = Val(false))
        is_in_outer = sign(new_dist)
        new_dist = abs(new_dist)
        dist = new_dist < dist ? new_dist : dist
    end
    dist = is_true(return_sqrt) ? sqrt(dist) : dist
    return is_in_outer * dist
end
function distance_to_polygon_multiple_curves(q, points, boundary_nodes)
    F = number_type(points)
    is_in_outer = -one(F)
    dist = typemax(F)
    nc = num_curves(boundary_nodes)
    for i in 1:nc
        bn = get_boundary_nodes(boundary_nodes, i)
        new_dist = distance_to_polygon_multiple_segments(
            q, points, bn, is_in_outer == one(F);
            return_sqrt = Val(false),
        )
        is_in_outer = sign(new_dist)
        new_dist = abs(new_dist)
        dist = new_dist < dist ? new_dist : dist
    end
    return is_in_outer * sqrt(dist)
end

function polygon_bounds(points, boundary_nodes, check_all_curves = Val(false))
    if has_multiple_curves(boundary_nodes)
        if !is_true(check_all_curves)
            return polygon_bounds_multiple_segments(points, get_boundary_nodes(boundary_nodes, 1)) # 1 is the outermost boundary, unless you have a multiple polygon. PolygonHierarchy would be better for this but too bad 
        else
            F = number_type(points)
            xmin, xmax, ymin, ymax = typemax(F), typemin(F), typemax(F), typemin(F)
            nc = num_curves(boundary_nodes)
            for i in 1:nc
                bn = get_boundary_nodes(boundary_nodes, i)
                xminᵢ, xmaxᵢ, yminᵢ, ymaxᵢ = polygon_bounds_multiple_segments(points, bn)
                xmin = min(xminᵢ, xmin)
                xmax = max(xmaxᵢ, xmax)
                ymin = min(yminᵢ, ymin)
                ymax = max(ymaxᵢ, ymax)
            end
            return xmin, xmax, ymin, ymax
        end
    elseif has_multiple_sections(boundary_nodes)
        return polygon_bounds_multiple_segments(points, boundary_nodes)
    else
        return polygon_bounds_single_segment(points, boundary_nodes)
    end
end
function polygon_bounds_single_segment(points, boundary_nodes)
    F = number_type(points)
    xmin, xmax, ymin, ymax = typemax(F), typemin(F), typemax(F), typemin(F)
    n_edge = num_boundary_edges(boundary_nodes)
    for i in 1:n_edge
        vᵢ = get_boundary_nodes(boundary_nodes, i)
        pᵢ = get_point(points, vᵢ)
        xᵢ, yᵢ = getxy(pᵢ)
        xmin = min(xᵢ, xmin)
        xmax = max(xᵢ, xmax)
        ymin = min(yᵢ, ymin)
        ymax = max(yᵢ, ymax)
    end
    return xmin, xmax, ymin, ymax
end
function polygon_bounds_multiple_segments(points, boundary_nodes)
    F = number_type(points)
    xmin, xmax, ymin, ymax = typemax(F), typemin(F), typemax(F), typemin(F)
    ns = num_sections(boundary_nodes)
    for i in 1:ns
        bn = get_boundary_nodes(boundary_nodes, i)
        _xmin, _xmax, _ymin, _ymax = polygon_bounds_single_segment(points, bn)
        xmin = xmin > _xmin ? _xmin : xmin
        xmax = xmax < _xmax ? _xmax : xmax
        ymin = ymin > _ymin ? _ymin : ymin
        ymax = ymax < _ymax ? _ymax : ymax
    end
    return xmin, xmax, ymin, ymax
end

function sort_convex_polygon!(vertices, points)
    cx, cy = mean_points(points, vertices)
    to_angle = p -> atan(gety(p) - cy, getx(p) - cx)
    vert_to_angle = v -> to_angle(get_point(points, v))
    sort!(vertices, by = vert_to_angle)
    return vertices
end