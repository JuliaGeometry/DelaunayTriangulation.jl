function triangle_area(p, q, r)
    A = AdaptivePredicates.orient2(getxy(p), getxy(q), getxy(r)) / 2
    return A
end

function triangle_area(ℓ₁²::Number, ℓ₂²::Number, ℓ₃²::Number)
    A² = (4ℓ₁² * ℓ₂² - (ℓ₁² + ℓ₂² - ℓ₃²)^2) / 16 # Herons' formula
    if A² < zero(A²)
        A² = let a = sqrt(ℓ₃²), b = sqrt(ℓ₂²), c = sqrt(ℓ₁²)
            (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)) / 16 # https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf
        end
    end
    A² < zero(A²) && return zero(A²)
    return sqrt(A²)
end

triangle_circumradius(A, ℓmin², ℓmed², ℓmax²) = sqrt(ℓmin² * ℓmed² * ℓmax²) / (4A)

triangle_perimeter(ℓmin::Number, ℓmed::Number, ℓmax::Number) = ℓmin + ℓmed + ℓmax

triangle_inradius(A, perimeter) = 2A / perimeter

triangle_aspect_ratio(inradius::Number, circumradius::Number) = inradius / circumradius

triangle_radius_edge_ratio(circumradius::Number, ℓmin::Number) = circumradius / ℓmin

function triangle_centroid(p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    qx′, qy′ = qx - px, qy - py
    rx′, ry′ = rx - px, ry - py
    cx′, cy′ = 2midpoint(qx′, rx′) / 3, 2midpoint(qy′, ry′) / 3
    cx, cy = cx′ + px, cy′ + py
    return cx, cy
end

function triangle_angles(p, q, r)
    A = triangle_area(p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    ax, by = px - qx, py - qy
    bx, ay = px - rx, py - ry
    dotab = ax * bx + ay * by
    θ₁ = if iszero(dotab)
        one(dotab) * π / 2
    else
        atan(2A / dotab)
    end
    if θ₁ < 0
        θ₁ += π
    end
    ax, by = qx - px, qy - py
    bx, ay = qx - rx, qy - ry
    dotab = ax * bx + ay * by
    θ₂ = if iszero(dotab)
        one(dotab) * π / 2
    else
        atan(2A / dotab)
    end
    if θ₂ < 0
        θ₂ += π
    end
    ax, by = rx - px, ry - py
    bx, ay = rx - qx, ry - qy
    dotab = ax * bx + ay * by
    θ₃ = if iszero(dotab)
        one(dotab) * π / 2
    else
        atan(2A / dotab)
    end
    if θ₃ < 0
        θ₃ += π
    end
    θ₁, θ₂, θ₃ = min_med_max(θ₁, θ₂, θ₃)
    return θ₁, θ₂, θ₃
end

function squared_triangle_lengths(p, q, r)
    ℓ₁², ℓ₂², ℓ₃², _ = squared_triangle_lengths_and_smallest_index(p, q, r)
    return ℓ₁², ℓ₂², ℓ₃²
end

function squared_triangle_lengths_and_smallest_index(p, q, r)
    p = getxy(p)
    q = getxy(q)
    r = getxy(r)
    ℓ₁² = dist_sqr(p, q)
    ℓ₂² = dist_sqr(q, r)
    ℓ₃² = dist_sqr(r, p)
    ℓmin², ℓmed², ℓmax² = min_med_max(ℓ₁², ℓ₂², ℓ₃²)
    ℓmin² == ℓ₁² && return ℓmin², ℓmed², ℓmax², 1
    ℓmin² == ℓ₂² && return ℓmin², ℓmed², ℓmax², 2
    return ℓmin², ℓmed², ℓmax², 3
end

function triangle_lengths(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    return sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
end

function triangle_circumcenter(p, q, r, A = triangle_area(p, q, r))
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    d11 = dist_sqr(p, r)
    d12 = py - ry
    d21 = dist_sqr(q, r)
    d22 = qy - ry
    ox = rx + (d11 * d22 - d12 * d21) / (4A)
    e11 = px - rx
    e12 = d11
    e21 = qx - rx
    e22 = d21
    oy = ry + (e11 * e22 - e12 * e21) / (4A)
    return (ox, oy)
end

function triangle_circumcenter(tri::Triangulation, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    return triangle_circumcenter(p, q, r)
end

function triangle_circumradius(p, q, r)
    ℓ₁², ℓ₂², ℓ₃² = squared_triangle_lengths(p, q, r)
    A = triangle_area(p, q, r)
    return triangle_circumradius(A, ℓ₁², ℓ₂², ℓ₃²)
end

function triangle_offcenter(p, q, r, c₁ = triangle_circumcenter(p, q, r), β = 1.0)
    ℓ₁², ℓ₂², _, idx = squared_triangle_lengths_and_smallest_index(p, q, r)
    ℓ₁ = sqrt(ℓ₁²)
    p, q, r = make_shortest_edge_first(p, q, r, idx)
    if ℓ₁² ≈ ℓ₂² # need to choose the edge out of the pair of shortest edges whose midpoint is furthest from c₁
        p, q, r = select_shortest_edge_for_offcenter(p, q, r, c₁, ℓ₁²)
    end
    h = distance_to_offcenter(β, ℓ₁)
    p = getxy(p)
    q = getxy(q)
    m = midpoint(p, q)
    c₁ = getxy(c₁)
    dist_to_c₁ = dist(m, c₁)
    c₁x, c₁y = c₁
    mx, my = m
    dirx, diry = (c₁x - mx) / dist_to_c₁, (c₁y - my) / dist_to_c₁
    ox = mx + h * dirx
    oy = my + h * diry
    return ox, oy
end

function distance_to_offcenter(β, ℓ)
    tβ = 2β
    Δ = sqrt(tβ^2 - 1)
    if β ≥ 1 / 2
        h = (tβ + Δ) / 2
    else
        h = 1 / (2Δ)
    end
    return h * ℓ
end

function make_shortest_edge_first(p, q, r, idx)
    if idx == 2
        return q, r, p
    elseif idx == 3
        return r, p, q
    else
        return p, q, r
    end
end

function select_shortest_edge_for_offcenter(p, q, r, c, ℓ²)
    p = getxy(p)
    q = getxy(q)
    r = getxy(r)
    pq = midpoint(p, q)
    qr = midpoint(q, r)
    rp = midpoint(r, p)
    ℓqr² = dist_sqr(r, q)
    ℓrp² = dist_sqr(r, p)
    ℓpqc² = dist_sqr(pq, c)
    if ℓqr² ≈ ℓ² && ℓrp² > ℓ² # shortest edges are pq and qr
        ℓqrc² = dist_sqr(qr, c)
        if ℓpqc² ≈ ℓqrc²
            if pq < qr
                return p, q, r
            else
                return q, r, p
            end
        elseif ℓpqc² > ℓqrc²
            return p, q, r
        else
            return q, r, p
        end
    elseif ℓrp² ≈ ℓ² && ℓqr² > ℓ² # shortest edges are pq and rp
        ℓrpc² = dist_sqr(rp, c)
        if ℓpqc² ≈ ℓrpc²
            if pq < rp
                return p, q, r
            else
                return r, p, q
            end
        elseif ℓpqc² > ℓrpc²
            return p, q, r
        else
            return r, p, q
        end
    end
    # The only other possibility is that the triangle is an equilateral triangle. This case doesn't matter since the circumcenter will be the 
    # same distance away from each midpoint anyway. To keep on with the idea of having the result be independent of the argument order, we
    # return the permutation that is lexicographically smallest.
    min_pt = min(p, q, r)
    min_idx = min_pt == p ? 1 : min_pt == q ? 2 : 3
    if min_idx == 1 # can't use min_med_max since we need to preserve orientation
        return p, q, r
    elseif min_idx == 2
        return q, r, p
    else
        return r, p, q
    end
end

function triangle_perimeter(p, q, r)
    ℓmin, ℓmed, ℓmax = triangle_lengths(p, q, r)
    return triangle_perimeter(ℓmin, ℓmed, ℓmax)
end

function triangle_inradius(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin, ℓmed, ℓmax = sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
    A = triangle_area(p, q, r)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    return triangle_inradius(A, perimeter)
end

function triangle_aspect_ratio(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin, ℓmed, ℓmax = sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
    A = triangle_area(p, q, r)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    inradius = triangle_inradius(A, perimeter)
    circumradius = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    return triangle_aspect_ratio(inradius, circumradius)
end

function triangle_radius_edge_ratio(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin = sqrt(ℓmin²)
    A = triangle_area(p, q, r)
    circumradius = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    return triangle_radius_edge_ratio(circumradius, ℓmin)
end

function triangle_edge_midpoints(p, q, r)
    return midpoint(p, q), midpoint(q, r), midpoint(r, p)
end

function min_med_max(a, b, c)
    b, c = minmax(b, c)
    a, c = minmax(a, c)
    a, b = minmax(a, b)
    return a, b, c
end

function triangle_sink(tri::Triangulation, T, prev_T = (integer_type(tri)(∅), integer_type(tri)(∅), integer_type(tri)(∅)); predicates::AbstractPredicateKernel = AdaptiveKernel())
    # TODO: This function would be faster if we just always search away from the largest angle.
    T = sort_triangle(T)
    c = triangle_circumcenter(tri, T)
    is_boundary_triangle(tri, T) && return c
    !has_ghost_triangles(tri) && any(e -> !edge_exists(tri, reverse_edge(e)), triangle_edges(T)) && return c
    flag = point_position_relative_to_triangle(predicates, tri, T, c)
    !is_outside(flag) && return c
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    #=
    TODO: This comment below might be outdated now that we use adaptive arithmetic for computing triangle areas.

    This function needs to be a bit slower than I'd like. Unfortunately, since we do not have a good 
    robust function for computing the circumcenter, it is possible for a circumcenter to be on an edge 
    but not recognisable as being on that triangle or on the adjoining triangle - no triangle contains 
    the circumcenter! As an example, consider the triangulation 
        triangulate_rectangle(0, 10, 0, 10, 14, 14)
    and the triangle (178, 165, 179). Its circumcenter is on the edge (178, 165), but 
    point_position_relative_to_circle returns `Outside`. If we consider `(165, 178, 164)` instead, 
    which adjoins the edge `(189, 165)`, then its circumcenter also returns `Outside`. To get around this, 
    we perform a secondary check to see if the triangle is obtuse. If it is, then we know that the
    circumcenter is outside of the triangle.
    =#
    _, _, θ₃ = triangle_angles(p, q, r)
    θ₃ ≤ π / 2 + ε(θ₃) && return c
    m = triangle_centroid(p, q, r)
    if !is_none(line_segment_intersection_type(predicates, p, q, m, c)) && !is_left(point_position_relative_to_line(predicates, p, q, c))
        next_T = (j, i, get_adjacent(tri, j, i))
    elseif !is_none(line_segment_intersection_type(predicates, q, r, m, c)) && !is_left(point_position_relative_to_line(predicates, q, r, c))
        next_T = (k, j, get_adjacent(tri, k, j))
    else # Must intersect the edge ki instead then 
        next_T = (i, k, get_adjacent(tri, i, k))
    end
    sort_triangle(next_T) == prev_T && return c
    return triangle_sink(tri, next_T, T; predicates)
end

function triangle_orthocenter(tri::Triangulation, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    a, b, c = get_weight(tri, i), get_weight(tri, j), get_weight(tri, k)
    return triangle_orthocenter(p, q, r, a, b, c)
end
function triangle_orthocenter(p, q, r, a, b, c)
    A = triangle_area(p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    d11 = dist_sqr(p, r) + c - a 
    d12 = py - ry 
    d21 = dist_sqr(q, r) + c - b 
    d22 = qy - ry 
    ox = rx + (d11 * d22 - d12 * d21) / (4A)
    e11 = px - rx
    e12 = d11
    e21 = qx - rx
    e22 = d21
    oy = ry + (e11 * e22 - e12 * e21) / (4A)
    return ox, oy
end

function triangle_orthoradius_squared(tri::Triangulation, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    a, b, c = get_weight(tri, i), get_weight(tri, j), get_weight(tri, k)
    return triangle_orthoradius_squared(p, q, r, a, b, c)
end
function triangle_orthoradius_squared(p, q, r, a, b, c)
    A = triangle_area(p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    d11 = dist_sqr(p, r) + c - a
    d21 = dist_sqr(q, r) + c - b
    d12 = py - ry
    d22 = qy - ry
    e11 = px - rx
    e21 = qx - rx
    e12 = d11
    e22 = d21
    t1 = d11 * d22 - d12 * d21
    t2 = e11 * e22 - e12 * e21
    return (t1^2 + t2^2)/(16A^2) - c
end

