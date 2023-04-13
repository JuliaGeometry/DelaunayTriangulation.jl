function compare_area(targets::RefinementTargets, T, A, p, q, r)
    flag = if targets.max_area isa Function
        targets.max_area(T, p, q, r, A)
    else
        A > targets.max_area
    end
    return flag 
end
function compare_ratio(targets::RefinementTargets, T, ρ, p, q, r)
    flag = if targets.max_radius_edge_ratio isa Function
        targets.max_radius_edge_ratio(T, p, q, r, ρ)
    else
        ρ > targets.max_radius_edge_ratio
    end
    return flag
end
function compare_points(targets::RefinementTargets, n)
    return n > targets.max_points
end
function compare_min_area(targets::RefinementTargets, A²)
    A² < zero(A²) && return true
    return sqrt(A²) ≤ targets.min_area
end

"""
    edge_is_seditious(tri::Triangulation, u, v, w, smallest_idx, ratio_flag, ρ)

Tests if the shortest edge of the triangle `(u, v, w)` (with radius-edge ratio `ρ`) is seditious, i.e. if it lies on two distinct segments, the two points lie on the same 
concentric circular shell of the common vertex of the two distinct segments, and the angle between the two distinct segments 
is small. The shortest edge is given by `smallest_idx`, with `smallest_idx == 1` mapping to `(u, v)`, `smallest_idx == 2`
mapping to `(v, w)`, and `smallest_idx == 3` mapping to `(w, u)`. 

!!! notes 

    We are not implementing the *precise* definition of seditious. To check that the points lie on the same 
    concentric circular shell, we just check if the points lie the same distance away from the common 
    vertex. This is sufficient as most points already lie on a concentric circular shell.

    Additionally, we do not actually compute the angle between the two segments. Instead, we just assume the 
    angle is small if `ratio_flag` is true. This is because computing the angle between two segments is
    just a bit expensive.
"""
function edge_is_seditious(tri::Triangulation, u, v, w, smallest_idx, ratio_flag, ρ)
    !ratio_flag && return false
    if smallest_idx == 1
        i, j = u, v
    elseif smallest_idx == 2
        i, j = v, w
    else
        i, j = w, u
    end
    shared_flag, common_vertex = edge_lies_on_two_distinct_segments(tri, i, j)
    !shared_flag && return false
    #ρ ≥ 5.0 && return true
    p, q, r = get_point(tri, i, j, common_vertex)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    ℓrp² = (px - rx)^2 + (py - ry)^2
    ℓrq² = (qx - rx)^2 + (qy - ry)^2
    if 0.999 < ℓrp² / ℓrq² < 1.001
        return true
    else
        return false
    end
end

"""
    assess_triangle_quality(tri::Triangulation, T, targets::RefinementTargets)

Assess the quality of a triangle `T` in a `Triangulation` `tri` with respect to the provided `RefinementTargets`. 

Returns `(ρ, flag)`, where `ρ` is the radius-edge ratio of the triangle and `flag` is `true` if the triangle should be refined.
"""
function assess_triangle_quality(tri::Triangulation, T, targets::RefinementTargets)
    u, v, w = indices(T)
    p, q, r = get_point(tri, u, v, w)
    ℓmin², ℓmed², ℓmax², idx = squared_triangle_lengths_and_smallest_index(p, q, r)
    A² = squared_triangle_area(ℓmin², ℓmed², ℓmax²)
    if compare_min_area(targets, A²)
        return typemax(number_type(tri)), false
    end
    A = sqrt(A²)
    r = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    ρ = r / sqrt(ℓmin²)
    area_flag = compare_area(targets, T, A, p, q, r)
    ratio_flag = compare_ratio(targets, T, ρ, p, q, r)
    seditious_flag = edge_is_seditious(tri, u, v, w, idx, ratio_flag, ρ)
    if seditious_flag
        Main.REFS[1] += 1 
    end
    return ρ, (area_flag || (ratio_flag && !seditious_flag))
end

"""
    assess_added_triangles!(tri::Triangulation, queue, events, targets)

Assess the quality of the triangles added to the triangulation `tri` after an insertion event,
as stored in `events` and add them to the `queue` if they should be refined according to the mesh 
targets defined in `targets`.
"""
function assess_added_triangles!(tri::Triangulation, queue, events, targets)
    for T in each_added_triangle(events)
        u, v, w = indices(T)
        if get_adjacent(tri, u, v) == w # It's not guaranteed that each triangle in events.added_triangles _actually_ still exists, edge flipping could have deleted it
            ρ, flag = assess_triangle_quality(tri, T, targets)
            flag && triangle_enqueue!(queue, T, ρ)
            for e in triangle_edges(T)
                if contains_constrained_edge(tri, e) && is_encroached(tri, e)
                    u, v = edge_indices(e)
                    p, q = get_point(tri, u, v)
                    px, py = getxy(p)
                    qx, qy = getxy(q)
                    ℓ² = (qx - px)^2 + (qy - py)^2
                    encroachment_enqueue!(queue, e, ℓ²)
                end
            end
        end
    end
    return nothing
end
