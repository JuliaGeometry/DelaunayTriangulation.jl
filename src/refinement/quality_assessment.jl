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

"""
    assess_triangle_quality(tri::Triangulation, T, targets::RefinementTargets)

Assess the quality of a triangle `T` in a `Triangulation` `tri` with respect to the provided `RefinementTargets`. 

Returns `(ρ, flag)`, where `ρ` is the radius-edge ratio of the triangle and `flag` is `true` if the triangle should be refined.
"""
function assess_triangle_quality(tri::Triangulation, T, targets::RefinementTargets)
    u, v, w = indices(T)
    p, q, r = get_point(tri, u, v, w)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    r = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    ρ = r / sqrt(ℓmin²)
    area_flag = compare_area(targets, T, A, p, q, r)
    ratio_flag = compare_ratio(targets, T, ρ, p, q, r)
    return ρ, (area_flag || ratio_flag)
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
