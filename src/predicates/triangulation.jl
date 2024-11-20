@inline @optarg1 DEFAULT_KERNEL function triangle_orientation(kernel::AbstractPredicateKernel, tri::Triangulation, i, j, k; ctr=get_predicate_diagnostics(tri))
    if is_exterior_ghost_triangle(tri, i, j, k)
        i, j, k = k, j, i
    end
    p, q, r = get_point(tri, i, j, k)
    return triangle_orientation(kernel, p, q, r; ctr)
end

@inline @optarg1 DEFAULT_KERNEL function triangle_orientation(kernel::AbstractPredicateKernel, tri::Triangulation, T; ctr=get_predicate_diagnostics(tri))
    i, j, k = triangle_vertices(T)
    return triangle_orientation(kernel, tri, i, j, k; ctr)
end

@inline @optarg1 DEFAULT_KERNEL function point_position_relative_to_line(kernel::AbstractPredicateKernel, tri::Triangulation, i, j, u; ctr=get_predicate_diagnostics(tri))
    a, b, p = get_point(tri, i, j, u)
    if is_exterior_ghost_edge(tri, i, j)
        return point_position_relative_to_line(kernel, b, a, p; ctr)
    else
        return point_position_relative_to_line(kernel, a, b, p; ctr)
    end
end

@inline @optarg1 DEFAULT_KERNEL function point_closest_to_line(kernel::AbstractPredicateKernel, tri::Triangulation, i, j, u, v; ctr=get_predicate_diagnostics(tri))
    a, b, p, q = get_point(tri, i, j, u, v)
    return point_closest_to_line(kernel, a, b, p, q; ctr)
end

@inline function point_position_on_line_segment(tri::Triangulation, i, j, u; ctr=get_predicate_diagnostics(tri))
    a, b, p = get_point(tri, i, j, u)
    return point_position_on_line_segment(a, b, p; ctr)
end

@inline @optarg1 DEFAULT_KERNEL function line_segment_intersection_type(kernel::AbstractPredicateKernel, tri::Triangulation, u, v, i, j; ctr=get_predicate_diagnostics(tri))
    p, q, a, b = get_point(tri, u, v, i, j)
    return line_segment_intersection_type(kernel, p, q, a, b; ctr)
end

@optarg1 DEFAULT_KERNEL function point_position_relative_to_triangle(kernel::AbstractPredicateKernel, tri::Triangulation, i, j, k, u; ctr=get_predicate_diagnostics(tri))
    if !is_exterior_ghost_triangle(tri, i, j, k)
        a, b, c, p = get_point(tri, i, j, k, u)
        return point_position_relative_to_triangle(kernel, a, b, c, p; ctr)
    else
        i, j, k = sort_triangle(i, j, k)
        a, b, c, p = get_point(tri, i, j, k, u) # a and b are solid vertices, but c is a ghost (mapped to a centroid in the opposite direction)
        edge_ab = point_position_relative_to_line(kernel, a, b, p; ctr)
        is_right(edge_ab) && return Outside
        if is_collinear(edge_ab)
            cert = point_position_on_line_segment(a, b, p; ctr)
            return (is_left(cert) || is_right(cert)) ? Outside : On
        end
        edge_bc = point_position_relative_to_line(kernel, c, b, p; ctr) # Flipped to match centroid location
        is_right(edge_bc) && return Outside
        #is_collinear(edge_bc) && return On # Don't need to check that it's not on the (c, b) part, since we already know we're to the left of (a, b) at this point
        edge_ca = point_position_relative_to_line(kernel, c, a, p; ctr)
        is_left(edge_ca) && return Outside
        #is_collinear(edge_ca) && return On
        return Inside

        # The collinear tests were deleted for the ghost edges. It doesn't really make much sense to see if a point is 
        # on the ghost edges. It's not like we can do anything with that information, and if we are using it then 
        # there's no point distinguishing between the two adjacent ghost triangles in that case.
    end
end

@inline @optarg1 DEFAULT_KERNEL function point_position_relative_to_triangle(kernel::AbstractPredicateKernel, tri::Triangulation, T, u; ctr=get_predicate_diagnostics(tri))
    i, j, k = triangle_vertices(T)
    return point_position_relative_to_triangle(kernel, tri, i, j, k, u; ctr)
end

@optarg1 DEFAULT_KERNEL function is_legal(kernel::AbstractPredicateKernel, tri::Triangulation, i, j; cache::PredicateCacheType=nothing, ctr=get_predicate_diagnostics(tri))
    if contains_segment(tri, i, j) ||
       is_boundary_edge(tri, j, i) || is_boundary_edge(tri, i, j) ||
       !edge_exists(tri, i, j) || !edge_exists(tri, j, i) ||
       is_ghost_edge(i, j)
        add_is_legal!(ctr)
        return Legal
    else
        k = get_adjacent(tri, i, j)
        ℓ = get_adjacent(tri, j, i)
        p, q, r, s = get_point(tri, i, j, k, ℓ)
        cert = is_legal(kernel, p, q, r, s; cache, ctr)
        return cert
    end
end

@inline @optarg1 DEFAULT_KERNEL function triangle_line_segment_intersection(kernel::AbstractPredicateKernel, tri::Triangulation, i, j, k, u, v; ctr=get_predicate_diagnostics(tri))
    p, q, r, a, b = get_point(tri, i, j, k, u, v)
    return triangle_line_segment_intersection(kernel, p, q, r, a, b; ctr)
end

@optarg1 DEFAULT_KERNEL function find_edge(kernel::AbstractPredicateKernel, tri::Triangulation, T, ℓ; ctr=get_predicate_diagnostics(tri))
    add_find_edge!(ctr)
    T′ = sort_triangle(T) # sort so that, if T is a ghost, then (i, j) is the solid edge
    i, j, k = triangle_vertices(T′)
    a, b, c = get_point(tri, i, j, k)
    E = edge_type(tri)
    if is_ghost_vertex(k) # must be on the solid edge
        return construct_edge(E, i, j)
    else
        is_collinear(point_position_relative_to_line(kernel, tri, a, b, ℓ; ctr)) && return construct_edge(E, i, j)
        is_collinear(point_position_relative_to_line(kernel, tri, b, c, ℓ; ctr)) && return construct_edge(E, j, k)
        return construct_edge(E, k, i)
    end
end

@optarg1 DEFAULT_KERNEL function point_position_relative_to_circumcircle(kernel::AbstractPredicateKernel, tri::Triangulation, i, j, k, ℓ; cache::PredicateCacheType=nothing, ctr=get_predicate_diagnostics(tri))
    add_point_position_relative_to_circumcircle!(ctr)
    u, v, w = sort_triangle(i, j, k)
    a, b, c, p = get_point(tri, u, v, w, ℓ)
    if is_ghost_vertex(w)
        cert = point_position_relative_to_oriented_outer_halfplane(kernel, a, b, p; ctr)
        if is_on(cert) && is_weighted(tri)
            u′, v′, w′ = replace_ghost_triangle_with_boundary_triangle(tri, (u, v, w))
            !edge_exists(w′) && return Inside # needed for the case tri = triangulate(get_points(triangulate_rectangle(0, 10, 0, 10, 3, 3)), weights = zeros(9), insertion_order = [9, 7, 6, 8, 4, 5, 3, 1, 2]) 
            sub_cert = point_position_relative_to_witness_plane(kernel, tri, u′, v′, w′, ℓ; cache=fix_orient3_cache(tri, cache), ctr)
            is_above(sub_cert) && return Outside
            return cert
        else
            return cert
        end
    elseif !is_weighted(tri)
        return point_position_relative_to_circle(kernel, a, b, c, p; cache=fix_incircle_cache(tri, cache), ctr)
    else
        cert = point_position_relative_to_witness_plane(kernel, tri, i, j, k, ℓ; cache=fix_orient3_cache(tri, cache), ctr)
        return is_above(cert) ? Outside :
               is_below(cert) ? Inside :
               On
    end
end

@inline @optarg1 DEFAULT_KERNEL function point_position_relative_to_circumcircle(kernel::AbstractPredicateKernel, tri::Triangulation, T, ℓ; cache::PredicateCacheType=nothing, ctr=get_predicate_diagnostics(tri))
    i, j, k = triangle_vertices(T)
    return point_position_relative_to_circumcircle(kernel, tri, i, j, k, ℓ; cache, ctr)
end

@optarg1 DEFAULT_KERNEL function point_position_relative_to_witness_plane(kernel::AbstractPredicateKernel, tri::Triangulation, i, j, k, ℓ; cache::PredicateCacheType=nothing, ctr=get_predicate_diagnostics(tri))
    add_point_position_relative_to_witness_plane!(ctr)
    p⁺ = get_lifted_point(tri, i)
    q⁺ = get_lifted_point(tri, j)
    r⁺ = get_lifted_point(tri, k)
    a⁺ = get_lifted_point(tri, ℓ)
    cert = orient_predicate(kernel, p⁺, q⁺, r⁺, a⁺; cache, ctr)
    return convert_certificate(cert, Above, On, Below)
end

@inline @optarg1 function test_visibility(kernel::AbstractPredicateKernel, tri::Triangulation, q, i, ctr=get_predicate_diagnostics(tri))
    add_test_visibility!(ctr)
    _, invisible_flag = find_triangle(tri, q; use_barriers=Val(true), k=i, concavity_protection=true, predicates=kernel, predicate_diagnostics=ctr, point_location_diagnostics=get_point_location_diagnostics(tri))
    if invisible_flag
        return Invisible
    else
        return Visible
    end
end

@optarg1 DEFAULT_KERNEL function test_visibility(kernel::AbstractPredicateKernel, tri::Triangulation, u, v, i; shift=0.0, attractor=get_point(tri, i), ctr=get_predicate_diagnostics(tri))
    add_test_visibility!(ctr)
    pu, pv = get_point(tri, u, v)
    pux, puy = getxy(pu)
    pvx, pvy = getxy(pv)
    qx, qy = getxy(attractor)
    ts = LinRange(0.00001, 0.99999, 10)
    for t in ts
        mx, my = pux + t * (pvx - pux), puy + t * (pvy - puy)
        m̃ = (mx + shift * (qx - mx), my + shift * (qy - my))
        cert = test_visibility(kernel, tri, m̃, i; ctr)
        is_visible(cert) && return Visible
    end
    return Invisible
end