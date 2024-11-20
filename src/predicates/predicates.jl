@inline @optarg1 DEFAULT_KERNEL function triangle_orientation(kernel::AbstractPredicateKernel, p, q, r; ctr)
    add_triangle_orientation!(ctr)
    cert = orient_predicate(kernel, p, q, r; ctr)
    return convert_certificate(cert, NegativelyOriented, Degenerate, PositivelyOriented)
end

@inline @optarg1 DEFAULT_KERNEL function point_position_relative_to_circle(kernel::AbstractPredicateKernel, a, b, c, p; cache::PredicateCacheType = nothing, ctr)
    add_point_position_relative_to_circle!(ctr)
    cert = incircle_predicate(kernel, a, b, c, p; cache, ctr)
    return convert_certificate(cert, Outside, On, Inside)
end

@inline @optarg1 DEFAULT_KERNEL function point_position_relative_to_line(kernel::AbstractPredicateKernel, a, b, p; ctr)
    add_point_position_relative_to_line!(ctr)
    cert = orient_predicate(kernel, a, b, p; ctr)
    return convert_certificate(cert, Right, Collinear, Left)
end

@inline @optarg1 DEFAULT_KERNEL function point_closest_to_line(kernel::AbstractPredicateKernel, a, b, p, q; ctr)
    add_point_closest_to_line!(ctr)
    cert = parallelorder_predicate(kernel, a, b, q, p; ctr)
    return convert_certificate(cert, Closer, Equidistant, Further)
end

@inline function point_position_on_line_segment(a, b, p; ctr)
    add_point_position_on_line_segment!(ctr)
    cert = sameside_predicate(a, b, p; ctr)
    converted_cert = convert_certificate(cert, On, Degenerate, Outside)
    if is_outside(converted_cert) # Do we have "a ---- b ---- p" or "p ---- a ---- b"?
        ap_cert = sameside_predicate(a, p, b; ctr)
        converted_ap_cert = convert_certificate(ap_cert, On, Degenerate, Outside)
        return is_on(converted_ap_cert) ? Right : Left
    end
    return converted_cert
end

@optarg1 DEFAULT_KERNEL function line_segment_intersection_type(kernel::AbstractPredicateKernel, p, q, a, b; ctr)
    add_line_segment_intersection_type!(ctr)
    cert = meet_predicate(kernel, p, q, a, b; ctr)
    converted_cert = convert_certificate(cert, None, Multiple, Single)
    if is_multiple(converted_cert)
        #=
        We need to check if the "multiple" really means one endpoint is touching the line, but 
        only one part of it. In particular, we are checking the scenarios (visualising just some of the possibilities):

         p --- a --- b --- q (Multiple)

               b 
               |
               |
         p --- a --- q       (On) 

               p 
               |
               |
         a --- q --- b       (On) 

        To check this, we can check if `(a, b)` is collinear with `p`, but not with `q`, as 
        this implies only a single intersection (if both are collinear, they are on the same 
        line, so there would be multiple). We need to also do the same for `(p, q)`.
        =#
        collinear_cert_1 = point_position_relative_to_line(kernel, a, b, p; ctr)
        collinear_cert_2 = point_position_relative_to_line(kernel, a, b, q; ctr)
        if (is_collinear(collinear_cert_1) && !is_collinear(collinear_cert_2)) || (!is_collinear(collinear_cert_1) && is_collinear(collinear_cert_2))
            return Touching
        end
        collinear_cert_3 = point_position_relative_to_line(kernel, p, q, a; ctr)
        collinear_cert_4 = point_position_relative_to_line(kernel, p, q, b; ctr)
        if (is_collinear(collinear_cert_3) && !is_collinear(collinear_cert_4)) || (!is_collinear(collinear_cert_3) && is_collinear(collinear_cert_4))
            return Touching
        end
    end
    return converted_cert
end

@optarg1 DEFAULT_KERNEL function point_position_relative_to_triangle(kernel::AbstractPredicateKernel, a, b, c, p; ctr)
    add_point_position_relative_to_triangle!(ctr)
    edge_ab = point_position_relative_to_line(kernel, a, b, p; ctr)
    edge_bc = point_position_relative_to_line(kernel, b, c, p; ctr)
    edge_ca = point_position_relative_to_line(kernel, c, a, p; ctr)
    all_edge_certs = (edge_ab, edge_bc, edge_ca)

    # In what follows, we are unrolling the loop "for (edge_cert, (a, b)) in zip((edge_ab, edge_bc, edge_ca), ((a, b), (b, c), (c, a)))"
    edge_cert, a, b = edge_ab, a, b
    if is_collinear(edge_cert) # Just being collinear isn't enough - what if it's collinear but off the edge?
        cert = point_position_on_line_segment(a, b, p; ctr)
        if is_left(cert) || is_right(cert)
            return Outside
        else
            return On
        end
    end

    edge_cert, a, b = edge_bc, b, c
    if is_collinear(edge_cert) # Just being collinear isn't enough - what if it's collinear but off the edge?
        cert = point_position_on_line_segment(a, b, p; ctr)
        if is_left(cert) || is_right(cert)
            return Outside
        else
            return On
        end
    end

    edge_cert, a, b = edge_ca, c, a
    if is_collinear(edge_cert) # Just being collinear isn't enough - what if it's collinear but off the edge?
        cert = point_position_on_line_segment(a, b, p; ctr)
        if is_left(cert) || is_right(cert)
            return Outside
        else
            return On
        end
    end

    all(is_left, all_edge_certs) && return Inside
    return Outside
end

@optarg1 DEFAULT_KERNEL function point_position_relative_to_oriented_outer_halfplane(kernel::AbstractPredicateKernel, a, b, p; ctr)
    add_point_position_relative_to_oriented_outer_halfplane!(ctr)
    in_open_halfplane = point_position_relative_to_line(kernel, a, b, p; ctr)
    if is_collinear(in_open_halfplane)
        is_on_boundary_edge = point_position_on_line_segment(a, b, p; ctr)
        if is_on(is_on_boundary_edge) || is_degenerate(is_on_boundary_edge)
            return On
        else
            return Outside
        end
    elseif is_left(in_open_halfplane)
        return Inside
    else
        return Outside
    end
end

@inline @optarg1 DEFAULT_KERNEL function is_legal(kernel::AbstractPredicateKernel, p, q, r, s; cache::PredicateCacheType = nothing, ctr)
    add_is_legal!(ctr)
    incirc = point_position_relative_to_circle(kernel, p, q, r, s; cache, ctr)
    if is_inside(incirc)
        return Illegal
    else
        return Legal
    end
end

@optarg1 DEFAULT_KERNEL function triangle_line_segment_intersection(kernel::AbstractPredicateKernel, p, q, r, a, b; ctr)
    add_triangle_line_segment_intersection!(ctr)

    ## Break into cases:
    # Case I: (a, b) is an edge of (p, q, r)
    a_eq = any(==(a), (p, q, r))
    b_eq = any(==(b), (p, q, r))
    if a_eq && b_eq
        return Touching
    end

    # Case II: a or b is one of the vertices of (p, q, r)
    if a_eq || b_eq
        a, b = a_eq ? (a, b) : (b, a) # rotate so that a is one of the vertices 
        p, q, r = choose_uvw(getxy(a) == getxy(p), getxy(a) == getxy(q), getxy(a) == getxy(r), p, q, r) # rotate so that a == p
        # First, find position of b relative to the edges pq and pr 
        pqb = point_position_relative_to_line(kernel, p, q, b; ctr)
        prb = point_position_relative_to_line(kernel, p, r, b; ctr)
        if is_right(pqb) || is_left(prb)
            # If b is right of pq or left of pr, it is away from the interior and thus outside
            return Outside
        elseif is_collinear(pqb) || is_collinear(prb)
            # If b is collinear with an edge, it could be on the edge or away from it
            s = is_collinear(pqb) ? q : r
            b_pos = point_position_on_line_segment(p, s, b; ctr)
            return is_left(b_pos) ? Outside : Touching
        else
            # If we are here, it means that b is between the rays through pq and pr. 
            # We just need to determine if b is inside or the triangle or outside. 
            # See that if b is inside, then it is left of qr, but if it is outside 
            # then it is right of qr. 
            qrb = point_position_relative_to_line(kernel, q, r, b; ctr)
            if is_right(qrb)
                return Multiple # Technically, one could call this Single, but Multiple is the most appropriate
            else
                return Inside # Also covers the is_on(qrb) case 
            end
        end
    end

    # We have now covered the cases where a or b, or both, are a vertex of (p, q, r).
    # The remaining degenerate case is the case where a or b, or both, are on an edge 
    # of (p, q, r). We can determine this possibility with the following certificates.
    a_tri = point_position_relative_to_triangle(kernel, p, q, r, a; ctr)
    b_tri = point_position_relative_to_triangle(kernel, p, q, r, b; ctr)

    # Case III: a and b are both on an edge of (p, q, r)
    if is_on(a_tri) && is_on(b_tri)
        # There are two options here, either a and b are both on the same edge,
        # or they are on different edges. In the latter case, (a, b) is entirely in 
        # (p, q, r). To determine this case, we just test against each edge.
        pqa = is_collinear(point_position_relative_to_line(kernel, p, q, a; ctr))
        qra = pqa ? false : is_collinear(point_position_relative_to_line(kernel, q, r, a; ctr)) # skip computation if we already know a result
        #rpa = (pqa || qra) ? false : is_collinear(point_position_relative_to_line(r, p, a))
        pqb = is_collinear(point_position_relative_to_line(kernel, p, q, b; ctr))
        qrb = pqb ? false : is_collinear(point_position_relative_to_line(kernel, q, r, b; ctr))
        #rpb = (pqb || qrb) ? false : is_collinear(point_position_relative_to_line(r, p, b))
        # Don't need rpa or rpb above, since checking the two edges automatically implies the result for the third edge
        a_edge = pqa ? (p, q) : (qra ? (q, r) : (r, p))
        b_edge = pqb ? (p, q) : (qrb ? (q, r) : (r, p))
        return a_edge == b_edge ? Touching : Inside
    end

    # Case IV: a or b are on an edge of (p, q, r)
    # For this case, the working is essentially the same as in Case II.
    if is_on(a_tri) || is_on(b_tri)
        # Rotate so that a is on (p, q)
        a, b = is_on(a_tri) ? (a, b) : (b, a)
        pqa = is_collinear(point_position_relative_to_line(kernel, p, q, a; ctr))
        qra = pqa ? false : is_collinear(point_position_relative_to_line(kernel, q, r, a; ctr))
        rpa = !(pqa || qra)
        p, q, r = choose_uvw(pqa, qra, rpa, p, q, r)
        # Similar to Case II, we start by considering the position of b relative to pq and qr 
        prb = point_position_relative_to_line(kernel, p, r, b; ctr)
        qrb = point_position_relative_to_line(kernel, q, r, b; ctr)
        # First, if we are left of right of pr or qr, respectively, we must have exited the triangle
        if is_left(prb) || is_right(qrb)
            return Multiple
        end
        # It is also possible that b is collinear with pq 
        pqb = point_position_relative_to_line(kernel, p, q, b; ctr)
        if is_collinear(pqb)
            #b_pos = point_position_on_line_segment(p, q, b)
            # We don't really need to know the position, we already know now that ab touches the triangle. 
            return Touching
            # Alternatively, b might be away from pq, meaning it is outside:
        elseif is_right(pqb)
            return Outside
        else
            # We now know that b is left of pqb, which from the above must imply that b is in the triangle. 
            return Inside
        end
    end

    # We can now assume that neither a nor b appear anywhere on (p, q, r). Just to make things a bit simpler, 
    # let us determine directly where a and b are relative to the triangle 
    a_pos = point_position_relative_to_triangle(kernel, p, q, r, a; ctr)
    b_pos = point_position_relative_to_triangle(kernel, p, q, r, b; ctr)

    # Case V: a and b are inside the triangle 
    if is_inside(a_pos) && is_inside(b_pos)
        return Inside
    elseif is_inside(a_pos) || is_inside(b_pos)
        # Case VI: Only one of a and b are inside the triangle 
        return Single
    else # is_outside(a_pos) && is_outside(b_pos)
        # Case VII: Both a and b are outside the triangle 
        # We consider the intersection of ab with each edge 
        pq_ab = line_segment_intersection_type(kernel, p, q, a, b; ctr)
        qr_ab = line_segment_intersection_type(kernel, q, r, a, b; ctr)
        rp_ab = line_segment_intersection_type(kernel, r, p, a, b; ctr)
        num_none = count(is_none, (pq_ab, qr_ab, rp_ab))
        num_touching = count(is_touching, (pq_ab, qr_ab, rp_ab))
        num_multiple = count(is_multiple, (pq_ab, qr_ab, rp_ab))
        # The possible cases:
        #   1. There are no intersections with any edge: num_none == 3
        #   2. ab only touches a single edge with an endpoint, but doesn't go through it. Fortunately, we tested for this case already. 
        #   3. It is possible that ab does not touch an edge with an endpoint, but its interior touches one of the points of the triangle, but doesn't go inside the triangle: num_touching == 2
        #   4. Another case is that ab is completely collinear with an edge. This happens when num_multiple == 1. This case is hard to distinguish from 3, so to handle it we need to do some more orientation tests.
        #   5. In all other cases, ab goes through the triangle and out the other side
        if num_none == 3
            return Outside
        elseif num_multiple == 1 || num_touching == 2
            # This is for cases 3 or 4 above. We need to see how ab line up with the vertices of (p, q, r)
            abp = point_position_relative_to_line(kernel, a, b, p; ctr)
            abq = point_position_relative_to_line(kernel, a, b, q; ctr)
            abr = point_position_relative_to_line(kernel, a, b, r; ctr)
            if count(is_collinear, (abp, abq, abr)) == 2
                return Touching # lines up with an edge 
            else
                return Outside # interior only just touches a single vertex 
            end
        else
            return Multiple
        end
    end
end

@inline @optarg1 DEFAULT_KERNEL function opposite_angle(kernel::AbstractPredicateKernel, p, q, r; ctr) # https://math.stackexchange.com/a/701682/861404 
    add_opposite_angle!(ctr)
    cert = angle_is_acute_predicate(kernel, p, q, r; ctr)
    return convert_certificate(cert, Obtuse, Right, Acute)
end

@inline @optarg1 DEFAULT_KERNEL function point_position_relative_to_diametral_circle(kernel::AbstractPredicateKernel, p, q, r; ctr)
    add_point_position_relative_to_diametral_circle!(ctr)
    d = opposite_angle(kernel, p, q, r; ctr)
    if is_obtuse(d)
        return Inside
    elseif is_right(d)
        return On
    else
        return Outside
    end
end

@optarg1 DEFAULT_KERNEL function point_position_relative_to_diametral_lens(p, q, r, lens_angle = 30.0; ctr)
    add_point_position_relative_to_diametral_lens!(ctr)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    ddot = inp(px - rx, py - ry, qx - rx, qy - ry) # (p - r) · (q - r)
    ddot > 0 && return Outside # no need to check - not obtuse
    # We could use an exact predicate for the above, but it would be too slow 
    # since we'd be essentially computing ddot twice for obtuse angles (since we need 
    # its square below).
    ddot² = ddot^2
    ℓpr² = dist_sqr(p, q)
    ℓqr² = dist_sqr(q, r)
    cosine_scale² = cosd(2lens_angle)^2 # ≥ 0 since lens_angle ∈ [0°, 45°]
    rhs = (ℓpr² * ℓqr²) * cosine_scale²
    if ddot² > rhs
        return Inside
    elseif ddot² == rhs
        return On
    else
        return Outside
    end
end