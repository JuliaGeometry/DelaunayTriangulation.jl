@doc """
    triangle_orientation(tri::Triangulation, i, j, k) -> Certificate 
    triangle_orientation(tri::Triangulation, T) -> Certificate
    triangle_orientation(p, q, r) -> Certificate

Computes the orientation of the triangle `T = (i, j, k)` with correspondig coordinates `(p, q, r)`. The returned value is a [`Certificate`](@ref), which is one of:

- `PositivelyOriented`: The triangle is positively oriented.
- `Degenerate`: The triangle is degenerate, meaning the coordinates are collinear.
- `NegativelyOriented`: The triangle is negatively oriented.
"""
triangle_orientation
function triangle_orientation(p, q, r)
    cert = orient_predicate(p, q, r)
    return convert_certificate(cert, Cert.NegativelyOriented, Cert.Degenerate,
        Cert.PositivelyOriented)
end
function triangle_orientation(tri::Triangulation, i, j, k)
    if is_exterior_ghost_triangle(tri, i, j, k)
        i, j, k = k, j, i
    end
    p, q, r = get_point(tri, i, j, k)
    return triangle_orientation(p, q, r)
end
triangle_orientation(tri::Triangulation, T) = triangle_orientation(tri, geti(T), getj(T), getk(T))

"""
    point_position_relative_to_circle(a, b, c, p) -> Certificate

Given a circle through the coordinates `(a, b, c)`, assumed to be positively oriented,
computes the position of `p` relative to the circle. Returns a [`Certificate`](@ref), which is one of:

- `Inside`: `p` is inside the circle.
- `On`: `p` is on the circle.
- `Outside`: `p` is outside the triangle.
"""
function point_position_relative_to_circle(a, b, c, p)
    cert = incircle_predicate(a, b, c, p)
    return convert_certificate(cert, Cert.Outside, Cert.On, Cert.Inside)
end

@doc """
    point_position_relative_to_line(a, b, p) -> Certificate
    point_position_relative_to_line(tri::Triangulation, i, j, u) -> Certificate

Tests the position of `p` (or the vertex `u` of `tri`) relative to the edge `(a, b)` (or the edge with vertices `(i, j)` of `tri`), returning 
a [`Certificate`](@ref) which is one of:

- `Left`: `p` is to the left of the line.
- `Collinear`: `p` is on the line.
- `Right`: `p` is to the right of the line.
"""
point_position_relative_to_line
function point_position_relative_to_line(a, b, p)
    cert = orient_predicate(a, b, p)
    return convert_certificate(cert, Cert.Right, Cert.Collinear, Cert.Left)
end
function point_position_relative_to_line(tri::Triangulation, i, j, u)
    a, b, p = get_point(tri, i, j, u)
    if is_exterior_ghost_edge(tri, i, j)
        return point_position_relative_to_line(b, a, p)
    else
        return point_position_relative_to_line(a, b, p)
    end
end

@doc """
    point_closest_to_line(a, b, p, q) -> Certificate
    point_closest_to_line(tri::Triangulation, i, j, u, v) -> Certificate

Given a line `ℓ` through `(a, b)` (or through the vertices `(i, j)`), tests if `p` (or the vertex `u`) is closer to `ℓ`
than `q` (or the vertex `v`), assuming that `p` and `q` are to the left of `ℓ`, returning a [`Certificate`](@ref) which is one of:

- `Closer`: `p` is closer to `ℓ`.
- `Further`: `q` is closer to `ℓ`.
- `Equidistant`: `p` and `q` are the same distance from `ℓ`.
"""
point_closest_to_line
function point_closest_to_line(a, b, p, q)
    cert = parallelorder_predicate(a, b, q, p)
    return convert_certificate(cert, Cert.Closer, Cert.Equidistant, Cert.Further)
end
function point_closest_to_line(tri::Triangulation, i, j, u, v)
    a, b, p, q = get_point(tri, i, j, u, v)
    return point_closest_to_line(a, b, p, q)
end

@doc """
    point_position_on_line_segment(a, b, p) -> Certificate
    point_position_on_line_segment(tri::Triangulation, i, j, u) -> Certificate

Given a point `p` (or vertex `u`) and the line segment `(a, b)` (or edge `(i, j)`), assuming `p`
to be collinear with `a` and `b`, computes the position of `p` relative to the line segment. The returned value is a 
[`Certificate`](@ref) which is one of:

- `On`: `p` is on the line segment, meaning between `a` and `b`.
- `Degenerate`: Either `p == a` or `p == b`, i.e. `p` is one of the endpoints.
- `Left`: `p` is off and to the left of the line segment.
- `Right`: `p` is off and to the right of the line segment.
"""
point_position_on_line_segment
function point_position_on_line_segment(a, b, p)
    cert = sameside_predicate(a, b, p)
    converted_cert = convert_certificate(cert, Cert.On, Cert.Degenerate, Cert.Outside)
    if is_outside(converted_cert) # Do we have "a ---- b ---- p" or "p ---- a ---- b"?
        ap_cert = sameside_predicate(a, p, b)
        converted_ap_cert = convert_certificate(ap_cert, Cert.On, Cert.Degenerate, Cert.Outside)
        return is_on(converted_ap_cert) ? Cert.Right : Cert.Left
    end
    return converted_cert
end
function point_position_on_line_segment(tri::Triangulation, i, j, u)
    a, b, p = get_point(tri, i, j, u)
    return point_position_on_line_segment(a, b, p)
end

@doc """
    line_segment_intersection_type(p, q, a, b) -> Certificate 
    line_segment_intersection_type(tri::Triangulation, u, v, i, j) -> Certificate

Given the coordinates `(p, q)` (or vertices `(u, v)`) and `(a, b)` (or vertices `(i, j)`) defining two line segments,
tests the number of intersections between the two segments. The returned value is a [`Certificate`](@ref), which is one of:

- `None`: The line segments do not meet at any points.
- `Multiple`: The closed line segments `[p, q]` and `[a, b]` meet in one or several points.
- `Single`: The open line segments `(p, q)` and `(a, b)` meet in a single point.
- `Touching`: One of the endpoints is on `[a, b]`, but there are no other intersections.
"""
line_segment_intersection_type
function line_segment_intersection_type(p, q, a, b)
    cert = meet_predicate(p, q, a, b)
    converted_cert = convert_certificate(cert, Cert.None, Cert.Multiple, Cert.Single)
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
        collinear_cert_1 = point_position_relative_to_line(a, b, p)
        collinear_cert_2 = point_position_relative_to_line(a, b, q)
        if (is_collinear(collinear_cert_1) && !is_collinear(collinear_cert_2)) || (!is_collinear(collinear_cert_1) && is_collinear(collinear_cert_2))
            return Cert.Touching
        end
        collinear_cert_3 = point_position_relative_to_line(p, q, a)
        collinear_cert_4 = point_position_relative_to_line(p, q, b)
        if (is_collinear(collinear_cert_3) && !is_collinear(collinear_cert_4)) || (!is_collinear(collinear_cert_3) && is_collinear(collinear_cert_4))
            return Cert.Touching
        end
    end
    return converted_cert
end
function line_segment_intersection_type(tri::Triangulation, u, v, i, j)
    p, q, a, b = get_point(tri, u, v, i, j)
    return line_segment_intersection_type(p, q, a, b)
end

@doc """
    point_position_relative_to_triangle(a, b, c, p) -> Certificate
    point_position_relative_to_triangle(tri::Triangulation, i, j, k, u) -> Certificate
    point_position_relative_to_triangle(tri::Triangulation, T, u) -> Certificate

Given a positively oriented triangle with coordinates `(a, b, c)` (or triangle `T = (i, j, k)` of `tri`), computes the
position of `p` (or vertex `u`) relative to the triangle. The returned value is a [`Certificate`](@ref), which is one of:

- `Outside`: `p` is outside of the triangle.
- `On`: `p` is on one of the edges.
- `Inside`: `p` is inside the triangle.

!!! note "Ghost triangles"

    If `T` is a ghost triangle, then it is not checked if the point is on any of the ghost edges,
"""
point_position_relative_to_triangle
function point_position_relative_to_triangle(a, b, c, p)
    edge_ab = point_position_relative_to_line(a, b, p)
    edge_bc = point_position_relative_to_line(b, c, p)
    edge_ca = point_position_relative_to_line(c, a, p)
    all_edge_certs = (edge_ab, edge_bc, edge_ca)
    for (edge_cert, (a, b)) in zip((edge_ab, edge_bc, edge_ca), ((a, b), (b, c), (c, a)))
        if is_collinear(edge_cert) # Just being collinear isn't enough - what if it's collinear but off the edge?
            cert = point_position_on_line_segment(a, b, p)
            if is_left(cert) || is_right(cert)
                return Cert.Outside
            else
                return Cert.On
            end
        end
    end
    all(is_left, all_edge_certs) && return Cert.Inside
    return Cert.Outside
end
function point_position_relative_to_triangle(tri::Triangulation, i, j, k, u)
    if !is_exterior_ghost_triangle(tri, i, j, k)
        a, b, c, p = get_point(tri, i, j, k, u)
        return point_position_relative_to_triangle(a, b, c, p)
    else
        i, j, k = sort_triangle(i, j, k)
        a, b, c, p = get_point(tri, i, j, k, u) # a and b are solid vertices, but c is a ghost (mapped to a centroid in the opposite direction)
        edge_ab = point_position_relative_to_line(a, b, p)
        is_right(edge_ab) && return Cert.Outside
        if is_collinear(edge_ab)
            cert = point_position_on_line_segment(a, b, p)
            return (is_left(cert) || is_right(cert)) ? Cert.Outside : Cert.On
        end
        edge_bc = point_position_relative_to_line(c, b, p) # Flipped to match centroid location
        is_right(edge_bc) && return Cert.Outside
        #is_collinear(edge_bc) && return Cert.On # Don't need to check that it's not on the (c, b) part, since we already know we're to the left of (a, b) at this point
        edge_ca = point_position_relative_to_line(c, a, p)
        is_left(edge_ca) && return Cert.Outside
        #is_collinear(edge_ca) && return Cert.On
        return Cert.Inside

        # The collinear tests were deleted for the ghost edges. It doesn't really make much sense to see if a point is 
        # on the ghost edges. It's not like we can do anything with that information, and if we are using it then 
        # there's no point distinguishing between the two adjacent ghost triangles in that case.
    end
end
point_position_relative_to_triangle(tri::Triangulation, T, u) = point_position_relative_to_triangle(tri, geti(T), getj(T), getk(T), u)

"""
    point_position_relative_to_oriented_outer_halfplane(a, b, p) -> Certificate 

Given an edge with coordinates `(a, b)` and a point `p`, 
tests the position of `p` relative to the oriented outer halfplane defined
by `(a, b)`. The returned value is a [`Certificate`](@ref), which is one of:

- `Outside`: `p` is outside of the oriented outer halfplane, meaning to the right of the line `(a, b)` or collinear with `a` and `b` but not on the line segment `(a, b)`.
- `On`: `p` is on the line segment `[a, b]`.
- `Inside`: `p` is inside of the oriented outer halfplane, meaning to the left of the line `(a, b)`.

# Extended help 
The oriented outer halfplane is the union of the open halfplane defined by the region to the left of the oriented line `(a, b)`, and the open line segment `(a, b)`.
"""
function point_position_relative_to_oriented_outer_halfplane(a, b, p)
    in_open_halfplane = point_position_relative_to_line(a, b, p)
    if is_collinear(in_open_halfplane)
        is_on_boundary_edge = point_position_on_line_segment(a, b, p)
        if is_on(is_on_boundary_edge) || is_degenerate(is_on_boundary_edge)
            return Cert.On
        else
            return Cert.Outside
        end
    elseif is_left(in_open_halfplane)
        return Cert.Inside
    else
        return Cert.Outside
    end
end

@doc """
    is_legal(tri::Triangulation, i, j) -> Certificate 
    is_legal(p, q, r, s) -> Certificate

Tests if the edge `(p, q)` (or the edge `(i, j)` of `tri`) is legal, where the edge `(p, q)`
is incident to two triangles `(p, q, r)` and `(q, p, s)`. In partiuclar, tests that `s` is not inside 
the triangle through `(p, q, r)`. The returned value is a [`Certificate`](@ref), which is one of:

- `Legal`: The edge `(p, q)` is legal.
- `Illegal`: The edge `(p, q)` is illegal.

If the edge `(i, j)` is a segment of `tri` or is a ghost edge, then the edge is always legal.
"""
is_legal
function is_legal(p, q, r, s)
    incirc = point_position_relative_to_circle(p, q, r, s)
    if is_inside(incirc)
        return Cert.Illegal
    else
        return Cert.Legal
    end
end

function is_legal(tri::Triangulation, i, j)
    (contains_segment(tri, i, j) ||
     is_boundary_edge(tri, j, i) ||
     is_boundary_edge(tri, i, j) ||
     !edge_exists(tri, i, j) ||
     !edge_exists(tri, j, i) ||
     is_ghost_edge(i, j)) && return Cert.Legal
    k = get_adjacent(tri, i, j)
    ℓ = get_adjacent(tri, j, i)
    p, q, r, s = get_point(tri, i, j, k, ℓ)
    cert = is_legal(p, q, r, s)
    return cert
end

@doc """
    triangle_line_segment_intersection(p, q, r, a, b) -> Certificate 
    triangle_line_segment_intersection(tri::Triangulation, i, j, k, u, v) -> Certificate

Classifies the intersection of the line segment `(a, b)` (or the edge `(u, v)` of `tri`) with the triangle `(p, q, r)` (or the triangle `(i, j, k)` of `tri`).
The returned value is a [`Certificate`](@ref), which is one of:

- `Inside`: `(a, b)` is entirely inside `(p, q, r)`.
- `Single`: `(a, b)` has one endpoint inside `(p, q, r)`, and the other is outside.
- `Outside`: `(a, b)` is entirely outside `(p, q, r)`.
- `Touching`: `(a, b)` is on `(p, q, r)`'s boundary, but not in its interior.
- `Multiple`: `(a, b)` passes entirely through `(p, q, r)`. This includes the case where a point is on the boundary of `(p, q, r)`.
"""
triangle_line_segment_intersection
function triangle_line_segment_intersection(p, q, r, a, b)
    ## Break into cases:
    # Case I: (a, b) is an edge of (p, q, r)
    a_eq = any(==(a), (p, q, r))
    b_eq = any(==(b), (p, q, r))
    if a_eq && b_eq
        return Cert.Touching
    end

    # Case II: a or b is one of the vertices of (p, q, r)
    if a_eq || b_eq
        a, b = a_eq ? (a, b) : (b, a) # rotate so that a is one of the vertices 
        p, q, r = choose_uvw(a == p, a == q, a == r, p, q, r) # rotate so that a == p
        # First, find position of b relative to the edges pq and pr 
        pqb = point_position_relative_to_line(p, q, b)
        prb = point_position_relative_to_line(p, r, b)
        if is_right(pqb) || is_left(prb)
            # If b is right of pq or left of pr, it is away from the interior and thus outside
            return Cert.Outside
        elseif is_collinear(pqb) || is_collinear(prb)
            # If b is collinear with an edge, it could be on the edge or away from it
            s = is_collinear(pqb) ? q : r
            b_pos = point_position_on_line_segment(p, s, b)
            return is_left(b_pos) ? Cert.Outside : Cert.Touching
        else
            # If we are here, it means that b is between the rays through pq and pr. 
            # We just need to determine if b is inside or the triangle or outside. 
            # See that if b is inside, then it is left of qr, but if it is outside 
            # then it is right of qr. 
            qrb = point_position_relative_to_line(q, r, b)
            if is_right(qrb)
                return Cert.Multiple # Technically, one could call this Single, but Multiple is the most appropriate
            else
                return Cert.Inside # Also covers the is_on(qrb) case 
            end
        end
    end

    # We have now covered the cases where a or b, or both, are a vertex of (p, q, r).
    # The remaining degenerate case is the case where a or b, or both, are on an edge 
    # of (p, q, r). We can determine this possibility with the following certificates.
    a_tri = point_position_relative_to_triangle(p, q, r, a)
    b_tri = point_position_relative_to_triangle(p, q, r, b)

    # Case III: a and b are both on an edge of (p, q, r)
    if is_on(a_tri) && is_on(b_tri)
        # There are two options here, either a and b are both on the same edge,
        # or they are on different edges. In the latter case, (a, b) is entirely in 
        # (p, q, r). To determine this case, we just test against each edge.
        pqa = is_collinear(point_position_relative_to_line(p, q, a))
        qra = pqa ? false : is_collinear(point_position_relative_to_line(q, r, a)) # skip computation if we already know a result
        #rpa = (pqa || qra) ? false : is_collinear(point_position_relative_to_line(r, p, a))
        pqb = is_collinear(point_position_relative_to_line(p, q, b))
        qrb = pqb ? false : is_collinear(point_position_relative_to_line(q, r, b))
        #rpb = (pqb || qrb) ? false : is_collinear(point_position_relative_to_line(r, p, b))
        # Don't need rpa or rpb above, since checking the two edges automatically implies the result for the third edge
        a_edge = pqa ? (p, q) : (qra ? (q, r) : (r, p))
        b_edge = pqb ? (p, q) : (qrb ? (q, r) : (r, p))
        return a_edge == b_edge ? Cert.Touching : Cert.Inside
    end

    # Case IV: a or b are on an edge of (p, q, r)
    # For this case, the working is essentially the same as in Case II.
    if is_on(a_tri) || is_on(b_tri)
        # Rotate so that a is on (p, q)
        a, b = is_on(a_tri) ? (a, b) : (b, a)
        pqa = is_collinear(point_position_relative_to_line(p, q, a))
        qra = pqa ? false : is_collinear(point_position_relative_to_line(q, r, a))
        rpa = !(pqa || qra)
        p, q, r = choose_uvw(pqa, qra, rpa, p, q, r)
        # Similar to Case II, we start by considering the position of b relative to pq and qr 
        prb = point_position_relative_to_line(p, r, b)
        qrb = point_position_relative_to_line(q, r, b)
        # First, if we are left of right of pr or qr, respectively, we must have exited the triangle
        if is_left(prb) || is_right(qrb)
            return Cert.Multiple
        end
        # It is also possible that b is collinear with pq 
        pqb = point_position_relative_to_line(p, q, b)
        if is_collinear(pqb)
            #b_pos = point_position_on_line_segment(p, q, b)
            # We don't really need to know the position, we already know now that ab touches the triangle. 
            return Cert.Touching
            # Alternatively, b might be away from pq, meaning it is outside:
        elseif is_right(pqb)
            return Cert.Outside
        else
            # We now know that b is left of pqb, which from the above must imply that b is in the triangle. 
            return Cert.Inside
        end
    end

    # We can now assume that neither a nor b appear anywhere on (p, q, r). Just to make things a bit simpler, 
    # let us determine directly where a and b are relative to the triangle 
    a_pos = point_position_relative_to_triangle(p, q, r, a)
    b_pos = point_position_relative_to_triangle(p, q, r, b)

    # Case V: a and b are inside the triangle 
    if is_inside(a_pos) && is_inside(b_pos)
        return Cert.Inside
    elseif is_inside(a_pos) || is_inside(b_pos)
        # Case VI: Only one of a and b are inside the triangle 
        return Cert.Single
    else # is_outside(a_pos) && is_outside(b_pos)
        # Case VII: Both a and b are outside the triangle 
        # We consider the intersection of ab with each edge 
        pq_ab = line_segment_intersection_type(p, q, a, b)
        qr_ab = line_segment_intersection_type(q, r, a, b)
        rp_ab = line_segment_intersection_type(r, p, a, b)
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
            return Cert.Outside
        elseif num_multiple == 1 || num_touching == 2
            # This is for cases 3 or 4 above. We need to see how ab line up with the vertices of (p, q, r)
            abp = point_position_relative_to_line(a, b, p)
            abq = point_position_relative_to_line(a, b, q)
            abr = point_position_relative_to_line(a, b, r)
            if count(is_collinear, (abp, abq, abr)) == 2
                return Cert.Touching # lines up with an edge 
            else
                return Cert.Outside # interior only just touches a single vertex 
            end
        else
            return Cert.Multiple
        end
    end
end
function triangle_line_segment_intersection(tri::Triangulation, i, j, k, u, v)
    p, q, r, a, b = get_point(tri, i, j, k, u, v)
    return triangle_line_segment_intersection(p, q, r, a, b)
end

"""
    find_edge(tri::Triangulation, T, ℓ) -> Edge 

Given a triangle `T = (i, j, k)` of `tri` and a vertex `ℓ` of `tri`, returns the edge of `T` that contains `ℓ`. 
If no such edge exists, the edge `(k, i)` is returned.
"""
function find_edge(tri::Triangulation, T, ℓ)
    T′ = sort_triangle(T) # sort so that, if T is a ghost, then (i, j) is the solid edge
    i, j, k = triangle_vertices(T′)
    a, b, c = get_point(tri, i, j, k)
    E = edge_type(tri)
    if is_ghost_vertex(k) # must be on the solid edge
        return construct_edge(E, i, j)
    else
        (is_collinear ∘ point_position_relative_to_line)(tri, a, b, ℓ) && return construct_edge(E, i, j)
        (is_collinear ∘ point_position_relative_to_line)(tri, b, c, ℓ) && return construct_edge(E, j, k)
        return construct_edge(E, k, i)
    end
end

@doc """
    point_position_relative_to_circumcircle(tri::Triangulation, i, j, k, ℓ) -> Certificate
    point_position_relative_to_circumcircle(tri::Triangulation, T, ℓ) -> Certificate

Tests the position of the vertex `ℓ` of `tri` relative to the circumcircle of the triangle `T = (i, j, k)`. The returned value is a [`Certificate`](@ref), which is one of:

- `Outside`: `ℓ` is outside of the circumcircle of `T`.
- `On`: `ℓ` is on the circumcircle of `T`.
- `Inside`: `ℓ` is inside the circumcircle of `T`.

!!! note "Ghost triangles"

    The circumcircle of a ghost triangle is defined as the oriented outer halfplane of the solid edge of the triangle. See [`point_position_relative_to_oriented_outer_halfplane`](@ref).

!!! note "Weighted triangulations"

    If `tri` is a weighted triangulation, then an orientation test is instead applied, testing the orientation of the lifted companions of each point to determine if 
    `ℓ` is above or below the witness plane relative to `(i, j, k)`. For ghost triangles, the same rule applies, although if the vertex is on the solid edge of the ghost triangle then,
    in addition to checking [`point_position_relative_to_oriented_outer_halfplane`](@ref), we must check if the new vertex is not submerged by the adjoining solid triangle.
"""
point_position_relative_to_circumcircle
function point_position_relative_to_circumcircle(tri::Triangulation, i, j, k, ℓ)
    u, v, w = sort_triangle(i, j, k)
    a, b, c, p = get_point(tri, u, v, w, ℓ)
    if is_ghost_vertex(w)
        cert = point_position_relative_to_oriented_outer_halfplane(a, b, p)
        if is_on(cert) && is_weighted(tri)
            u′, v′, w′ = replace_ghost_triangle_with_boundary_triangle(tri, (u, v, w))
            sub_cert = point_position_relative_to_witness_plane(tri, u′, v′, w′, ℓ)
            is_above(sub_cert) && return Cert.Outside
            return cert
        else
            return cert
        end
    elseif !is_weighted(tri)
        return point_position_relative_to_circle(a, b, c, p)
    else
        cert = point_position_relative_to_witness_plane(tri, i, j, k, ℓ)
        return is_above(cert) ? Cert.Outside :
               is_below(cert) ? Cert.Inside :
               Cert.On
    end
end
point_position_relative_to_circumcircle(tri::Triangulation, T, ℓ) = point_position_relative_to_circumcircle(tri, geti(T), getj(T), getk(T), ℓ)

"""
    point_position_relative_to_witness_plane(tri::Triangulation, i, j, k, ℓ) -> Certificate

Given a positively oriented triangle `T = (i, j, k)` of `tri` and a vertex `ℓ` of `tri`, returns the position of `ℓ` relative to the witness plane of `T`. The returned value is a [`Certificate`](@ref), which is one of:

- `Above`: `ℓ` is above the witness plane of `T`.
- `On`: `ℓ` is on the witness plane of `T`.
- `Below`: `ℓ` is below the witness plane of `T`.

# Extended help 
The witness plane of a triangle is defined as the plane through the triangle `(i⁺, j⁺, k⁺)` where, for example, `pᵢ⁺ = (x, y, x^2 + y^2 - wᵢ)` is the lifted companion of `i` 
and `(x, y)` are the coordinates of the `i`th vertex. Moreover, by the orientation of `ℓ` relative to this witness plane we are referring to `ℓ⁺`'s position, not the plane point `ℓ`.
"""
function point_position_relative_to_witness_plane(tri::Triangulation, i, j, k, ℓ)
    p⁺ = get_lifted_point(tri, i)
    q⁺ = get_lifted_point(tri, j)
    r⁺ = get_lifted_point(tri, k)
    a⁺ = get_lifted_point(tri, ℓ)
    cert = orient_predicate(p⁺, q⁺, r⁺, a⁺)
    return convert_certificate(cert, Cert.Above, Cert.On, Cert.Below)
end

ExactPredicates.Codegen.@genpredicate function angle_is_acute(p::2, q::2, r::2)
    pr = p - r
    qr = q - r
    ExactPredicates.Codegen.group!(pr...)
    ExactPredicates.Codegen.group!(qr...)
    return pr[1] * qr[1] + pr[2] * qr[2]
end

"""
    opposite_angle(p, q, r) -> Certificate

Tests the angle opposite to the edge `(p, q)` in the triangle `(p, q, r)`, meaning `∠prq`. The returned value is a 
[`Certificate`](@ref), which is one of:

- `Obtuse`: The angle opposite to `(p, q)` is obtuse.
- `Right`: The angle opposite to `(p, q)` is a right angle.
- `Acute`: The angle opposite to `(p, q)` is acute.

!!! note "Angle between two vectors"

    This function computes `(p - r) ⋅ (q - r)`. If you want the angle between vectors `a = pq` and `b = pr`,
    then you should use `opposite_angle(r, q, p) = (r - p) ⋅ (q - p)`.
"""
function opposite_angle(p, q, r) # https://math.stackexchange.com/a/701682/861404 
    p, q, r = _getxy(p), _getxy(q), _getxy(r)
    cert = angle_is_acute(p, q, r)
    return convert_certificate(cert, Certificate.Obtuse, Certificate.Right, Certificate.Acute)
end

"""
    point_position_relative_to_diametral_circle(p, q, r) -> Certificate

Given an edge `(p, q)` and a point `r`, returns the position of `r` relative to the diametral circle of `(p, q)`. The returned value 
is a [`Certificate`](@ref), which is one of:

- `Inside`: `r` is inside the diametral circle.
- `On`: `r` is on the diametral circle.
- `Outside`: `r` is outside the diametral circle.

# Extended help
The check is done by noting that a point is in the diametral circle if, and only if, the angle at `r` is obtuse.
"""
function point_position_relative_to_diametral_circle(p, q, r)
    d = opposite_angle(p, q, r)
    if is_obtuse(d)
        return Certificate.Inside
    elseif is_right(d)
        return Certificate.On
    else
        return Certificate.Outside
    end
end

"""
    point_position_relative_to_diametral_lens(p, q, r, lens_angle=30.0) -> Certificate 

Given an edge `(p, q)` and a point `r`, returns the position of `r` relative to the diametral lens of `(p, q)` with lens angle `lens_angle` (in degrees). The returned value 
is a [`Certificate`](@ref), which is one of:
    
- `Inside`: `r` is inside the diametral lens.
- `On`: `r` is on the diametral lens.
- `Outside`: `r` is outside the diametral lens.

!!! warning 

    This function assumes that the lens angle is at most 45°.
    
# Extended help 
The diametral lens is slightly similar to a diametral circle. Let us first define the standard definition of a diametral lens, where `lens_angle = 30°`, and then we define its 
generalisation. The standard definition was introduced by Shewchuk in 2002 in [this paper](https://doi.org/10.1016/S0925-7721(01)00047-5), and the generalisation was originally 
described in Section 3.4 of Shewchuk's 1997 PhD thesis [here](https://www.cs.cmu.edu/~quake-papers/delaunay-refinement.pdf) (as was the standard definition, of course).

Standard definition: Let the segment be `s ≔ (p, q)` and let `ℓ` be the perpendicular bisector of `s`. Define two circles `C⁺` and `C⁻` whose centers lie left and right 
of `s`, respectively, both the same distance away from the midpoint `m = (p + q)/2`. By placing each disk a distance `|s|/(2√3)` away from `m` and extending the disks so that 
they both touch `p` and `q`, meaning they each have radius `|s|/√3`, we obtain disks that touch `p` and `q` as well as the center of the other disk. The intersection of the two 
disks defines the diametral lens. The lens angle associated with this lens is `30°`, as we could draw lines between `p` and `q` and the poles of the disks to form isosceles triangles 
on each side, giving angles of `30°` at `p` and `q`.

Generalised definition: The generalisation of the above definition aims to generalise the lens angle to some angle `θ`. In particular, let us draw a line from `p` towards some point 
`u` which is left of `s` and on the perpendicular bisector, where the line is at an angle `θ`. This defines a triplet of points `(p, q, u)` from which we can define a circle, and similarly for a point `v` right of `s`.
The intersection of these two circles is the diametral lens with angle `θ`. We can also figure out how far `u` is along the perpendicular bisector, i.e. how far away it is from `(p + q)/2`,
using basic trigonometry. Let `m = (p +  q)/2`, and consider the triangle `△pmu`. The side lengths of this triangle are `|pm| = |s|/2`, `|mu| ≔ b`, and `|up| ≔ c`. Let us first compute `c`. 
We have `cos(θ) = |pm|/|up| = |s|/(2c)`, and so `c = |s|/(2cos(θ))`. So, `sin(θ) = |mu|/|up| = b/c`, which gives `b = csin(θ) = |s|sin(θ)/(2cos(θ)) = |s|tan(θ)/2`. So, `u` is a distance `|s|tan(θ)/2`
from the midpoint. Notice that in the case `θ = 30°`, `tan(θ) = √3/3`, giving `b = |s|√3/6 = |s|/(2√3)`, which is the same as the standard definition.

To test if a point `r` is inside the diametral lens with lens angle `θ°`, we simply have to check the angle 
at that point, i.e. the angle at `r` in the triangle `pqr`. If this angle is greater than `180° - 2θ°`, then `r` is inside of the lens. This result 
comes from [Shewchuk (2002)](https://doi.org/10.1016/S0925-7721(01)00047-5). Note that this is the same as a diametral circle in the case `θ° = 45°`. 
"""
function point_position_relative_to_diametral_lens(p, q, r, lens_angle=30.0)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    rx, ry = _getxy(r)
    ddot = (px - rx) * (qx - rx) + (py - ry) * (qy - ry) # (p - r) · (q - r)
    ddot > 0 && return Certificate.Outside # no need to check - not obtuse
    ddot² = ddot^2
    ℓpr² = dist_sqr((px, py), (rx, ry))
    ℓqr² = dist_sqr((qx, qy), (rx, ry))
    cosine_scale² = cosd(2lens_angle)^2 # ≥ 0 since lens_angle ∈ [0°, 45°]
    rhs = (ℓpr² * ℓqr²) * cosine_scale²
    if ddot² > rhs
        return Certificate.Inside
    elseif ddot² == rhs
        return Certificate.On
    else
        return Certificate.Outside
    end
end

"""
    test_visibility(tri::Triangulation, q, i) -> Certificate 

Tests if the vertex `i` and the point `q` can see each other. Here, visibility means that the line segment joining 
the two does not intersect any segments.

# Arguments 
- `tri`: The [`Triangulation`](@ref).
- `1`: The point from which we are testing visibility.
- `i`: The vertex we are testing visibility of.

# Outputs 
- `cert`: A [`Certificate`](@ref). This will be `Visible` if `i` is visible from `q`, and `Invisible` otherwise.
"""
function test_visibility(tri::Triangulation, q, i)
    V, invisible_flag = jump_and_march(tri, q; use_barriers=Val(true), k=i, concavity_protection=true)
    if invisible_flag
        return Certificate.Invisible
    else
        return Certificate.Visible
    end
end

"""
    test_visibility(tri::Triangulation, u, v, i; shift=0.0, attractor=get_point(tri,i)) -> Certificate 

Tests if the edge `(u, v)` and the point `i` can see each other. Here, visibility means that any point in the interior 
of `(u, v)` can see `i`. To test this, we only check `10` points equally spaced between `u` and `v`, excluding `u` and `v`.

# Arguments 
- `tri`: The [`Triangulation`](@ref).
- `u`: The first vertex of the edge.
- `v`: The second vertex of the edge.
- `i`: The vertex we are testing visibility of.

# Keyword Arguments 
- `shift=0.0`: The amount by which to shift each point on the edge towards `attractor`, i.e. if `p` is a point on the edge, then `p .+ shift .* (attractor - p)` is the point used to test visibility rather than `p` itself.

# Outputs
- `cert`: A [`Certificate`](@ref). This will be `Visible` if `i` is visible from `(u, v)`, and `Invisible` otherwise.
"""
function test_visibility(tri::Triangulation, u, v, i; shift=0.0, attractor=get_point(tri,i))
    pu, pv = get_point(tri, u, v)
    pux, puy = getxy(pu)
    pvx, pvy = getxy(pv)
    qx, qy = getxy(attractor)
    ts = LinRange(0.00001, 0.99999, 10)
    for t in ts
        mx, my = pux + t * (pvx - pux), puy + t * (pvy - puy)
        m̃ = (mx + shift * (qx - mx), my + shift * (qy - my))
        cert = test_visibility(tri, m̃, i)
        is_visible(cert) && return Certificate.Visible
    end
    return Cert.Invisible
end