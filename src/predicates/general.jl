"""
    orient_predicate(p, q, r)

Returns `ExactPredicates.orient(p, q, r)`.
"""
orient_predicate(p, q, r) = orient(getxy(p), getxy(q), getxy(r))

"""
    incircle_predicate(a, b, c, p)

Returns `ExactPredicates.incircle(a, b, c, p)`.
"""
incircle_predicate(a, b, c, p) = incircle(getxy(a), getxy(b), getxy(c), getxy(p))

"""
    sameside_predicate(a, b, p)

Returns `ExactPredicates.sameside(p, a, b)` (but we redefine it here).

(The difference in the argument order is to match the convention that the 
main point being tested is the last argument.)
"""
function sameside_predicate(a, b, p)
    _p = getxy(p)
    _a = getxy(a)
    _b = getxy(b)
    if _a < _p && _b < _p || _a > _p && _b > _p
        return 1
    elseif _a < _p && _b > _p || _a > _p && _b < _p
        return -1
    else
        return 0
    end
end

opposite_signs(x, y) = xor(x, y) == -2 # also from ExactPredicates.jl 

"""
    meet_predicate(p, q, a, b) 

Returns `ExactPredicates.meet(p, q, a, b)`  (but we redefine it here).
"""
function meet_predicate(p, q, a, b)
    pqa = orient_predicate(p, q, a)
    pqb = orient_predicate(p, q, b)
    abp = orient_predicate(a, b, p)
    abq = orient_predicate(a, b, q)
    if opposite_signs(pqa, pqb) && opposite_signs(abp, abq)
        return 1
    elseif (pqa ≠ 0 && pqa == pqb) || (abq ≠ 0 && abp == abq)
        return -1
    elseif pqa == 0 && pqb == 0
        if sameside_predicate(a, b, p) == 1 &&
           sameside_predicate(a, b, q) == 1 &&
           sameside_predicate(p, q, a) == 1 &&
           sameside_predicate(p, q, b) == 1
            return -1
        else
            return 0
        end
    else
        return 0
    end
end

"""
    triangle_orientation(p, q, r)

Given a triangle with coordinates `(p, q, r)`, computes its orientation, returning:

- `Certificate.PositivelyOriented`: The triangle is positively oriented.
- `Certificate.Degenerate`: The triangle is degenerate, meaning the coordinates are collinear. 
- `Certificate.NegativelyOriented`: The triangle is negatively oriented.

See also [`orient_predicate`](@ref).
"""
@inline function triangle_orientation(p, q, r)
    cert = orient_predicate(p, q, r)
    return convert_certificate(cert, Cert.NegativelyOriented, Cert.Degenerate,
        Cert.PositivelyOriented)
end

"""
    point_position_relative_to_circle(a, b, c, p)

Given a circle through the coordinates `(a, b, c)`, assumed to be positively oriented, 
computes the position of `p` relative to the circle. In particular, returns:

- `Certificate.Inside`: `p` is inside the circle.
- `Certificate.On`: `p` is on the circle. 
- `Certificate.Outside`: `p` is outside the triangle.

See also [`incircle_predicate`](@ref).
"""
function point_position_relative_to_circle(a, b, c, p)
    cert = incircle_predicate(a, b, c, p)
    return convert_certificate(cert, Cert.Outside, Cert.On, Cert.Inside)
end

"""
    point_position_relative_to_line(a, b, p)

Given a point `p` and the oriented line `(a, b)`, computes the position 
of `p` relative to the line, returning:

- `Certificate.Left`: `p` is to the left of the line. 
- `Certificate.Collinear`: `p` is on the line.
- `Certificate.Right`: `p` is to the right of the line. 

See also [`orient_predicate`](@ref).
"""
function point_position_relative_to_line(a, b, p)
    cert = orient_predicate(a, b, p)
    return convert_certificate(cert, Cert.Right, Cert.Collinear, Cert.Left)
end

"""
    point_position_on_line_segment(a, b, p)

Given a point `p` and the line segment `(a, b)`, assuming `p` 
to be collinear with `a` and `b`, computes the position of `p`
relative to the line segment. In particular, returns:

- `Certificate.On`: `p` is on the line segment, meaning between `a` and `b`.
- `Certificate.Degenerate`: Either `p == a` or `p == b`, i.e. `p` is one of the endpoints. 
- `Certificate.Left`: `p` is off and to the left of the line segment.
- `Certificate.Right`: `p` is off and to the right of the line segment.

See also [`sameside_predicate`](@ref).
"""
function point_position_on_line_segment(a, b, p)
    cert = sameside_predicate(a, b, p)
    converted_cert = convert_certificate(cert, Cert.On, Cert.Degenerate, Cert.Outside)
    if is_outside(converted_cert) # Do we have "a ---- b ---- p" or "p ---- a ---- b"?
        ap_cert = sameside_predicate(a, p, b)
        converted_ap_cert = convert_certificate(ap_cert, Cert.On, Cert.Degenerate,
            Cert.Outside)
        return is_on(converted_ap_cert) ? Cert.Right : Cert.Left
    end
    return converted_cert
end

"""
    line_segment_intersection_type(p, q, a, b)

Given the coordinates `(p, q)` and `(a, b)` defining two line segments, 
tests the number of intersections between the two segments. In particular, 
we return:

- `Certificate.None`: The line segments do not meet at any points. 
- `Certificate.Multiple`: The closed line segments `[p, q]` and `[a, b]` meet in one or several points. 
- `Certificate.Single`: The open line segments `(p, q)` and `(a, b)` meet in a single point. 
- `Certificate.Touching`: One of the endpoints is on `[a, b]`, but there are no other intersections.

See also [`meet_predicate`](@ref).
"""
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
        if (is_collinear(collinear_cert_1) && !is_collinear(collinear_cert_2)) ||
           (!is_collinear(collinear_cert_1) && is_collinear(collinear_cert_2))
            return Cert.Touching
        end
        collinear_cert_3 = point_position_relative_to_line(p, q, a)
        collinear_cert_4 = point_position_relative_to_line(p, q, b)
        if (is_collinear(collinear_cert_3) && !is_collinear(collinear_cert_4)) ||
           (!is_collinear(collinear_cert_3) && is_collinear(collinear_cert_4))
            return Cert.Touching
        end
    end
    return converted_cert
end

"""
    point_position_relative_to_triangle(a, b, c, p)

Given a positively oriented triangle with coordinates `(a, b, c)`, computes the 
position of `p` relative to the triangle. In particular, returns: 

- `Certificate.Outside`: `p` is outside of the triangle. 
- `Certificate.On`: `p` is on one of the edges. 
- `Certificate.Inside`: `p` is inside the triangle.
"""
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

"""
    point_position_relative_to_oriented_outer_halfplane(a, b, p)

Given an edge with coordinates `(a, b)` and a point `p`, 
tests the position of `p` relative to the oriented outer halfplane defined 
by `(a, b)`. The oriented outer halfplane is the union of the open halfplane 
defined by the region to the left of the oriented line `(a, b)`, and the 
open line segment `(a, b)`. The returned values are:

- `Cert.Outside`: `p` is outside of the oriented outer halfplane, meaning to the right of the line `(a, b)` or collinear with `a` and `b` but not on the line segment `(a, b)`.
- `Cert.On`: `p` is on the open line segment `(a, b)`.
- `Cert.Inside`: `p` is inside of the oriented outer halfplane, meaning to the left of the line `(a, b)`.
"""
function point_position_relative_to_oriented_outer_halfplane(a, b, p)
    in_open_halfplane = point_position_relative_to_line(a, b, p)
    if is_collinear(in_open_halfplane)
        is_on_boundary_edge = point_position_on_line_segment(a, b, p)
        if is_on(is_on_boundary_edge)
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

"""
    is_legal(p, q, r, s)

Given an edge `pq`, incident to two triangles `pqr` and `qps`, tests 
if the edge `pq` is legal, i.e. if `s` is not inside the triangle through 
`p`, `q`, and `r`.
"""
function is_legal(p, q, r, s)
    incirc = point_position_relative_to_circle(p, q, r, s)
    if is_inside(incirc)
        return Cert.Illegal 
    else 
        return Cert.Legal 
    end
end