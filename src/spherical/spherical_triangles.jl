"""
    spherical_circumcenter(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint) -> SphericalPoint

Computes the circumcenter of the spherical triangle defined by the points `p`, `q`, and `r`. 

See https://brsr.github.io/2021/05/02/spherical-triangle-centers.html
"""
function spherical_circumcenter(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint)
    num = cross(p, q) .+ cross(q, r) .+ cross(r, p)
    den = _hypot(num...)
    return SphericalPoint(num ./ den)
end

"""
    spherical_triangle_area(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint) -> Number

Computes the area of the spherical triangle defined by the points `p`, `q`, and `r`.

See https://brsr.github.io/2023/11/04/spherical-areal.html.
"""
function spherical_triangle_area(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint)
    num = abs(dot(p, SphericalPoint(cross(q, r))))
    den = 1 + dot(p, q) + dot(q, r) + dot(r, p)
    rat = num / den
    A = 2atan(num, den)
    return A
end

"""
    spherical_triangle_circumradius(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint) -> Number

Computes the circumradius of the spherical triangle defined by the points `p`, `q`, and `r`.

See https://brsr.github.io/2021/05/02/spherical-triangle-centers.html.
"""
function spherical_triangle_circumradius(p, q, r)
    num = scalar_triple_product(p, q, r)
    _den = cross(p, q) .+ cross(q, r) .+ cross(r, p)
    den = _hypot(_den...)
    return acos(clamp(num / den, -1, 1))
end

"""
    spherical_distance(p::SphericalPoint, q::SphericalPoint) -> Number

Computes the spherical distance between the points `p` and `q`.
"""
function spherical_distance(p::SphericalPoint, q::SphericalPoint)
    # https://brsr.github.io/2021/05/01/vector-spherical-geometry.html
    _px, _py, _pz = __getxyz(p)
    _qx, _qy, _qz = __getxyz(q)
    dx, dy, dz = _px - _qx, _py - _qy, _pz - _qz
    ex, ey, ez = _px + _qx, _py + _qy, _pz + _qz
    num = _hypot(dx, dy, dz)
    den = _hypot(ex, ey, ez)
    return 2atan(num / den)
end

"""
    spherical_angles(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint) -> NTuple{3, Number}

Computes the spherical angles of the spherical triangle defined by the points `p`, `q`, and `r`, returned 
as a tuple `(α, β, γ)` for the interior angles at `p`, `q`, and `r`, respectively.

See https://brsr.github.io/2021/05/01/vector-spherical-geometry.html.
"""
function spherical_angles(p, q, r)
    return (
        _spherical_angle(p, q, r),
        _spherical_angle(q, r, p),
        _spherical_angle(r, p, q)
    )
end

function _spherical_angle(p, q, r)
    # Computes the interior angle at p
    y = scalar_triple_product(p, q, r)
    dot_qr = dot(q, r)
    dot_pq = dot(p, q)
    dot_rp = dot(r, p)
    x = dot_qr - dot_pq * dot_rp
    return atan(y, x)
end

"""
    spherical_triangle_centroid(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint) -> SphericalPoint

Computes the centroid of the spherical triangle defined by the points `p`, `q`, and `r`.
This is the point that separates the triangle into three equal-area regions.

See https://brsr.github.io/2023/11/04/spherical-areal.html.
"""
function spherical_triangle_centroid(p, q, r)
    A, B, C = spherical_angles(p, q, r)
    S = (A + B + C) / 2
    Ω = (2S - π) / 3
    x = csc(A - Ω)
    y = csc(B - Ω)
    z = csc(C - Ω)
    _px, _py, _pz = __getxyz(p)
    _qx, _qy, _qz = __getxyz(q)
    _rx, _ry, _rz = __getxyz(r)
    sa = _hypot(cross(q, r)...)
    sb = _hypot(cross(r, p)...)
    sc = _hypot(cross(p, q)...)
    Px = _px * sa * x + _qx * sb * y + _rx * sc * z
    Py = _py * sa * x + _qy * sb * y + _ry * sc * z
    Pz = _pz * sa * x + _qz * sb * y + _rz * sc * z
    scale = _hypot(Px, Py, Pz)
    return SphericalPoint(Px / scale, Py / scale, Pz / scale)
end

"""
    point_position_relative_to_spherical_triangle(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint, s::SphericalPoint) -> Certificate 

Tests if `s` is in the spherical triangle defined by the points `p`, `q`, and `r`. Returns a `Certificate` object, which is

- `Certificate.Inside`: if `s` is in the triangle,
- `Certificate.On`: if `s` is on the edge of the triangle, or
- `Certificate.Outside`: if `s` is outside the triangle.
"""
function point_position_relative_to_spherical_triangle(p::SphericalPoint, q::SphericalPoint, r::SphericalPoint, s::SphericalPoint)
    x1, y1, z1 = __getxyz(p)
    x2, y2, z2 = __getxyz(q)
    x3, y3, z3 = __getxyz(r)
    x, y, z = __getxyz(s)
    D = x * y1 * z2 - x * y2 * z1 - x1 * y * z2 + x1 * y2 * z + x2 * y * z1 - x2 * y1 * z - x * y1 * z3 + x * y3 * z1 + x1 * y * z3 - x1 * y3 * z - x3 * y * z1 + x3 * y1 * z + x * y2 * z3 - x * y3 * z2 - x2 * y * z3 + x2 * y3 * z + x3 * y * z2 - x3 * y2 * z
    a = (x * y2 * z3 - x * y3 * z2 - x2 * y * z3 + x2 * y3 * z + x3 * y * z2 - x3 * y2 * z) / D
    b = -(x * y1 * z3 - x * y3 * z1 - x1 * y * z3 + x1 * y3 * z + x3 * y * z1 - x3 * y1 * z) / D
    c = (D + x * y1 * z3 - x * y3 * z1 - x1 * y * z3 + x1 * y3 * z + x3 * y * z1 - x3 * y1 * z - x * y2 * z3 + x * y3 * z2 + x2 * y * z3 - x2 * y3 * z - x3 * y * z2 + x3 * y2 * z) / D
    lambda = (x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1) / D
    if a > 0 && b > 0 && c > 0 && lambda > 0
        return Certificate.Inside
    elseif a == 0 || b == 0 || c == 0
        return Certificate.On
    else
        return Certificate.Outside
    end
end