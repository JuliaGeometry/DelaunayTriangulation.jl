# Some formulas in this file come from https://perso.uclouvain.be/jean-francois.remacle/LMECA2170/robnotes.pdf.
# Unfortunately, I can't use the robust forms of the formulas because we don't have implementations for them (only the versions that take the sign of the orient predicate).

"""
    IndividualTriangleStatistics{T}

Struct storing statistics of a single triangle.

# Fields 
- `area::T`: The area of the triangle.
- `lengths::NTuple{3,T}`: The lengths of the edges of the triangle, given in sorted order. 
- `circumcenter::NTuple{2,T}`: The circumcenter of the triangle.
- `circumradius::T`: The circumradius of the triangle.
- `angles::NTuple{3, T}`: The angles of the triangle, given in sorted order.
- `radius_edge_ratio::T`: The ratio of the circumradius to the shortest edge length.
- `edge_midpoints::NTuple{3,NTuple{2,T}}`: The midpoints of the edges of the triangle.
- `aspect_ratio::T`: The ratio of the inradius to the circumradius.
- `inradius::T`: The inradius of the triangle.
- `perimeter::T`: The perimeter of the triangle.
- `centroid::NTuple{2,T}`: The centroid of the triangle.
- `offcenter::NTuple{2,T}`: The offcenter of the triangle with radius-edge ratio cutoff `β=1`. See [this paper](https://doi.org/10.1016/j.comgeo.2008.06.002).
- `sink::NTuple{2,T}`: The sink of the triangle relative to the parent triangulation. See [this paper](https://doi.org/10.1145/378583.378644).

# Constructors 
The constructor is 

    IndividualTriangleStatistics(p, q, r, sink = (NaN, NaN))

where `p`, `q`, and `r` are the coordinates of the triangle given in 
counter-clockwise order. `sink` is the triangle's sink. This must be provided 
separately since it is only computed relative to a triangulation, and so requires 
vertices rather than coordinates; see [`triangle_sink`](@ref).

# Extended help 
The relevant functions used for computing these statistics are 

- [`squared_triangle_lengths`](@ref)
- [`triangle_area`](@ref)
- [`triangle_circumcenter`](@ref)
- [`triangle_circumradius`](@ref)
- [`triangle_radius_edge_ratio`](@ref)
- [`triangle_edge_midpoints`](@ref)
- [`triangle_perimeter`](@ref)
- [`triangle_inradius`](@ref)
- [`triangle_aspect_ratio`](@ref)
- [`triangle_centroid`](@ref)
- [`triangle_angles`](@ref)
- [`triangle_offcenter`](@ref)
- [`triangle_sink`](@ref).
"""
struct IndividualTriangleStatistics{T}
    area::T
    lengths::NTuple{3, T}
    circumcenter::NTuple{2, T}
    circumradius::T
    angles::NTuple{3, T}
    radius_edge_ratio::T
    edge_midpoints::NTuple{3, NTuple{2, T}}
    aspect_ratio::T
    inradius::T
    perimeter::T
    centroid::NTuple{2, T}
    offcenter::NTuple{2, T}
    sink::NTuple{2, T}
end
function IndividualTriangleStatistics(p, q, r, sink = (NaN, NaN))
    F = number_type(p)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin, ℓmed, ℓmax = sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    circumcenter = triangle_circumcenter(p, q, r, A)
    circumradius = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    radius_edge_ratio = triangle_radius_edge_ratio(circumradius, ℓmin)
    edge_midpoints = triangle_edge_midpoints(p, q, r)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    inradius = triangle_inradius(A, perimeter)
    aspect_ratio = triangle_aspect_ratio(inradius, circumradius)
    centroid = triangle_centroid(p, q, r)
    angles = triangle_angles(p, q, r)
    offcenter = triangle_offcenter(p, q, r, circumcenter)
    sx, sy = getxy(sink)
    return IndividualTriangleStatistics(
        F(A),
        F.((ℓmin, ℓmed, ℓmax)),
        F.(circumcenter),
        F(circumradius),
        F.(angles),
        F(radius_edge_ratio),
        (F.(edge_midpoints[1]), F.(edge_midpoints[2]), F.(edge_midpoints[3])),
        F(aspect_ratio),
        F(inradius),
        F(perimeter),
        F.(centroid),
        F.(offcenter),
        (F(sx), F(sy)),
    )
end

## We could simplify some of this even further, e.g. keeping some more things as their square temporarily.
## Is it really worth the effort, though?

"""
    triangle_area(ℓ₁², ℓ₂², ℓ₃²) -> Number

Compute the area of a triangle given the squares of its edge lengths. If there are precision issues that cause the area to be negative, then the area is set to zero.

See also [`squared_triangle_area`](@ref).
"""
function triangle_area(ℓ₁²::Number, ℓ₂²::Number, ℓ₃²::Number)
    A² = squared_triangle_area(ℓ₁², ℓ₂², ℓ₃²)
    A = A² < 0.0 ? zero(A²) : sqrt(A²) # needed e.g. if A² is like -1.084020e-19
    return A
end

@doc raw"""
    squared_triangle_area(ℓ₁², ℓ₂², ℓ₃²) -> Number

Compute the squared area of a triangle given the squares of its edge lengths. Heron's formula is used, so that the squared area is 

```math 
A^2 = \dfrac{1}{16}\left[4\ell_1^2\ell_2^2 - \left(\ell_1^2 + \ell_2^2 - \ell_3^2\right)^2\right]..
```

See also [`squared_triangle_area_v2`](@ref).
"""
squared_triangle_area(ℓ₁²::Number, ℓ₂²::Number, ℓ₃²::Number) = (4ℓ₁² * ℓ₂² - (ℓ₁² + ℓ₂² - ℓ₃²)^2) / 16 # Heron's formula

@doc raw"""
    squared_triangle_area_v2(ℓ₁², ℓ₂², ℓ₃²) -> Number

Compute the squared area of a triangle given the squares of its edge lengths, given in sorted order so that `ℓ₁² ≤ ℓ₂² ≤ ℓ₃²`. This is a more numerically stable version of [`squared_triangle_area`](@ref) using 
the formula from [Kahan (2014)](https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf):

```math 
A^2 = \dfrac{1}{16}\left\{\left[\ell_3 + \left(\ell_2 + \ell_1\right)\right]\left[\ell_1 - \left(\ell_3 - \ell_2\right)\right]\left[\ell_1 + \left(\ell_3 - \ell_2\right)\right]\left[\ell_3 + \left(\ell_2 - \ell_1\right)\right]\right\}.
```
"""
squared_triangle_area_v2(ℓ₁²::Number, ℓ₂²::Number, ℓ₃²::Number) =
    let a = sqrt(ℓ₃²), b = sqrt(ℓ₂²), c = sqrt(ℓ₁²)
    return (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)) / 16 # https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf
end

@doc raw"""
    triangle_circumradius(A, ℓmin², ℓmed², ℓmax²) -> Number

Computes the circumradius of a triangle with area `A` and squared edge lengths `ℓmin² ≤ ℓmed² ≤ ℓmax²`. The circumradius is given by 

```math
r = \dfrac{\ell_{\min}\ell_{\text{med}}\ell_{\max}}{4A}.
```
"""
triangle_circumradius(A, ℓmin², ℓmed², ℓmax²) = sqrt(ℓmin² * ℓmed² * ℓmax²) / (4A)

@doc raw"""
    triangle_perimeter(ℓmin::Number, ℓmed::Number, ℓmax::Number) -> Number 

Computes the perimeter of a triangle with edge lengths `ℓmin ≤ ℓmed ≤ ℓmax`. The perimeter is given by 

```math
P = \ell_{\min} + \ell_{\text{med}} + \ell_{\max}.
```
"""
triangle_perimeter(ℓmin::Number, ℓmed::Number, ℓmax::Number) = ℓmin + ℓmed + ℓmax

@doc raw"""
    triangle_inradius(A, perimeter) -> Number

Computes the inradius of a triangle with area `A` and perimeter `perimeter`. The inradius is given by

```math
r_i = \dfrac{2A}{P},
```
where ``P`` is the `perimeter`.
"""
triangle_inradius(A, perimeter) = 2A / perimeter

@doc raw"""
    triangle_aspect_ratio(inradius::Number, circumradius::Number) -> Number 

Computes the aspect ratio of a triangle with inradius `inradius` and circumradius `circumradius`. The aspect ratio is given by

```math
\tau = \dfrac{r_i}{r},
```
where ``r_i`` is the inradius and ``r`` is the circumradius.
"""
triangle_aspect_ratio(inradius::Number, circumradius::Number) = inradius / circumradius

@doc raw"""
    triangle_radius_edge_ratio(circumradius::Number, ℓmin::Number) -> Number

Computes the radius-edge ratio of a triangle with circumradius `circumradius` and minimum edge length `ℓmin`, given by 

```math 
\rho = \dfrac{r}{\ell_{\min}},
```
where ``r`` is the circumradius and ``\ell_{\min}`` is the shortest edge length.
"""
triangle_radius_edge_ratio(circumradius::Number, ℓmin::Number) = circumradius / ℓmin

@doc raw"""
    triangle_centroid(p, q, r) -> (Number, Number)

Computes the centroid of a triangle with vertices `p`, `q`, and `r`, given by

```math
c = \dfrac{p + q + r}{3}.
```
"""
triangle_centroid(p, q, r) = ((getx(p) + getx(q) + getx(r)) / 3, (gety(p) + gety(q) + gety(r)) / 3)

@doc raw"""
    triangle_angles(p, q, r) -> (Number, Number, Number)

Computes the angles of a triangle with vertices `p`, `q`, and `r`. The formula for, say, the angle at `p` is given by 

```math
\theta_1 = \arctan\left(\dfrac{2A}{\left(p - q\right)\cdot\left(p - r\right)}\right),
```
where `A` is the area of the triangle. The angles are returned in sorted order.
"""
function triangle_angles(p, q, r)
    ℓ₁², ℓ₂², ℓ₃² = squared_triangle_lengths(p, q, r)
    A = triangle_area(ℓ₁², ℓ₂², ℓ₃²)
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

"""
    squared_triangle_area(p, q, r) -> Number

Computes the squared area of the triangle with coordinates `p`, `q`, `r`. Initially, `squared_triangle_area` is used for this, unless the squared area 
is found to be negative due to precision issues, in which case [`squared_triangle_area_v2`](@ref) is used instead.

!!! note "Precision"

    All coordinates are converted into Float64, but the returned area is converted back into the original precision.
"""
function squared_triangle_area(_p, _q, _r)
    p, q, r = _getxy(_p), _getxy(_q), _getxy(_r)
    ℓ₁², ℓ₂², ℓ₃² = squared_triangle_lengths(p, q, r)
    A² = squared_triangle_area(ℓ₁², ℓ₂², ℓ₃²)
    if A² ≤ zero(A²)
        A² = squared_triangle_area_v2(ℓ₁², ℓ₂², ℓ₃²)
    end
    if A² ≤ zero(A²)
        A² = zero(A²)
    end
    return number_type(_p)(A²)
end

"""
    triangle_area(p, q, r) -> Number

Computes the area of the triangle with coordinates `p`, `q`, `r`.

!!! note "Precision"

    All coordinates are converted into Float64, but the returned area is converted back into the original precision.
"""
function triangle_area(p, q, r)
    A² = squared_triangle_area(p, q, r)
    return A² < zero(A²) ? zero(A²) : sqrt(A²)
end

"""
    squared_triangle_lengths(p, q, r) -> (Number, Number, Number)

Computes the squared lengths of the edges of the triangle with coordinates `p`, `q`, `r`. The squared lengths are returned in sorted order.
"""
function squared_triangle_lengths(p, q, r)
    ℓ₁², ℓ₂², ℓ₃², _ = squared_triangle_lengths_and_smallest_index(p, q, r)
    return ℓ₁², ℓ₂², ℓ₃²
end

"""
    squared_triangle_lengths_and_smallest_index(p, q, r) -> (Number, Number, Number, Integer)

Computes the squared lengths of the edges of the triangle with coordinates `p`, `q`, `r`. The squared lengths are returned in sorted order, and the index of the shortest edge is returned as well. Here, 
the index refers to which edge in the order `(p, q)`, `(q, r)`, `(q, p)`.
"""
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

"""
    triangle_lengths(p, q, r) -> (Number, Number, Number)

Computes the lengths of the edges of the triangle with coordinates `p`, `q`, `r`. The lengths are returned in sorted order.
"""
function triangle_lengths(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    return sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
end

@doc raw"""
    triangle_circumcenter(p, q, r, A=triangle_area(p, q, r)) -> (Number, Number)

Computes the circumcenter of the triangle with coordinates `(p, q, r)`. The circumcenter is given by 

```math 
c_x = r_x + \dfrac{d_{11}d_{22} - d_{12}d_{21}}{4A}, \quad c_y = r_y + \dfrac{e_{11}e_{22} - e_{12}e_{21}}{4A},
```
where ``d_{11} = \|p - r\|_2^2``, ``d_{12} = p_y - r_y``, ``d_{21} = \|q - r\|_2^2``, ``d_{22} = q_y - r_y``, ``e_{11} = p_x - r_x``
``e_{12} = d_{11}``, ``e_{21} = q_x - r_x``, and ``e_{22} = d_{21}``.

!!! note "Precision"

    All coordinates are converted into Float64, but the returned area is converted back into the original precision.
"""
function triangle_circumcenter(_p, _q, _r, _A = triangle_area(_p, _q, _r))
    p, q, r = _getxy(_p), _getxy(_q), _getxy(_r)
    A = Float64(_A)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    d11 = (px - rx)^2 + (py - ry)^2
    d12 = py - ry
    d21 = (qx - rx)^2 + (qy - ry)^2
    d22 = qy - ry
    ox = rx + (d11 * d22 - d12 * d21) / (4A)
    e11 = px - rx
    e12 = d11
    e21 = qx - rx
    e22 = d21
    oy = ry + (e11 * e22 - e12 * e21) / (4A)
    F = number_type(_p)
    return (F(ox), F(oy))
end

"""
    triangle_circumcenter(tri::Triangulation, T) -> (Number, Number)

Computes the circumcenter of the triangle `T` in the triangulation `tri`.
"""
function triangle_circumcenter(tri::Triangulation, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    return triangle_circumcenter(p, q, r)
end

"""
    triangle_circumradius(p, q, r) -> Number

Computes the circumradius of the triangle with coordinates `(p, q, r)`. 
"""
function triangle_circumradius(p, q, r)
    ℓ₁², ℓ₂², ℓ₃² = squared_triangle_lengths(p, q, r)
    A = triangle_area(ℓ₁², ℓ₂², ℓ₃²)
    return triangle_circumradius(A, ℓ₁², ℓ₂², ℓ₃²)
end

"""
    triangle_offcenter(p, q, r, c₁=triangle_circumcenter(p, q, r), β=1.0) -> (Number, Number)

Computes the off-center of the triangle `(p, q, r)`.

# Arguments 
- `p`, `q`, `r`: The coordinates of the triangle, given in counter-clockwise order.
- `c₁=triangle_circumcenter(p, q, r)`: The circumcenter of the triangle.
- `β=1.0`: The radius-edge ratio cutoff. 

# Output 
- `cx`: The x-coordinate of the off-center.
- `cy`: The y-coordinate of the off-center.

!!! danger "Difference in definitions" 

    In the original [this paper](https://doi.org/10.1016/j.comgeo.2008.06.002), the off-center is defined to instead 
    be the circumcenter if it the triangle `pqc₁` has radius-edge ratio less than `β`. Here, we just let the off-center
    be the point `c` so that `pqc` has radius-edge ratio of exactly `β`.
"""
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

"""
    distance_to_offcenter(β, ℓ) -> Number 

Given a triangle with shortest edge length `ℓ`, computes the distance from the edge
to the offcenter of the triangle with radius-edge ratio cutoff `β`.
"""
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

"""
    make_shortest_edge_first(p, q, r, idx) -> (NTuple{2, Number}, NTuple{2, Number}, NTuple{2, Number})

Given a triangle `(p, q, r)`, rotate it (preserving orientation) so that the shortest edge is first. 
The argument `idx` gives the index of the shortest edge, where `idx == 1` means `(p, q)`,
`idx == 2` means `(q, r)`, and `idx == 3` means `(r, p)`.
"""
function make_shortest_edge_first(p, q, r, idx)
    if idx == 2
        return q, r, p
    elseif idx == 3
        return r, p, q
    else
        return p, q, r
    end
end

"""
    select_shortest_edge_for_offcenter(p, q, r, c, ℓ²) -> (NTuple{2, Number}, NTuple{2, Number}, NTuple{2, Number})

Given a triangle `(p, q, r)` with more than one edge attaining the shortest length, selects the appropriate shortest edge for [`triangle_offcenter`](@ref).

# Arguments
- `p`: The first vertex of the triangle.
- `q`: The second vertex of the triangle.
- `r`: The third vertex of the triangle.
- `c`: The circumcenter of the triangle.
- `ℓ²`: The squared shortest edge length.

The arguments should be so that `(p, q, r)` is positively oriented and `ℓ² = |p - q|²` is the squared shortest edge length.

# Outputs 
- `u`: The first vertex of the rotated triangle. 
- `v`: The second vertex of the rotated triangle.
- `w`: The third vertex of the rotated triangle.

These outputs `(u, v, w)` are a permutation of `(p, q, r)` (maintaining positive orientation) such that `|m - c₁|` is maximised over all other shortest edges, 
where `m = (u + v)/2`. If there is no unique maximiser, then the output is the permutation that is lexicographically smallest (i.e., sorted by x and then by y).
"""
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

"""
    triangle_perimeter(p, q, r) -> Number

Computes the perimeter of the triangle with coordinates `(p, q, r)`.
"""
function triangle_perimeter(p, q, r)
    ℓmin, ℓmed, ℓmax = triangle_lengths(p, q, r)
    return triangle_perimeter(ℓmin, ℓmed, ℓmax)
end

"""
    triangle_inradius(p, q, r) -> Number

Computes the inradius of the triangle with coordinates `(p, q, r)`.
"""
function triangle_inradius(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin, ℓmed, ℓmax = sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    return triangle_inradius(A, perimeter)
end

"""
    triangle_aspect_ratio(p, q, r) -> Number

Computes the aspect ratio of the triangle with coordinates `(p, q, r)`.
"""
function triangle_aspect_ratio(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin, ℓmed, ℓmax = sqrt(ℓmin²), sqrt(ℓmed²), sqrt(ℓmax²)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    perimeter = triangle_perimeter(ℓmin, ℓmed, ℓmax)
    inradius = triangle_inradius(A, perimeter)
    circumradius = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    return triangle_aspect_ratio(inradius, circumradius)
end

"""
    triangle_radius_edge_ratio(p, q, r) -> Number

Computes the radius-edge ratio of the triangle with coordinates `(p, q, r)`.
"""
function triangle_radius_edge_ratio(p, q, r)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    ℓmin = sqrt(ℓmin²)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    circumradius = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    return triangle_radius_edge_ratio(circumradius, ℓmin)
end

"""
    triangle_edge_midpoints(p, q, r) -> (Number, Number), (Number, Number), (Number, Number)

Computes the midpoints of the edges of the triangle with coordinates `(p, q, r)`.
"""
function triangle_edge_midpoints(p, q, r)
    return midpoint(p, q), midpoint(q, r), midpoint(r, p)
end

"""
    min_med_max(a, b, c) -> (Number, Number, Number)

Returns the arguments in sorted order.
"""
function min_med_max(a, b, c)
    b, c = minmax(b, c)
    a, c = minmax(a, c)
    a, b = minmax(a, b)
    return a, b, c
end

"""
    triangle_sink(p, q, r, tri::Triangulation) -> (Number, Number)

Computes the sink of each triangle in `tri`. See [this paper](https://doi.org/10.1145/378583.378644) for more information.

# Extended help 
Sinks were introduced in [this paper](https://doi.org/10.1145/378583.378644). For a given triangle `T`, the sink of 
`T` is defined as follows:

1. If `c`, the circumcenter of `T`, is in the interior of `T`, then the sink of `T` is `T`.
2. If `T` is a boundary triangle, then the sink of `T` is `T`.
3. If neither 1 or 2, then the sink is defined as the sink of the triangle `V`, where `V` is the triangle adjoining the edge of `T`
   which intersects the line `mc`, where `m` is the centroid of `T`.

In cases where the triangulation has holes, this definition can lead to loops. In such a case, we just pick one of the triangles 
in the loop as the sink triangle.
"""
function triangle_sink(tri::Triangulation, T, prev_T = construct_triangle(triangle_type(tri), integer_type(tri)(∅), integer_type(tri)(∅), integer_type(tri)(∅)))
    # TODO: This function would be faster if we just always search away from the largest angle.
    T = sort_triangle(T)
    c = triangle_circumcenter(tri, T)
    is_boundary_triangle(tri, T) && return c
    !has_ghost_triangles(tri) && any(e -> !edge_exists(tri, reverse_edge(e)), triangle_edges(T)) && return c
    flag = point_position_relative_to_triangle(tri, T, c)
    !is_outside(flag) && return c
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(tri, i, j, k)
    #=
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
    θ₃ ≤ π / 2 + ε && return c
    m = triangle_centroid(p, q, r)
    if !is_none(line_segment_intersection_type(p, q, m, c)) && !is_left(point_position_relative_to_line(p, q, c))
        next_T = construct_triangle(triangle_type(tri), j, i, get_adjacent(tri, j, i))
    elseif !is_none(line_segment_intersection_type(q, r, m, c)) && !is_left(point_position_relative_to_line(q, r, c))
        next_T = construct_triangle(triangle_type(tri), k, j, get_adjacent(tri, k, j))
    else # Must intersect the edge ki instead then 
        next_T = construct_triangle(triangle_type(tri), i, k, get_adjacent(tri, i, k))
    end
    sort_triangle(next_T) == prev_T && return c
    return triangle_sink(tri, next_T, T)
end
