"""
    CircularArc <: AbstractParametricCurve

Curve for representing a circular arc, parametrised over `0 ≤ t ≤ 1`. This curve can be evaluated 
using `circular_arc(t)` and returns a tuple `(x, y)` of the coordinates of the point on the curve at `t`.

# Fields 
- `center::NTuple{2,Float64}`: The center of the arc.
- `radius::Float64`: The radius of the arc.
- `start_angle::Float64`: The angle of the initial point of the arc, in radians.
- `sector_angle::Float64`: The angle of the sector of the arc, in radians. This is given by `end_angle - start_angle`, where `end_angle` is the angle at `last`, and so might be negative for negatively oriented arcs.
- `first::NTuple{2,Float64}`: The first point of the arc.
- `last::NTuple{2,Float64}`: The last point of the arc.
- `pqr::NTuple{3, NTuple{2, Float64}}`: Three points on the circle through the arc. This is needed for [`point_position_relative_to_curve`](@ref).

!!! warning "Orientation"

    The angles `start_angle` and `end_angle` should be setup such that `start_angle > end_angle` implies a positively oriented arc, 
    and `start_angle < end_angle` implies a negatively oriented arc. Moreover, they must be in `[0°, 2π°)`.

# Constructor 
You can construct a `CircularArc` using 

    CircularArc(first, last, center; positive=true)

It is up to you to ensure that `first` and `last` are equidistant from `center` - the radius used will be the 
distance between `center` and `first`. The `positive` keyword argument is used to determine if the 
arc is positively oriented or negatively oriented.
"""
struct CircularArc <: AbstractParametricCurve
    center::NTuple{2,Float64}
    radius::Float64
    start_angle::Float64
    sector_angle::Float64
    first::NTuple{2,Float64}
    last::NTuple{2,Float64}
    pqr::NTuple{3,NTuple{2,Float64}}
end
function Base.:(==)(c₁::CircularArc, c₂::CircularArc)
    c₁.center ≠ c₂.center && return false
    c₁.radius ≠ c₂.radius && return false
    c₁.start_angle ≠ c₂.start_angle && return false
    c₁.sector_angle ≠ c₂.sector_angle && return false
    return true
end

Base.copy(c::CircularArc) = c

function CircularArc(p, q, c; positive=true)
    px, py = getxy(p)
    qx, qy = getxy(q)
    cx, cy = getxy(c)
    r = dist((px, py), (cx, cy))
    θ₀ = mod(atan(py - cy, px - cx), 2π)
    if p == q
        θ₁ = θ₀
    else
        θ₁ = mod(atan(qy - cy, qx - cx), 2π)
    end
    θ₀, θ₁ = adjust_θ(θ₀, θ₁, positive)
    sector_angle = θ₁ - θ₀
    p′ = (cx + r, cy)
    q′ = (cx, cy + r)
    r′ = (cx - r, cy)
    if positive
        pqr = (p′, q′, r′)
    else
        pqr = (r′, q′, p′)
    end
    return CircularArc(c, r, θ₀, sector_angle, p, q, pqr)
end
function (c::CircularArc)(t)
    if iszero(t)
        return c.first
    elseif isone(t)
        return c.last
    else
        θ₀, Δθ = c.start_angle, c.sector_angle
        θ = Δθ * t + θ₀
        sθ, cθ = sincos(θ)
        cx, cy = getxy(c.center)
        x = cx + c.radius * cθ
        y = cy + c.radius * sθ
        return (x, y)
    end
end

function differentiate(c::CircularArc, t)
    θ₀, Δθ = c.start_angle, c.sector_angle
    θ = Δθ * t + θ₀
    sθ, cθ = sincos(θ)
    x = -c.radius * sθ
    y = c.radius * cθ
    return (x * Δθ, y * Δθ)
end

function twice_differentiate(c::CircularArc, t)
    θ₀, Δθ = c.start_angle, c.sector_angle
    θ = Δθ * t + θ₀
    sθ, cθ = sincos(θ)
    x = -c.radius * cθ
    y = -c.radius * sθ
    return (x * Δθ^2, y * Δθ^2)
end

function point_position_relative_to_curve(kernel::AbstractPredicateKernel, c::CircularArc, p)
    a, b, c = c.pqr
    cert = point_position_relative_to_circle(kernel, a, b, c, p)
    if is_outside(cert)
        return Cert.Right
    elseif is_inside(cert)
        return Cert.Left
    else
        return Cert.On
    end
end

arc_length(c::CircularArc) = c.radius * abs(c.sector_angle)
arc_length(c::CircularArc, t₁, t₂) = c.radius * abs(c.sector_angle) * (t₂ - t₁)

curvature(c::CircularArc, t) = sign(c.sector_angle) / c.radius

total_variation(c::CircularArc) = abs(c.sector_angle)
function total_variation(c::CircularArc, t₁, t₂)
    Δθ = c.sector_angle
    return abs(Δθ) * (t₂ - t₁)
end

get_equidistant_split(c::CircularArc, t₁, t₂) = midpoint(t₁, t₂)
get_equivariation_split(c::CircularArc, t₁, t₂) =
    let t = midpoint(t₁, t₂)
        (t, total_variation(c, t₁, t))
    end

function get_inverse(c::CircularArc, p)
    if p == c.first
        return 0.0
    elseif p == c.last
        return 1.0
    end
    px, py = getxy(p)
    cx, cy = getxy(c.center)
    r = c.radius
    cθ = (px - cx) / r
    sθ = (py - cy) / r
    θ = atan(sθ, cθ)
    Δθ, θ₀ = c.sector_angle, c.start_angle
    t = (θ - θ₀) / Δθ
    while t < 0
        t += 2π / abs(Δθ)
    end
    while t > 1
        t -= 2π / abs(Δθ)
    end
    return t
end

function _reverse(c::CircularArc)
    return CircularArc(
        c.center,
        c.radius,
        c.start_angle + c.sector_angle,
        -c.sector_angle,
        c.last,
        c.first,
        reverse(c.pqr)
    )
end