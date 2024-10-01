"""
    EllipticalArc <: AbstractParametricCurve

Curve for representing an elliptical arc, parametrised over `0 ≤ t ≤ 1`. This curve can be evaluated
using `elliptical_arc(t)` and returns a tuple `(x, y)` of the coordinates of the point on the curve at `t`.

# Fields
- `center::NTuple{2,Float64}`: The center of the ellipse.
- `horz_radius::Float64`: The horizontal radius of the ellipse. 
- `vert_radius::Float64`: The vertical radius of the ellipse.
- `rotation_scales::NTuple{2,Float64}`: If `θ` is the angle of rotation of the ellipse, then this is `(sin(θ), cos(θ))`.
- `start_angle::Float64`: The angle of the initial point of the arc measured from `center`, in radians. This angle is measured from the center prior to rotating the ellipse.
- `sector_angle::Float64`: The angle of the sector of the arc, in radians. This is given by `end_angle - start_angle`, where `end_angle` is the angle at `last`, and so might be negative for negatively oriented arcs.
- `first::NTuple{2,Float64}`: The first point of the arc.
- `last::NTuple{2,Float64}`: The last point of the arc.

# Constructor
You can construct an `EllipticalArc` using 

    EllipticalArc(first, last, center, major_radius, minor_radius, rotation; positive=true)

where `rotation` is the angle of rotation of the ellipse, in degrees. The `positive` keyword argument is used to determine if the 
arc is positively oriented or negatively oriented.
"""
struct EllipticalArc <: AbstractParametricCurve
    center::NTuple{2,Float64}
    horz_radius::Float64
    vert_radius::Float64
    rotation_scales::NTuple{2,Float64}
    start_angle::Float64
    sector_angle::Float64
    first::NTuple{2,Float64}
    last::NTuple{2,Float64}
end
function Base.:(==)(e₁::EllipticalArc, e₂::EllipticalArc)
    e₁.center ≠ e₂.center && return false
    e₁.horz_radius ≠ e₂.horz_radius && return false
    e₁.vert_radius ≠ e₂.vert_radius && return false
    e₁.rotation_scales ≠ e₂.rotation_scales && return false
    e₁.start_angle ≠ e₂.start_angle && return false
    e₁.sector_angle ≠ e₂.sector_angle && return false
    return true
end

function EllipticalArc(p, q, c, α, β, θ°; positive=true)
    px, py = getxy(p)
    qx, qy = getxy(q)
    cx, cy = getxy(c)
    θ = deg2rad(θ°)
    sθ, cθ = sincos(θ)
    start_cost = inv(α) * (cθ * (px - cx) + sθ * (py - cy))
    start_sint = inv(β) * (-sθ * (px - cx) + cθ * (py - cy))
    start_angle = mod(atan(start_sint, start_cost), 2π)
    if p == q
        end_angle = start_angle
    else
        end_cost = inv(α) * (cθ * (qx - cx) + sθ * (qy - cy))
        end_sint = inv(β) * (-sθ * (qx - cx) + cθ * (qy - cy))
        end_angle = mod(atan(end_sint, end_cost), 2π)
    end
    start_angle, end_angle = adjust_θ(start_angle, end_angle, positive)
    sector_angle = end_angle - start_angle
    return EllipticalArc(c, α, β, (sθ, cθ), start_angle, sector_angle, p, q)
end
function (e::EllipticalArc)(t)
    if iszero(t)
        return e.first
    elseif isone(t)
        return e.last
    else
        c, α, β, (sinθ, cosθ), θ₀, Δθ = e.center, e.horz_radius, e.vert_radius, e.rotation_scales, e.start_angle, e.sector_angle
        t′ = Δθ * t + θ₀
        st, ct = sincos(t′)
        cx, cy = getxy(c)
        x = cx + α * ct * cosθ - β * st * sinθ
        y = cy + α * ct * sinθ + β * st * cosθ
        return (x, y)
    end
end

function differentiate(e::EllipticalArc, t)
    α, β, (sinθ, cosθ), θ₀, Δθ = e.horz_radius, e.vert_radius, e.rotation_scales, e.start_angle, e.sector_angle
    t′ = Δθ * t + θ₀
    st, ct = sincos(t′)
    x = -α * st * cosθ - β * ct * sinθ
    y = -α * st * sinθ + β * ct * cosθ
    return (x * Δθ, y * Δθ)
end

function twice_differentiate(e::EllipticalArc, t)
    α, β, (sinθ, cosθ), θ₀, Δθ = e.horz_radius, e.vert_radius, e.rotation_scales, e.start_angle, e.sector_angle
    t′ = Δθ * t + θ₀
    st, ct = sincos(t′)
    x = -α * ct * cosθ + β * st * sinθ
    y = -α * ct * sinθ - β * st * cosθ
    return (x * Δθ^2, y * Δθ^2)
end

function curvature(e::EllipticalArc, t)
    α, β, θ₀, Δθ = e.horz_radius, e.vert_radius, e.start_angle, e.sector_angle
    t′ = Δθ * t + θ₀
    st, ct = sincos(t′)
    return sign(Δθ) * α * β / (α^2 * st^2 + β^2 * ct^2)^(3 / 2)
end

function total_variation(e::EllipticalArc, t₁, t₂)
    if e.first == e.last && (t₁ == 0 && t₂ == 1)
        return 2π
    end
    T₁ = differentiate(e, t₁)
    T₂ = differentiate(e, t₂)
    if e.sector_angle > 0
        θ = angle_between(T₂, T₁)
    else
        θ = angle_between(T₁, T₂)
    end
    return θ
end

function point_position_relative_to_curve(kernel::AbstractPredicateKernel, e::EllipticalArc, p)
    x, y = getxy(p)
    c, α, β, (sinθ, cosθ), Δθ = e.center, e.horz_radius, e.vert_radius, e.rotation_scales, e.sector_angle
    cx, cy = getxy(c)
    x′ = x - cx
    y′ = y - cy
    x′′ = x′ * cosθ + y′ * sinθ
    y′′ = -x′ * sinθ + y′ * cosθ
    x′′′ = x′′ / α
    y′′′ = y′′ / β
    positive = Δθ > 0
    if positive
        a, b, c = (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)
    else
        a, b, c = (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)
    end
    cert = point_position_relative_to_circle(kernel, a, b, c, (x′′′, y′′′))
    if is_outside(cert)
        return Cert.Right
    elseif is_inside(cert)
        return Cert.Left
    else
        return Cert.On
    end
end

function get_inverse(e::EllipticalArc, p)
    if p == e.first
        return 0.0
    elseif p == e.last
        return 1.0
    end
    px, py = getxy(p)
    c, α, β, (sinθ, cosθ), θ₀, Δθ = e.center, e.horz_radius, e.vert_radius, e.rotation_scales, e.start_angle, e.sector_angle
    cx, cy = getxy(c)
    px′ = px - cx
    py′ = py - cy
    ct′ = inv(α) * (px′ * cosθ + py′ * sinθ)
    st′ = inv(β) * (-px′ * sinθ + py′ * cosθ)
    t′ = mod(atan(st′, ct′), 2π)
    t = (t′ - θ₀) / Δθ
    while t < 0
        t += 2π / abs(Δθ)
    end
    while t > 1
        t -= 2π / abs(Δθ)
    end
    return t
end

function _reverse(c::EllipticalArc)
    return EllipticalArc(
        c.center, 
        c.horz_radius,
        c.vert_radius,
        c.rotation_scales,
        c.start_angle + c.sector_angle,
        -c.sector_angle,
        c.last,
        c.first
    )
end