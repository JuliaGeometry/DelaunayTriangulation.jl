"""
    LineSegment <: AbstractParametricCurve

Curve for representing a line segment, parametrised over `0 ≤ t ≤ 1`. This curve can be using 
`line_segment(t)` and returns a tuple `(x, y)` of the coordinates of the point on the curve at `t`.

# Fields
- `first::NTuple{2,Float64}`: The first point of the line segment.
- `last::NTuple{2,Float64}`: The last point of the line segment.
- `length::Float64`: The length of the line segment.

# Constructor
You can construct a `LineSegment` using 

    LineSegment(first, last)
"""
struct LineSegment <: AbstractParametricCurve # line segment
    first::NTuple{2,Float64}
    last::NTuple{2,Float64}
    length::Float64
end
Base.:(==)(L₁::LineSegment, L₂::LineSegment) = L₁.first == L₂.first && L₁.last == L₂.last

function LineSegment(p₀, p₁)
    return LineSegment(p₀, p₁, dist(p₀, p₁))
end
function (L::LineSegment)(t)
    if iszero(t)
        return L.first
    elseif isone(t)
        return L.last
    else
        x₀, y₀ = getxy(L.first)
        x₁, y₁ = getxy(L.last)
        x = x₀ + t * (x₁ - x₀)
        y = y₀ + t * (y₁ - y₀)
        return (x, y)
    end
end

function differentiate(L::LineSegment, t)
    x₀, y₀ = getxy(L.first)
    x₁, y₁ = getxy(L.last)
    return (x₁ - x₀, y₁ - y₀)
end

twice_differentiate(::LineSegment, t) = (0.0, 0.0)

curvature(::LineSegment, t) = 0.0

total_variation(::LineSegment) = 0.0
total_variation(::LineSegment, t₁, t₂) = 0.0

"""
    point_position_relative_to_curve([kernel::AbstractPredicateKernel=AdaptiveKernel(),] L::LineSegment, p) -> Certificate

Returns the position of `p` relative to `L`, returning a [`Certificate`](@ref):

- `Left`: `p` is to the left of `L`.
- `Right`: `p` is to the right of `L`.
- `On`: `p` is on `L`.

See also [`point_position_relative_to_line`](@ref).

The `kernel` argument determines how this result is computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.
"""
function point_position_relative_to_curve(kernel::AbstractPredicateKernel, L::LineSegment, p)
    cert = point_position_relative_to_line(kernel, L.first, L.last, p)
    if is_collinear(cert)
        return Cert.On
    else
        return cert
    end
end

arc_length(L::LineSegment) = L.length
arc_length(L::LineSegment, t₁, t₂) = L.length * (t₂ - t₁)

get_equidistant_split(L::LineSegment, t₁, t₂) = midpoint(t₁, t₂)
get_equivariation_split(L::LineSegment, t₁, t₂) = (midpoint(t₁, t₂), 0.0)

function get_inverse(L::LineSegment, p)
    if p == L.first
        return 0.0
    elseif p == L.last
        return 1.0
    end
    px, py = getxy(p)
    x₀, y₀ = getxy(L.first)
    x₁, y₁ = getxy(L.last)
    if iszero(x₁ - x₀)
        return (py - y₀) / (y₁ - y₀)
    else
        return (px - x₀) / (x₁ - x₀)
    end
end

"""
    angle_between(L₁::LineSegment, L₂::LineSegment) -> Float64

Returns the angle between `L₁` and `L₂`, assuming that `L₁.last == L₂.first` (this is not checked). For consistency with 
If the segments are part of some domain, then the line segments should be oriented so that the interior is to the left of both segments.
"""
function angle_between(L₁::LineSegment, L₂::LineSegment)
    T₁ = differentiate(L₁, 1.0)
    T₂ = differentiate(L₂, 0.0)
    T₁x, T₁y = getxy(T₁)
    T₁′ = (-T₁x, -T₁y)
    θ = angle_between(T₁′, T₂)
    return θ
end

function get_circle_intersection(L::LineSegment, t₁, t₂, r)
    ℓ = L.length
    if iszero(t₁)
        t = r / ℓ
    elseif isone(t₁)
        t = 1 - r / ℓ
    else
        p, q = L(t₁), L(t₂)
        Ls = LineSegment(p, q)
        return get_circle_intersection(Ls, 0.0, 1.0, r)
    end
    return t, L(t)
end

function _reverse(L::LineSegment)
    return LineSegment(L.last, L.first, L.length)
end

is_linear(::LineSegment) = true