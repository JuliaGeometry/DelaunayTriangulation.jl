"""
    BezierCurve <: AbstractParametricCurve

Curve for representing a Bezier curve, parametrised over `0 ≤ t ≤ 1`. This curve can be evaluated
using `bezier_curve(t)` and returns a tuple `(x, y)` of the coordinates of the point on the curve at `t`.

A good reference on Bezier curves is [this](https://pomax.github.io/bezierinfo/).

See also [`BSpline`](@ref) and [`CatmullRomSpline`](@ref).

!!! danger "Loops"

    This curve is only tested on loop-free curves (and closed curves that otherwise have no self-intersections). It is not guaranteed to work on curves with loops, especially for finding the nearest point on the curve to a given point.

!!! danger "Interpolation"

    Remember that Bezier curves are not interpolation curves. They only go through the first and last control points, but not the intermediate ones. If you want an interpolation curve, use [`CatmullRomSpline`](@ref).

# Fields
- `control_points::Vector{NTuple{2,Float64}}`: The control points of the Bezier curve. The curve goes through the first and last control points, but not the intermediate ones.
- `cache::Vector{NTuple{2,Float64}}`: A cache of the points on the curve. This is used to speed up evaluation of the curve using de Casteljau's algorithm. 
- `lookup_table::Vector{NTuple{2,Float64}}`: A lookup table for the Bezier curve, used for finding the point on the curve closest to a given point. The `i`th entry of the lookup table
   corresponds to the `t`-value `i / (length(lookup_table) - 1)`.
- `orientation_markers::Vector{Float64}`: The orientation markers of the curve. These are defined so that the orientation of the curve is monotone between any two consecutive markers. The first and last markers are always `0` and `1`, respectively. See [`orientation_markers`](@ref).

!!! warning "Concurrency"

    The cache is not thread-safe, and so you should not evaluate this curve in parallel.

# Constructor
You can construct a `BezierCurve` using 

    BezierCurve(control_points::Vector{NTuple{2,Float64}}; lookup_steps=5000, kwargs...)

The keyword argument `lookup_steps=100` controls how many time points in `[0, 1]` are used for the lookup table. The `kwargs...` are keyword arguments passed to [`orientation_markers`](@ref).
"""
struct BezierCurve <: AbstractParametricCurve
    control_points::Vector{NTuple{2,Float64}}
    cache::Vector{NTuple{2,Float64}}
    lookup_table::Vector{NTuple{2,Float64}}
    orientation_markers::Vector{Float64}
end
function Base.:(==)(b₁::BezierCurve, b₂::BezierCurve)
    b₁.control_points ≠ b₂.control_points && return false
    return true
end

function Base.copy(b::BezierCurve)
    return BezierCurve(
        copy(b.control_points),
        copy(b.cache),
        copy(b.lookup_table),
        copy(b.orientation_markers)
    )
end

function BezierCurve(control_points::Vector{NTuple{2,Float64}}; lookup_steps=5000, kwargs...)
    cache = similar(control_points) # will be copyto! later
    lookup_table = similar(control_points, lookup_steps)
    markers = Float64[]
    spl = BezierCurve(control_points, cache, lookup_table, markers)
    for i in 1:lookup_steps
        t = (i - 1) / (lookup_steps - 1)
        spl.lookup_table[i] = spl(t)
    end
    markers = orientation_markers(spl; kwargs...)
    resize!(spl.orientation_markers, length(markers))
    copyto!(spl.orientation_markers, markers)
    return spl
end

function (b::BezierCurve)(t)::NTuple{2,Float64}
    return _eval_bezier_curve(b.control_points, b.cache, t)
end
function de_casteljau!(control_points, t)
    if iszero(t)
        return control_points[begin]
    elseif isone(t)
        return control_points[end]
    else
        n = length(control_points) - 1
        for j in 1:n
            for k in 1:(n-j+1)
                xₖ, yₖ = getxy(control_points[k])
                xₖ₊₁, yₖ₊₁ = getxy(control_points[k+1])
                x = (one(t) - t) * xₖ + t * xₖ₊₁
                y = (one(t) - t) * yₖ + t * yₖ₊₁
                control_points[k] = (x, y)
            end
        end
        return control_points[begin]
    end
end
function _eval_bezier_curve(control_points, cache, t) # de Casteljau's algorithm
    copyto!(cache, control_points)
    return de_casteljau!(cache, t)
end

function differentiate(b::BezierCurve, t)
    copyto!(b.cache, b.control_points)
    n = length(b.control_points) - 1
    for i in 1:n
        xᵢ, yᵢ = getxy(b.control_points[i])
        xᵢ₊₁, yᵢ₊₁ = getxy(b.control_points[i+1])
        b.cache[i] = (n * (xᵢ₊₁ - xᵢ), n * (yᵢ₊₁ - yᵢ))
    end
    return @views de_casteljau!(b.cache[begin:(end-1)], t)
end

function twice_differentiate(b::BezierCurve, t)
    copyto!(b.cache, b.control_points)
    n = length(b.control_points) - 1
    if n == 1
        return (0.0, 0.0)
    end
    for i in 1:(n-1)
        xᵢ, yᵢ = getxy(b.control_points[i])
        xᵢ₊₁, yᵢ₊₁ = getxy(b.control_points[i+1])
        xᵢ₊₂, yᵢ₊₂ = getxy(b.control_points[i+2])
        # To compute the second derivative control points, we note that e.g. 
        # for points [A, B, C, D], the first derivative gives [3(B - A), 3(C - B), 3(D - C)].
        # Let these new points be [A', B', C']. Thus, differentiating again, we obtain 
        # [2(B' - A'), 2(C - B')].
        # So, the ith point is given by 
        # Qᵢ′′ = (p-1) * (Qᵢ₊₁′ - Qᵢ′),
        # where 
        # Qᵢ′ = p * (Pᵢ₊₁ - Pᵢ).
        # Thus, Qᵢ′′ = p(p-1) * (Pᵢ₊₂ - 2Pᵢ₊₁ + Pᵢ).
        scale = n * (n - 1)
        b.cache[i] = (scale * (xᵢ₊₂ - 2xᵢ₊₁ + xᵢ), scale * (yᵢ₊₂ - 2yᵢ₊₁ + yᵢ))
    end
    return @views de_casteljau!(b.cache[begin:(end-2)], t)
end

function thrice_differentiate(b::BezierCurve, t)
    copyto!(b.cache, b.control_points)
    n = length(b.control_points) - 1
    if n ≤ 2
        return (0.0, 0.0)
    end
    for i in 1:(n-2)
        xᵢ, yᵢ = getxy(b.control_points[i])
        xᵢ₊₁, yᵢ₊₁ = getxy(b.control_points[i+1])
        xᵢ₊₂, yᵢ₊₂ = getxy(b.control_points[i+2])
        xᵢ₊₃, yᵢ₊₃ = getxy(b.control_points[i+3])
        # We know that Qᵢ′ = p(Pᵢ₊₁ - Pᵢ), where Qᵢ′ is the ith control point for the first derivative,
        # p is the degree of b, and Pᵢ is the ith control point of b. Thus,
        # Qᵢ′′ = (p-1)(Qᵢ₊₁′ - Qᵢ′) = p(p-1)(Pᵢ₊₂ - 2Pᵢ₊₁ + Pᵢ), and then 
        # Qᵢ′′′ = (p-2)(Qᵢ₊₁′′ - Qᵢ′′) = p(p-1)(p-2)(Pᵢ₊₃ - 3Pᵢ₊₂ + 3Pᵢ₊₁ - Pᵢ).
        scale = n * (n - 1) * (n - 2)
        b.cache[i] = (scale * (xᵢ₊₃ - 3xᵢ₊₂ + 3xᵢ₊₁ - xᵢ), scale * (yᵢ₊₃ - 3yᵢ₊₂ + 3yᵢ₊₁ - yᵢ))
    end
    return @views de_casteljau!(b.cache[begin:(end-3)], t)
end

total_variation(b::BezierCurve, t₁, t₂) = marked_total_variation(b, t₁, t₂)

has_lookup_table(b::BezierCurve) = true

function _reverse(c::BezierCurve)
    return BezierCurve(
        reverse(c.control_points),
        reverse(c.cache),
        reverse(c.lookup_table),
        1 .- reverse(c.orientation_markers)
    )
end