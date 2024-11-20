struct CatmullRomSplineSegment <: AbstractParametricCurve
    a::NTuple{2,Float64}
    b::NTuple{2,Float64}
    c::NTuple{2,Float64}
    d::NTuple{2,Float64}
    p₁::NTuple{2,Float64}
    p₂::NTuple{2,Float64}
end

Base.copy(c::CatmullRomSplineSegment) = c

function (c::CatmullRomSplineSegment)(t)
    if iszero(t)
        return c.p₁
    elseif isone(t)
        return c.p₂
    else
        ax, ay = getxy(c.a)
        bx, by = getxy(c.b)
        cx, cy = getxy(c.c)
        dx, dy = getxy(c.d)
        cx = evalpoly(t, (dx, cx, bx, ax))
        cy = evalpoly(t, (dy, cy, by, ay))
        return (cx, cy)
    end
end

function catmull_rom_spline_segment(p₀, p₁, p₂, p₃, α, τ)
    x₀, y₀ = getxy(p₀)
    x₁, y₁ = getxy(p₁)
    x₂, y₂ = getxy(p₂)
    x₃, y₃ = getxy(p₃)
    t₀₁ = dist(p₀, p₁)^α
    t₁₂ = dist(p₁, p₂)^α
    t₂₃ = dist(p₂, p₃)^α
    τ′ = one(τ) - τ
    if iszero(τ′)
        m₁x, m₁y, m₂x, m₂y = zero(τ), zero(τ), zero(τ), zero(τ)
    else
        m₁x = τ′ * (x₂ - x₁ + t₁₂ * ((x₁ - x₀) / t₀₁ - (x₂ - x₀) / (t₀₁ + t₁₂)))
        m₁y = τ′ * (y₂ - y₁ + t₁₂ * ((y₁ - y₀) / t₀₁ - (y₂ - y₀) / (t₀₁ + t₁₂)))
        m₂x = τ′ * (x₂ - x₁ + t₁₂ * ((x₃ - x₂) / t₂₃ - (x₃ - x₁) / (t₁₂ + t₂₃)))
        m₂y = τ′ * (y₂ - y₁ + t₁₂ * ((y₃ - y₂) / t₂₃ - (y₃ - y₁) / (t₁₂ + t₂₃)))
    end
    ax = 2(x₁ - x₂) + m₁x + m₂x
    ay = 2(y₁ - y₂) + m₁y + m₂y
    bx = -3(x₁ - x₂) - 2m₁x - m₂x
    by = -3(y₁ - y₂) - 2m₁y - m₂y
    cx = m₁x
    cy = m₁y
    dx = x₁
    dy = y₁
    a = (ax, ay)
    b = (bx, by)
    c = (cx, cy)
    d = (dx, dy)
    return CatmullRomSplineSegment(a, b, c, d, p₁, p₂)
end

function differentiate(c::CatmullRomSplineSegment, t)
    ax, ay = getxy(c.a)
    bx, by = getxy(c.b)
    cx, cy = getxy(c.c)
    a′x, a′y = 3ax, 3ay
    b′x, b′y = 2bx, 2by
    x = evalpoly(t, (cx, b′x, a′x))
    y = evalpoly(t, (cy, b′y, a′y))
    return (x, y)
end

function twice_differentiate(c::CatmullRomSplineSegment, t)
    ax, ay = getxy(c.a)
    bx, by = getxy(c.b)
    a′′x, a′′y = 6ax, 6ay
    b′′x, b′′y = 2bx, 2by
    x = evalpoly(t, (b′′x, a′′x))
    y = evalpoly(t, (b′′y, a′′y))
    return (x, y)
end

function thrice_differentiate(c::CatmullRomSplineSegment, t)
    ax, ay = getxy(c.a)
    a′′′x, a′′′y = 6ax, 6ay
    return (a′′′x, a′′′y)
end

function _reverse(c::CatmullRomSplineSegment)
    #=
    If q(t) = at³ + bt² + ct + d, then
       q(1-t) = a′t³ + b′t² + c′t + d,
    where 
       a′ = -a 
       b′ = 3a + b 
       c′ = -3a - 2b - c 
       d′ = a + b + c + d.
    =#
    a, b, c, d, p₁, p₂ = c.a, c.b, c.c, c.d, c.p₁, c.p₂
    return CatmullRomSplineSegment(
        .-a,
        3 .* a .+ b,
        -3 .* a .- 2 .* b .- c,
        a .+ b .+ c .+ d,
        p₂,
        p₁
    )
end

struct CatmullRomSpline <: AbstractParametricCurve
    control_points::Vector{NTuple{2,Float64}}
    knots::Vector{Float64}
    lookup_table::Vector{NTuple{2,Float64}}
    alpha::Float64
    tension::Float64
    left::NTuple{2,Float64}
    right::NTuple{2,Float64}
    lengths::Vector{Float64}
    segments::Vector{CatmullRomSplineSegment}
    orientation_markers::Vector{Float64}
end

function Base.copy(c::CatmullRomSpline)
    return CatmullRomSpline(
        copy(c.control_points),
        copy(c.knots),
        copy(c.lookup_table),
        c.alpha,
        c.tension,
        c.left,
        c.right,
        copy(c.lengths),
        copy(c.segments),
        copy(c.orientation_markers)
    )
end

function _reverse(c::CatmullRomSpline)
    return CatmullRomSpline(
        reverse(c.control_points),
        1 .- reverse(c.knots),
        reverse(c.lookup_table),
        c.alpha,
        c.tension,
        c.right,
        c.left,
        reverse(c.lengths),
        reverse(reverse.(c.segments)),
        1 .- reverse(c.orientation_markers)
    )
end

function Base.:(==)(spl1::CatmullRomSpline, spl2::CatmullRomSpline)
    spl1.control_points ≠ spl2.control_points && return false
    spl1.knots ≠ spl2.knots && return false
    spl1.alpha ≠ spl2.alpha && return false
    spl1.tension ≠ spl2.tension && return false
    return true
end

is_interpolating(spl::CatmullRomSpline) = true

function CatmullRomSpline(control_points; _alpha=1 / 2, _tension=0.0, lookup_steps=5000, kwargs...)
    alpha = _alpha
    tension = _tension
    @assert length(control_points) ≥ 4 lazy"Catmull-Rom splines require at least 4 control points, got $(length(control_points))."
    nc = length(control_points)
    @assert 0 ≤ alpha ≤ 1 lazy"Alpha must be in [0, 1], got $alpha."
    @assert 0 ≤ tension ≤ 1 lazy"Tension must be in [0, 1], got $tension."
    knots = zeros(nc)
    for i in 2:nc
        knots[i] = knots[i-1] + dist(control_points[i-1], control_points[i])^alpha
    end
    left = extend_left_control_point(control_points)
    right = extend_right_control_point(control_points)
    scale = knots[end]
    knots ./= scale
    knots[end] = 1.0
    lookup_table = similar(control_points, lookup_steps)
    lengths = zeros(nc - 1)
    markers = Float64[]
    segments = Vector{CatmullRomSplineSegment}(undef, nc - 1)
    spl = CatmullRomSpline(control_points, knots, lookup_table, alpha, tension, left, right, lengths, segments, markers)
    for i in 1:(nc-1)
        spl.segments[i] = _get_segment(spl, i)
    end
    for i in 1:lookup_steps
        t = (i - 1) / (lookup_steps - 1)
        spl.lookup_table[i] = spl(t)
    end
    for i in 1:(nc-1)
        segment = get_segment(spl, i)
        spl.lengths[i] = arc_length(segment, 0.0, 1.0)
    end
    markers = orientation_markers(spl; kwargs...)
    resize!(spl.orientation_markers, length(markers))
    copyto!(spl.orientation_markers, markers)
    return spl
end

function extend_left_control_point(control_points)
    is_closed = control_points[begin] == control_points[end]
    if is_closed
        return control_points[end-1]
    else
        c₁, c₂, c₃, c₄ = control_points[begin], control_points[begin+1], control_points[begin+2], control_points[begin+3]
        x₁, x₂ = getx(c₁), getx(c₂)
        reverse_flag = x₁ == x₂
        if reverse_flag
            c₁, c₂, c₃, c₄ = reverse(getxy(c₁)), reverse(getxy(c₂)), reverse(getxy(c₃)), reverse(getxy(c₄))
            x₁, x₂ = getx(c₁), getx(c₂)
        end
        x = 2x₁ - x₂ # reflection
        y = thiele4(c₁, c₂, c₃, c₄, x)
        if !reverse_flag
            return (x, y)
        else
            return (y, x)
        end
    end
end
function extend_right_control_point(control_points)
    is_closed = control_points[begin] == control_points[end]
    if is_closed
        return control_points[begin+1]
    else
        cₙ₋₃, cₙ₋₂, cₙ₋₁, cₙ = control_points[end-3], control_points[end-2], control_points[end-1], control_points[end]
        xₙ₋₁, xₙ = getx(cₙ₋₁), getx(cₙ)
        reverse_flag = xₙ₋₁ == xₙ
        if reverse_flag
            cₙ₋₃, cₙ₋₂, cₙ₋₁, cₙ = reverse(getxy(cₙ₋₃)), reverse(getxy(cₙ₋₂)), reverse(getxy(cₙ₋₁)), reverse(getxy(cₙ))
            xₙ₋₁, xₙ = getx(cₙ₋₁), getx(cₙ)
        end
        x = 2xₙ - xₙ₋₁
        y = thiele4(cₙ₋₃, cₙ₋₂, cₙ₋₁, cₙ, x)
        if !reverse_flag
            return (x, y)
        else
            return (y, x)
        end
    end
end

# See https://github.com/JeffreySarnoff/CatmullRom.jl/tree/49e6536c184dfc4200b980f89aa55bc5cf357b82/src/fewpoints
function thiele4(p₁, p₂, p₃, p₄, x)
    x₁, y₁ = getxy(p₁)
    x₂, y₂ = getxy(p₂)
    x₃, y₃ = getxy(p₃)
    x₄, y₄ = getxy(p₄)
    y = thiele4(x₁, x₂, x₃, x₄, y₁, y₂, y₃, y₄, x)
    if !isfinite(y)
        if x ≤ x₃
            y = thiele3(p₁, p₂, p₃, x)
        else
            y = thiele3(p₂, p₃, p₄, x)
        end
    end
    return y
end
function thiele4(x₁, x₂, x₃, x₄, y₁, y₂, y₃, y₄, x)
    t₂ = (x₂ - x₃) / (y₂ - y₃)
    t₁ = (x₁ - x₂) / (y₁ - y₂)
    t₄ = -(x₃ - x₄) / (y₃ - y₄) + t₂
    t₃ = (x₁ - x₃) / (t₁ - t₂)
    t₄ = -(x₂ - x₄) / t₄ + y₂ - y₃ + t₃
    t₂ = -(x₁ - x₄) / t₄ + t₁ - t₂
    t₂ = (x₃ - x) / t₂ - y₁ + y₂ + t₃
    t₁ = -(x₂ - x) / t₂ + t₁
    y = -(x₁ - x) / t₁ + y₁
    return y
end
function thiele3(p₁, p₂, p₃, x)
    x₁, y₁ = getxy(p₁)
    x₂, y₂ = getxy(p₂)
    x₃, y₃ = getxy(p₃)
    y = thiele3(x₁, x₂, x₃, y₁, y₂, y₃, x)
    if !isfinite(y)
        y = quadratic_interp(p₁, p₂, p₃, x)
    end
    return y
end
function thiele3(x₁, x₂, x₃, y₁, y₂, y₃, x)
    t₁ = (x₁ - x₂) / (y₁ - y₂)
    t₂ = -(x₂ - x₃) / (y₂ - y₃) + t₁
    t₂ = (x₁ - x₃) / t₂ - y₁ + y₂
    t₁ = -(x₂ - x) / t₂ + t₁
    y = -(x₁ - x) / t₁ + y₁
    return y
end
function quadratic_interp(p₁, p₂, p₃, x)
    x₁, y₁ = getxy(p₁)
    x₂, y₂ = getxy(p₂)
    x₃, y₃ = getxy(p₃)
    y = quadratic_interp(x₁, x₂, x₃, y₁, y₂, y₃, x)
    return y
end
function quadratic_interp(x₁, x₂, x₃, y₁, y₂, y₃, x)
    t₁ = x₂ - x₁
    t₂ = x₁ - x₃
    t₃ = x₃ - x₂
    t₄ = t₂ * y₂
    t₅ = x₃^2
    t₆ = x₂^2
    t₇ = x₁^2
    s = -inv(t₁ * t₂ * t₃)
    a = (y₁ * t₃ + y₃ * t₁ + t₄) * x
    b = t₅ * (y₁ - y₂)
    c = t₆ * (y₃ - y₁)
    d = t₇ * (y₂ - y₃)
    q = t₆ * (y₁ * x₃ - y₃ * x₁)
    r = t₄ * x₁ * x₃
    n = (-y₁ * t₅ + y₃ * t₇) * x₂
    α = a - b - c - d
    β = r - n - q
    y = s * evalpoly(x, (β, α))
    return y
end

function (c::CatmullRomSpline)(t)
    if iszero(t)
        return c.control_points[begin]
    elseif isone(t)
        return c.control_points[end]
    else
        segment, i = get_segment(c, t)
        t′ = map_t_to_segment(c, i, t)
        return segment(t′)
    end
end

function map_t_to_segment(c::CatmullRomSpline, i, t)
    tᵢ, tᵢ₊₁ = c.knots[i], c.knots[i+1]
    t′ = (t - tᵢ) / (tᵢ₊₁ - tᵢ)
    return t′
end

function differentiate(c::CatmullRomSpline, t)
    segment, i = get_segment(c, t)
    tᵢ, tᵢ₊₁ = c.knots[i], c.knots[i+1]
    t′ = map_t_to_segment(c, i, t)
    ∂x, ∂y = getxy(differentiate(segment, t′))
    scale = inv(tᵢ₊₁ - tᵢ)
    return (scale * ∂x, scale * ∂y)
end

function twice_differentiate(c::CatmullRomSpline, t)
    segment, i = get_segment(c, t)
    tᵢ, tᵢ₊₁ = c.knots[i], c.knots[i+1]
    t′ = map_t_to_segment(c, i, t)
    ∂x, ∂y = getxy(twice_differentiate(segment, t′))
    scale = inv(tᵢ₊₁ - tᵢ)
    return (scale^2 * ∂x, scale^2 * ∂y)
end

function thrice_differentiate(c::CatmullRomSpline, t)
    segment, i = get_segment(c, t)
    tᵢ, tᵢ₊₁ = c.knots[i], c.knots[i+1]
    t′ = map_t_to_segment(c, i, t)
    ∂x, ∂y = getxy(thrice_differentiate(segment, t′))
    scale = inv(tᵢ₊₁ - tᵢ)
    return (scale^3 * ∂x, scale^3 * ∂y)
end

function get_segment(c::CatmullRomSpline, t)
    i = min(lastindex(c.knots) - 1, searchsortedlast(c.knots, t)) # avoid issues at t = 1
    return get_segment(c, i), i
end
function _get_segment(c::CatmullRomSpline, i::Int)
    if i == firstindex(c.control_points)
        pᵢ₋₁ = c.left
    else
        pᵢ₋₁ = c.control_points[i-1]
    end
    if i == lastindex(c.control_points) - 1
        pᵢ₊₂ = c.right
    else
        pᵢ₊₂ = c.control_points[i+2]
    end
    pᵢ, pᵢ₊₁ = c.control_points[i], c.control_points[i+1]
    segment = catmull_rom_spline_segment(pᵢ₋₁, pᵢ, pᵢ₊₁, pᵢ₊₂, c.alpha, c.tension)
    return segment
end
function get_segment(c::CatmullRomSpline, i::Int)
    return c.segments[i]
end

function arc_length(c::CatmullRomSpline)
    return sum(c.lengths)
end
function arc_length(c::CatmullRomSpline, t₁, t₂)
    segment₁, i₁ = get_segment(c, t₁)
    segment₂, i₂ = get_segment(c, t₂)
    if i₁ == i₂ # same segment - just integrate directly
        t₁′ = map_t_to_segment(c, i₁, t₁)
        t₂′ = map_t_to_segment(c, i₂, t₂)
        return arc_length(segment₁, t₁′, t₂′)
    elseif i₁ == i₂ - 1 # adjacent segments - integrate on each segment 
        t₁′ = map_t_to_segment(c, i₁, t₁)
        t₂′ = map_t_to_segment(c, i₂, t₂)
        s₁ = arc_length(segment₁, t₁′, 1.0)
        s₂ = arc_length(segment₂, 0.0, t₂′)
        return s₁ + s₂
    else # at least one complete segment separates the outer segments 
        s = 0.0
        for i in (i₁+1):(i₂-1)
            s += c.lengths[i]
        end
        t₁′ = map_t_to_segment(c, i₁, t₁)
        t₂′ = map_t_to_segment(c, i₂, t₂)
        s₁ = arc_length(segment₁, t₁′, 1.0)
        s₂ = arc_length(segment₂, 0.0, t₂′)
        return s + s₁ + s₂
    end
end

total_variation(c::CatmullRomSpline, t₁, t₂) = marked_total_variation(c, t₁, t₂)

has_lookup_table(c::CatmullRomSpline) = true