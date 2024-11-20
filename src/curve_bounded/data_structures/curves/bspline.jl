struct BSpline <: AbstractParametricCurve
    control_points::Vector{NTuple{2,Float64}}
    knots::Vector{Int}
    cache::Vector{NTuple{2,Float64}}
    lookup_table::Vector{NTuple{2,Float64}}
    orientation_markers::Vector{Float64}
end
function Base.:(==)(b₁::BSpline, b₂::BSpline)
    b₁.control_points ≠ b₂.control_points && return false
    b₁.knots ≠ b₂.knots && return false
    return true
end

function Base.copy(b::BSpline)
    return BSpline(
        copy(b.control_points),
        copy(b.knots),
        copy(b.cache),
        copy(b.lookup_table),
        copy(b.orientation_markers)
    )
end

function BSpline(control_points::Vector{NTuple{2,Float64}}; degree=3, lookup_steps=5000, kwargs...)
    nc = length(control_points)
    @assert degree ≥ 1 lazy"Degree must be at least 1, got $degree."
    @assert degree ≤ nc - 1 lazy"Degree must be at most n - 1 = $(nc - 1), where n is the number of control points, got $degree."
    order = degree + 1
    cache = similar(control_points)
    knots = zeros(nc + order)
    for i in eachindex(knots)
        if i ≤ order
            knots[i] = 0
        elseif i < nc + 1
            knots[i] = knots[i-1] + 1
        else
            knots[i] = knots[nc] + 1
        end
    end
    lookup_table = similar(control_points, lookup_steps)
    markers = Float64[]
    spl = BSpline(control_points, knots, cache, lookup_table, markers)
    for i in 1:lookup_steps
        t = (i - 1) / (lookup_steps - 1)
        spl.lookup_table[i] = spl(t)
    end
    markers = orientation_markers(spl; kwargs...)
    resize!(spl.orientation_markers, length(markers))
    copyto!(spl.orientation_markers, markers)
    return spl
end

function (b::BSpline)(t)::NTuple{2,Float64}
    return _eval_bspline(b.control_points, b.knots, b.cache, t)
end

function de_boor!(control_points, knots, t)
    if iszero(t)
        return control_points[begin]
    elseif isone(t)
        return control_points[end]
    end
    nc = length(control_points)
    nk = length(knots)
    order = nk - nc
    domain = (order, nc + 1) # nc + 1 = nk - degree (order = degree + 1)
    a, b = knots[domain[1]], knots[domain[2]]
    t = a + t * (b - a)
    s = @views searchsortedfirst(knots[domain[1]:domain[2]], t) + domain[1] - 2
    for L in 1:order
        for i in s:-1:(s-order+L+1)
            numerator = t - knots[i]
            denominator = knots[i+order-L] - knots[i]
            α = numerator / denominator
            α′ = 1 - α
            xᵢ₋₁, yᵢ₋₁ = getxy(control_points[i-1])
            xᵢ, yᵢ = getxy(control_points[i])
            control_points[i] = (α′ * xᵢ₋₁ + α * xᵢ, α′ * yᵢ₋₁ + α * yᵢ)
        end
    end
    return control_points[s]
end
function _eval_bspline(control_points, knots, cache, t) # de Boor's algorithm
    if iszero(t)
        return control_points[begin]
    elseif isone(t)
        return control_points[end]
    end
    copyto!(cache, control_points)
    return de_boor!(cache, knots, t)
end

function differentiate(b::BSpline, t)
    copyto!(b.cache, b.control_points)
    nc = length(b.control_points)
    nk = length(b.knots)
    degree = nk - nc - 1
    for i in 1:(nc-1)
        xᵢ, yᵢ = getxy(b.control_points[i])
        xᵢ₊₁, yᵢ₊₁ = getxy(b.control_points[i+1])
        scale = degree / (b.knots[i+degree+1] - b.knots[i+1])
        b.cache[i] = (scale * (xᵢ₊₁ - xᵢ), scale * (yᵢ₊₁ - yᵢ))
    end
    deriv = @views de_boor!(b.cache[begin:(end-1)], b.knots[(begin+1):(end-1)], t)
    # Need to scale, since the formula used assumes that the knots are all in [0, 1]
    range = b.knots[end] - b.knots[begin]
    return (deriv[1] * range, deriv[2] * range)
end

function twice_differentiate(b::BSpline, t)
    copyto!(b.cache, b.control_points)
    nc = length(b.control_points)
    nk = length(b.knots)
    degree = nk - nc - 1
    if degree == 1
        return (0.0, 0.0)
    end
    for i in 1:(nc-2)
        xᵢ, yᵢ = getxy(b.control_points[i])
        xᵢ₊₁, yᵢ₊₁ = getxy(b.control_points[i+1])
        xᵢ₊₂, yᵢ₊₂ = getxy(b.control_points[i+2])
        scale1 = degree / (b.knots[i+degree+1] - b.knots[i+1])
        scale2 = degree / (b.knots[i+degree+2] - b.knots[i+2])
        scale3 = (degree - 1) / (b.knots[i+degree+1] - b.knots[i+2]) # different shifts between the knots are knots[begin+1:end-1]
        Qᵢ′x, Qᵢ′y = (scale1 * (xᵢ₊₁ - xᵢ), scale1 * (yᵢ₊₁ - yᵢ))
        Qᵢ₊₁′x, Qᵢ₊₁′y = (scale2 * (xᵢ₊₂ - xᵢ₊₁), scale2 * (yᵢ₊₂ - yᵢ₊₁))
        b.cache[i] = (scale3 * (Qᵢ₊₁′x - Qᵢ′x), scale3 * (Qᵢ₊₁′y - Qᵢ′y))
    end
    deriv = @views de_boor!(b.cache[begin:(end-2)], b.knots[(begin+2):(end-2)], t)
    range = (b.knots[end] - b.knots[begin])^2
    return (deriv[1] * range, deriv[2] * range)
end

function thrice_differentiate(b::BSpline, t) # yes there is a way to evaluate (B, B', B'', B''') all in one pass. just haven't implemented it
    copyto!(b.cache, b.control_points)
    nc = length(b.control_points)
    nk = length(b.knots)
    degree = nk - nc - 1
    if degree ≤ 2
        return (0.0, 0.0)
    end
    for i in 1:(nc-3)
        xᵢ, yᵢ = getxy(b.control_points[i])
        xᵢ₊₁, yᵢ₊₁ = getxy(b.control_points[i+1])
        xᵢ₊₂, yᵢ₊₂ = getxy(b.control_points[i+2])
        xᵢ₊₃, yᵢ₊₃ = getxy(b.control_points[i+3])
        scale1 = degree / (b.knots[i+degree+1] - b.knots[i+1])
        scale2 = degree / (b.knots[i+degree+2] - b.knots[i+2])
        scale3 = degree / (b.knots[i+degree+3] - b.knots[i+3])
        scale4 = (degree - 1) / (b.knots[i+degree+1] - b.knots[i+2])
        scale5 = (degree - 1) / (b.knots[i+degree+2] - b.knots[i+3])
        scale6 = (degree - 2) / (b.knots[i+degree+1] - b.knots[i+3])
        Qᵢ′x, Qᵢ′y = (scale1 * (xᵢ₊₁ - xᵢ), scale1 * (yᵢ₊₁ - yᵢ))
        Qᵢ₊₁′x, Qᵢ₊₁′y = (scale2 * (xᵢ₊₂ - xᵢ₊₁), scale2 * (yᵢ₊₂ - yᵢ₊₁))
        Qᵢ₊₂′x, Qᵢ₊₂′y = (scale3 * (xᵢ₊₃ - xᵢ₊₂), scale3 * (yᵢ₊₃ - yᵢ₊₂))
        Qᵢ′′x, Qᵢ′′y = (scale4 * (Qᵢ₊₁′x - Qᵢ′x), scale4 * (Qᵢ₊₁′y - Qᵢ′y))
        Qᵢ₊₁′′x, Qᵢ₊₁′′y = (scale5 * (Qᵢ₊₂′x - Qᵢ₊₁′x), scale5 * (Qᵢ₊₂′y - Qᵢ₊₁′y))
        b.cache[i] = (scale6 * (Qᵢ₊₁′′x - Qᵢ′′x), scale6 * (Qᵢ₊₁′′y - Qᵢ′′y))
    end
    deriv = @views de_boor!(b.cache[begin:(end-3)], b.knots[(begin+3):(end-3)], t)
    range = (b.knots[end] - b.knots[begin])^3
    return (deriv[1] * range, deriv[2] * range)
end

total_variation(b::BSpline, t₁, t₂) = marked_total_variation(b, t₁, t₂)

has_lookup_table(b::BSpline) = true

function _reverse(c::BSpline)
    return BSpline(
        reverse(c.control_points),
        c.knots, 
        reverse(c.cache),
        reverse(c.lookup_table),
        1 .- reverse(c.orientation_markers)
    )
end