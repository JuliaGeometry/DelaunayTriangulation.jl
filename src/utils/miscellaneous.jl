@inline number_type(x) = number_type(typeof(x))
@inline number_type(::Type{T}) where {T<:AbstractArray} = number_type(eltype(T))
@inline number_type(::Type{<:NTuple{N,T}}) where {N,T} = number_type(T)
@inline number_type(::Type{<:NTuple{0}}) = Any
@inline number_type(::Type{Tuple{}}) = Any
@inline number_type(::Type{T}) where {T} = T

@inline is_true(b::Bool) = b
@inline is_true(b::Val{true}) = true
@inline is_true(b::Val{false}) = false

@inline function choose_uvw(e1, e2, _, u, v, w)
    e1 && return (u, v, w)
    e2 && return (v, w, u)
    return (w, u, v)
end

function get_ordinal_suffix(i)
    j = i % 10
    k = i % 100
    if j == 1 && k ≠ 11
        return "st"
    elseif j == 2 && k ≠ 12
        return "nd"
    elseif j == 3 && k ≠ 13
        return "rd"
    else
        return "th"
    end
end


@inline function adjust_θ(θ₁, θ₂, positive)
    if positive
        if θ₁ < θ₂
            return (θ₁, θ₂)
        elseif θ₁ > θ₂
            return (θ₁, θ₂ + 2π)
        else
            return (θ₁, θ₂ + 2π)
        end
    else
        if θ₁ < θ₂
            return (θ₁, θ₂ - 2π)
        elseif θ₁ > θ₂
            return (θ₁, θ₂)
        else
            return (θ₁, θ₂ - 2π)
        end
    end
end

function uniquetol(A::Vector{Float64}; tol=1.0e-12)
    isempty(A) && return Float64[]
    M = max(abs(A[begin]), abs(A[end])) # assuming A is sorted
    intol = let M = M, tol = tol
        (x, y) -> abs(x - y) ≤ M * tol
    end
    Auniq = [A[1]]
    for i in 2:lastindex(A)
        a = A[i]
        if !intol(a, Auniq[end])
            push!(Auniq, a)
        end
    end
    return Auniq
end

@inline _to_val(v::V) where {V} = Val(v)::Val{v}
@inline _to_val(v::Val{B}) where {B} = v

function angle_between(p, q)
    px, py = getxy(p)
    qx, qy = getxy(q)
    a = px * qx + py * qy
    b = px * -qy + py * qx
    return mod2pi(atan(b, a))
end

function project_onto_line(p, q, r)
    # Projects r onto the line through p and q. This is taken from the function 
    # two_point_interpolate from https://github.com/DanielVandH/NaturalNeighbours.jl/blob/13b807d4a0f733719d67dfd120a67d7c7cc8ce0a/src/interpolation/extrapolation.jl
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    ℓ² = dist_sqr(p, q)
    t = (rx - px) * (qx - px) + (ry - py) * (qy - py)
    t /= ℓ²
    cx = px + t * (qx - px)
    cy = py + t * (qy - py)
    return cx, cy
end

function iterated_neighbourhood(tri::Triangulation, i, d)
    I = integer_type(tri)
    neighbours = Set{I}()
    sizehint!(neighbours, ceil(I, 6^(d / 2)))
    return iterated_neighbourhood!(neighbours, tri, i, d)
end
function iterated_neighbourhood!(neighbours, tri, i, d)
    empty!(neighbours)
    i_neighbours = get_neighbours(tri, i)
    I = integer_type(tri)
    for j in i_neighbours
        if !is_ghost_vertex(j)
            push!(neighbours, j)
        end
    end
    for _ in 2:d
        new_neighbours = Set{I}() # don't want to mutate the iterator while iterating
        for j in neighbours
            for k in get_neighbours(tri, j)
                if k ≠ i && !is_ghost_vertex(k)
                    push!(new_neighbours, k)
                end
            end
        end
        union!(neighbours, new_neighbours)
    end
    return neighbours
end
