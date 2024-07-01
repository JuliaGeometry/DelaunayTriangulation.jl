"""
    RepresentativeCoordinates{IntegerType, NumberType}

A mutable struct for representing the coordinates of a representative point of polygon or a set of points. 

# Fields 
- `x::NumberType`: The x-coordinate of the representative point.
- `y::NumberType`: The y-coordinate of the representative point.
- `n::IntegerType`: The number of points represented by the representative point.
"""
mutable struct RepresentativeCoordinates{I,T}
    x::T
    y::T
    n::I
end
function RepresentativeCoordinates{I,T}() where {I,T}
    return RepresentativeCoordinates{I,T}(zero(T), zero(T), zero(I))
end
function Base.:(==)(p::RepresentativeCoordinates, q::RepresentativeCoordinates)
    getx(p) ≠ getx(q) && return false
    gety(p) ≠ gety(q) && return false
    getn(p) ≠ getn(q) && return false
    return true
end
function Base.convert(::Type{RepresentativeCoordinates{I,T}}, c::RepresentativeCoordinates) where {I,T}
    x = getx(c)
    y = gety(c)
    n = getn(c)
    return RepresentativeCoordinates(T(x), T(y), I(n))
end
function Base.convert(::Type{RepresentativeCoordinates{I,T}}, c::RepresentativeCoordinates{I,T}) where {I,T}
    return c
end

"""
    getx(c::RepresentativeCoordinates) -> Number

Returns the x-coordinate of `c`.
"""
getx(c::RepresentativeCoordinates) = c.x

"""
    gety(c::RepresentativeCoordinates) -> Number

Returns the y-coordinate of `c`.
"""
gety(c::RepresentativeCoordinates) = c.y

"""
    getn(c::RepresentativeCoordinates) -> Integer

Returns the number of points represented by `c`.
"""
getn(c::RepresentativeCoordinates) = c.n

"""
    reset!(c::RepresentativeCoordinates)

Resets the coordinates of `c` to zero.
"""
function reset!(c::RepresentativeCoordinates{I,T}) where {I,T}
    c.x = zero(T)
    c.y = zero(T)
    c.n = zero(I)
    return c
end

"""
    add_point!(c::RepresentativeCoordinates, p)

Treating `c` as an arithmetic average, updates the coordinates of `c` to include `p`.
"""
function add_point!(c::RepresentativeCoordinates, p)
    n = getn(c)
    c.x = 1 / (n + 1) * (n * c.x + getx(p))
    c.y = 1 / (n + 1) * (n * c.y + gety(p))
    c.n += 1
    return c
end

"""
    delete_point!(c::RepresentativeCoordinates, p)

Treating `c` as an arithmetic average, updates the coordinates of `c` to exclude `p`.
"""
function delete_point!(c::RepresentativeCoordinates, p)
    n = getn(c)
    c.x = 1 / (n - 1) * (n * c.x - getx(p))
    c.y = 1 / (n - 1) * (n * c.y - gety(p))
    c.n -= 1
    return nothing
end

"""
    compute_centroid!(c::RepresentativeCoordinates, points)

Computes the centroid of `points` and stores the result in `c`.
"""
function compute_centroid!(c::RepresentativeCoordinates, points)
    reset!(c)
    for p in each_point(points)
        add_point!(c, p)
    end
    return c
end