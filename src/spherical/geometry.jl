"""
    abstract type AbstractSurface{T}

An abstract type for surfaces.
"""
abstract type AbstractSurface{T} end

"""
    struct UnitSphere{T} <: AbstractSurface{T}

A unit sphere centered at the origin. The type parameter `T` is the number type.
The constructor `UnitSphere()` creates a `UnitSphere{Float64}`.
"""
struct UnitSphere{T} <: AbstractSurface{T} end
UnitSphere() = UnitSphere{Float64}()

function Random.eltype(::Type{UnitSphere{T}}) where {T}
    return SphericalPoint{NTuple{3,T}}
end
function rand(rng::AbstractRNG, ::SamplerTrivial{UnitSphere{T}}) where {T}
    x, y, z = ntuple(_ -> randn(rng, T), Val(3))
    r⁻¹ = 1 / _hypot(x, y, z)
    return SphericalPoint((x * r⁻¹, y * r⁻¹, z * r⁻¹))
end

"""
    struct SphericalPoint{P}

Wrapper for representing a point on the surface of a [`UnitSphere`](@ref). The type parameter `P` is type of the wrapped point.

# Fields
- `p::P`: The point on the sphere.

# Constructors
- `SphericalPoint(x, y, z)`: Creates a `SphericalPoint` from the Cartesian coordinates `x`, `y`, and `z`.
- `SphericalPoint(p::SphericalPoint)`: Returns `p`.
- `SphericalPoint(latitude, longitude)`: Creates a `SphericalPoint` from the latitude and longitude. The resulting point is still on the unit sphere. The angles are assumed to be in degrees.
"""
struct SphericalPoint{P}
    p::P
end
SphericalPoint(x, y, z) = SphericalPoint((x, y, z))
SphericalPoint(p::SphericalPoint) = p
function SphericalPoint(latitude, longitude)
    lat, lon = latitude, longitude
    slat, clat = sincosd(lat)
    slon, clon = sincosd(lon)
    x = clat * clon
    y = clat * slon
    z = slat
    return SphericalPoint(x, y, z)
end
Base.show(io::IO, p::SphericalPoint) = print(io, "SphericalPoint", __getxyz(p), "")


Base.:(==)(p::SphericalPoint, q::SphericalPoint) = __getxyz(p) == __getxyz(q)

"""
    getp(p::SphericalPoint{P}) -> P

Returns the point `p` wrapped by `SphericalPoint`.
"""
getp(p::SphericalPoint) = p.p

"""
    __getxyz(p::SphericalPoint) -> Tuple{Number, Number, Number}

Returns the Cartesian coordinates of the point `p` wrapped by `SphericalPoint`.
"""
__getxyz(p::SphericalPoint) = getxyz(getp(p))

"""
    projectx(p::SphericalPoint) -> Number

Computes the `x`-coordinate of the stereographic projection of `p`, i.e. returns ``x/(1-z)``, where ``p = (x, y, z)``.
"""
function projectx(p::SphericalPoint)
    x, _, z = __getxyz(p)
    return x / (1 - z)
end

"""
    projecty(p::SphericalPoint) -> Number

Computes the `y`-coordinate of the stereographic projection of `p`, i.e. returns ``y/(1-z)``, where ``p = (x, y, z)``.
"""
function projecty(p::SphericalPoint)
    _, y, z = __getxyz(p)
    return y / (1 - z)
end

"""
    project(p::SphericalPoint) -> Tuple{Number, Number}

Computes the stereographic projection of `p`, i.e. returns the `Tuple` `(projectx(p), projecty(p))`.
"""
project(p::SphericalPoint) = (projectx(p), projecty(p))

"""
    invproject(p::Tuple{Number, Number}) -> SphericalPoint

Computes the inverse stereographic projection of `p`, i.e. returns the point on the sphere corresponding to `p`. The 
wrapped point is an `NTuple{3, Number}`.
"""
function invproject(p)
    X, Y = getxy(p)
    if isnan(X) && isnan(Y)
        return north_pole(number_type(p))
    end
    normsq = norm_sqr(p)
    den = 1 + normsq
    p = (2X / den, 2Y / den, (normsq - 1) / den)
    return SphericalPoint(p)
end

#=
"""
    getx(p::SphericalPoint) -> Number

Returns the `x`-coordinate of the stereographic projection of `p`. Note that this is 
*not* the `x`-coordinate of the point `p` on the sphere.
"""
=#
getx(p::SphericalPoint) = projectx(p)

#=
"""
    gety(p::SphericalPoint) -> Number

Returns the `y`-coordinate of the stereographic projection of `p`. Note that this is
*not* the `y`-coordinate of the point `p` on the sphere.
"""
=#
gety(p::SphericalPoint) = projecty(p)

#=
"""
    number_type(::Type{SphericalPoint{P}}) -> T

Returns the number type of the point `P` wrapped by `SphericalPoint`.
"""
=#
number_type(::Type{SphericalPoint{P}}) where {P} = number_type(P)

"""
    struct SphericalPoints{P}

Wrapper for representing a collection of points on the surface of a [`UnitSphere`](@ref). The type parameter `P` is the type of the wrapped points.
"""
struct SphericalPoints{T,P} <: AbstractVector{SphericalPoint{T}}
    points::P
    function SphericalPoints(points::AbstractVector{<:SphericalPoint{P}}) where {P}
        return new{P,typeof(points)}(points)
    end
end
SphericalPoints(p::SphericalPoints) = p

Base.eachindex(points::SphericalPoints) = Base.eachindex(points.points)
Base.iterate(points::SphericalPoints, state...) = Base.iterate(points.points, state...)
Base.size(points::SphericalPoints) = size(points.points)
Base.length(points::SphericalPoints) = length(points.points)
Base.getindex(points::SphericalPoints, i) = SphericalPoint(points.points[i])
number_type(::Type{SphericalPoints{P}}) where {P} = number_type(P)
set_point!(points::SphericalPoints, i, x, y) = setindex!(points.points, invproject((x, y)), i)
Base.pop!(points::SphericalPoints) = pop!(points.points)
push_point!(points::SphericalPoints, x, y) = push!(points.points, invproject((x, y)))
is_planar(points::SphericalPoints) = true # hack to avoid warnings

north_pole(::Type{T}) where {T} = SphericalPoint((zero(T), zero(T), one(T)))

function _getindex(points::SphericalPoints, i)
    if is_ghost_vertex(i)
        return north_pole(number_type(points))
    else 
        return getindex(points, i)
    end
end