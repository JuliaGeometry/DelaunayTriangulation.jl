"""
    getx(p::P) where {P}

Given a point `p`, returns the `x`-coordinate. The 
only methods currently defined are

    getx(p::NTuple{N,T}) where {N,T}
    getx(p::AbstractVector)

You can extend this function as you need.
"""
function getx end
function getx(::F) where {F}
    return error("The getx function has not been defined for the type $F.")
end
getx(p::NTuple{N,T}) where {N,T} = p[1]
getx(p::AbstractVector) = p[1]

"""
    gety(p::P) where {P}

Given a point `p`, returns the `y`-coordinate. The 
only methods currently defined are

    gety(p::NTuple{N,T}) where {N,T}
    gety(p::AbstractVector)

You can extend this function as you need.
"""
function gety end
function gety(::F) where {F}
    return error("The gety function has not been defined for the type $F.")
end
gety(p::NTuple{N,T}) where {N,T} = p[2]
gety(p::AbstractVector) = p[2]

"""
    getxy(p)

Given a point `p`, returns `(getx(p), gety(p))`.
"""
getxy(p) = (getx(p), gety(p))

@inline _getx(p) = Float64(getx(p))
@inline _gety(p) = Float64(gety(p))
@inline _getxy(p) = (_getx(p), _gety(p))

"""
    getpoint(pts::P, i)

Given a collection of points `pts`, returns a `Tuple` 
of the `x` and `y` coordinates of the `i`th point in 
the collection. The methods currently defined are 

    getpoint(pts::AbstractVecOrTuple, i)
    getpoint(pts::AbstractMatrix, i)

You can extend this function as you need. 

It is assumed that whenever `i` is not an integer, `i` is meant to be 
a point, so `(getx(i), gety(i))` would be returned in that case. This 
makes it easier to use some predicates without having to know the index 
of the point, simply passing the point directly.
"""
function getpoint end
function getpoint(::P, ::Integer) where {P}
    return error("The getpoint function has not been defined for the type $P.")
end
function getpoint(pts::Union{AbstractVector,Tuple}, i::Integer)
    pt = pts[i]
    x, y = getx(pt), gety(pt)
    return (x, y)
end
function getpoint(pts::AbstractMatrix, i::Integer)
    pt = @views pts[:, i]
    x, y = getx(pt), gety(pt)
    return (x, y)
end
getpoint(pts, p) = (getx(p), gety(p)) # so that we can mix points and vertices

"""
    get_point(pts::P, i...)

Given a collection of points `pts`, returns the points 
corresponding to the indices in `i...`. This simply 
calls [`getpoint`](@ref) - you do not need to 
extend this method. 
"""
function get_point end
get_point(pts::P, i) where {P} = getpoint(pts, i)
function get_point(pts::P, i::Vararg{Any,N}) where {P,N}
    return ntuple(j -> get_point(pts, i[j]), Val(N))
end

"""
    get_point(pts::P, representative_point_list, boundary_map, i...)

Given points `pts`, a boundary map from [`construct_boundary_map`](@ref) that 
maps boundary indices to their corresponding curve, and some indices `i...`, returns 
the `i`th point from `pts`. If `i` is just a single integer, 
the result is a `Tuple` of the coordinates, but if there are multiple 
integers provided then a `Tuple` of `Tuples` of these coordinates 
is returned, with the `j`th `Tuple` the coordinates for the `i[j]`th 
point. 

In the single `i` case, if it is not an integer then we simply return 
`(getx(i), gety(i))`, assuming that it is supposed to represent a single point.

If `is_boundary_index(i)`, then instead of returning the `i`th point, 
the representative point for the boundary curve corresponding to the 
boundary index `i` is returned, using the `representative_point_list` argument.
See also [`get_representative_point_coordinates`](@ref).
"""
function get_point(pts::P, representative_point_list::AbstractDict, boundary_map::AbstractDict, i::I) where {P,I<:Integer}
    if !is_boundary_index(i)
        return get_point(pts, i)
    elseif i == I(DefaultAdjacentValue)
        throw(BoundsError(pts, I(DefaultAdjacentValue)))
    else
        curve_index = get_curve_index(boundary_map, i)
        return get_representative_point_coordinates(representative_point_list, curve_index)
    end
end
get_point(::P, ::AbstractDict, ::AbstractDict, i) where {P} = (getx(i), gety(i))
function get_point(pts::P, representative_point_list::AbstractDict, boundary_map::AbstractDict, i::Vararg{Any,N}) where {P,N}
    return ntuple(j -> get_point(pts, representative_point_list, boundary_map, i[j]), N)
end

"""
    each_point_index(pts::P) where {P}

Given a collection of points `pts`, returns an iterator 
over the indices of the collection. The methods currently 
defined are 

    each_point_index(pts::AbstractVecOrTuple)
    each_point_index(pts::AbstractMatrix) 

with the first returning `eachindex(pts)` and the second 
returning `axes(pts, 2)`. You can extend this function 
as you need.
"""
function each_point_index end
function each_point_index(::P) where {P}
    return error("The each_point_index function has not been defined for the type $P.")
end
each_point_index(pts::Union{AbstractVector,Tuple}) = eachindex(pts)
each_point_index(pts::AbstractMatrix) = axes(pts, 2)

"""
    each_point(pts::P) where {p}

For a given collection of points `p`, returns an iterator that 
goes over each point in the collection. The methods currently 
defined are 

    each_point(pts::AbstractVecOrTuple)
    each_point(pts::AbstractMatrix)

with the first method simply returning `pts`, and the second returning 
`eachcol(pts)`. You can extend this function as you need.
"""
function each_point end
function each_point(::P) where {P}
    return error("The each_point function has not been defined for the type $P.")
end
each_point(pts::Union{AbstractVector,Tuple}) = pts
each_point(pts::AbstractMatrix) = eachcol(pts)

"""
    num_points(pts)

Returns the number of points in `pts`. The methods currently defined are 

    num_points(pts::AbstractVecOrTuple)
    num_points(pts::AbstractMatrix)

with the first returning `length(pts)` and the second returning
`size(pts, 2)`. You can extend this function as you need.
"""
function num_points end
function num_points(::P) where {P}
    return error("The num_points function has not been defined for the type $P.")
end
num_points(pts::Union{AbstractVector,Tuple}) = length(pts)
num_points(pts::AbstractMatrix) = size(pts, 2)

"""
    points_are_unique(pts)

Returns `true` if `pts` has no duplicate points, and `false` otherwise.
"""
function points_are_unique(pts)
    n = num_points(pts)
    p = first(each_point(pts))
    seen = Set{typeof(p)}()
    for p in each_point(pts)
        p ∈ seen && return false
        push!(seen, p)
    end
    return n == length(seen)
end

"""
    lexicographic_order(pts)

Returns a set of indices `idx` that gives the lexicographic ordering 
of the set of points `pts`, i.e. sorting by `x` and then sorting points 
with duplicate `x`-coordinates by `y`. The implementation is simply 

    lexicographic_order(pts) = (sortperm ∘ collect ∘ each_point)(pts)

which you might want to specialise for an easier representation of your 
points `pts`.
"""
function lexicographic_order end
lexicographic_order(pts) = (sortperm ∘ collect ∘ each_point)(pts)

"""
    push_point!(pts, x, y)

Pushes the point `(x, y)` into `pts`. The only methods currently 
defined are  

    push_point!(pts::AbstractVector{T}, x, y) where {F,T<:NTuple{2,F}} = push!(pts, (F(x), F(y)))
    push_point!(pts::AbstractVector{T}, x, y) where {F<:Number,T<:AbstractVector{F}} = push!(pts, F[x, y])

You can extend this function as needed. We also provide the method 

    push_point!(pts, p) = push_point!(pts, getx(p), gety(p))

which you can extend if you have a point type `p` that has `getx` and `gety`.
"""
function push_point! end
function push_point!(::P, ::Any, ::Any) where {P}
    return error("The push_point! function has not been defined for the collection type $P.")
end
push_point!(pts::AbstractVector{T}, x, y) where {F,T<:NTuple{2,F}} = push!(pts, (F(x), F(y)))
push_point!(pts::AbstractVector{T}, x, y) where {F<:Number,T<:AbstractVector{F}} = push!(pts, F[x, y])
push_point!(pts, p) = push_point!(pts, getx(p), gety(p))
push_point!(pts::AbstractMatrix{T}, x, y) where {T<:Number} = append!(pts, (x, y)) # ElasticArrays 

"""
    pop_point!(pts)

Pops the last point from `pts`. The only method currently defined is

    pop_point!(pts::AbstractVector) = pop!(pts)

You can extend this function as needed.
"""
function pop_point! end
function pop_point!(::P) where {P}
    return error("The pop_point! function has not been defined for the collection type $P.")
end
pop_point!(pts::AbstractVector) = pop!(pts)
pop_point!(pts::AbstractMatrix) = resize!(pts, (2, size(pts, 2) - 1)) # ElasticArrays

"""
    mean_points(points, vertices = each_point_index(points))

Returns the mean of the points in `points` indexed by `vertices`, given as a `Tuple` of the form `(mean_x, mean_y)`.
"""
function mean_points(points, vertices = each_point_index(points))
    F = number_type(points)
    cx = zero(F)
    cy = zero(F)
    n = 0
    for v in vertices 
        p = get_point(points, v)
        px, py = getxy(p)
        cx += px
        cy += py
        n += 1
    end
    return (cx / n, cy / n)
end

"""
    set_point!(points, i, x, y)

Sets the point at index `i` in `points` to `(x, y)`. The only methods currently
defined are 

    set_point!(points::AbstractVector{T}, i, x, y) = (points[i] = (x, y))
    set_point!(points::AbstractMatrix{T}, i, x, y) where {T} = (points[1, i] = x; points[2, i] = y)

You can extend this function as needed. We also define 

    set_point!(points, i, p) = set_point!(points, i, getx(p), gety(p))
"""
function set_point! end
function set_point!(points::AbstractVector, i, x, y)
    points[i] = (x, y) # also works for eltype(points) <: Tuple, StaticArrays, etc.
end
function set_point!(points::AbstractVector{T}, i, x, y) where {T<:Vector}
    points[i] = [x, y]
end
function set_point!(points::AbstractMatrix{T}, i, x, y) where {T}
    points[1, i] = x
    points[2, i] = y
end
set_point!(points, i, p) = set_point!(points, i, getx(p), gety(p))
