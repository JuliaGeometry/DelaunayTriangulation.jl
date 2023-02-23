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

"""
    getpoint(pts::P, i)

Given a collection of points `pts`, returns a `Tuple` 
of the `x` and `y` coordinates of the `i`th point in 
the collection. The methods currently defined are 

    getpoint(pts::AbstractVector, i)
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
function getpoint(pts::AbstractVector, i::Integer)
    pt = pts[i]
    x, y = getx(pt), gety(pt)
    return (x, y)
end
function getpoint(pts::AbstractMatrix, i::Integer)
    pt = @views pts[:, i]
    x, y = getx(pt), gety(pt)
    return (x, y)
end
getpoint(pts, p) = (getx(p), gety(p))

"""
    get_point(pts::P, i...)

Given a collection of points `pts`, returns the points 
corresponding to the indices in `i...`. This simply 
calls [`getpoint`](@ref) - you do not need to 
extend this method. 
"""
function get_point end
function get_point(pts::P, i) where {P}
    return getpoint(pts, i)
end
#function get_point(pts::P, i1, i2::Vararg{Any, N}) where {P,I,N} # need to split into two args to disambiguate with the method below
#    return (get_point(pts, i1), ntuple(j -> get_point(pts, i2[j]), Val(N))...)
#end
function get_point(pts::P, i::Vararg{Any,N}) where {P,N}
    return ntuple(j -> get_point(pts, i[j]), Val(N))
end

"""
    get_point(pts::P, boundary_map, i...)

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
the centroid for the boundary curve corresponding to the 
boundary index `i` is returned. See also [`get_representative_point_coordinates`](@ref).
"""
function get_point(pts::P, boundary_map::AbstractDict, i::I) where {P,I<:Integer}
    if !is_boundary_index(i)
        return get_point(pts, i)
    elseif i == I(DefaultAdjacentValue)
        throw(BoundsError(pts, I(DefaultAdjacentValue)))
    else
        F = number_type(pts)
        curve_index = get_curve_index(boundary_map, i)
        return get_representative_point_coordinates(curve_index, F)
    end
end
get_point(::P, ::AbstractDict, i) where {P} = (getx(i), gety(i))
function get_point(pts::P, boundary_map::AbstractDict, i::Vararg{Any,N}) where {P,N}
    return ntuple(j -> get_point(pts, boundary_map, i[j]), Val(N))
end

"""
    each_point_index(pts::P) where {P}

Given a collection of points `pts`, returns an iterator 
over the indices of the collection. The methods currently 
defined are 

    each_point_index(pts::AbstractVector)
    each_point_index(pts::AbstractMatrix) 

with the first returning `eachindex(pts)` and the second 
returning `axes(pts, 2)`. You can extend this function 
as you need.
"""
function each_point_index end
function each_point_index(::P) where {P}
    return error("The each_point_index function has not been defined for the type $P.")
end
each_point_index(pts::AbstractVector) = eachindex(pts)
each_point_index(pts::AbstractMatrix) = axes(pts, 2)

"""
    each_point(pts::P) where {p}

For a given collection of points `p`, returns an iterator that 
goes over each point in the collection. The methods currently 
defined are 

    each_point(pts::AbstractVector)
    each_point(pts::AbstractMatrix)

with the first method simply returning `pts`, and the second returning 
`eachcol(pts)`. You can extend this function as you need.
"""
function each_point end
function each_point(::P) where {P}
    return error("The each_point function has not been defined for the type $P.")
end
each_point(pts::AbstractVector) = pts
each_point(pts::AbstractMatrix) = eachcol(pts)

"""
    num_points(pts)

Returns the number of points in `pts`.
"""
function num_points end
function num_points(::P) where {P}
    return error("The num_points function has not been defined for the type $P.")
end
num_points(pts::AbstractVector) = length(pts)
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