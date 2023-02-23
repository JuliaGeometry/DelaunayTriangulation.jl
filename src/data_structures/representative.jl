# The code in this file used to be implemented based on centroids of points, but this is problematic
# for non-convex polygons since the centroid could even be outside of the domain. This was realised 
# after seeing https://blog.jochentopf.com/2022-11-10-finding-representative-points-for-polygons.html,
# so now our implementation is inspired from this. The idea is to take the "visual center" of a polygon 
# (https://blog.mapbox.com/a-new-algorithm-for-finding-a-visual-center-of-a-polygon-7c77e6492fbc),
# essentially the point at the center of the largest inscribed circle of the region. This 
# implementation is in the file polygon_utils.jl.
#
# While the Delaunay triangulations are built, though, the region is convex, so the centroid view is 
# still fine - we compute a proper representative point only after a triangulation is complete, 
# meaning point location can still work fine.
#
# It is still possible for interior ghost triangles to leave their domain though (see the ghost triangles 
# in the "simple_geometry" function from helper_functions.jl) - I'm not sure how to fix this yet. Use the 
# representative point as an initial point and then calibrate the ghost vertex accordingly? Is this even 
# always possible? Maybe multiple ghost vertices are needed, which creates even more issues. Will 
# probably need to have a fallback to the brute force method. 
# 
# Note that, again, these issues do no matter when the triangulation is actually being built - it only 
# matters when using point location after the fact.
#

"""
    mutable struct RepresentativeCoordinates{I,T}

Struct representing the representative point of a region.

See also [`RepresentativePointList`](@ref).

# Fields 
- `x`: The `x`-coordinate. 
- `y`: The `y`-coordinate. 
- `n`: The number of points defining the centroid, if that is what is being used.

# Usage 

See `reset!`, `add_point!`, `delete_point!`, and `compute_representative_point!`.
"""
mutable struct RepresentativeCoordinates{I,T}
    x::T
    y::T
    n::I
end
function RepresentativeCoordinates{I,T}() where {I,T}
    return RepresentativeCoordinates{I,T}(T(0.0), T(0.0), I(0))
end

getx(c::RepresentativeCoordinates) = c.x
gety(c::RepresentativeCoordinates) = c.y
getn(c::RepresentativeCoordinates) = c.n

function Base.convert(::Type{RepresentativeCoordinates{I,T}},
                      c::RepresentativeCoordinates) where {I,T}
    x = getx(c)
    y = gety(c)
    n = getn(c)
    return RepresentativeCoordinates(T(x), T(y), I(n))
end
function Base.convert(::Type{RepresentativeCoordinates{I,T}},
                      c::RepresentativeCoordinates{I,T}) where {I,T}
    return c
end

function reset!(c::RepresentativeCoordinates{I,T}) where {I,T}
    c.x = zero(T)
    c.y = zero(T)
    c.n = zero(T)
    return nothing
end

function add_point!(c::RepresentativeCoordinates, p)
    n = getn(c)
    c.x = 1 / (n + 1) * (n * c.x + getx(p))
    c.y = 1 / (n + 1) * (n * c.y + gety(p))
    c.n += 1
    return nothing
end

function delete_point!(c::RepresentativeCoordinates, p)
    n = getn(c)
    c.x = 1 / (n - 1) * (n * c.x - getx(p))
    c.y = 1 / (n - 1) * (n * c.y - gety(p))
    c.n -= 1
    return nothing
end

function compute_centroid!(c::RepresentativeCoordinates, pts)
    reset!(c)
    for p in each_point(pts)
        add_point!(c, p)
    end
    return nothing
end

"""
    const RepresentativePointList = Dict{Int64,RepresentativeCoordinates{Int64, Float64}}()

This is a list of points for a set of representative points corresponding to multiple 
curves.

See also `update_centroid_after_addition`, `update_centroid_after_deletion`,
`reset_centroids!`, `empty_centroids!`, [`get_representative_point_coordinates`](@ref),
and `new_centroid`.

!!! warning
    While the function [`compute_representative_points!`](@ref) computes 
    an appropriate visual center of the polygon represented by the curves, i.e. by joining 
    points, the update functions like `update_centroid_after_addition` 
    and `update_centroid_after_deletion` treat the centroid as if it is a 
    mean of all the points. If you need to use the actual visual centers of the 
    polygons, you will need to reuse [`compute_representative_points!`](@ref) again (as is 
    done at the end of [`triangulate`](@ref)).
"""
const RepresentativePointList = Dict{Int64,RepresentativeCoordinates{Int64,Float64}}()

function reset_representative_points!()
    for c in values(RepresentativePointList)
        reset!(c)
    end
    return nothing
end

function empty_representative_points!()
    empty!(RepresentativePointList)
    return nothing
end

function update_centroid_after_addition!(i::I, p) where {I}
    centroid = get!(RepresentativeCoordinates{Int64,Float64}, RepresentativePointList, i)
    add_point!(centroid, p)
    return nothing
end

function update_centroid_after_deletion!(i::I, p) where {I}
    centroid = RepresentativePointList[i]
    delete_point!(centroid, p)
    return nothing
end

"""
    get_representative_point_coordinates(i::I, ::Type{T}) where {I, T}

Given an index `i` and some number type `T`, returns the coordinates 
corresponding to the centroid in [`RepresentativePointList`](@ref) at the key `i`.
The returned value is a `Tuple` `(T(x), T(y))` of the coordinates.
"""
function get_representative_point_coordinates(i::I, ::Type{T}) where {I,T}
    centroid = RepresentativePointList[i]
    x = getx(centroid)
    y = gety(centroid)
    return (T(x), T(y))
end

function new_representative_point(i::I, ::Type{T}) where {I,T}
    RepresentativePointList[i] = RepresentativeCoordinates(zero(T), zero(T), zero(I))
    return nothing
end

"""
    compute_representative_points!(points, boundary_nodes; precision = 1.0)

Given a list of `points` and a list of `boundary_nodes`, computes visual centers for 
the polygons represented by these curves using [`pole_of_inaccessibility`](@ref). The 
keyword argument `precision` is the precision used in [`pole_of_inaccessibility`](@ref).

See also [`RepresentativePointList`](@ref).

!!! warning
    While this function computes an appropriate visual center of the polygon represented
    by the curves, i.e. by joining points, the update functions like `update_centroid_after_addition` 
    and `update_centroid_after_deletion` treat the centroid as if it is a 
    mean of all the points. If you need to use the actual visual centers of the 
    polygons, you will need to reuse this function again (as is 
    done at the end of [`triangulate`](@ref)).
"""
function compute_representative_points!(points, boundary_nodes;
                                        precision=one(number_type(points)))
    empty_representative_points!()
    if has_multiple_curves(boundary_nodes)
        nc = num_curves(boundary_nodes)
        x, y = pole_of_inaccessibility(points, boundary_nodes; precision) # this is separate so that we make sure the center is not in any of the holes
        RepresentativePointList[1] = RepresentativeCoordinates(x, y, 0)
        for i in 2:nc
            bn = get_boundary_nodes(boundary_nodes, i)
            x, y = pole_of_inaccessibility(points, bn; precision)
            RepresentativePointList[i] = RepresentativeCoordinates(x, y, 0)
        end
    else
        x, y = pole_of_inaccessibility(points, boundary_nodes; precision)
        RepresentativePointList[1] = RepresentativeCoordinates(x, y, 0)
    end
    return nothing
end
