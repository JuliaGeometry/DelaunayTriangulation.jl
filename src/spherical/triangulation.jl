"""
    SphericalTriangulation = Triangulation{<:SphericalPoints}

A type alias for a triangulation of points on the surface of a [`UnitSphere`](@ref).

For this triangulation, note that while ghost triangles are present, they should be interpreted as solid triangles - 
they are the triangles that help patch up the north pole. 

Note that the ghost vertex does not map to the north pole due to technical reasons. Use `north_pole` to obtain 
the north pole.

If you want to replace all these triangles with solid triangles, use [`solidify!`](@ref).
"""
const SphericalTriangulation = Triangulation{<:SphericalPoints}

north_pole(tri::SphericalTriangulation) = north_pole(number_type(tri))

function triangle_circumcenter(tri::SphericalTriangulation, T, _=nothing) # the unused argument is the triangle area
    points = get_points(tri)
    u, v, w = triangle_vertices(T)
    p, q, r = _getindex(points, u), _getindex(points, v), _getindex(points, w) # can't do get_point since DelaunayTriangulation projects it
    sc = spherical_circumcenter(p, q, r)
    return project(sc) # Voronoi vertices are assumed to be NTuple{2}
end

function triangle_area(tri::SphericalTriangulation, T)
    T = sort_triangle(T)
    points = get_points(tri)
    u, v, w = triangle_vertices(T)
    p, q, r = _getindex(points, u), _getindex(points, v), _getindex(points, w) # can't do get_point since DelaunayTriangulation projects it
    sc = spherical_triangle_area(p, q, r)
    return sc
end

function triangle_circumradius(tri::SphericalTriangulation, T, _=nothing) # the unused argument is the triangle area
    points = get_points(tri)
    u, v, w = triangle_vertices(T)
    p, q, r = _getindex(points, u), _getindex(points, v), _getindex(points, w) # can't do get_point since DelaunayTriangulation projects it
    sc = spherical_triangle_circumradius(p, q, r)
    return sc
end

function triangle_centroid(tri::SphericalTriangulation, T)
    points = get_points(tri)
    u, v, w = triangle_vertices(T)
    p, q, r = _getindex(points, u), _getindex(points, v), _getindex(points, w) # can't do get_point since DelaunayTriangulation projects it
    sc = spherical_triangle_centroid(p, q, r)
    return project(sc) 
end

function _check_args!(kwargs)
    if haskey(kwargs, :weights)
        throw(ArgumentError("Spherical triangulations do not support weighted points."))
    end
    if haskey(kwargs, :boundary_nodes)
        throw(ArgumentError("Spherical triangulations do not support boundary nodes."))
    end
    if haskey(kwargs, :segments)
        throw(ArgumentError("Spherical triangulations do not support segments."))
    end
    delete!(kwargs, :delete_ghosts)
    delete!(kwargs, :recompute_representative_points)
    return kwargs
end

"""
    spherical_triangulate(points::AbstractVector{<:SphericalPoint}; kwargs...) -> SphericalTriangulation

Computes the spherical Delaunay triangulation of the points `points` on the unit sphere. The keyword arguments are the same as
those for `triangulate`, except that `weights`, `boundary_nodes`, and `segments` are not supported.

See the docstring for [`SphericalTriangulation`](@ref) for more information.
"""
function spherical_triangulate(points::AbstractVector{<:SphericalPoint}; rng::AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel(), kwargs...)
    wpoints = SphericalPoints(points)
    kwargs = _check_args!(Dict(kwargs))
    T = number_type(wpoints)
    np = north_pole(T)
    idx = findfirst(==(np), each_point(wpoints))
    _skip_points = isnothing(idx) ? () : idx
    skip_points = haskey(kwargs, :skip_points) ? union(_skip_points, pop!(kwargs, :skip_points)) : _skip_points
    tri = triangulate(wpoints; skip_points, kwargs..., delete_ghosts=false, recompute_representative_points=false, rng, predicates)
    # lock_convex_hull!(tri; rng, predicates)

    #=
    empty_representative_points!(tri)
    I = integer_type(tri)
    rep = new_representative_point!(tri, one(I))
    rep[1] = RepresentativeCoordinates(T(NaN), T(NaN), zero(I))
    =# # This causes issues with point location :-(
    return tri::SphericalTriangulation
end

"""
    solidify!(tri::SphericalTriangulation)

Replaces all the ghost triangles representing triangles covering the patch around the north pole with 
solid triangles. The north pole is added as a point to the triangulation. Note that, after applying this 
function, the projected triangulation will no longer be planar as the new solid triangles will intersect 
other triangles. 

It is recommended that you just handle the ghost triangles directly instead of using this function.
"""
function solidify!(tri::SphericalTriangulation)
    np = north_pole(tri)
    idx = findfirst(==(np), each_point(tri))
    if isnothing(idx)
        push_point!(tri, np)
        n = num_points(tri)
    else 
        n = idx
    end
    for T in each_ghost_triangle(tri)
        u, v, w = triangle_vertices(sort_triangle(T))
        u′, v′, w′ = sort_triangle(u, v, w)
        add_triangle!(tri, v′, u′, n)
    end
    delete_ghost_triangles!(tri)
    delete_ghost_vertices_from_graph!(tri)
    return tri
end