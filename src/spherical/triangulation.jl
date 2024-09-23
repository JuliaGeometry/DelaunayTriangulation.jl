"""
    SphericalTriangulation = Triangulation{<:SphericalPoints}

A type alias for a triangulation of points on the surface of a [`UnitSphere`](@ref).

For this triangulation, note that while ghost triangles are present, they should be interpreted as solid triangles - 
they are the triangles that help patch up the north pole. 

Note that the ghost vertex does not map to the north pole due to technical reasons. Use `north_pole` to obtain 
the north pole.

If you want to replace all these triangles with solid triangles, use [`solidify!`](@ref). This is not recommended.

Note also that the triangulations, and e.g. point location, are not as robust as the planar case. You may 
experience numerical issues near the north pole for example.
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

function point_position_relative_to_spherical_triangle(tri::SphericalTriangulation, V, s)
    u, v, w = triangle_vertices(V)
    points = get_points(tri)
    p, q, r = _getindex(points, u), _getindex(points, v), _getindex(points, w)
    return point_position_relative_to_spherical_triangle(p, q, r, s)
end

"""
    find_triangle(tri::SphericalTriangulation, q; check_sphere=true, kwargs...) 

Finds the spherical triangle in `tri` that contains the point `q`. The keyword arguments are the same as those for 
the usual [`find_triangle`], except that `check_sphere` can be used to check if the found triangle actually contains `q` 
(since the check is done in the stereographic projection). If it's not, it will be found by searching the triangles near the found 
triangle or using an exhaustive search.
"""
function find_triangle(tri::SphericalTriangulation, s; check_sphere = true, kwargs...)
    V = @invoke find_triangle(tri::Triangulation, s; kwargs..., check_sphere = false, use_barriers = Val(false))
    check_sphere || return V
    cert = point_position_relative_to_spherical_triangle(tri, V, s)
    if !is_outside(cert)
        return V 
    else 
        T = triangle_type(tri)
        u, v, w = triangle_vertices(V)
        ℓuv = get_adjacent(tri, v, u)
        ℓvw = get_adjacent(tri, w, v)
        ℓwu = get_adjacent(tri, u, w)
        V1 = construct_triangle(T, v, u, ℓuv)
        V2 = construct_triangle(T, w, v, ℓvw)
        V3 = construct_triangle(T, u, w, ℓwu)
        for _V in (V1, V2, V3)
            cert = point_position_relative_to_spherical_triangle(tri, _V, s)
            !is_outside(cert) && return _V
        end
        # Still not found. Find the nearest vertex and then look at its neighbours
        points = get_points(tri)
        p, q, r = _getindex(points, u), _getindex(points, v), _getindex(points, w)
        d1 = spherical_distance(p, s)
        d2 = spherical_distance(q, s)
        d3 = spherical_distance(r, s)
        d = min(d1, d2, d3)
        minidx = d == d1 ? u : d == d2 ? v : w
        for (i, j) in each_edge(get_adjacent2vertex(tri, minidx))
            _V = construct_triangle(T, i, j, minidx)
            cert = point_position_relative_to_spherical_triangle(tri, _V, s)
            !is_outside(cert) && return _V
        end
        # Still not found. Do an exhaustive search 
        for _V in each_triangle(tri)
            cert = point_position_relative_to_spherical_triangle(tri, _V, s)
            !is_outside(cert) && return _V
        end
        throw(PointNotFoundError(tri, s))
    end
end