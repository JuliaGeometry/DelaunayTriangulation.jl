"""
    move_generator_to_centroid!(points, vorn::VoronoiTessellation, generator) -> Number 

Moves the generator `generator` to the centroid of its Voronoi cell. Returns the distance moved.

# Arguments 
- `points`: The underlying point set. This is a `deepcopy` of the points of the underlying triangulation. 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `generator`: The generator to move.

# Outputs 
- `dist`: The distance moved.
"""
function move_generator_to_centroid!(points, vorn::VoronoiTessellation, generator)
    c = get_centroid(vorn, generator)
    cx, cy = _getxy(c)
    p = get_generator(vorn, generator)
    δ = dist(c, p)
    set_point!(points, generator, cx, cy)
    return δ
end

"""
    default_displacement_tolerance(vorn::VoronoiTessellation) -> Number

Returns the default displacement tolerance for the centroidal smoothing algorithm. The default is given by 
`max_extent / 1e4`, where `max_extent = max(width, height)`, where `width` and `height` are the width and height 
of the bounding box of the underlying point set.

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).

# Outputs
- `tol`: The default displacement tolerance.
"""
function default_displacement_tolerance(vorn::VoronoiTessellation)
    xmin, xmax, ymin, ymax = polygon_bounds(vorn)
    max_extent = max(xmax - xmin, ymax - ymin)
    return max_extent / 1e4
end

"""
    _centroidal_smooth_itr(vorn::VoronoiTessellation, set_of_boundary_nodes, points, rng; kwargs...) -> (VoronoiTessellation, Number)

Performs a single iteration of the centroidal smoothing algorithm. 

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `set_of_boundary_nodes`: The set of boundary nodes in the underlying triangulation.
- `points`: The underlying point set. This is a `deepcopy` of the points of the underlying triangulation.
- `rng`: The random number generator.

# Keyword Arguments
- `kwargs...`: Extra keyword arguments passed to [`retriangulate`](@ref).

# Outputs 
- `vorn`: The updated [`VoronoiTessellation`](@ref).
- `max_dist`: The maximum distance moved by any generator.
"""
function _centroidal_smooth_itr(vorn::VoronoiTessellation, set_of_boundary_nodes, points, rng; kwargs...)
    F = number_type(vorn)
    max_dist = zero(F)
    for i in each_generator(vorn)
        if i ∉ set_of_boundary_nodes
            dist = move_generator_to_centroid!(points, vorn, i)
            max_dist = max(max_dist, dist)::F
        end
    end
    _tri = retriangulate(get_triangulation(vorn), points; rng, kwargs...)
    vorn = voronoi(_tri, clip=true)
    return vorn, max_dist
end

"""
    centroidal_smooth(vorn::VoronoiTessellation; maxiters=1000, tol=default_displacement_tolerance(vorn), rng=Random.default_rng(), kwargs...) -> VoronoiTessellation

Smooths `vorn` into a centroidal tessellation so that the new tessellation is of a set of generators whose associated Voronoi polygon is that polygon's centroid.

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).

# Keyword Arguments 
- `maxiters=1000`: The maximum number of iterations.
- `tol=default_displacement_tolerance(vorn)`: The displacement tolerance. See [`default_displacement_tolerance`](@ref) for the default. 
- `rng=Random.default_rng()`: The random number generator.
- `kwargs...`: Extra keyword arguments passed to [`retriangulate`](@ref).

# Outputs 
- `vorn`: The updated [`VoronoiTessellation`](@ref). This is not done in-place.

# Extended help 
The algorithm is simple. We iteratively smooth the generators, moving them to the centroid of their associated Voronoi polygon for the current tessellation, 
continuing until the maximum distance moved of any generator is less than `tol`. Boundary generators are not moved.
"""
function centroidal_smooth(vorn::VoronoiTessellation; maxiters=1000, tol=default_displacement_tolerance(vorn), rng=Random.default_rng(), kwargs...) 
    iter = 0
    F = number_type(vorn)
    max_dist = typemax(F)
    tri = get_triangulation(vorn)
    has_ghost = has_ghost_triangles(tri)
    !has_ghost && add_ghost_triangles!(tri)
    has_bnds = has_boundary_nodes(tri)
    !has_bnds && lock_convex_hull!(tri)
    set_of_boundary_nodes = get_all_boundary_nodes(tri)
    points = (deepcopy ∘ get_points)(tri)
    while iter < maxiters && max_dist > tol
        vorn, max_dist = _centroidal_smooth_itr(vorn, set_of_boundary_nodes, points, rng; kwargs...)
        iter += 1
    end
    !has_bnds && unlock_convex_hull!(tri)
    !has_ghost && delete_ghost_triangles!(tri)
    return vorn
end