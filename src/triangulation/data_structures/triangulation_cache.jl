struct TriangulationCache{T, M, I, S, F, IC, OC, IS}
    triangulation::T
    triangulation_2::T
    marked_vertices::M
    interior_segments_on_hull::I
    surrounding_polygon::S
    fan_triangles::F 
    incircle_cache::IC
    orient3_cache::OC
    insphere_cache::IS
end

const EMPTY_TRIANGULATION_CACHE = TriangulationCache(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)

@inline function Base.copy(cache::C)::C where {C<:TriangulationCache}
    if cache === EMPTY_TRIANGULATION_CACHE
        return cache
    else
        return _copy_cache(cache)
    end
end

function _copy_cache(cache::TriangulationCache; weights = copy(get_weights(get_triangulation(cache))))
    # Doesn't actually copy the values in the cache. Just generates a new cache.
    tri = get_triangulation(cache)
    points = copy(get_points(tri))
    I = integer_type(tri)
    E = edge_type(tri)
    Es = edges_type(tri)
    new_cache = TriangulationCache(points, I, E, Es, weights, Val(true))
    return new_cache
end

@inline get_triangulation(cache::TriangulationCache) = cache.triangulation
@inline get_triangulation_2(cache::TriangulationCache) = cache.triangulation_2
@inline get_marked_vertices(cache::TriangulationCache) = cache.marked_vertices
@inline get_interior_segments_on_hull(cache::TriangulationCache) = cache.interior_segments_on_hull
@inline get_surrounding_polygon(cache::TriangulationCache) = cache.surrounding_polygon
@inline get_fan_triangles(cache::TriangulationCache) = cache.fan_triangles
@inline get_incircle_cache(cache::TriangulationCache) = cache.incircle_cache
@inline get_orient3_cache(cache::TriangulationCache) = cache.orient3_cache
@inline get_insphere_cache(cache::TriangulationCache) = cache.insphere_cache

function Base.empty!(cache::TriangulationCache)
    tri = get_triangulation(cache)
    tri_2 = get_triangulation_2(cache)
    marked_vertices = get_marked_vertices(cache)
    interior_segments_on_hull = get_interior_segments_on_hull(cache)
    surrounding_polygon = get_surrounding_polygon(cache)
    fan_triangles = get_fan_triangles(cache)
    empty_unconstrained_triangulation(tri)
    empty_unconstrained_triangulation(tri_2)
    empty!(marked_vertices)
    empty!(interior_segments_on_hull)
    empty!(surrounding_polygon)
    empty!(fan_triangles)
    # No need to empty the predicate caches 
    return cache
end