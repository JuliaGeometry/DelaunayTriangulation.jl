"""
    TriangulationCache{T,M,I,S}

A cache to be used as a field in [`Triangulation`](@ref).

# Fields 
- `triangulation::T`: The cached triangulation. This will only refer to an unconstrained triangulation,
   meaning it cannot have any segments or boundary nodes. It will contain the `weights`. This 
   is used for constrained triangulations.
- `triangulation_2::T`: An extra cached triangulation. This is needed for retriangulating fans for constrained triangulations.
- `marked_vertices::M`: Marked vertices cache for use in constrained triangulations.
- `interior_segments_on_hull::I`: Interior segments in the triangulation that also appear on the 
   convex hull of `tri`. This is needed for [`lock_convex_hull!`](@ref) in case the convex hull also contains 
   interior segments.
- `surrounding_polygon::S`: The polygon surrounding the triangulation. This is needed for [`delete_point!`](@ref).
- `fan_triangles::F`: Triangles in a fan. This is needed for sorting fans for constrained triangulations.

!!! note "Caches of caches"

    The triangulation cache itself does not have a cache. Instead, it stores a 
    `TriangulationCache(nothing)`.

!!! danger "Aliasing"

    The `points` of the cache's `triangulation` will be aliased to the `points` of the parent triangulation.
"""
struct TriangulationCache{T, M, I, S, F}
    triangulation::T
    triangulation_2::T
    marked_vertices::M
    interior_segments_on_hull::I
    surrounding_polygon::S
    fan_triangles::F
end
function Base.show(io::IO, ::MIME"text/plain", cache::TriangulationCache)
    if isnothing(get_triangulation(cache)) # all fields are nothing, or none are 
        return print(io, "TriangulationCache with no storage.")
    else
        return print(io, "TriangulationCache with storage.")
    end
end


"""
    get_triangulation(cache::TriangulationCache) -> Triangulation

Returns the [`Triangulation`](@ref) stored in `cache`.
"""
get_triangulation(cache::TriangulationCache) = cache.triangulation


"""
    get_triangulation_2(cache::TriangulationCache) -> Triangulation

Returns the second [`Triangulation`](@ref) stored in `cache`.
"""
get_triangulation_2(cache::TriangulationCache) = cache.triangulation_2


"""
    get_marked_vertices(cache::TriangulationCache) -> Vector{Vertex}

Returns the marked vertices stored in `cache`.
"""
get_marked_vertices(cache::TriangulationCache) = cache.marked_vertices


"""
    get_interior_segments_on_hull(cache::TriangulationCache) -> Set{Edge}

Returns the interior segments on the convex hull of the triangulation stored in `cache`.
"""
get_interior_segments_on_hull(cache::TriangulationCache) = cache.interior_segments_on_hull


"""
    get_surrounding_polygon(cache::TriangulationCache) -> Vector{Vertex}

Returns the polygon surrounding the triangulation stored in `cache`.
"""
get_surrounding_polygon(cache::TriangulationCache) = cache.surrounding_polygon


"""
    get_fan_triangles(cache::TriangulationCache) -> Triangles 

Returns the triangles in a fan stored in `cache`.
"""
get_fan_triangles(cache::TriangulationCache) = cache.fan_triangles


"""
    empty!(cache::TriangulationCache) 

Empties the cache by emptying the triangulation stored in it.
"""
function Base.empty!(cache::TriangulationCache)
    tri = get_triangulation(cache)
    tri_2 = get_triangulation_2(cache)
    empty_unconstrained_triangulation!(tri)
    empty_unconstrained_triangulation!(tri_2)
    marked_vertices = get_marked_vertices(cache)
    empty!(marked_vertices)
    interior_segments_on_hull = get_interior_segments_on_hull(cache)
    empty!(interior_segments_on_hull)
    surrounding_polygon = get_surrounding_polygon(cache)
    empty!(surrounding_polygon)
    fan_triangles = get_fan_triangles(cache)
    empty!(fan_triangles)
    return cache
end


"""
    empty_unconstrained_triangulation!(tri::Triangulation)

Empties the triangulation `tri` by emptying its fields. Only works for unconstrained triangulaions.
"""
function empty_unconstrained_triangulation!(tri::Triangulation)
    triangles = get_triangles(tri)
    adjacent = get_adjacent(tri)
    adjacent2vertex = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    convex_hull = get_convex_hull(tri)
    representative_points = get_representative_point_list(tri)
    empty!(triangles)
    empty!(adjacent)
    empty!(adjacent2vertex)
    empty!(graph)
    empty!(convex_hull)
    empty!(representative_points)
    return tri
end