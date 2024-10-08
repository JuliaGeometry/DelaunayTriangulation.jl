using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StructEquality
using Setfield
using AdaptivePredicates
using ..DelaunayTriangulation: add_weight!, get_weight, get_weights

@struct_equal DT.TriangulationCache

cache = DT.TriangulationCache()
@test isnothing(DT.get_triangulation(cache))
@test isnothing(DT.get_triangulation_2(cache))
@test isnothing(DT.get_marked_vertices(cache))
@test isnothing(DT.get_interior_segments_on_hull(cache))
@test isnothing(DT.get_fan_triangles(cache))
@test isnothing(DT.get_incircle_cache(cache))
@test isnothing(DT.get_orient3_cache(cache))
@test isnothing(DT.get_insphere_cache(cache))
@test isnothing(DT.get_surrounding_polygon(cache))
@test typeof(cache) == DT.EmptyTriangulationCache

tri = triangulate(rand(2, 50); weights=rand(50))
cache = DT.get_cache(tri)
@test get_weights(DT.get_triangulation(DT.get_cache(tri))) === get_weights(tri)
DT.unconstrained_triangulation!(DT.get_triangulation(cache))
@test !isempty(DT.get_triangles(DT.get_triangulation(cache)))
@test !isempty(DT.get_adjacent(DT.get_adjacent(DT.get_triangulation(cache))))
@test !isempty(DT.get_adjacent2vertex(DT.get_adjacent2vertex(DT.get_triangulation(cache))))
@test !isempty(DT.get_vertices(DT.get_graph(DT.get_triangulation(cache))))
@test !isempty(DT.get_edges(DT.get_graph(DT.get_triangulation(cache))))
@test !isempty(DT.get_neighbours(DT.get_graph(DT.get_triangulation(cache))))
@test !isempty(DT.get_vertices(DT.get_convex_hull(DT.get_triangulation(cache))))
@test get_weights(DT.get_triangulation_2(DT.get_cache(tri))) === get_weights(tri)
DT.unconstrained_triangulation!(DT.get_triangulation_2(cache))
@test !isempty(DT.get_triangles(DT.get_triangulation_2(cache)))
@test !isempty(DT.get_adjacent(DT.get_adjacent(DT.get_triangulation_2(cache))))
@test !isempty(DT.get_adjacent2vertex(DT.get_adjacent2vertex(DT.get_triangulation_2(cache))))
@test !isempty(DT.get_vertices(DT.get_graph(DT.get_triangulation_2(cache))))
@test !isempty(DT.get_edges(DT.get_graph(DT.get_triangulation_2(cache))))
@test !isempty(DT.get_neighbours(DT.get_graph(DT.get_triangulation_2(cache))))
@test !isempty(DT.get_vertices(DT.get_convex_hull(DT.get_triangulation_2(cache))))
buf = IOBuffer()
show(buf, MIME"text/plain"(), cache)
str = take!(buf)
@test String(str) == "TriangulationCache with storage."
@test isempty(cache.interior_segments_on_hull)
push!(cache.interior_segments_on_hull, (5, 7))
@test cache.interior_segments_on_hull == DT.get_interior_segments_on_hull(cache) == Set([(5, 7)])
@test !isempty(cache.interior_segments_on_hull)
empty!(cache)
@test isempty(cache.interior_segments_on_hull)
@test sprint(show, MIME"text/plain"(), cache.triangulation.cache) == "TriangulationCache with no storage."
@test isempty(DT.get_triangles(DT.get_triangulation(cache)))
@test isempty(DT.get_adjacent(DT.get_adjacent(DT.get_triangulation(cache))))
@test isempty(DT.get_adjacent2vertex(DT.get_adjacent2vertex(DT.get_triangulation(cache))))
@test isempty(DT.get_vertices(DT.get_graph(DT.get_triangulation(cache))))
@test isempty(DT.get_edges(DT.get_graph(DT.get_triangulation(cache))))
@test isempty(DT.get_neighbours(DT.get_graph(DT.get_triangulation(cache))))
@test isempty(DT.get_vertices(DT.get_convex_hull(DT.get_triangulation(cache))))
@test isempty(DT.get_triangles(DT.get_triangulation_2(cache)))
@test isempty(DT.get_adjacent(DT.get_adjacent(DT.get_triangulation_2(cache))))
@test isempty(DT.get_adjacent2vertex(DT.get_adjacent2vertex(DT.get_triangulation_2(cache))))
@test isempty(DT.get_vertices(DT.get_graph(DT.get_triangulation_2(cache))))
@test isempty(DT.get_edges(DT.get_graph(DT.get_triangulation_2(cache))))
@test isempty(DT.get_neighbours(DT.get_graph(DT.get_triangulation_2(cache))))
@test isempty(DT.get_vertices(DT.get_convex_hull(DT.get_triangulation_2(cache))))
marked_vertices = DT.get_marked_vertices(cache)
@test isempty(marked_vertices)
push!(marked_vertices, 5)
@test marked_vertices == DT.get_marked_vertices(cache) == [5]
@test !isempty(marked_vertices)
empty!(cache)
@test isempty(marked_vertices)
fan_triangles = DT.get_fan_triangles(cache)
@test isempty(fan_triangles)
push!(fan_triangles, (5, 7, 9))
@test fan_triangles == DT.get_fan_triangles(cache) == Set([(5, 7, 9)])
@test !isempty(fan_triangles)
empty!(cache)
@test isempty(fan_triangles)
surrounding_polygon = DT.get_surrounding_polygon(cache)
@test isempty(surrounding_polygon)
push!(surrounding_polygon, 5)
@test surrounding_polygon == DT.get_surrounding_polygon(cache) == [5]
@test !isempty(surrounding_polygon)
empty!(cache)
@test isempty(surrounding_polygon)
@test length.(DT.get_incircle_cache(cache)) == length.(AdaptivePredicates.incircleadapt_cache(Float64))
@test length.(DT.get_orient3_cache(cache)) == length.(AdaptivePredicates.orient3adapt_cache(Float64))
@test length.(DT.get_insphere_cache(cache)) == length.(AdaptivePredicates.insphereexact_cache(Float64))

@testset "copy/deepcopy" begin
    cache = DT.TriangulationCache()
    @test copy(cache) === cache

    tri = triangulate(rand(2, 50))
    cache = DT.get_cache(tri)
    cache2 = copy(cache)
    @inferred copy(cache2)
    @test typeof(cache2) == typeof(cache) && !(cache2 === cache)

    tri = triangulate(rand(2, 50); weights=rand(50))
    cache = DT.get_cache(tri)
    cache2 = copy(cache)
    @inferred copy(cache2)
    @test typeof(cache2) == typeof(cache) && !(cache2 === cache)
end