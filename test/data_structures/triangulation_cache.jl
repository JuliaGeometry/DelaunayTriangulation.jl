using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StructEquality
using Setfield
using ..DelaunayTriangulation: add_weight!, get_weight, get_weights


@struct_equal DT.TriangulationCache


tri = triangulate(rand(2, 50); weights = DT.ZeroWeight())
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
