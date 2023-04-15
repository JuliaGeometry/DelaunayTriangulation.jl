using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using StableRNGs

include("../helper_functions.jl")

rng = StableRNG(8888)

@testset "Adding a point" begin
    pts = rand(rng, 2, 500)
    tri = DT.initialise_bowyer_watson(pts)
    for i in setdiff(each_point_index(pts), get_vertices(tri))
        add_point!(tri, i)
        @test validate_triangulation(tri)
    end
    convex_hull!(tri; reconstruct=false)
    delete_ghost_triangles!(tri)
    clear_empty_features!(tri)
    _tri = triangulate(pts; rng)
    clear_empty_features!(_tri)
    @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
    @test get_adjacent(tri) == get_adjacent(_tri)
    @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
    @test get_graph(tri) == get_graph(_tri)
    @test get_convex_hull(tri) == get_convex_hull(_tri)
end

@testset "Returning the triangle" begin
    tri = triangulate([rand(2) for _ in 1:50])
    V = add_point!(tri, 0.5, 0.3)
    @test DT.is_inside(DT.point_position_relative_to_triangle(tri, V..., (0.5, 0.3)))
end