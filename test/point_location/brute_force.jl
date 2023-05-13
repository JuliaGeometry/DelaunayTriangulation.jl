using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StatsBase

include("../helper_functions.jl")

global tri, label_map, index_map = simple_geometry()

@testset "Finding points in ghost triangles" begin
    add_ghost_triangles!(tri)
    DT.compute_representative_points!(tri)
    rep = DT.get_representative_point_list(tri)
    rep[1].x = 10.0
    rep[1].y = 10.0
    _pts = tri.points[[12, 11, 10, 9]]
    rep[2].x = mean([8.0, 8.0, 4.0, 4.0])
    rep[2].y = mean([16.0, 6.0, 6.0, 16.0])
    _pts = tri.points[[18, 17, 16, 15, 14, 13]]
    rep[3].x = mean([18.0, 18.0, 14.0, 12.0, 14.0, 14.0])
    rep[3].y = mean([12.0, 6.0, 2.0, 4.0, 6.0, 10.0])
    pts = get_points(tri)
    c1 = (5.3197, 26.51)
    d1 = (15.0, 25.0)
    e1 = (7.579498, 13.29)
    f1 = (9.16, 12.49)
    g1 = (14.358, 0.86)
    h1 = (14.5282416299238, 8.0932516857679)
    push!(pts, c1, d1, e1, f1, g1, h1)
    h1_i = length(pts)
    g1_i = length(pts) - 1
    f1_i = length(pts) - 2
    e1_i = length(pts) - 3
    d1_i = length(pts) - 4
    c1_i = length(pts) - 5
    @test DT.compare_triangles(DT.brute_force_search(tri, c1_i),
        (index_map["g"], index_map["f"], DT.BoundaryIndex))
    @test DT.compare_triangles(DT.brute_force_search(tri, d1_i),
        (index_map["f"], index_map["e"], DT.BoundaryIndex))
    @test DT.compare_triangles(DT.brute_force_search(tri, e1_i),
        (index_map["k"], index_map["ℓ"], DT.BoundaryIndex - 1))
    @test DT.compare_triangles(DT.brute_force_search(tri, f1_i),
        (index_map["ℓ"], index_map["k"], index_map["w"]))
    @test DT.compare_triangles(DT.brute_force_search(tri, g1_i),
        (index_map["b"], index_map["c"], index_map["p"]))
    @test DT.compare_triangles(DT.brute_force_search(tri, h1_i),
        (index_map["m"], index_map["n"], DT.BoundaryIndex - 3))
    @inferred DT.brute_force_search(tri, c1_i)
end

@testset "Finding points in each triangle" begin
    pts = get_points(tri)
    rep = DT.get_representative_point_list(tri)
    c1 = (5.3197, 26.51)
    d1 = (15.0, 25.0)
    e1 = (7.579498, 13.29)
    f1 = (9.16, 12.49)
    g1 = (14.358, 0.86)
    h1 = (14.5282416299238, 8.0932516857679)
    push!(pts, c1, d1, e1, f1, g1, h1)
    h1_i = length(pts)
    g1_i = length(pts) - 1
    f1_i = length(pts) - 2
    e1_i = length(pts) - 3
    d1_i = length(pts) - 4
    c1_i = length(pts) - 5
    for T in each_triangle(tri.triangles)
        for i in indices(T)
            if i ≠ DT.BoundaryIndex
                p = get_point(pts, rep, tri.boundary_map, i)
                V = DT.brute_force_search(tri, p)
                @inferred DT.brute_force_search(tri, p)
                if !DT.is_ghost_triangle(T)
                    @test i ∈ V
                else
                    @test i ∈ V || i - 1 ∈ V || i + 1 ∈ V
                end
            end
        end
    end
end