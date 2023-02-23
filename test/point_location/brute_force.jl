using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StatsBase

include("../helper_functions.jl")

tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)
DT.compute_representative_points!(tri)
DT.RepresentativePointList[1].x = 10.0
DT.RepresentativePointList[1].y = 10.0
_pts = tri.points[[12, 11, 10, 9]]
DT.RepresentativePointList[2].x = mean([8.0, 8.0, 4.0, 4.0])
DT.RepresentativePointList[2].y = mean([16.0, 6.0, 6.0, 16.0])
_pts = tri.points[[18, 17, 16, 15, 14, 13]]
DT.RepresentativePointList[3].x = mean([18.0, 18.0, 14.0, 12.0, 14.0, 14.0])
DT.RepresentativePointList[3].y = mean([12.0, 6.0, 2.0, 4.0, 6.0, 10.0])
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
@test DT.compare_triangles(DT.brute_force_search(tri.triangles, c1_i, pts,
                                                 tri.boundary_map),
                           (index_map["g"], index_map["f"], DT.BoundaryIndex))
@test DT.compare_triangles(DT.brute_force_search(tri.triangles, d1_i, pts,
                                                 tri.boundary_map),
                           (index_map["f"], index_map["e"], DT.BoundaryIndex))
@test DT.compare_triangles(DT.brute_force_search(tri.triangles, e1_i, pts,
                                                 tri.boundary_map),
                           (index_map["k"], index_map["ℓ"], DT.BoundaryIndex - 1))
@test DT.compare_triangles(DT.brute_force_search(tri.triangles, f1_i, pts,
                                                 tri.boundary_map),
                           (index_map["ℓ"], index_map["k"], index_map["w"]))
@test DT.compare_triangles(DT.brute_force_search(tri.triangles, g1_i, pts,
                                                 tri.boundary_map),
                           (index_map["b"], index_map["c"], index_map["p"]))
@test DT.compare_triangles(DT.brute_force_search(tri.triangles, h1_i, pts,
                                                 tri.boundary_map),
                           (index_map["m"], index_map["n"], DT.BoundaryIndex - 3))

@inferred DT.brute_force_search(tri.triangles, c1_i, pts, tri.boundary_map)

for T in each_triangle(tri.triangles)
    for i in indices(T)
        if i ≠ DT.BoundaryIndex
            p = get_point(pts, tri.boundary_map, i)
            V = DT.brute_force_search(tri.triangles, p, pts, tri.boundary_map)
            @inferred DT.brute_force_search(tri.triangles, p, pts, tri.boundary_map)
            if !DT.is_ghost_triangle(T)
                @test i ∈ V
            else
                @test i ∈ V || i - 1 ∈ V || i + 1 ∈ V
            end
        end
    end
end
