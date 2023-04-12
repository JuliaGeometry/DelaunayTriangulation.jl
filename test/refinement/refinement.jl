using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using StableRNGs

include("../helper_functions.jl")

Random.seed!(1919191)
p1 = (0.0, 0.0)
p2 = (1.0, 0.0)
p3 = (0.0, 1.0)
p4 = (1.0, 1.0)
pts = [p1, p2, p3, p4]
tri = triangulate(pts; delete_ghosts=false)
refine!(tri; max_area=0.0001)
stats = statistics(tri)
@test DT.get_smallest_angle(stats) ≥ deg2rad(30.0)
@test DT.get_largest_area(stats) ≤ 0.0001
@test !DT.is_constrained(tri)
@test DT.convex_hull(tri).indices == DT.convex_hull(tri.points).indices
@test validate_triangulation(tri)
validate_statistics(tri)

Random.seed!(1919191)
p1 = (0.0, 0.0)
p2 = (1.0, 0.0)
p3 = (0.0, 1.0)
p4 = (1.0, 1.0)
pts = [p1, p2, p3, p4]
tri = triangulate(pts; delete_ghosts=false)
add_edge!(tri, 1, 4)
refine!(tri; max_area=0.01)
