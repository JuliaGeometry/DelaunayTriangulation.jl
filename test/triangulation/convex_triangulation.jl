using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using Test
using DataStructures
using CairoMakie

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

S = [8, 17, 21, 28, 35, 37]
k, next, prev, seen, shuffled_indices = DT.prepare_convex_triangulation_vectors(S)
@test k == 6
@test next == zeros(Int64, 6)
@test prev == zeros(Int64, 6)
@test seen == Set{Int64}()
@test shuffled_indices == [1, 2, 3, 4, 5, 6]
push!(seen, 2)
DT.reset_convex_triangulation_vectors!(next, prev, seen, k)
@test seen == Set{Int64}()
@test next == [2, 3, 4, 5, 6, 1]
@test prev == [6, 1, 2, 3, 4, 5]

a = [4.0, 8.0]
b = [7.0, 7.0]
c = [8.0, 3.0]
d = [5.0, 1.0]
e = [2.0, 4.0]
f = [1.0, 7.0]
pts = [f,e,d,c,b,a]  # counterclockwise
S = [1, 2, 3, 4, 5, 6]
tri = triangulate_convex(pts, S)
triplot(tri;recompute_centers=false)



pts = rand(2, 50)
tri = triangulate(pts)
ch = get_convex_hull_indices(tri)
_tri = triangulate_convex(pts, ch)