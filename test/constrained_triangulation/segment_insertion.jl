using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using StableRNGs
include("../helper_functions.jl")

#=
We test constrained Delaunay triangulations by 
comparing results to MATLAB. For example, 
for fixed_shewchuk_example_constrained,
we use the MATLAB code for generating 
a triangulation with a constrained edge 
(2, 7):

    P = [
        0.0 0.0;
        0.0 1.0;
        0.0 2.5;
        2.0 0.0;
        6.0 0.0;
        8.0 0.0;
        8.0 0.5;
        7.5 1.0;
        4.0 1.0;
        4.0 2.5;
        8.0 2.5
    ];
    C = [2 7];
    tri = delaunayTriangulation(P, C);
    tri.ConnectivityList
    ans =

     1     4     2
     8    11    10
     8     7    11
     2     9     3
     3     9    10
     9     8    10
     9     7     8
     2     7     9
     7     5     6
     7     2     5
     5     2     4
=#


@testset "vertex_is_closer_than_neighbours" begin
    tri = fixed_shewchuk_example_constrained()
    @test !DT.vertex_is_closer_than_neighbours(tri, 2, 7, 10, 8, 9)
    @test !DT.vertex_is_closer_than_neighbours(tri, 2, 7, 10, 9, 8)
    @test DT.vertex_is_closer_than_neighbours(tri, 2, 7, 9, 8, 3)
    @test DT.vertex_is_closer_than_neighbours(tri, 2, 7, 9, 3, 8)
    @test !DT.vertex_is_closer_than_neighbours(tri, 2, 7, 8, 9, 11)
    @test DT.vertex_is_closer_than_neighbours(tri, 2, 7, 3, 10, 11)
    @test !DT.vertex_is_closer_than_neighbours(tri, 7, 2, 1, 4, 5)
    @test !DT.vertex_is_closer_than_neighbours(tri, 7, 2, 4, 1, 5)
    @test DT.vertex_is_closer_than_neighbours(tri, 7, 2, 5, 1, 4)
    @test DT.vertex_is_closer_than_neighbours(tri, 7, 2, 6, 5, 4)
end


rng = StableRNG(191919)
tri = fixed_shewchuk_example_constrained()
e = (2, 7)
T, C, L, R = DT.locate_intersecting_triangles(tri, e)
DT.delete_intersected_triangles!(tri, T)
points = get_points(tri)
_tri = Triangulation(points)
V = L
v = V[begin]
u = V[end]
prev, next, shuffled_indices = DT.prepare_vertex_linked_list(V)
@test prev[2:end-1] == [1, 2, 3, 4, 5]
@test V[prev[2:end-1]] == [7, 8, 10, 9, 10]
@test next[2:end-1] == [3, 4, 5, 6, 7]
@test V[next[2:end-1]] == [10, 9, 10, 3, 2]
@test shuffled_indices[2:end-1] == [2, 3, 4, 5, 6]
DT.delete_polygon_vertices_in_random_order!(_tri, V, shuffled_indices, prev, next, u, v, rng)
add_triangle!(tri, V[begin], V[shuffled_indices[2]], V[end]; protect_boundary=true, update_ghost_edges=false)

tri = fixed_shewchuk_example_constrained()


fig, ax, sc = triplot(tri)
let vert = each_solid_vertex(tri)
    text!(ax, collect(get_point(tri, vert...)); text=string.(vert))
end
lines!(ax, [get_point(tri, 2, 7)...], color=:blue, linestyle=:dash)
fig

e = (2, 7)
T, C, L, R = DT.locate_intersecting_triangles(tri, e)
DT.delete_intersected_triangles!(tri, T)
_tri_1 = DT.triangulate_cavity_cdt(points, L)
_tri_2 = DT.triangulate_cavity_cdt(points, R)

fig, ax, sc = triplot(tri)
triplot!(ax, _tri_1, strokecolor=:blue)
triplot!(ax, _tri_2, strokecolor=:red)
fig

DT.add_new_triangles!(tri, _tri_1, _tri_2)

points = get_points(tri)
V = L
_tri = Triangulation(points)
v = V[begin]
u = V[end]
prev, next, shuffled_indices = DT.prepare_vertex_linked_list(V)
DT.delete_polygon_vertices_in_random_order!(_tri, V, shuffled_indices, prev, next, u, v)

kp = [1 4 2
    8 11 10
    8 7 11
    2 9 3
    3 9 10
    9 8 10
    9 7 8
    2 7 9
    7 5 6
    7 2 5
    5 2 4]'
kp = collect(eachcol(kp))

