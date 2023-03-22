using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
import SimpleGraphs: relabel, UndirectedGraph
using DataStructures
using StableRNGs

include("../test_setup.jl")
save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

include("../helper_functions.jl")

tri = example_with_special_corners()
fig, ax, sc = triplot(tri)
let vert = each_solid_vertex(tri)
    text!(ax, collect(get_point(tri, vert...)); text=string.(vert))
end

rng = StableRNG(292929)

point = 16
delete_point!(tri, 16; rng)
validate_triangulation(tri)
_tri = Triangulation(get_points(tri))
true_T = [
    10 9 11
    11 9 15
    11 15 12
    10 11 12
    10 12 13
    12 18 13
    13 18 5
    18 4 5
    18 17 4
    14 17 18
    12 14 18
    12 15 14
    14 15 17
    9 8 15
    15 8 7
    15 7 17
    17 6 4
    17 7 6
    6 7 3
    4 6 3
    4 3 1
    1 3 2
    5 4 1
]
for T in eachrow(true_T)
    add_triangle!(_tri, T; update_ghost_edges=true)
end
clear_empty_features!(_tri)
clear_empty_features!(tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
_tri = triangulate(get_points(tri); skip_points=16, delete_ghosts=false)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 15; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 14; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 13; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 17; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 2; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17, 2), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)




################
#####
#########
########
#
####
####
######### REMINDER: Maybe for the case of deleting a boundary point, just connect the neighbouring points 
######### on the boundary since those will be connected, and then legalise the two other edges?


point = 13
S = DT.get_surrounding_polygon(tri, point)
neighbouring_edges = DT.get_adjacent2vertex(tri, point)
for uv in each_edge(neighbouring_edges)
    @show DT.get_adjacent(tri, 1, 5)
    u = initial(uv)
    v = terminal(uv)
    @show (u, v, point)
    delete_triangle!(tri, point, u, v; protect_boundary=true)
end
next, prev, k, S, shuffled_indices = DT.prepare_convex_triangulation_vectors(S)
DT.reset_convex_triangulation_vectors!(next, prev, shuffled_indices, k, rng)
DT.delete_vertices_in_random_order!(next, prev, shuffled_indices, tri, S, k, rng)
u, v, w = DT.index_shuffled_linked_list(S, next, prev, shuffled_indices, 1)
add_triangle!(tri, u, v, w; protect_boundary=true)
@show DT.get_adjacent(tri, 1, 5)
set_S = Set(S)
for i in 4:k
    @show DT.get_adjacent(tri, 1, 5), i
    u, v, w = DT.index_shuffled_linked_list(S, next, prev, shuffled_indices, i)
    DT.add_point_convex_triangulation!(tri, u, v, w, set_S, false, true)
end
DT.delete_adjacent2vertex!(tri, point)
DT.delete_vertex!(tri, point)


DT.fix_edges_after_deletion!(tri,S,-1,Val(true))



push!(S, S[begin])
i = 2
u = S[i]
v = S[i+1]
get_adjacent(tri, u, v)
get_adjacent(tri, v, u) - get_adjacent(_tri, v, u)






push!(S, S[begin])
i = 3
u = S[i]
v = S[i+1]
w = get_adjacent(tri, v, u)



i = 5
u, v, w = DT.index_shuffled_linked_list(S, next, prev, shuffled_indices, i)
x = get_adjacent(tri, w, v)











delete_point!(tri, 16; rng)
validate_triangulation(tri)
_tri = Triangulation(get_points(tri))
true_T = [
    10 9 11
    11 9 15
    11 15 12
    10 11 12
    10 12 13
    12 18 13
    13 18 5
    18 4 5
    18 17 4
    14 17 18
    12 14 18
    12 15 14
    14 15 17
    9 8 15
    15 8 7
    15 7 17
    17 6 4
    17 7 6
    6 7 3
    4 6 3
    4 3 1
    1 3 2
    5 4 1
]
for T in eachrow(true_T)
    add_triangle!(_tri, T; update_ghost_edges=true)
end
clear_empty_features!(_tri)
clear_empty_features!(tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
_tri = triangulate(get_points(tri); skip_points=16, delete_ghosts=false)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 15)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16), delete_ghosts=false)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)
