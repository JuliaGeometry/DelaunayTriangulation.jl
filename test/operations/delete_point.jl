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
fig, ax, sc = triplot(_tri)
let vert = each_solid_vertex(_tri)
    text!(ax, collect(get_point(_tri, vert...)); text=string.(vert))
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

delete_point!(tri, 9; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17, 2, 9), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 7; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17, 2, 9, 7), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 10; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17, 2, 9, 7, 10), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 8; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17, 2, 9, 7, 10, 8), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 11; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17, 2, 9, 7, 10, 8, 11), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 18; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17, 2, 9, 7, 10, 8, 11, 18), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 5; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17, 2, 9, 7, 10, 8, 11, 18, 5), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)

delete_point!(tri, 12; rng)
validate_triangulation(tri)
_tri = triangulate(get_points(tri); skip_points=(15, 16, 14, 13, 17, 2, 9, 7, 10, 8, 11, 18, 5, 12), delete_ghosts=false)
clear_empty_features!(tri)
clear_empty_features!(_tri)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
@test get_adjacent(tri) == get_adjacent(_tri)
@test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
@test get_graph(tri) == get_graph(_tri)
convex_hull!(tri)
@test get_convex_hull(tri) == get_convex_hull(_tri)


point = 5
nn = DT.num_neighbours(tri, point)
@show nn
is_bnd, boundary_index = DT.is_boundary_node(tri, point)
check_orientation = Val(is_bnd)
if is_bnd
    left_bnd = DT.get_left_boundary_node(tri, point, boundary_index)
    right_bnd = DT.get_right_boundary_node(tri, point, boundary_index)
else
    I = integer_type(tri)
    left_bnd = I(DT.DefaultAdjacentValue)
    right_bnd = I(DT.DefaultAdjacentValue)
end
S = DT.get_surrounding_polygon(tri, point)
neighbouring_edges = get_adjacent2vertex(tri, point)
for uv in each_edge(neighbouring_edges)
    u = initial(uv)
    v = terminal(uv)
    delete_triangle!(tri, point, u, v; protect_boundary=true)
end
next, prev, k, S, shuffled_indices = DT.prepare_convex_triangulation_vectors(S)
DT.triangulate_convex!(tri, S; rng, update_ghost_edges=false, check_orientation)
DT.delete_adjacent2vertex!(tri, point)
DT.delete_vertex!(tri, point)
DT.fix_edges_after_deletion!(tri, S)
DT.fix_ghost_edges_after_deletion!(tri, boundary_index, check_orientation, left_bnd, right_bnd)


S 
for i in firstindex(S):(lastindex(S)-1)
    u = S[i]
    v = S[i+1]
    @show get_adjacent(tri,u, v)
end

## CONSIDER INTERSECTION OF THE BOUNDARY NEIGHBOURS WITH S, FIND GHOST TRIANGLES THAT EXIST AND THEN DETERMINE IF THEIR ADJACENCIES ARE WRONG???