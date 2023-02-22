using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random

include("../helper_functions.jl")

pts = rand(2, 500)
tri = DT.initialise_bowyer_watson(pts)
for i in setdiff(each_point_index(pts), get_vertices(tri))
    add_point!(tri, i)
    validate_triangulation(tri)
end
convex_hull!(tri; reconstruct=false)
delete_ghost_triangles!(tri)
clear_empty_features!(tri)
_tri = triangulate(pts)
@test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri)) &&
      get_adjacent(tri) == get_adjacent(_tri) &&
      get_adjacent2vertex(tri) == get_adjacent2vertex(_tri) &&
      get_graph(tri) == get_graph(_tri) &&
      get_convex_hull(tri) == get_convex_hull(_tri)
