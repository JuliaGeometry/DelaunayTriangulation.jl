using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
include("../helper_functions.jl")

tri = shewchuk_example_constrained()

fig, ax, sc = triplot(tri)
let vert = each_solid_vertex(tri)
    text!(ax, collect(get_point(tri, vert...)); text=string.(vert))
end
lines!(ax, [get_point(tri, 2, 7)...], color=:blue, linestyle=:dash)
fig

intersecting_triangles = DT.find_all_intersecting_triangles_with_segment(
    (2, 7),
    get_points(tri),
    get_adjacent(tri),
    get_adjacent2vertex(tri),
    get_graph(tri),
    get_boundary_index_ranges(tri),
    get_boundary_map(tri),
    NTuple{3,Int64},
    Vector{NTuple{3,Int64}}
)
@test all(T->DT.is_positively_oriented(DT.triangle_orientation(tri,T)), intersecting_triangles)

