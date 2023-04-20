using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using ColorSchemes

A = (-1.0, 7.0)
B = (4.0, 4.0)
C = (-2.0, -1.0)
D = (-1.0, 3.0)
E = (3.0, -1.0)
F = (1.0, 4.0)
G = (-3.0, 5.0)
pts = [A, B, C, D, E, F, G]
tri = triangulate(pts; delete_ghosts=false, randomise=false)
add_edge!(tri, 6, 5)

fig, ax, sc = triplot(tri)
let vert = each_solid_vertex(tri)
    text!(ax, collect(get_point(tri, vert...)); text=string.(vert))
end
fig

vorn = voronoi(tri)
xmin, xmax, ymin, ymax = DT.polygon_bounds(vorn, 0.5)
cmap = Makie.cgrad(:jet)
colors = get_cell_colors(vorn, cmap)
fig, ax, sc = voronoiplot(vorn, cell_color=colors)
triplot!(ax, tri)
xlims!(ax, xmin, xmax)
ylims!(ax, ymin, ymax)
fig



push!(vorn.unbounded_cells, 1, 7, 3, 5, 2)
coords = DT.get_cell_coordinates(vorn, 4)
@test collect.(coords) ≈ collect.([
    (1.1666666666666665, 1.1666666666666667)
    (-0.75, 5.0)
    (-1.0, 5.0)
    (-4.3, 1.7000000000000002)
    (0.5, 0.5)
    (1.1666666666666665, 1.1666666666666667)
])
coords = DT.get_cell_coordinates(vorn, 6)
@test collect.(coords) ≈ collect.([
    (2.5, 1.7000000000000002)
    (2.5, 7.166666666666667)
    (-0.75, 5.0)
    (1.1666666666666665, 1.1666666666666667)
    (2.5, 1.7000000000000002)
])
coords = DT.get_cell_coordinates(vorn, 2, bbox, bbox_order)


cmap = Makie.cgrad(:hsv)
colors = get_cell_colors(vorn, cmap)
fig, ax, sc = voronoiplot(vorn, cell_color=colors, colormap=cmap)
xlims!(ax, -8, 8)
ylims!(ax, -8, 8)
fig

tri = triangulate(rand(2, 50))
vorn = voronoi(tri)
cmap = Makie.cgrad(:matter)
colors = get_cell_colors(vorn, cmap)
fig, ax, sc = voronoiplot(vorn, cell_color=colors, colormap=cmap, markersize=6)
xlims!(ax,-1,2)
ylims!(ax,-1,2)
fig
