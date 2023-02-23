using CairoMakie
using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

include("./helper_functions.jl")

x, y = complicated_geometry()
tri = generate_mesh(x, y, 0.1; convert_result=true, add_ghost_triangles=true)

fig1 = Figure()
ax = Axis(fig1[1, 1])
triplot!(ax, tri)
save("$save_path/test_figure_1.png", fig1)

fig2 = Figure()
ax = Axis(fig2[1, 1])
triplot!(ax, tri;
         show_ghost_edges=true,
         show_all_points=true,
         markersize=6)
xlims!(ax, -0.4, 2.4)
ylims!(ax, -0.4, 6.4)
save("$save_path/test_figure_2.png", fig2)

fig3 = Figure()
ax = Axis(fig3[1, 1])
triplot!(ax, tri;
         show_ghost_edges=true,
         show_all_points=true,
         markersize=6,
         point_color=:blue,
         strokecolor=:red,
         ghost_edge_color=:darkgreen)
xlims!(ax, -0.4, 2.4)
ylims!(ax, -0.4, 6.4)
save("$save_path/test_figure_3.png", fig3)

fig4 = Figure()
ax = Axis(fig4[1, 1])
triplot!(ax, tri;
         show_ghost_edges=true,
         show_all_points=true,
         markersize=6,
         point_color=:blue,
         strokecolor=:red,
         ghost_edge_color=:darkgreen)
xlims!(ax, -0.4, 2.4)
ylims!(ax, -0.4, 6.4)
save("$save_path/test_figure_4.png", fig4)

fig5 = Figure()
ax = Axis(fig5[1, 1])
triplot!(ax, tri;
         show_ghost_edges=true,
         show_all_points=true,
         markersize=6,
         point_color=:blue,
         strokecolor=:red,
         ghost_edge_color=:darkgreen,
         triangle_color=(:blue, 0.3),
         recompute_centroid=false,
         strokewidth=0.0)
xlims!(ax, -0.4, 2.4)
ylims!(ax, -0.4, 6.4)
save("$save_path/test_figure_5.png", fig5)

pts = get_points(tri)
T = get_triangles(tri)
bn = get_boundary_nodes(tri)
fig6 = Figure()
ax = Axis(fig6[1, 1])
triplot!(ax, pts, T, bn, get_convex_hull(tri);
         show_ghost_edges=true,
         show_all_points=true,
         markersize=6,
         point_color=:blue,
         strokecolor=:red,
         ghost_edge_color=:darkgreen,
         triangle_color=(:blue, 0.3),
         recompute_centroid=false,
         strokewidth=0.0)
xlims!(ax, -0.4, 2.4)
ylims!(ax, -0.4, 6.4)
save("$save_path/test_figure_6.png", fig6)

fig7 = Figure()
ax = Axis(fig7[1, 1])
triplot!(ax, pts, T, bn, get_convex_hull(tri);
         show_ghost_edges=true,
         show_all_points=true,
         markersize=6,
         point_color=:blue,
         strokecolor=:red,
         ghost_edge_color=:darkgreen,
         triangle_color=(:blue, 0.3),
         recompute_centroid=false,
         plot_convex_hull=true,
         convex_hull_linewidth=7,
         convex_hull_linestyle=:dash,
         convex_hull_color=:red,
         strokewidth=0.0)
xlims!(ax, -0.4, 2.4)
ylims!(ax, -0.4, 6.4)
save("$save_path/test_figure_7.png", fig7)
