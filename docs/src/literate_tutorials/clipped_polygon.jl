# # Clipped Voronoi Tessellations 
# ## Clipping to a Generic Convex Polygon
#

# In this tutorial we show how to clip a Voronoi tessellation to 
# more generic convex polygons (non-convex polygons are not currently supported) than just a convex hull or a rectangle. This is
# done by using the `clip_polygon` keyword argument in `voronoi`.

# We start by clipping the tessellation to a rectangle, showing an alternative
# to the [previous tutorial](clipped_rectangle.md). To start, we load in the packages 
# we need and generate some data.
using DelaunayTriangulation
using DelaunayTriangulation: EllipticalArc
using CairoMakie
using StableRNGs
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

rng = StableRNG(123)
points = randn(rng, 2, 50)
tri = triangulate(points; rng)
vorn = voronoi(tri)

# To define the polygon, we define the points and vertices just as we would, for example,
# the boundary of a triangulation.
xmin, xmax, ymin, ymax = -1/2, 1/2, -1.0, 1.0
clip_points = ((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax))
clip_vertices = (1, 2, 3, 4, 1)
clip_polygon = (clip_points, clip_vertices)

# Now we simply pass the polygon into `voronoi`.
clipped_vorn = voronoi(tri, clip = true, clip_polygon = clip_polygon)

# Now let's look at the results.
fig = Figure()
ax1 = Axis(fig[1, 1], title = "Unclipped", width = 600, height = 400)
ax2 = Axis(fig[1, 2], title = "Clipped", width = 600, height = 400)
voronoiplot!(ax1, vorn, show_generators = false, color = :white, strokewidth = 4)
xlims!(ax1, -2, 2)
ylims!(ax1, -2, 2)
lines!(ax1, [clip_points..., clip_points[begin]], color = :black, linewidth = 4, linestyle = :dash)
voronoiplot!(ax2, clipped_vorn, show_generators = false, color = :white, strokewidth = 4)
xlims!(ax2, -2, 2)
ylims!(ax2, -2, 2)
resize_to_layout!(fig)
fig
@test_reference joinpath(fig_path, "generic_clipped_voronoi_ex_1.png") fig #src

# We can clip to any convex polygon that we want to. For example, below we clip to an elliptical boundary.
rng = StableRNG(123333)
points = randn(rng, 2, 50)
tri = triangulate(points; rng)
vorn = voronoi(tri)
ellip = EllipticalArc((1/2, 0.0), (1/2, 0.0), (0.0, 0.0), 1/2, 1.0, 0.0)
t = LinRange(0, 1, 50)
clip_points = ellip.(t)
clip_vertices = [1:(length(clip_points)-1); 1]
clip_polygon = (clip_points, clip_vertices)
clipped_vorn = voronoi(tri, clip = true, clip_polygon = clip_polygon)
fig = Figure()
ax1 = Axis(fig[1, 1], title = "Unclipped", width = 600, height = 400)
ax2 = Axis(fig[1, 2], title = "Clipped", width = 600, height = 400)
voronoiplot!(ax1, vorn, show_generators = false, color = :white, strokewidth = 4)
xlims!(ax1, -2, 2)
ylims!(ax1, -2, 2)
lines!(ax1, [clip_points..., clip_points[begin]], color = :black, linewidth = 4, linestyle = :dash)
voronoiplot!(ax2, clipped_vorn, show_generators = false, color = :white, strokewidth = 4)
xlims!(ax2, -2, 2)
ylims!(ax2, -2, 2)
resize_to_layout!(fig)
fig
@test_reference joinpath(fig_path, "generic_clipped_voronoi_ex_2.png") fig #src