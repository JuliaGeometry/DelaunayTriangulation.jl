# # Convex Hulls 
# 

# In this tutorial, we show how we can compute convex hulls. 
# The main method for this is to compute them directly from the 
# triangulation, noting that for an unconstrained triangulation the boundary 
# is the convex hull. We also provide a method for computing the convex hull 
# from the point set directly, using the [monotone chain algorithm](https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain).

# Let us first demonstrate how we compute convex hulls from a triangulation. 
# This is only possible from triangulations without a constrained boundary 
using DelaunayTriangulation
using CairoMakie
using StableRNGs
using Test #src
using ReferenceTests #src
fig_path = joinpath(@__DIR__, "../figures") #src

# We construct a triangulation of a point set.
rng = StableRNG(123)
points = randn(rng, 2, 250)
tri = triangulate(points; rng)

# To get the convex hull, just use [`get_convex_hull`](@ref):
get_convex_hull(tri)
@test get_convex_hull(tri) == tri.convex_hull #src

# You can also obtain the vertices directly from 
# [`get_convex_hull_vertices`](@ref):
get_convex_hull_vertices(tri)
@test get_convex_hull_vertices(tri) == tri.convex_hull.vertices #src

# The vertices refer to indices of points in `points`, and they are 
# given in counter-clockwise order. To obtain the convex hull directly 
# from the points without constructing the triangulation, use [`convex_hull`](@ref):
ch = convex_hull(points)
ch_points = [get_point(tri, i) for i in DelaunayTriangulation.get_vertices(ch)]
fig, ax, sc = lines(ch_points, color=:red, linewidth=4)
scatter!(ax, points)
fig
@test_reference joinpath(fig_path, "convex_hull_ex_1.png") fig #src
