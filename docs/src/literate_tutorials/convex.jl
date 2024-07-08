# # Triangulating Convex Polygons 
# 

# In this tutorial, we show how we can triangulate convex 
# polygons. The function [`triangulate_convex`](@ref) is used for this.
# Let us start with a simple example.
using DelaunayTriangulation
using CairoMakie
using ReferenceTests #src
using Test #src
using Preferences #src
fig_path = joinpath(@__DIR__, "../figures") #src

points = [
    (10.0, 12.0), (7.0, 11.0), (8.0, 6.0),
    (10.0, 3.0), (14.0, 5.0), (15.0, 10.0),
    (13.0, 12.0)
]
S = 1:7
tri = triangulate_convex(points, 1:7)

#-
fig, ax, sc = triplot(tri)
fig 
load_preference(DelaunayTriangulation, "USE_EXACTPREDICATES", true) && @test_reference joinpath(fig_path, "convex_ex_1.png") fig #src

# This `tri` is our triangulation of the convex polygon. 
# The first input is the set of points, and `S` defines 
# the vertices to take from these points and their order, 
# which must be provided in a counter-clockwise order.
# Note that points does not have to contain only the points 
# of the polygon, since `S` will be used to define the points 
# needed.

# Let us give a larger example. For simplicity, we triangulate a  
# a discretised circle. 
θ = LinRange(0, 2π, 5000) |> collect 
pop!(θ)
x = cos.(θ)
y = sin.(θ)
points = tuple.(x, y)
S = 1:4999 # can also be [1:4999; 1], if you want the array to be circular 
tri = triangulate_convex(points, S)

#-
fig, ax, sc = triplot(tri)
fig 
load_preference(DelaunayTriangulation, "USE_EXACTPREDICATES", true) && @test_reference joinpath(fig_path, "convex_ex_2.png") fig #src

# Here is a comparison of the time it takes to triangulate this 
# using `triangulate_convex` or `triangulate`.
using BenchmarkTools 
@benchmark triangulate_convex($points, $S) 

#-
@benchmark triangulate($points) 

# For the smaller example that we started with above, `triangulate_convex` is also 
# faster, although not by much (≈15.10 μs versus ≈10.7 μs).