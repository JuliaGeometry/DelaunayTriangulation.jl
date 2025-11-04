# # Triangulation Operations 
# ## Edge Splitting 

# Sometimes you have a point `r` that you know is on, or is very 
# close to being on, an edge `(i, j)`. In this tutorial, we show how 
# the [`split_edge!`](@ref) function can be used for putting a point on this 
# edge. First, let us consider the following triangulation. 
using DelaunayTriangulation
using CairoMakie
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

points = [
    (0.0, 0.0), (0.0, 4.0), (2.0, 3.0), (-2.0, 3.0),
    (-2.0, 7.0), (3.0, 6.0), (2.0, -2.0), (-4.0, 1.0),
    (1.0, 5.0),
]
p = (0.0, 3.0)
tri = triangulate(points)
fig, ax, sc = triplot(tri)
scatter!(ax, [p], markersize = 14)
fig
@test_reference joinpath(fig_path, "triangulation_operations_17.png") fig #src

# We want to add the blue point onto the edge shown, which is `(1, 2)`. To do this, 
# we can use the function `split_edge!`.
push!(points, p)
r = length(points)
i, j = 1, 2
split_edge!(tri, i, j, r)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "triangulation_operations_18.png") fig #src

# Notice that this has only split the edge in one direction. This is because the 
# edges in this case are treated as being oriented. To split the edge in the other
# direction, we simply swap the indices.
split_edge!(tri, j, i, r)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "triangulation_operations_19.png") fig #src

# If you also want to restore the Delaunay property of the triangulation 
# following this splitting, you need to use [`legalise_edge!`](@ref). In this 
# example, though, there are no illegal edges. If there were, we would use
k = get_adjacent(tri, i, r) # get_adjacent(tri, i, j) before split_edge!(tri, i, j)
legalise_edge!(tri, j, k, r)
legalise_edge!(tri, k, i, r)
k = get_adjacent(tri, j, r) # get_adjacent(tri, j, i) before split_edge!(tri, j, i)
legalise_edge!(tri, i, k, r)
legalise_edge!(tri, k, j, r)

# These steps, in particular the steps of splitting both sides of the edge and then 
# legalising, are also implemented in [`DelaunayTriangulation.complete_split_edge_and_legalise!`](@ref).
