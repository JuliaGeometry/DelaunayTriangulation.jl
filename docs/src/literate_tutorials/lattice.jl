# # Triangulating Rectangular Regions 
# 

# In this tutorial, we show how you can easily 
# triangulate rectangular regions of the form 
# $[a, b] \times [c, d]$. Rather than using `triangulate`,
# you can use `triangulate_rectangle` for this purpose. To start, 
# we give a simple example 
using DelaunayTriangulation
using CairoMakie
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

a, b, c, d = 0.0, 2.0, 0.0, 10.0
nx, ny = 10, 25
tri = triangulate_rectangle(a, b, c, d, nx, ny)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "rectangle_ex_1.png") fig #src

# This can be much faster than if we just construct the points in the 
# lattice manually and `triangulate` those. Here's a comparison of the times. 
using BenchmarkTools
points = get_points(tri)
@benchmark triangulate($points; randomise=$false) # randomise=false because points are already in lattice order, i.e. spatially sorted

#-
@benchmark triangulate_rectangle($a, $b, $c, $d, $nx, $ny)

# This difference would be more pronounced for larger `nx, ny`.

# Note that the output of `triangulate_rectangle` treats the boundary 
# as a constrained boundary: 
get_boundary_nodes(tri)
@test !isempty(get_boundary_nodes(tri)) #src
@test DelaunayTriangulation.has_boundary_nodes(tri) #src

# This boundary is split into four separate segments, one for each 
# side of the rectangle. If you would prefer to keep the boundary as one 
# contiguous segment, use `single_boundary=true`. Moreover, note that 
# this `tri` has ghost triangles:
tri
@test DelaunayTriangulation.has_ghost_triangles(tri) #src

# You can opt into not having these by using `add_ghost_triangles=false`.
tri = triangulate_rectangle(a, b, c, d, nx, ny; single_boundary=true, add_ghost_triangles=false)
tri

#-
get_boundary_nodes(tri)

#-
DelaunayTriangulation.has_ghost_triangles(tri)