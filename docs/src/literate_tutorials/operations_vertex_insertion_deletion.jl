# # Triangulation Operations 
# ## Vertex Insertion and Deletion 

# This tutorial demonstrates how to insert and delete vertices from a 
# triangulation while maintaining the Delaunay property of the 
# triangulation. First, load the packages we need:
using DelaunayTriangulation
using CairoMakie
using StableRNGs
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

# Let us now define our initial triangulation.
points = [(0.0, 0.0), (2.0, 0.0), (1.0, 2.0)]
tri = triangulate(points)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "triangulation_operations_1.png") fig #src

# Note that we use a structure for `points` that is mutable so that 
# points can be pushed into it. 

# We now want to insert a vertex into the triangulation. To do this, 
# we use [`add_point!`](@ref) and provide the coordinates of the new point. There 
# is a method which uses the index of the vertex rather than the coordinates, 
# but we don't use that here as the points to be added are not already in 
# `points`. Here, we add a point at `(1.0, 0.5)`.
add_point!(tri, 1.0, 0.5)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "triangulation_operations_2.png") fig #src

# This is still a valid Delaunay triangulation, unsurprisingly due to the 
# small number of points. We can also add points that are outside of the triangulation:
add_point!(tri, 0.0, 1.0)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "triangulation_operations_3.png") fig #src

# One important thing to note here is that, if not for the ghost triangles inside `tri`,
# adding a point outside of the triangulation would not work. Here is an example of this failing.
delete_ghost_triangles!(tri)
try #hide
    add_point!(tri, 2.0, 1.5)
catch e #hide
println(e) #hide
end #hide

# This is a `BoundsError`, because the triangulation has had to 
# walk into the boundary but failed to find any ghost triangles. 
# Anytime you need to 
# perform some operation including a point outside of the boundary, 
# you need to be sure that you have ghost triangles, which you can query 
# using [`DelaunayTriangulation.has_ghost_triangles`](@ref). 
DelaunayTriangulation.has_ghost_triangles(tri)

#-
add_ghost_triangles!(tri)

#-
DelaunayTriangulation.has_ghost_triangles(tri)

# Another issue is that the convex hull is not updated as we add (or delete) 
# points for performance reasons:
get_convex_hull_vertices(tri)

# If we do want to fix the convex hull, we can use [`convex_hull!(tri)`](@ref).
convex_hull!(tri)
fig, ax, sc = triplot(tri, show_convex_hull = true)
fig
@test_reference joinpath(fig_path, "triangulation_operations_4.png") fig #src

# We now have the same triangulation that we would have had if we had done 
# `triangulate` on this set of points originally. To now push this further, 
# let's add in a bunch of random points. 
rng = StableRNG(123)
for _ in 1:1000
    new_point = 2rand(rng, 2)
    add_point!(tri, new_point)
end
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "triangulation_operations_5.png") fig #src

# Let us now demonstrate how to delete points. To do this, we use [`delete_point!`](@ref).
# This function takes vertices rather than coordinates, identifying the corresponding point 
# in `points` by its index. Let us demonstrate this operation by 
# deleting all points within a distance of `1/2` around `(1.0, 1.0)`.
vertices_to_delete = Iterators.filter(each_solid_vertex(tri)) do i
    p = get_point(tri, i)
    r2 = (getx(p) - 1.0)^2 + (gety(p) - 1.0)^2
    return r2 < 1 / 2^2
end
for i in vertices_to_delete
    delete_point!(tri, i)
end
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "triangulation_operations_6.png") fig #src

# Note that in this situation, `points` still contains those points that we have now deleted. 
# This is the reason to be careful about using, say, [`DelaunayTriangulation.each_point`](@ref) rather than [`each_solid_vertex`](@ref).
# This triangulation is also still Delaunay. 
