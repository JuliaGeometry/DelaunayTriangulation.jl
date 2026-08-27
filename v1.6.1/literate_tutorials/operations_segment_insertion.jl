# # Triangulation Operations 
# ## Segment Insertion 

# This tutorial shows how we can add segments into a triangulation. First, load the packages we need:
using DelaunayTriangulation
using CairoMakie
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

# Let us now define our initial triangulation.
points = [
    (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0),
    (0.9, 0.9), (0.5, 0.5), (0.2, 0.5), (0.5, 0.8),
]
tri = triangulate(points)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "triangulation_operations_8.png") fig #src

# To add a segment, we use [`add_segment!`](@ref), providing the vertices for the points 
# that a segment should be added between. Let us add a segment between `(0.0, 0.0)`
# and `(1.0, 1.0)`, which corresponds to vertices `1` and `3`.
add_segment!(tri, 1, 3)
fig, ax, sc = triplot(tri, show_constrained_edges = true)
fig
@test_reference joinpath(fig_path, "triangulation_operations_9.png") fig #src

# Of course, this changed nothing since the segment was already there. We do note,
# though, that if we look at the constrained edges
get_interior_segments(tri)

# then we notice that the segment `(1, 3)` was converted into the segments `(1, 6)`, `(5, 3)`,
# and `(5, 6)`. This is because the segment `(1, 3)` crossed through other vertices, 
# and so the algorithm automatically breaks down the segments into a sequence of connected collinear 
# segments. 

# Now we add a segment that was not already there. 
add_segment!(tri, 1, 8)
fig, ax, sc = triplot(tri, show_constrained_edges = true)
fig
@test_reference joinpath(fig_path, "triangulation_operations_10.png") fig #src

# Currently, the segments that you add must not intersect at an angle (they can be 
# collinear with other edges as we have demonstrated above). To see what happens 
# if we do this:
add_segment!(tri, 8, 2)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "triangulation_operations_11.png") fig #src

# The other constrained edge was partially removed.
