# # Clipped Voronoi Tessellations 
# ## Clipping to a Rectangle
# 

# In the previous tutorial, we demonstrated how we can clip to
# the convex hull of the point set. However, it is often useful to clip
# to a rectangle, for example if you want to clip to a region of interest
# in a simulation. We do not yet support this within `voronoi` itself, 
# but we provide the function `get_polygon_coordinates` for this (this is what 
# `voronoiplot` uses to plot inside a bounding box).

# Let us now demonstrate. First, we construct a tessellation of 
# some example point set.
using DelaunayTriangulation 
using CairoMakie
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src
A = (-3.0, 7.0)
B = (1.0, 6.0)
C = (-1.0, 3.0)
D = (-2.0, 4.0)
E = (3.0, -2.0)
F = (5.0, 5.0)
G = (-4.0, -3.0)
H = (3.0, 8.0)
points = [A, B, C, D, E, F, G, H]
tri = triangulate(points)
vorn = voronoi(tri)

# Let us show the tessellation, and the rectangle we want to clip 
# the tessellation to. 
fig, ax, sc = voronoiplot(vorn)
a, b, c, d = -2.0, 3.0, 0.0, 7.0
lines!(ax, [(a,c),(b,c),(b,d),(a,d),(a,c)], color = :black, linewidth = 4)
fig
@test_reference joinpath(fig_path, "voronoi_ex_3.png") fig #src

# To apply this clipping, we need to provide a bounding box of the form 
# `(xmin, xmax, ymin, ymax)`. Here, we will use 
bounding_box = (a, b, c, d)

# You can obtain some reasonable defaults for this bounding box using 
# DelaunayTriangulation.polygon_bounds(vorn); see the docstring for more information. 
# The coordinates for each polygon clipped to this box can be obtained as follows.
clipped_coords = Vector{Vector{NTuple{2,Float64}}}(undef, num_polygons(vorn))
for i in each_polygon_index(vorn)
    clipped_coords[i] = get_polygon_coordinates(vorn, i, bounding_box)
end
clipped_coords 

# Now let's plot these. 
filter!(!isempty, clipped_coords) #hide
fig, ax, sc = poly(clipped_coords, color = :white, strokewidth = 4)
@test_reference joinpath(fig_path, "voronoi_ex_4.png") fig #src

# As we can see, the polygons have been clipped to the rectangle.
# Note that if you just want this for plotting, you can also call `voronoiplot` with the 
# `bounding_box` keyword argument.