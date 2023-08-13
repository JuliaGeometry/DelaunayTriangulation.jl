# # Clipped Voronoi Tessellations 
# ## Clipping to the Convex Hull 
# 

# One issue that may arise when dealing with Voronoi tessellations is the 
# presence of unbounded polygons occurring on the boundary. One way to deal with this 
# is to clip polygons to the convex hull of the tessellation. (Arbitrary clipping boundaries 
# are on the to-do list, but they are not yet implemented - see [Issue #48](https://github.com/DanielVandH/DelaunayTriangulation.jl/issues/48).)

# In the example below, we clip the tessellation to the convex hull of the point set by setting 
# the second argument of `voronoi` to `true`.
