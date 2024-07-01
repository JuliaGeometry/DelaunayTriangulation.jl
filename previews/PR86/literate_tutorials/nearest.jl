# # Nearest Neighbour Queries
# 

# One useful feature of Voronoi tessellations is that they can be used to
# obtain nearest neighbours, since by definition a Voronoi tile contains all 
# points in the plane closest to the associated generator. This implies that 
# finding a nearest neighbour is the same as a point location problem, meaning 
# given a point `p` find the Voronoi tile `P` containing it. Here we give an example 
# of how we can use triangulations or tessellations to find the nearest neighbour 
# in the point set to a given point. First, we load in the packages we need.
using DelaunayTriangulation
using CairoMakie
using ReferenceTests #src 
using Test #src 

# Now we define the tessellation we will use for this example. The white points 
# shown are the points that we will query.
points = [
    (-3.0, 7.0), (2.0, 6.0), (0.0, 3.0),
    (0.0, 0.0), (-5.0, 5.0), (-3.0, 1.0),
    (2.0, -3.0), (5.0, 5.0), (-4.0, 3.0)
]
tri = triangulate(points, delete_ghosts=false)
vorn = voronoi(tri)
p, q = (-2.0, 7.5), (0.0, 4.0)
fig, ax, sc = voronoiplot(vorn, markersize=14)
scatter!(ax,[p,q],color=:white,strokecolor=:black,strokewidth=2,markersize=14)
fig

# To get the nearest neighbour of a point, we use `get_nearest_neighbour`. 
np = get_nearest_neighbour(vorn, p)
@test np == 1 #src 

#-
nq = get_nearest_neighbour(vorn, q)
@test nq == 3 #src 

# We see that the nearest point in `points` to `p` is the first point, and to `q` it is 
# the third point. We note that we could have also performed this query without constructing `vorn` directly,
# instead using `tri`:
np_tri = get_nearest_neighbour(tri, p)
@test np_tri == 1 #src

#-
nq_tri = get_nearest_neighbour(tri, q)
@test nq_tri == 3 #src

# Both methods lead to the same results because they use the same algorithm.