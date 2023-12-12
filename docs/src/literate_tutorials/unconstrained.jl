# # Unconstrained Triangulations
# 

# Our first example considers constructing unconstrained 
# Delaunay triangulations. To start, let us load in the 
# packages we will need. 
using DelaunayTriangulation
using CairoMakie # for plotting
using StableRNGs # for reproducibility
using LinearAlgebra # used for computing norms later
using StatsBase #src
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

# We consider just triangulating a random set of points. First, generating
# the points:
rng = StableRNG(123)
points = rand(rng, 2, 500) # just do rand(2, 500) if you are not concerned about the RNG

# We now triangulate these points by using `triangulate`. We pass the `rng` 
# as a keyword argument, but again if you are not concerned about the RNG (or 
# set the seed using `Random.seed!`) then you can ignore this.
tri = triangulate(points; rng=rng)

#-
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "unconstrained_ex_1.png") fig #src

# In this plot, the red line shows the convex hull. This `tri` object 
# is our triangulation, and we can interact with it in many ways. 

# ## Iterating over vertices

# For example, we can iterate over the points in the triangulation using `each_solid_vertex`. 
# Here we compute the centroid of the point cloud:
function compute_centroid(tri)
    s = [0.0, 0.0]
    for i in each_solid_vertex(tri)
        p = get_point(tri, i)
        s .+= p
    end
    s ./= num_solid_vertices(tri)
    return s
end
s = compute_centroid(tri)
@test s ≈ [sum(points[1, :]) / 500, sum(points[2, :]) / 500] atol = 1e-10 #src

# We need to use the `solid` identifier because triangulations are made up of both _solid_ 
# and _ghost_ vertices/edges/triangles, for reasons described in the [manual](../manual/ghost_triangles.md).
# If we just use `each_vertex(tri)`, then we also find a vertex `-1` that corresponds to the boundary. For example, 
# the points on the boundary can be obtained using:
get_neighbours(tri, -1)
@test (sort ∘ collect ∘ get_neighbours)(tri, -1) == sort(unique(get_convex_hull_indices(tri))) #src

# One reason to be careful with this especially is that `get_point(tri, -1)` does actually correspond to a coordinate,
get_point(tri, -1)
@test get_point(tri, -1) == DelaunayTriangulation.get_representative_point_coordinates(tri, 1) #src

# (This is the pole of inaccessibility for the domain; see [here](../tutorials/pole_of_inaccessibility.md) for more details.)
# You can use `each_vertex` or `each_ghost_vertex` to consider all types of vertices or only the ghost vertices.
# If you just want the vertices, use `get_vertices(tri)`, which will also include the boundary vertex.  
get_vertices(tri)
@test get_vertices(tri) == tri.graph.graph.V #src

# ## Iterating over edges 

# We can also iterate over the edges of the triangulation using `each_solid_edge`, or 
# `each_edge` for both solid and ghost edges and `each_ghost_edge` for only the ghost edges. 
# To give an example, here we compute the average length of an edge.
function compute_mean_edge_length(tri)
    ℓ = 0.0
    for e in each_solid_edge(tri)
        u, v = edge_indices(e)
        p, q = get_point(tri, u, v)
        ℓ += norm(p .- q)
    end
    ℓ /= num_solid_edges(tri)
    return ℓ
end
ℓ = compute_mean_edge_length(tri)
@test ℓ ≈ 0.056456938080335355 #src

# By default, the triangulation has no ghost edges, so `each_edge` and `each_solid_edge` are the same.
# You can also use the accessor `get_edges` to get all the edges, including the ghost edges if there were any.
get_edges(tri)
@test get_edges(tri) == tri.graph.graph.E #src

# Note also that the edges are all given as unordered, so the set of edges only includes 
# one of `(i, j)` and `(j, i)` for each edge `(i, j)`.

# ## Iterating over triangles 

# Similarly, we can iterate over the triangles using `each_triangle`, `each_ghost_triangle`, 
# or `each_triangle`. By default, ghost triangles are not included in the output. 
# Here we compute the area of the domain by getting the area of each triangle.
area(p, q, r) = 0.5 * ((getx(q) - getx(p)) * (gety(r) - gety(p)) - (gety(q) - gety(p)) * (getx(r) - getx(p)))
function compute_triangulation_area(tri)
    A = 0.0
    for T in each_solid_triangle(tri)
        i, j, k = indices(T)
        p, q, r = get_point(tri, i, j, k)
        A += area(p, q, r)
    end
    return A
end
A = compute_triangulation_area(tri)
@test A ≈ get_total_area(tri) #src

# (You can compute areas like this using `get_total_area(tri)`.) You can access the set of 
# `triangles` using `get_triangles(tri)`:
get_triangles(tri)
@test get_triangles(tri) == tri.triangles #src

# The triangles are all positively oriented, meaning the triangles are given such that the 
# corresponding points are traversed in counter-clockwise order. 

# ## Querying neighbour information

# We can query neighbour information at several levels. 

# ### Points 
# For a given point, there are two type of neighbours: The neighbouring vertices, 
# and the neighbouring triangles. The neighbours can be obtained using `get_neighbours`.
# For example, the set of vertices that share an edge with the fifth vertex is:
get_neighbours(tri, 5)
@test get_neighbours(tri, 5) == Set([93, 117, 301, 214]) #src

# The set of triangles that share an edge with the fifth vertex is obtained 
# using `get_adjacent2vertex`. This returns a set of edges `(v, w)` such that,
# for a given vertex `u`, `(u, v, w)` is a positively oriented triangle in the 
# triangulation. For example,
get_adjacent2vertex(tri, 5)
@test get_adjacent2vertex(tri, 5) == Set(((93, 117), (117, 301), (301, 214), (214, 93))) #src

# means that the triangles that contain `5` as a vertex are `(5, 93, 117)`, `(5, 117, 301)`,
# `(5, 301, 214)`, and `(5, 214, 93)`. We can verify this: 
filter(T -> 5 ∈ indices(T), get_triangles(tri))

# These queries can also be applied to the ghost vertices, in which information about the 
# boundary is provided.
get_neighbours(tri, -1)

#-
get_adjacent2vertex(tri, -1)

# ### Edges 
# For a given edge `(u, v)`, the relevant neighbours are the vertices that are next to it 
# so that a triangle is formed. We can find the vertex `w` such that `(u, v, w)` is a positively 
# oriented triangle in the triangulation using `get_adjacent(tri, u, v)`. For example, 
get_adjacent(tri, 163, 365)

# means that `(163, 365, 295)` is a positively oriented triangle, as we can verify:
DelaunayTriangulation.triangle_orientation(tri, 163, 365, 295)
@test DelaunayTriangulation.is_positively_oriented((DelaunayTriangulation.triangle_orientation(tri, 163, 365, 295))) #src

# (The representation of this predicate using a `Certificate` is described in more detail 
# in the [manual](../manual/predicates.md).) The other triangle adjoining the unordered 
# edge `(u, v)`, meaning the oriented edge `(v, u)`, is obtained similarly:
get_adjacent(tri, 365, 163)
@test DelaunayTriangulation.is_positively_oriented((DelaunayTriangulation.triangle_orientation(tri, 365, 163, get_adjacent(tri, 365, 163)))) #src

# If an edge `(u, v)` is on the boundary, oriented so that there is no solid vertex `w` 
# such that `(u, v, w)` is a triangle in the triangulation, then `get_adjacent(tri, u, v)`
# returns the boundary vertex. For example,
get_adjacent(tri, 398, 258)
@test get_adjacent(tri, 398, 258) == -1 #src

# means that `(398, 258)` is a boundary edge and `(398, 258, -1)` is a ghost triangle. 
# You can test for this case using `is_boundary_edge`:
DelaunayTriangulation.is_boundary_edge(tri, 398, 258)
@test DelaunayTriangulation.is_boundary_edge(tri, 398, 258) #src