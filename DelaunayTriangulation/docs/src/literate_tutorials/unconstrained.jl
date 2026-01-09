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

# We now triangulate these points by using [`triangulate`](@ref). We pass the `rng` 
# as a keyword argument, but again if you are not concerned about the RNG (or 
# set the seed using `Random.seed!`) then you can ignore this.
tri = triangulate(points; rng = rng)

#-
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "unconstrained_ex_1.png") fig #src

# This `tri` object is our [`Triangulation`](@ref), and we can interact with it in many ways. 

# ## Iterating over vertices
# For example, we can iterate over the points in the triangulation using [`each_solid_vertex`](@ref). 
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
@test s ≈ [sum(points[1, :]) / 500, sum(points[2, :]) / 500] atol = 1.0e-10 #src

# We need to use the `solid` identifier because triangulations are made up of both _solid_ 
# and _ghost_ vertices/edges/triangles, for reasons described in the [manual](../manual/ghost_triangles.md).
# If we just use [`each_vertex(tri)`](@ref each_vertex), then we also find a vertex `-1` that corresponds to the boundary. For example, 
# the points on the boundary can be obtained using:
get_neighbours(tri, -1)
@test (sort ∘ collect ∘ get_neighbours)(tri, -1) == sort(unique(get_convex_hull_vertices(tri))) #src

# One reason to be careful with this especially is that `get_point(tri, -1)` does actually correspond to a coordinate,
get_point(tri, -1)
@test get_point(tri, -1) == DelaunayTriangulation.get_representative_point_coordinates(tri, 1) #src

# (This is the pole of inaccessibility for the domain; see [here](../tutorials/pole_of_inaccessibility.md) for more details.)
# You can use [`each_vertex`](@ref) or [`each_ghost_vertex`](@ref) to consider all types of vertices or only the ghost vertices.
# If you just want the vertices, use `each_vertex(tri)`, which will also include the ghost vertex.  
each_vertex(tri)
@test each_vertex(tri) == tri.graph.vertices #src

# ## Iterating over edges 
# We can also iterate over the edges of the triangulation using [`each_solid_edge`](@ref), or 
# [`each_edge`](@ref) for both solid and ghost edges and [`each_ghost_edge`](@ref) for only the ghost edges. 
# To give an example, here we compute the average length of an edge.
function compute_mean_edge_length(tri)
    ℓ = 0.0
    for e in each_solid_edge(tri)
        u, v = edge_vertices(e)
        p, q = get_point(tri, u, v)
        ℓ += norm(p .- q)
    end
    ℓ /= num_solid_edges(tri)
    return ℓ
end
ℓ = compute_mean_edge_length(tri)
@test ℓ ≈ 0.056456938080335355 #src

# By default, the triangulation has ghost edges, so [`each_edge`](@ref) and [`each_solid_edge`](@ref) are not the same.
each_edge(tri)
@test each_edge(tri) == tri.graph.edges #src
@test each_edge(tri) != each_solid_edge(tri) #src

# Note also that the edges are all given as unordered, so the set of edges only includes 
# one of `(i, j)` and `(j, i)` for each edge `(i, j)`.

# ## Iterating over triangles 
# Similarly, we can iterate over the triangles using [`each_solid_triangle`](@ref), [`each_ghost_triangle`](@ref), 
# or [`each_triangle`](@ref). By default, ghost triangles are included in the output. 
# Here we compute the area of the domain by getting the area of each triangle.
area(p, q, r) = 0.5 * ((getx(q) - getx(p)) * (gety(r) - gety(p)) - (gety(q) - gety(p)) * (getx(r) - getx(p)))
function compute_triangulation_area(tri)
    A = 0.0
    for T in each_solid_triangle(tri)
        i, j, k = triangle_vertices(T)
        p, q, r = get_point(tri, i, j, k)
        A += area(p, q, r)
    end
    return A
end
A = compute_triangulation_area(tri)
@test A ≈ get_area(tri) #src

# (You can compute areas like this using [`get_area(tri)`](@ref get_area).) You can access the set of 
# `triangles` using [`get_triangles(tri)`](@ref get_triangles):
get_triangles(tri)
@test get_triangles(tri) == tri.triangles #src

# The triangles are all positively oriented, meaning the triangles are given such that the 
# corresponding points are traversed in counter-clockwise order. 

# ## Querying neighbour information
# We can query neighbour information at several levels. 

# ### Points 
# For a given point, there are two type of neighbours: The neighbouring vertices, 
# and the neighbouring triangles. The neighbours can be obtained using [`get_neighbours`](@ref).
# For example, the set of vertices that share an edge with the fifth vertex is:
get_neighbours(tri, 5)
@test get_neighbours(tri, 5) == Set([93, 117, 301, 214]) #src

# The set of triangles that share an edge with the fifth vertex is obtained 
# using [`get_adjacent2vertex`](@ref). This returns a set of edges `(v, w)` such that,
# for a given vertex `u`, `(u, v, w)` is a positively oriented triangle in the 
# triangulation. For example,
get_adjacent2vertex(tri, 5)
@test get_adjacent2vertex(tri, 5) == Set(((93, 117), (117, 301), (301, 214), (214, 93))) #src

# means that the triangles that contain `5` as a vertex are `(5, 93, 117)`, `(5, 117, 301)`,
# `(5, 301, 214)`, and `(5, 214, 93)`. We can verify this: 
filter(T -> 5 ∈ triangle_vertices(T), get_triangles(tri))

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

# (The representation of this predicate using a [`DelaunayTriangulation.Certificate`](@ref) is described in more detail 
# in the [manual](../manual/predicates.md).) The other triangle adjoining the unordered 
# edge `(u, v)`, meaning the oriented edge `(v, u)`, is obtained similarly:
get_adjacent(tri, 365, 163)
@test DelaunayTriangulation.is_positively_oriented((DelaunayTriangulation.triangle_orientation(tri, 365, 163, get_adjacent(tri, 365, 163)))) #src

# If an edge `(u, v)` is on the boundary, oriented so that there is no solid vertex `w` 
# such that `(u, v, w)` is a triangle in the triangulation, then `get_adjacent(tri, u, v)`
# returns the ghost vertex. For example,
get_adjacent(tri, 398, 258)
@test get_adjacent(tri, 398, 258) == -1 #src

# means that `(398, 258)` is a boundary edge and `(398, 258, -1)` is a ghost triangle. 
# You can test for this case using [`DelaunayTriangulation.is_boundary_edge`](@ref):
DelaunayTriangulation.is_boundary_edge(tri, 258, 398)
@test DelaunayTriangulation.is_boundary_edge(tri, 258, 398) #src
