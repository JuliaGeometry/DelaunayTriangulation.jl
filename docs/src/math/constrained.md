# Constrained Delaunay Triangulations 

Now we introduce constrained Delaunay triangulations. Constrained Delaunay triangulations are similar to Delaunay triangulations $\mathcal D\mathcal T(\mathcal P)$, except with the additional requirement that a provided set of segments $\mathcal S$ are constrained to be part of the triangulation, giving a triangulation denoted $\mathcal D\mathcal T(\mathcal P, \mathcal S)$.

An important concept for defining constrained Delaunay triangulations is that of _visibility_. We say that two points $p$ and $q$ are visible to each other if the open line segment $pq$ does not intersect any other segment in $\mathcal D\mathcal T(\mathcal P, \mathcal S)$. We see that a triangle $T \in \mathcal D\mathcal T(\mathcal P, \mathcal S)$ is _constrained Delaunay_ if its open circumcircle contains no point in $\mathcal P$ that is visible from any point inside $T$. Using this definition, we can thus define the _constrained Delaunay triangulation_ of $\mathcal P$ with respect to $\mathcal S$ as the triangulation in which every triangle is constrained Delaunay. There is a corresponding constrained Delaunay lemma just as for unconstrained Delaunay triangulations, but we do not list it here.

## Incremental Insertion Algorithm

The algorithm for constructing constrained Delaunay triangulations is incremental and reasonably simply to describe. To compute $\mathcal D\mathcal T(\mathcal P, \mathcal S)$, we:

1. Compute $\mathcal D\mathcal T(\mathcal P)$, the unconstrained Delaunay triangulation.
2. Insert each segment $s \in \mathcal S$ into $\mathcal D\mathcal T(\mathcal P)$ one by one.

Of course, this second step is where all the complication lies. Most of the discussion in this section will be about the problem of inserting segments into a triangulation.

### Inserting a Segment into a Triangulation

Let's now discuss how we can add a segment into a triangulation. To understand the algorithm, let's first look at an example. Consider the figure below.

```@setup segins1
using DelaunayTriangulation 
using StableRNGs
using CairoMakie
a = (0.0, 0.0)
b = (0.0, 1.0)
c = (0.0, 2.5)
d = (2.0, 0.0)
e = (6.0, 0.0)
f = (8.0, 0.0)
g = (8.0, 0.5)
h = (7.5, 1.0)
i = (4.0, 1.0)
j = (4.0, 2.5)
k = (8.0, 2.5)
pts = [a, b, c, d, e, f, g, h, i, j, k]
rng = StableRNG(213)
tri = triangulate(pts; rng)
e = (2, 7)
fig, ax, sc = triplot(tri)
lines!(ax, [get_point(tri, e...)...], color = :blue, linewidth = 4)
history = DelaunayTriangulation.PointLocationHistory{NTuple{3,Int},NTuple{2,Int},Int}()
q = get_point(tri, 7)
find_triangle(tri, q; k = 2, history, store_history = Val(true))
lines!(ax, [get_point(tri, 2, history.left_vertices..., 7)...], color = :red, linewidth = 4)
lines!(ax, [get_point(tri, 2, history.right_vertices..., 7)...], color = :green, linewidth = 4)
for T in history.triangles
    i, j, k = triangle_vertices(T)
    pp, qq, rr = get_point(tri, i, j, k)
    poly!(ax, [pp, qq, rr], color = (:magenta, 0.3))
end
```

```@example segins1 
fig # hide
```

To develop an algorithm, we need to notice one important thing from this figure: Since the blue segment will occlude visibility between any points on either side of the segment, the blue segment effectively divides, locally, the triangulation into two parts that no longer interact with each other. In the figure above, this means that any changes to the triangles bounded between the red curve and the blue segment will not interact those in the region bounded between the green curve and the bkue segment. This is a key observation that will allow us to develop an algorithm for inserting segments: we can handle the two sides separately.

We need to also understand what we are showing in this figure. The highlighted triangles show all triangles intersected by the blue segment. The red boundary shows the chain of vertices intersected by the blue segment above the segment. There is a clear problem with this boundary though: there is a dangling edge, caused by the segment intersecting through two triangles that share an edge, except they are not hit right after each other. The main problem with this is that the object defined by the union of the blue segment and the red boundary is not technically a polygon. We treat it as if it were a polygon, though, by imagingin an edge walking around the boundary of this object and splitting the dangling vertex into two copies of itself, so that the ant walking around the boundary essentially sees two different vertices. (It is possible for a vertex to be repeated three or more times, but this is much more rare.) The green boundary is defined similarly.

What can we do now with this information? We have now established that the blue segment defines two polygonal cavities $\mathcal C_1$ and $\mathcal C_2$ that are adjacent to each other, sharing only the blue segment, and any changes to the triangulation within $\mathcal C_1$ and $\mathcal C_2$ will have no affect on the other cavity. We thus devise the following algorithm:

1. Given a segment $e_{ij}$ to be inserted, find the triangles intersected by the segment.
2. Delete all triangles intersected by the segment from the triangulation and add $e_{ij}$.
3. Identify the two polygonal cavities $\mathcal C_1$ and $\mathcal C_2$ on each side of $e_{ij}$, taking care for any danling edges by repeating the vertices as needed.
4. Retriangulate $\mathcal C_1$ and $\mathcal C_2$ separately.
5. Add the triangles from the triangulated polygonal cavities into the original triangulation.

Once step five is complete, the triangulation has now successfully added $e_{ij}$ into the constrained triangulation. We are of course skipping over some key facts, like demonstrating why there are no other changes to the triangulation away from $\mathcal C_1$ and $\mathcal C_2$; for these details, see the original paper by [Shewchuk and Brown (2015)](https://www.sciencedirect.com/science/article/pii/S0925772115000322).

To actually implement this algorithm, there are some key details that need to be worked out:

1. How can we find all triangles intersected by $e_{ij}$?
2. How can we identify the polygonal cavities $\mathcal C_1$ and $\mathcal C_2$?
3. How can we efficiently retriangulate $\mathcal C_1$ and $\mathcal C_2$?

Let's address these one at a time.

### Finding Triangles Intersected by a Segment 

The problem is: Given a segment $e_{ij}$, find the set of all triangles intersected by the segment. This is exactly the problem that is solved by our point location algorithm - remember that our algorithm marches along all triangles from $p_i$ to $p_j$ until it finds the triangle containing $p_j$. Thus, by keeping track of all triangles visited during a point location algorithm jumping from $p_i$ to $p_j$, we can find all the triangles intersected by the segment.

### Finding the Polygonal Cavities

Now we need to consider the problem of finding the two polygonal cavities $\mathcal C_1$ and $\mathcal C_2$. These can also be found by the point location algorithm. In addition to keeping track of all triangles visited, whenever an edge $e_{k\ell}$ is stepped over we keep track of the vertices $k$ and $\ell$ and put them into the list of vertices for either $\mathcal C_1$ and $\mathcal C_2$, depending on which side of the segment $e_{ij}$ they are on. There are some important details here to consider for e.g. collinear edges, but we do not address this here.

### Efficiency Triangulating the Polygonal Cavities

Now we must address the problem of retriangulating the polygonal cavities. This issue is the same for both $\mathcal C_1$ and $\mathcal C_2$, so we discuss this problem in general for a given polygonal cavity $\mathcal C$ which could be either $\mathcal C_1$ or $\mathcal C_2$. While we could just simply triangulate $\mathcal C$ using a Bowyer-Watson algorithm, there is a better approach. It turns out that Chew's algorithm for triangulating convex polygons, as described in the [convex triangulation section](convex.md), also works for the polygonal cavity $\mathcal C$. This is what we use for triangulating $\mathcal C$. 

### Putting Everything Together

Now understanding all these details, the algorithm for inserting a segment into a triangulation follows the same procedure as outlined above. We first find all triangles intersected by the segment, then find the polygonal cavities $\mathcal C_1$ and $\mathcal C_2$, and then retriangulate these polygonal cavities using Chew's algorithm. Once this is done, we have successfully inserted the segment into the triangulation.

## Boundaries 

One other important part of a constrained Delaunay triangulation is the enforcement of boundaries and holes in a triangulation. This part of the algorithm comes after all segments have been computed. The problem is as follows: Given a set of boundaries defined by segments, delete all triangles that fall outside of the boundaries. The approach we take for this follows the description given by Shewchuk in his [Triangle paper](https://link.springer.com/chapter/10.1007/BFb0014497):

1. First, we find all points that will be deleted. To do this, we go across each boundary edge (from the user-provided boundary, not necessarily of the triangulation) and identify the triangle adjacent to the edge away from the interior. The vertex of this triangle away from the boundary is used to find further points, checking all neighbours of this vertex and so on, deleting all vertices that are not in the boundary.
2. Using the set of points identified for deletion, we then identify all triangles to delete. This is done by looking all triangles that have a vertex in the set of points to delete, using the `Adjacent2Vertex` map to do so. We also need to assess the boundary edges once again. Safety is important here since, for example, a triangle might be comprised of vertices that all lie on the boundary but the triangle's interior is outside of the boundary. So, we check each triangle and ensure that those we've identified for deletion have a centroid that appear outside of the boundary. If this is not the case, we remove it from the set of triangles marked for deletion.
3. Now having all the triangles to delete, we delete them all from the triangulation.

With these steps complete, we will have our complete triangulation. To make these steps clear, the figure below shows an example.

```@setup segins1
a = (-3.0, 4.0)
b = (3.0, 4.0)
c = (2.0, 2.0)
d = (3.0, 0.0)
e = (0.0, 1.0)
f = (2.0, -2.0)
g = (0.0, -1.0)
h = (-1.84, -2.04)
i = (-3.0, 0.0)
j = (-2.0, 3.0)
k = (-2.0, 0.0)
l = (1.0, 3.0)
m = (-1.0, 2.0)
n = (-5.0, 4.0)
o = (-3.0, 6.0)
p = (5.0, 6.0)
q = (7.0, 3.0)
r = (-4.82, 1.0)
points = [o, p, n, q, m, r]
boundary_nodes = [[[h, g, f, e, d, c, b, a, i, h]], [[j, l, k, j]]]
boundary_nodes, points = convert_boundary_points_to_indices(boundary_nodes; existing_points = points)
tri = triangulate(points; boundary_nodes, delete_holes = false)
fig = Figure()
ax = Axis(fig[1, 1], width = 300, height = 300, title = "Triangulation", titlealign = :left)
triplot!(ax, tri, show_points = true)
lines!(ax, [h, g, f, e, d, c, b, a, i, h, (NaN, NaN), j, l, k, j], color = :red, linewidth = 4)
points_to_process = DelaunayTriangulation.find_all_points_to_delete(tri)
ax2 = Axis(fig[1, 2], width = 300, height = 300, title = "Points for deletion", titlealign = :left)
triplot!(ax2, tri, show_points = true)
lines!(ax2, [h, g, f, e, d, c, b, a, i, h, (NaN, NaN), j, l, k, j], color = :red, linewidth = 4)
scatter!(ax2, [get_point(tri, filter(!DelaunayTriangulation.is_ghost_vertex, points_to_process)...)...], color = :blue, markersize = 10)
triangles_to_delete = DelaunayTriangulation.find_all_triangles_to_delete(tri, points_to_process)
ax3 = Axis(fig[1, 3], width = 300, height = 300, title = "Triangles for deletion", titlealign = :left)
triplot!(ax3, tri, show_points = true)
lines!(ax3, [h, g, f, e, d, c, b, a, i, h, (NaN, NaN), j, l, k, j], color = :red, linewidth = 4)
for T in triangles_to_delete
    DelaunayTriangulation.is_ghost_triangle(T) && continue
    i, j, k = triangle_vertices(T)
    pp, qq, rr = get_point(tri, i, j, k)
    poly!(ax3, [pp, qq, rr], color = (:magenta, 0.3))
end
ax4 = Axis(fig[1, 4], width = 300, height = 300, title = "Final triangulation", titlealign = :left)
DelaunayTriangulation.delete_all_exterior_triangles!(tri, triangles_to_delete)
triplot!(ax4, tri, show_points = true)
resize_to_layout!(fig)
```

```@example segins1
fig # hide
```

## Adding Points into a Constrained Delaunay Triangulation

Adding points into a constrained Delaunay triangulation is similar to the unconstrained case. The Bowyer-Watson algorithm is again used, except that, when we are excavating the cavity for the new point, we avoid stepping across any segments. If a point is inserted onto a segment, we dig cavities separately on each side of the segment.