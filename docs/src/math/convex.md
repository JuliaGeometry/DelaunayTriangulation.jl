# Triangulating Convex Polygons

Here we describe an algorithm used for triangulating convex polygons, called Chew's algorithm after Chew's 1990 paper [_Building Voronoi diagrams for convex polygons in linear expected time_](https://dl.acm.org/doi/10.5555/867875). 

Suppose we have a counter-clockwise sequence of vertices $\mathcal S$ defining a convex polygon, and generate a random permutation of $\mathcal S$ denoted $\pi(\mathcal S)$ that defines the insertion order. The idea of the algorithm is to, starting from a triangle constructed by the first three vertices of $\pi(\mathcal S)$, add the remaining points in one at a time using the Bowyer-Watson algorithm, but leveraging the convexity of $\mathcal S$ to greatly simplify the point location step.

We can determine how to avoid the point location step by considering an example. Using the Bowyer-Watson algorithm, let's look at how a convex polygon gets triangulated.

```@setup convpoly
using DelaunayTriangulation
using StableRNGs 
using CairoMakie 
a = (-3.0, 4.0)
b = (-5.0, 3.0)
c = (-4.0, -2.0)
d = (-2.0, -1.0)
e = (0.0, 1.0)
f = (0.0, 3.0)
points = [a, b, c, d, e, f]
tri = triangulate(points; skip_points = [1, 3, 6])
fig = Figure()
ax = Axis(fig[1, 1], width = 300, height = 300)
triplot!(ax, tri)
scatter!(ax, points, color = :black)
scatter!(ax, [a], color = :red, markersize = 16)
add_point!(tri, 1)
ax = Axis(fig[1, 2], width = 300, height = 300)
triplot!(ax, tri)
scatter!(ax, points, color = :black)
scatter!(ax, [f], color = :red, markersize = 16)
add_point!(tri, 6)
ax = Axis(fig[1, 3], width = 300, height = 300)
triplot!(ax, tri)
scatter!(ax, points, color = :black)
scatter!(ax, [c], color = :red, markersize = 16)
add_point!(tri, 3)
ax = Axis(fig[1, 4], width = 300, height = 300)
triplot!(ax, tri)
scatter!(ax, points, color = :black)
resize_to_layout!(fig)
```

```@example convpoly
fig # hide
```

See that at each stage the vertex $v_u$ to be added, shown in red, it lies outside of the triangulation, and only a single edge $e_{vw}$ separates $v_u$ from the triangulation. Thus, we can identify that the point location step amounts to finding this edge $e_{vw}$ so that we inserting $v_u$ into the triangulation can be done by retriangulating the cavity formed by the union of the triangles $T_{uvw}$ and the triangles containing $u$ in their circumcircles.

We now need to find an efficient way to find this edge $e_{vw}$. Imagine running the Bowyer-Watson algorithm in reverse, meaning removing vertices one at a time in the reverse order of $\pi(\mathcal S)$. When we remove $v_u$, this is the same as connecting its neighbours $v_v$ and $v_w$ with an edge $e_{vw}$, which is exactly the edge needed for our point location problem. Thus, if we keep track of the polygon vertices and their neighbours, we can easily find the edge $e_{vw}$ in constant time. 

Using this insight, we can now present Chew's algorithm:

1. Write $\mathcal S = \{v_1, \ldots, v_n\}$ and obtain some random permutation $\pi(\mathcal S)$ of $\mathcal S$, representing the permutation $\pi(\mathcal S)$ as a permutation of $\{1, 2, \ldots, n\}$ chosen uniformly at random.
2. Construct a circularly- and doubly-linked list of the vertices in $\mathcal S$ by defining $\mathcal S_{\text{next}} = \{2, 3, \ldots, n, 1\}$ and $\mathcal S_{\text{prev}} = \{n, 1, 2, \ldots, n-1\}$.
3. For $i = n,n-1,\ldots,4$: Delete $v_{\pi(\mathcal S)[i]}$ from the list by setting $\mathcal S_{\text{next}}[\mathcal S_{\text{prev}}[\pi[i]]] = \mathcal S_{\text{next}}[\pi[i]]$ and $\mathcal S_{\text{prev}}[\mathcal S_{\text{next}}[\pi[i]]] = \mathcal S_{\text{prev}}[\pi[i]]$, where $\pi[i] \equiv \mathcal \pi(\mathcal S)[i]$.
4. Initialise the triangulation by adding $T_{v_{\pi[1]}v_{\pi[2]}v_{\pi[3]}}$ to the triangulation.
5. For $i = 4, \ldots, k$: Add $v_{\pi[i]}$ into the triangulation using the Bowyer-Watson algorithm, noting that the polygonal cavity can be evacuated starting from the edge $e_{\mathcal S_{\text{next}}[\pi[i]]\mathcal S_{\text{prev}}[\pi[i]]}$.

