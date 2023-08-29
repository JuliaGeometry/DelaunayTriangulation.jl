```@meta
EditURL = "https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/operations_segment_insertion.jl"
```

# Triangulation Operations
## Edge Insertion

This tutorial shows how we can add edges into a triangulation. First, load the packages we need:

````@example operations_segment_insertion
using DelaunayTriangulation
using CairoMakie
````

Let us now define our initial triangulation.

````@example operations_segment_insertion
points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0),
(0.9, 0.9), (0.5, 0.5), (0.2, 0.5), (0.5, 0.8)]
tri = triangulate(points)
fig, ax, sc = triplot(tri)
fig
````

To add an edge, we use `add_edge!`, providing the indices for the points
that an edge should be added between. Let us add an edge between `(0.0, 0.0)`
and `(1.0, 1.0)`, which corresponds to vertices `1` and `3`.

````@example operations_segment_insertion
add_edge!(tri, 1, 3)
fig, ax, sc = triplot(tri, show_constrained_edges = true)
fig
````

Of course, this changed nothing since the edge was already there. We do note,
though, that if we look at the constrained edges

````@example operations_segment_insertion
get_constrained_edges(tri)
````

then we notice that the edge `(1, 3)` was converted into the edges `(1, 6)`, `(5, 3)`,
and `(5, 6)`. This is because the edge `(1, 3)` crossed through other vertices,
and so the algorithm automatically breaks down the edges into a sequence of connected collinear
edges.

Now we add an edge that was not already there.

````@example operations_segment_insertion
add_edge!(tri, 1, 8)
fig, ax, sc = triplot(tri, show_constrained_edges = true)
fig
````

Currently, the edges that you add must not intersect at an angle (they can be
collinear with other edges as we have demonstrated above). To see what happens
if we do this:

````@example operations_segment_insertion
add_edge!(tri, 8, 2)
fig, ax, sc = triplot(tri)
fig
````

The other constrained edge was partially removed. Eventually, we would like to
make it so that the intersection is automatically detected and a point
is added at the intersection to break up the inserted edge.
## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/new-docs/docs/src/literate_tutorials/operations_segment_insertion.jl).

```julia
using DelaunayTriangulation
using CairoMakie

points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0),
(0.9, 0.9), (0.5, 0.5), (0.2, 0.5), (0.5, 0.8)]
tri = triangulate(points)
fig, ax, sc = triplot(tri)
fig

add_edge!(tri, 1, 3)
fig, ax, sc = triplot(tri, show_constrained_edges = true)
fig

get_constrained_edges(tri)

add_edge!(tri, 1, 8)
fig, ax, sc = triplot(tri, show_constrained_edges = true)
fig

add_edge!(tri, 8, 2)
fig, ax, sc = triplot(tri)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

