```@meta
EditURL = "https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/operations_split_edge.jl"
```

# Triangulation Operations
## Edge Splitting

Sometimes you have a point `r` that you know is on, or is very
close to being on, an edge `(i, j)`. In this tutorial, we show how
the `split_edge!` function can be used for putting a point on this
edge. First, let us consider the following triangulation.

````@example operations_split_edge
using DelaunayTriangulation
using CairoMakie

points = [(0.0, 0.0), (0.0, 4.0), (2.0, 3.0), (-2.0, 3.0),
    (-2.0, 7.0), (3.0, 6.0), (2.0, -2.0), (-4.0, 1.0),
    (1.0, 5.0)]
p = (0.0, 3.0)
tri = triangulate(points)
fig, ax, sc = triplot(tri)
scatter!(ax, [p], markersize=14)
fig
````

We want to add the blue point onto the edge shown, which is `(1, 2)`. To do this,
we can use the function `split_edge!`.

````@example operations_split_edge
push!(points, p)
r = length(points)
i, j = 1, 2
split_edge!(tri, i, j, r)
fig, ax, sc = triplot(tri)
fig
````

Notice that this has only split the edge in one direction. This is because the
edges in this case are treated as being oriented. To split the edge in the other
direction, we simply swap the indices.

````@example operations_split_edge
split_edge!(tri, j, i, r)
fig, ax, sc = triplot(tri)
fig
````

If you also want to restore the Delaunay property of the triangulation
following this splitting, you need to use `legalise_edge!`. In this
example, though, there are no illegal edges. If there were, we would use

````@example operations_split_edge
k = get_adjacent(tri, i, r) # get_adjacent(tri, i, j) before split_edge!(tri, i, j)
legalise_edge!(tri, j, k, r)
legalise_edge!(tri, k, i, r)
k = get_adjacent(tri, j, r) # get_adjacent(tri, j, i) before split_edge!(tri, j, i)
legalise_edge!(tri, i, k, r)
legalise_edge!(tri, k, j, r)
````

These steps, in particular the steps of splitting both sides of the edge and then
legalising, are also implemented in `DelaunayTriangulation.complete_split_edge_and_legalise!`.
## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/new-docs/docs/src/literate_tutorials/operations_split_edge.jl).

```julia
using DelaunayTriangulation
using CairoMakie

points = [(0.0, 0.0), (0.0, 4.0), (2.0, 3.0), (-2.0, 3.0),
    (-2.0, 7.0), (3.0, 6.0), (2.0, -2.0), (-4.0, 1.0),
    (1.0, 5.0)]
p = (0.0, 3.0)
tri = triangulate(points)
fig, ax, sc = triplot(tri)
scatter!(ax, [p], markersize=14)
fig

push!(points, p)
r = length(points)
i, j = 1, 2
split_edge!(tri, i, j, r)
fig, ax, sc = triplot(tri)
fig

split_edge!(tri, j, i, r)
fig, ax, sc = triplot(tri)
fig

k = get_adjacent(tri, i, r) # get_adjacent(tri, i, j) before split_edge!(tri, i, j)
legalise_edge!(tri, j, k, r)
legalise_edge!(tri, k, i, r)
k = get_adjacent(tri, j, r) # get_adjacent(tri, j, i) before split_edge!(tri, j, i)
legalise_edge!(tri, i, k, r)
legalise_edge!(tri, k, j, r)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

