```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/operations_split_triangle.jl"
```

# Triangulation Operations
## Triangle Splitting

As we briefly discussed in the [`legalise_edge!` tutorial](operations_legalise_edge.md),
it is sometimes useful to put a point inside of an existing triangle, and then to
connect all the vertices of that triangle to that new point, thereby splitting
that triangle into three new triangles. This is called *triangle splitting*.

Let us give an example.

````@example operations_split_triangle
using DelaunayTriangulation
using CairoMakie

points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
p = (0.2, 0.5)
tri = triangulate(points)
fig, ax, sc = triplot(tri)
scatter!(ax, [p], markersize = 14)
fig
````

The blue point shows the point we want to add into the triangulation
using [`split_triangle!`](@ref). To use this, we provide (1) the index of the point
in `points` and (2) the triangle that the point is in. The index of the point
will be `4` after pushing `p` into `points`, and in this simple example
the triangle that `p` is in is `(1, 2, 3)`.

````@example operations_split_triangle
push!(points, p)
r = length(points)
i, j, k = 1, 2, 3
split_triangle!(tri, i, j, k, r)
fig, ax, sc = triplot(tri)
fig
````

See the [`legalise_edge!` tutorial](operations_legalise_edge.md) for more discussion
about restoring the Delaunay property of the triangulation after using
`split_triangle!`.

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/operations_split_triangle.jl).

```julia
using DelaunayTriangulation
using CairoMakie

points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
p = (0.2, 0.5)
tri = triangulate(points)
fig, ax, sc = triplot(tri)
scatter!(ax, [p], markersize = 14)
fig

push!(points, p)
r = length(points)
i, j, k = 1, 2, 3
split_triangle!(tri, i, j, k, r)
fig, ax, sc = triplot(tri)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

