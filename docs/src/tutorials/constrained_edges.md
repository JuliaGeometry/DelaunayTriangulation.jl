```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/constrained_edges.jl"
```

# Constrained Triangulations
## Constrained Segments

In this tutorial, we introduce constrained triangulations,
starting with the simple case of only having constrained segments,
meaning edges that are forced to be in the final triangulation.
To start, let us load in the packages we will need.

````@example constrained_edges
using DelaunayTriangulation
using CairoMakie
````

We consider triangulating the following set of points:

````@example constrained_edges
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
````

To now define the segments, we define:

````@example constrained_edges
C = Set([(2, 1), (2, 11), (2, 7), (2, 5)])
````

With this notation, each `Tuple` is an individual
edge to include in the triangulation, with
`(i, j)` meaning the edge connecting the points
`pts[i]` and `pts[j]` together. Let us now
make the triangulation, comparing it to its
unconstrained counterpart.

````@example constrained_edges
tri = triangulate(pts)
````

````@example constrained_edges
cons_tri = triangulate(pts; segments=C)
````

````@example constrained_edges
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel="x", ylabel=L"y",
    title="(a): Unconstrained", titlealign=:left,
    width=300, height=300)
ax2 = Axis(fig[1, 2], xlabel="x", ylabel=L"y",
    title="(b): Unconstrained", titlealign=:left,
    width=300, height=300)
triplot!(ax1, tri)
triplot!(ax2, cons_tri, show_constrained_edges = true)
resize_to_layout!(fig)
fig
````

As you can see, the constrained edges in magenta
have now been included in the triangulation in (b),
whereas in (a) most were previously not included.

You can view the constrained edges by using

````@example constrained_edges
get_interior_segments(cons_tri)
````

There is also a function [`get_all_segments`](@ref),
which in this case is the same as [`get_interior_segments`](@ref),
but in the case of a triangulation with constrained
boundaries, it will also include the boundary segments
whereas `get_interior_segments` will not; this is
demonstrated in the later tutorials.
## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/constrained_edges.jl).

```julia
using DelaunayTriangulation
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

C = Set([(2, 1), (2, 11), (2, 7), (2, 5)])

tri = triangulate(pts)

cons_tri = triangulate(pts; segments=C)

fig = Figure()
ax1 = Axis(fig[1, 1], xlabel="x", ylabel=L"y",
    title="(a): Unconstrained", titlealign=:left,
    width=300, height=300)
ax2 = Axis(fig[1, 2], xlabel="x", ylabel=L"y",
    title="(b): Unconstrained", titlealign=:left,
    width=300, height=300)
triplot!(ax1, tri)
triplot!(ax2, cons_tri, show_constrained_edges = true)
resize_to_layout!(fig)
fig

get_interior_segments(cons_tri)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

