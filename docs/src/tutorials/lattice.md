```@meta
EditURL = "https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/lattice.jl"
```

# Triangulating Rectangular Regions

In this tutorial, we show how you can easily
triangulate rectangular regions of the form
$[a, b] \times [c, d]$. Rather than using `triangulate`,
you can use `triangulate_rectangle` for this purpose. To start,
we give a simple example

````@example lattice
using DelaunayTriangulation
using CairoMakie

a, b, c, d = 0.0, 2.0, 0.0, 10.0
nx, ny = 10, 25
tri = triangulate_rectangle(a, b, c, d, nx, ny)
fig, ax, sc = triplot(tri)
fig
````

This can be much faster than if we just construct the points in the
lattice manually and `triangulate` those. Here's a comparison of the times.

````@example lattice
using BenchmarkTools
points = get_points(tri)
@benchmark triangulate($points; randomise=$false) # randomise=false because points are already in lattice order, i.e. spatially sorted
````

````@example lattice
@benchmark triangulate_rectangle($a, $b, $c, $d, $nx, $ny)
````

This difference would be more pronounced for larger `nx, ny`.

Note that the output of `triangulate_rectangle` treats the boundary
as a constrained boundary:

````@example lattice
get_boundary_nodes(tri)
````

This boundary is split into four separate segments, one for each
side of the rectangle. If you would prefer to keep the boundary as one
contiguous segment, use `single_boundary=true`. Moreover, note that
this `tri` has ghost triangles:

````@example lattice
tri
````

You can opt into not having these by using `add_ghost_triangles=false`.

````@example lattice
tri = triangulate_rectangle(a, b, c, d, nx, ny; single_boundary=true, add_ghost_triangles=false)
tri
````

````@example lattice
get_boundary_nodes(tri)
````

````@example lattice
DelaunayTriangulation.has_ghost_triangles(tri)
````

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/new-docs/docs/src/literate_tutorials/lattice.jl).

```julia
using DelaunayTriangulation
using CairoMakie

a, b, c, d = 0.0, 2.0, 0.0, 10.0
nx, ny = 10, 25
tri = triangulate_rectangle(a, b, c, d, nx, ny)
fig, ax, sc = triplot(tri)
fig

using BenchmarkTools
points = get_points(tri)
@benchmark triangulate($points; randomise=$false) # randomise=false because points are already in lattice order, i.e. spatially sorted

@benchmark triangulate_rectangle($a, $b, $c, $d, $nx, $ny)

get_boundary_nodes(tri)

tri

tri = triangulate_rectangle(a, b, c, d, nx, ny; single_boundary=true, add_ghost_triangles=false)
tri

get_boundary_nodes(tri)

DelaunayTriangulation.has_ghost_triangles(tri)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

