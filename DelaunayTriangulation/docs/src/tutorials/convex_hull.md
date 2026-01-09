```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/convex_hull.jl"
```

# Convex Hulls

In this tutorial, we show how we can compute convex hulls.
The main method for this is to compute them directly from the
triangulation, noting that for an unconstrained triangulation the boundary
is the convex hull. We also provide a method for computing the convex hull
from the point set directly, using the [monotone chain algorithm](https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain).

Let us first demonstrate how we compute convex hulls from a triangulation.
This is only possible from triangulations without a constrained boundary

````@example convex_hull
using DelaunayTriangulation
using CairoMakie
using StableRNGs
````

We construct a triangulation of a point set.

````@example convex_hull
rng = StableRNG(123)
points = randn(rng, 2, 250)
tri = triangulate(points; rng)
````

To get the convex hull, just use [`get_convex_hull`](@ref):

````@example convex_hull
get_convex_hull(tri)
````

You can also obtain the vertices directly from
[`get_convex_hull_vertices`](@ref):

````@example convex_hull
get_convex_hull_vertices(tri)
````

The vertices refer to indices of points in `points`, and they are
given in counter-clockwise order. To obtain the convex hull directly
from the points without constructing the triangulation, use [`convex_hull`](@ref):

````@example convex_hull
ch = convex_hull(points)
ch_points = [get_point(tri, i) for i in DelaunayTriangulation.get_vertices(ch)]
fig, ax, sc = lines(ch_points, color = :red, linewidth = 4)
scatter!(ax, points)
fig
````

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/convex_hull.jl).

```julia
using DelaunayTriangulation
using CairoMakie
using StableRNGs

rng = StableRNG(123)
points = randn(rng, 2, 250)
tri = triangulate(points; rng)

get_convex_hull(tri)

get_convex_hull_vertices(tri)

ch = convex_hull(points)
ch_points = [get_point(tri, i) for i in DelaunayTriangulation.get_vertices(ch)]
fig, ax, sc = lines(ch_points, color = :red, linewidth = 4)
scatter!(ax, points)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

