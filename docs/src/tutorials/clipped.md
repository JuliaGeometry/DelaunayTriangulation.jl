```@meta
EditURL = "https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/clipped.jl"
```

# Clipped Voronoi Tessellations
## Clipping to the Convex Hull

One issue that may arise when dealing with Voronoi tessellations is the
presence of unbounded polygons occurring on the boundary. One way to deal with this
is to clip polygons to the convex hull of the tessellation. (Arbitrary clipping boundaries
are on the to-do list, but they are not yet implemented - see [Issue #48](https://github.com/DanielVandH/DelaunayTriangulation.jl/issues/48).)

In the example below, we clip the tessellation to the convex hull of the point set by setting
the second argument of `voronoi` to `true`.

````@example clipped
using DelaunayTriangulation
using CairoMakie
using StableRNGs
using ReferenceTests #src

rng = StableRNG(123)
points = randn(rng, 2, 50)
tri = triangulate(points; rng)
vorn = voronoi(tri)
````

````@example clipped
clipped_vorn = voronoi(tri, true)
````

Note that the clipping has put more polygon vertices in. We compare
the clipped tessellations below.

````@example clipped
fig = Figure()
ax1 = Axis(fig[1, 1], title="Unclipped", width=600, height=400)
ax2 = Axis(fig[1, 2], title="Clipped", width=600, height=400)
voronoiplot!(ax1, vorn, show_generators=false, colormap=:matter, strokewidth=4)
voronoiplot!(ax2, clipped_vorn, show_generators=false, colormap=:matter, strokewidth=4)
resize_to_layout!(fig)
fig
````

As you can see, the unbounded polygons, and any polygons that included points
outside of the convex hull, have now been clipped to the convex hull.
## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/new-docs/docs/src/literate_tutorials/clipped.jl).

```julia
using DelaunayTriangulation
using CairoMakie
using StableRNGs
using ReferenceTests #src

rng = StableRNG(123)
points = randn(rng, 2, 50)
tri = triangulate(points; rng)
vorn = voronoi(tri)

clipped_vorn = voronoi(tri, true)

fig = Figure()
ax1 = Axis(fig[1, 1], title="Unclipped", width=600, height=400)
ax2 = Axis(fig[1, 2], title="Clipped", width=600, height=400)
voronoiplot!(ax1, vorn, show_generators=false, colormap=:matter, strokewidth=4)
voronoiplot!(ax2, clipped_vorn, show_generators=false, colormap=:matter, strokewidth=4)
resize_to_layout!(fig)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

