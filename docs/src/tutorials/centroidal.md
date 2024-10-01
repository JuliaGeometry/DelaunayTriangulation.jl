```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/centroidal.jl"
```

# Centroidal Voronoi Tessellations
In this tutorial, we demonstrate how we can compute
centroidal Voronoi tessellations. These are tessellations
for which each generator is approximately the centroid of its
associated polygon, except for boundary generators which are fixed
in place. This method is only applicable to clipped tessellations.

We give a simple example. First, we compute the clipped tessellation
of a point set.

````@example centroidal
using DelaunayTriangulation
using CairoMakie
using StableRNGs
rng = StableRNG(123)
points = 25randn(rng, 2, 500)
tri = triangulate(points; rng)
vorn = voronoi(tri, clip = true)
````

To now compute the centroidal tessellation, use [`centroidal_smooth`](@ref). (
If you want to straight from a triangulation to a centroidal tessellation, you
can also just use `smooth_vorn = voronoi(tri, clip = true, smooth = true)`.)

````@example centroidal
smooth_vorn = centroidal_smooth(vorn; rng)
````

Let us now compare the two tessellations.

````@example centroidal
fig = Figure()
ax1 = Axis(fig[1, 1], title = "Original", width = 600, height = 400)
ax2 = Axis(fig[1, 2], title = "Smoothed", width = 600, height = 400)
voronoiplot!(ax1, vorn, colormap = :matter, strokewidth = 2)
voronoiplot!(ax2, smooth_vorn, colormap = :matter, strokewidth = 2)
resize_to_layout!(fig)
fig
````

As you can see, the tiles are all reasonably uniform, and their generators
do look to be near the centroid of their corresponding tile. Note that
this function `centroidal_smooth` is iterative, and you can control the iteration limits
and the tolerance of the iterations (based on the maximum displacement of any generator at
each iteration) using the `maxiters` and `tol` keyword arguments.

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/centroidal.jl).

```julia
using DelaunayTriangulation
using CairoMakie
using StableRNGs
rng = StableRNG(123)
points = 25randn(rng, 2, 500)
tri = triangulate(points; rng)
vorn = voronoi(tri, clip = true)

smooth_vorn = centroidal_smooth(vorn; rng)

fig = Figure()
ax1 = Axis(fig[1, 1], title = "Original", width = 600, height = 400)
ax2 = Axis(fig[1, 2], title = "Smoothed", width = 600, height = 400)
voronoiplot!(ax1, vorn, colormap = :matter, strokewidth = 2)
voronoiplot!(ax2, smooth_vorn, colormap = :matter, strokewidth = 2)
resize_to_layout!(fig)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

