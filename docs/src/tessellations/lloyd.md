```@meta
CurrentModule = DelaunayTriangulation
```

# Centroidal Voronoi Tessellation

Now we show how to compute a centroidal Voronoi tessellation, i.e. a Voronoi tessellation where the generators are at the centroids of the polygons. This method only works on clipped tessellations, and generators on the boundary are fixed in place. The function that handles this computation is `centroidal_smooth`:

```@docs 
centroidal_smooth
```

The algorithm for computing the tessellation is very simple, making use of Lloyd's algorithm:

1. First, compute the tessellation.
2. Compute the centroids of each polygon and move the generator to the centroid.
3. If the maximum displacement of a generator is below some tolerance, or if the maximum number of iterations has been reached, exit. Otherwise, go back to step 1 with the new generators.

## Example

Let's now give an example.

```julia
pts = 25randn(2, 500)
tri = triangulate(pts)
vorn = voronoi(tri, true)
smooth_vorn = centroidal_smooth(vorn)

cmap = Makie.cgrad(:jet)
colors = get_polygon_colors(vorn, cmap)
fig = Figure()
ax = Axis(fig[1, 1], aspect=1)
voronoiplot!(ax, vorn, strokecolor=:red, strokewidth=0.2, polygon_color=colorsmarkersize=4)
xlims!(ax, -100, 100)
ylims!(ax, -100, 100)

ax = Axis(fig[1, 2], aspect=1)
voronoiplot!(ax, smooth_vorn, strokecolor=:red, strokewidth=0.2polygon_color=colors, markersize=4)
xlims!(ax, -100, 100)
ylims!(ax, -100, 100)
```

```@raw html
<figure>
    <img src='../figs/lloyd.png', alt='Centroidal Voronoi Tessellation'><br>
</figure>
```

See that we indeed have smoothed out the points significantly. Out of interest, here's the difference between the triangulations.

```@raw html
<figure>
    <img src='../figs/lloyd_tri.png', alt='Centroidal Voronoi Tessellation Underlying Triangulation'><br>
</figure>
```