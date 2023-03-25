```@meta
CurrentModule = DelaunayTriangulation
```

# Convex Polygons 

We have an algorithm available for computing the Delaunay triangulation of a convex polygon:

```@docs 
triangulate_convex 
```

The function takes as input some set of points $\mathcal P$ and a corresponding sequence of indices $\{v_1, \ldots, v_n\}$, in counter-clockwise order, defining some convex polygon and returns the polygon's Delaunay triangulation. Note that the keyword arguments here differ from `triangulate`. Here is an example of `triangulate_convex` in action.

```julia
using DelaunayTriangulation 
using CairoMakie
p1 = [10.0, 12.0]
p2 = [7.0, 11.0]
p3 = [8.0, 6.0]
p4 = [10.0, 3.0]
p5 = [14.0, 5.0]
p6 = [15.0, 10.0]
p7 = [13.0, 12.0]
pts = [p1, p2, p3, p4, p5, p6, p7]
S = collect(1:7)
tri = triangulate_convex(pts, S)
fig, ax, sc = triplot(tri; plot_convex_hull=false)
```

```@raw html
<figure>
    <img src='../figs/convex_triangulation_example.png', alt='Triangulation'><br>
</figure>
```