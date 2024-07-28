```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/clipped_rectangle.jl"
```

# Clipped Voronoi Tessellations
## Clipping to a Rectangle

In the previous tutorial, we demonstrated how we can clip to
the convex hull of the point set. However, it is often useful to clip
to a rectangle, for example if you want to clip to a region of interest
in a simulation. We do not yet support this within [`voronoi`](@ref) itself,
but we provide the function [`get_polygon_coordinates`](@ref) for this (this is what
`voronoiplot` uses to plot inside a bounding box).

Let us now demonstrate. First, we construct a tessellation of
some example point set.

````@example clipped_rectangle
using DelaunayTriangulation
using CairoMakie
A = (-3.0, 7.0)
B = (1.0, 6.0)
C = (-1.0, 3.0)
D = (-2.0, 4.0)
E = (3.0, -2.0)
F = (5.0, 5.0)
G = (-4.0, -3.0)
H = (3.0, 8.0)
points = [A, B, C, D, E, F, G, H]
tri = triangulate(points)
vorn = voronoi(tri)
````

Let us show the tessellation, and the rectangle we want to clip
the tessellation to.

````@example clipped_rectangle
fig, ax, sc = voronoiplot(vorn)
a, b, c, d = -2.0, 3.0, 0.0, 7.0
lines!(ax, [(a, c), (b, c), (b, d), (a, d), (a, c)], color = :black, linewidth = 4)
fig
````

To apply this clipping, we need to provide a bounding box of the form
`(xmin, xmax, ymin, ymax)`. Here, we will use

````@example clipped_rectangle
bounding_box = (a, b, c, d)
````

You can obtain some reasonable defaults for this bounding box using
[DelaunayTriangulation.polygon_bounds(vorn)](@ref polygon_bounds).
The coordinates for each polygon clipped to this box can be obtained as follows.

````@example clipped_rectangle
clipped_coords = Vector{Vector{NTuple{2, Float64}}}(undef, num_polygons(vorn))
for i in each_polygon_index(vorn)
    clipped_coords[i] = get_polygon_coordinates(vorn, i, bounding_box)
end
clipped_coords
````

Now let's plot these.

````@example clipped_rectangle
fig, ax, sc = poly(clipped_coords, color = :white, strokewidth = 4)
fig
````

As we can see, the polygons have been clipped to the rectangle.
Note that if you just want this for plotting, you can also call `voronoiplot` with the
`bounding_box` keyword argument.

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/clipped_rectangle.jl).

```julia
using DelaunayTriangulation
using CairoMakie
A = (-3.0, 7.0)
B = (1.0, 6.0)
C = (-1.0, 3.0)
D = (-2.0, 4.0)
E = (3.0, -2.0)
F = (5.0, 5.0)
G = (-4.0, -3.0)
H = (3.0, 8.0)
points = [A, B, C, D, E, F, G, H]
tri = triangulate(points)
vorn = voronoi(tri)

fig, ax, sc = voronoiplot(vorn)
a, b, c, d = -2.0, 3.0, 0.0, 7.0
lines!(ax, [(a, c), (b, c), (b, d), (a, d), (a, c)], color = :black, linewidth = 4)
fig

bounding_box = (a, b, c, d)

clipped_coords = Vector{Vector{NTuple{2, Float64}}}(undef, num_polygons(vorn))
for i in each_polygon_index(vorn)
    clipped_coords[i] = get_polygon_coordinates(vorn, i, bounding_box)
end
clipped_coords

fig, ax, sc = poly(clipped_coords, color = :white, strokewidth = 4)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

