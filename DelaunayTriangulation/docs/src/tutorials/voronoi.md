```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/voronoi.jl"
```

# Voronoi Tessellations

In this tutorial, we demonstrate how we can construct Voronoi
tessellations and work with them. Voronoi tessellations are built from
a dual Delaunay triangulation using [`voronoi`](@ref). To start, let us
load in the packages.

````@example voronoi
using DelaunayTriangulation
using CairoMakie
using StableRNGs
````

We build the tessellation by constructing the triangulation,
and then passing that triangulation into `voronoi`.

````@example voronoi
points = [
    (-3.0, 7.0), (1.0, 6.0), (-1.0, 3.0),
    (-2.0, 4.0), (3.0, -2.0), (5.0, 5.0),
    (-4.0, -3.0), (3.0, 8.0),
]
rng = StableRNG(123)
tri = triangulate(points; rng)
vorn = voronoi(tri)
````

To visualise the tessellation, you can use `voronoiplot`. Here,
we also compare the tessellation with its dual triangulation.

````@example voronoi
fig, ax, sc = voronoiplot(vorn, markersize = 13, colormap = :matter, strokecolor = :white, strokewidth = 5)
triplot!(ax, tri)
fig
````

The polygons each correspond to a *generator*, which is the black point
inside it coming from `points`, i.e. the vertices of the triangulation.
The polygons are all convex. Note also that the unbounded polygons,
those coming from generators on the convex hull of the point set,
are clipped in the above plot but in reality they go on to infinity.
In the [clipping to rectangular regions tutorial](clipped_rectangle.md),
we discuss how these polygons are clipped to a rectangle in more detail.

Let us now demonstrate in more detail how we can work with `vorn`.

## Iterating over generators
The generators are stored in a `Dict`, mapping vertices to coordinates;
this is a `Dict` rather than a `Vector` because the original triangulation
may not contain all the points, noting that [`get_generators(vorn)`](@ref DelaunayTriangulation.get_generators) is just a repackaged
version of `get_points(tri)`. (A separate field is needed so that clipped and centroidal
tessellations can add new generators without affecting the points in `tri`.)
These generators can be accessed using `get_generators(vorn)`:

````@example voronoi
DelaunayTriangulation.get_generators(vorn)
````

It is preferred that you use [`each_generator(vorn)`](@ref each_generator), though, which is an iterator
over the generators:

````@example voronoi
each_generator(vorn)
````

Note that this is just the keys of the above `Dict`. To access specific generators, you use [`get_generator`](@ref). For example,

````@example voronoi
get_generator(vorn, 3)
````

To give an example, here is how we can compute the average
generator position.

````@example voronoi
function average_generator(vorn)
    cx, cy = 0.0, 0.0
    for i in each_generator(vorn)
        p = get_generator(vorn, i)
        px, py = getxy(p)
        cx += px
        cy += py
    end
    n = num_polygons(vorn) # same as DelaunayTriangulation.num_generators(vorn)
    cx /= n
    cy /= n
    return cx, cy
end
cx, cy = average_generator(vorn)
````

## Iterating over polygon and polygon vertices
You can also look at the individual polygons, and all their vertices. These polygons
and their vertices are stored as below:

````@example voronoi
DelaunayTriangulation.get_polygons(vorn)
````

````@example voronoi
DelaunayTriangulation.get_polygon_points(vorn)
````

You should not work with these fields directly, though. If you want to
look at a specific polygon, you should use [`get_polygon`](@ref). For example,

````@example voronoi
get_polygon(vorn, 1)
````

This is a `Vector` of the vertices of the polygon, in counter-clockwise order,
and such that the first and last vertices are the same. The vertices refer to points in
the `polygon_points` field, which you could then obtain using [`get_polygon_point`](@ref). For example,
the first vertex corresponds to the coordinates:

````@example voronoi
get_polygon_point(vorn, 1)
````

In `get_polygon(vorn, 1)`, notice that there are two negative indices. These negative indices
correspond to vertices out at infinity; their actual values do not matter, just that they
are negative. Thus, this polygon is actually an unbounded polygon. This can be checked in two ways:

````@example voronoi
1 ∈ DelaunayTriangulation.get_unbounded_polygons(vorn)
````

````@example voronoi
get_area(vorn, 1) == Inf
````

There are several ways that you can iterate over the polygons. If you want each
polygon, then you can use [`each_polygon`](@ref):

````@example voronoi
each_polygon(vorn)
````

This is an iterator over all the polygon vertices, but you do not get the
associated polygon index. If you want an iterator that is over the polygon indices,
you use [`each_polygon_index`](@ref):

````@example voronoi
each_polygon_index(vorn)
````

If you did want to have both the indices and the vertices together, you can use `zip`
on these two iterators. This would be the same as iterating over the internal `polygons`
field, but this could be subject to change in the future. To get an iterator over the polygon
vertices rather than caring about a specific polygon, you use [`each_polygon_vertex`](@ref):

````@example voronoi
each_polygon_vertex(vorn)
````

This is just an iterator of the indices to pass into `get_polygon_point`. You can also
query the number of polygons and polygon vertices as follows:

````@example voronoi
num_polygons(vorn)
````

````@example voronoi
num_polygon_vertices(vorn)
````

To give an example of how we might use these iterators, here we compute the
area of all polygons.

````@example voronoi
function get_polygon_area(vorn, i)
    i ∈ DelaunayTriangulation.get_unbounded_polygons(vorn) && return Inf
    area = 0.0
    vertices = get_polygon(vorn, i)
    vⱼ = vertices[begin]
    pⱼ = get_polygon_point(vorn, vⱼ)
    xⱼ, yⱼ = getxy(pⱼ)
    for j in (firstindex(vertices) + 1):lastindex(vertices) # same as 2:length(vertices)
        vⱼ₊₁ = vertices[j]
        pⱼ₊₁ = get_polygon_point(vorn, vⱼ₊₁)
        xⱼ₊₁, yⱼ₊₁ = getxy(pⱼ₊₁)
        area += xⱼ * yⱼ₊₁ - xⱼ₊₁ * yⱼ
        vⱼ, pⱼ, xⱼ, yⱼ = vⱼ₊₁, pⱼ₊₁, xⱼ₊₁, yⱼ₊₁
    end
    return area / 2
end
function get_polygon_areas(vorn)
    areas = zeros(num_polygons(vorn))
    for i in each_polygon_index(vorn)
        areas[i] = get_polygon_area(vorn, i)
    end
    return areas
end
vorn_areas = get_polygon_areas(vorn)
````

Note that this could have also been obtained using [`get_area`](@ref):

````@example voronoi
function direct_polygon_areas(vorn)
    areas = zeros(num_polygons(vorn))
    for i in each_polygon_index(vorn)
        areas[i] = get_area(vorn, i)
    end
    return areas
end
vorn_areas ≈ direct_polygon_areas(vorn)
````

Moreover, note that the following is false:

````@example voronoi
vorn_areas ≈ [get_area(vorn, i) for i in each_polygon_index(vorn)]
````

because `each_polygon_index` does not return the polygons in a sorted order.

Before we move on, we emphasise that it is not guaranteed that the values of the vertices
for a given polygon are the same if you were to recompute the tessellation. The actual coordinates
will all be the same, but they just might correspond to different vertex values. The indices
of the polygons will be the same, though, as they are derived from the point indices.

## Getting polygons adjacent to an edge
Given an edge in the tessellation, you can use [`get_adjacent`](@ref) to get the polygon
that it is a part of (taking care of order). For example,

````@example voronoi
get_adjacent(vorn, 1, 8)
````

means that the edge `(1, 8)` belongs to the third polygon, as we can easily verify:

````@example voronoi
get_polygon(vorn, 3)
````

see that `(1, 8)` at the end. The order is important here, since

````@example voronoi
get_adjacent(vorn, 8, 1)
````

means that the edge `(8, 1)` belongs to the seventh polygon:

````@example voronoi
get_polygon(vorn, 7)
````

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/voronoi.jl).

```julia
using DelaunayTriangulation
using CairoMakie
using StableRNGs

points = [
    (-3.0, 7.0), (1.0, 6.0), (-1.0, 3.0),
    (-2.0, 4.0), (3.0, -2.0), (5.0, 5.0),
    (-4.0, -3.0), (3.0, 8.0),
]
rng = StableRNG(123)
tri = triangulate(points; rng)
vorn = voronoi(tri)

fig, ax, sc = voronoiplot(vorn, markersize = 13, colormap = :matter, strokecolor = :white, strokewidth = 5)
triplot!(ax, tri)
fig

DelaunayTriangulation.get_generators(vorn)

each_generator(vorn)

get_generator(vorn, 3)

function average_generator(vorn)
    cx, cy = 0.0, 0.0
    for i in each_generator(vorn)
        p = get_generator(vorn, i)
        px, py = getxy(p)
        cx += px
        cy += py
    end
    n = num_polygons(vorn) # same as DelaunayTriangulation.num_generators(vorn)
    cx /= n
    cy /= n
    return cx, cy
end
cx, cy = average_generator(vorn)

DelaunayTriangulation.get_polygons(vorn)

DelaunayTriangulation.get_polygon_points(vorn)

get_polygon(vorn, 1)

get_polygon_point(vorn, 1)

1 ∈ DelaunayTriangulation.get_unbounded_polygons(vorn)

get_area(vorn, 1) == Inf

each_polygon(vorn)

each_polygon_index(vorn)

each_polygon_vertex(vorn)

num_polygons(vorn)

num_polygon_vertices(vorn)

function get_polygon_area(vorn, i)
    i ∈ DelaunayTriangulation.get_unbounded_polygons(vorn) && return Inf
    area = 0.0
    vertices = get_polygon(vorn, i)
    vⱼ = vertices[begin]
    pⱼ = get_polygon_point(vorn, vⱼ)
    xⱼ, yⱼ = getxy(pⱼ)
    for j in (firstindex(vertices) + 1):lastindex(vertices) # same as 2:length(vertices)
        vⱼ₊₁ = vertices[j]
        pⱼ₊₁ = get_polygon_point(vorn, vⱼ₊₁)
        xⱼ₊₁, yⱼ₊₁ = getxy(pⱼ₊₁)
        area += xⱼ * yⱼ₊₁ - xⱼ₊₁ * yⱼ
        vⱼ, pⱼ, xⱼ, yⱼ = vⱼ₊₁, pⱼ₊₁, xⱼ₊₁, yⱼ₊₁
    end
    return area / 2
end
function get_polygon_areas(vorn)
    areas = zeros(num_polygons(vorn))
    for i in each_polygon_index(vorn)
        areas[i] = get_polygon_area(vorn, i)
    end
    return areas
end
vorn_areas = get_polygon_areas(vorn)

function direct_polygon_areas(vorn)
    areas = zeros(num_polygons(vorn))
    for i in each_polygon_index(vorn)
        areas[i] = get_area(vorn, i)
    end
    return areas
end
vorn_areas ≈ direct_polygon_areas(vorn)

vorn_areas ≈ [get_area(vorn, i) for i in each_polygon_index(vorn)]

get_adjacent(vorn, 1, 8)

get_polygon(vorn, 3)

get_adjacent(vorn, 8, 1)

get_polygon(vorn, 7)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

