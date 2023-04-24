```@meta
CurrentModule = DelaunayTriangulation
```

# Voronoi Tessellations 

We provide support for constructing Voronoi tessellations, the dual graph to a Delaunay triangulation (strictly true only for unconstrained triangulations). The main function for this is `voronoi`:

```@docs 
voronoi
```

The algorithm is really quite simple:

1. Construct the Delaunay triangulation.
2. For each point in the triangulation, connect the circumcenters of all triangles with that point as a vertex.

(Omitting, of course, some details like duplicate circumcenters from cocircular susets -- as with many things in computational geometry, anything is "simple" until you think about reality.) Note that extra care is taken for the unbounded edges, which we simply represent using negative indices. The negative indices we use have no real meaning, but they are unique so that we can uniquely identify for any unbounded ray the associated ghost triangle (shown by example later).

## Example

Let us give an example. The first step in any tessellation is the construction of the triangulation.

```julia 
using DelaunayTriangulation, CairoMakie
pts = rand(2, 25)
tri = triangulate(pts)
```

To now construct the tessellation, use `voronoi`.

```julia
vorn = voronoi(tri)
```

We can visualise the diagram using `voronoiplot`. This function uses `get_polygon_coordinates` to get a finite representation of the unbounded edges, so we need to use `xlims` and `ylims` to get a reasonable visualisation of the unbounded edges; a better way would be to somehow use observables to track the axis so that the edges truly are infinite, like in `hlines` and `vlines`, but I don't know how to do that / haven't got around to figuring out the details yet (without depending on Makie.jl). We also use `get_polygon_colors` to get a nice spread of colours to use.

```julia
cmap = Makie.cgrad(:jet)
colors = get_polygon_colors(vorn, cmap)
fig, ax, sc = voronoiplot(vorn, strokecolor=:red, polygon_color=colors)
triplot!(ax, tri, strokewidth=0.2, strokecolor=(:black, 0.4))
xlims!(ax, -1 / 2, 3 / 2)
ylims!(ax, 0, 3 / 2)`
```

```@raw html
<figure>
    <img src='../figs/unbounded.png', alt='Voronoi Tessellation'><br>
</figure>
```

(See also `polygon_bounds` to get appropriate bounds programatically.) The function that handles chopping to a box is `get_polygon_coordinates`, which primarily relies on `intersection_of_ray_with_boundary`:

```@docs 
intersection_of_ray_with_boundary 
```

We will see later how to clip these polygons to the convex hull show in the red dashed line.

Let us now actually show the data structure. There are several fields:

```julia-repl
julia> vorn.
adjacent                  circumcenter_to_triangle  generators                polygons                  triangulation
boundary_polygons         cocircular_circumcenters  polygon_points            triangle_to_circumcenter  unbounded_polygon
```

The `adjacent` field is similar to the one we deal with for triangulations, except edges now map to their associated polygon. For example, we have

```julia-repl
julia> get_polygon(vorn, 1)
7-element Vector{Int64}:
 23
 17
 -8
 -6
 32
 34
 23

julia> get_adjacent(vorn, 23, 17) == get_adjacent(vorn, 17, -8) == 1
true
```

The `circumcenter_to_triangle` field allows us to go from the index of a circumcenter to its associated triangle. For instance, if we look at the vertices of the first polygon above, we may want to know where this vertex `34` came from:

```julia-repl
julia> DelaunayTriangulation.get_circumcenter_to_triangle(vorn, 34)
(10, 4, 1)
```

Thus, it came from the triangle `(10, 4, 1)`. Note that the triangles used in this structure are always sorted so that the minimum index is last (while preserving the counter-clockwise orientation).

The `generators` field is a repackaged version of `tri.points`:

```julia-repl
julia> DelaunayTriangulation.get_generators(vorn)
Dict{Int64, Tuple{Float64, Float64}} with 25 entries:
  5  => (0.684537, 0.618674)
  16 => (0.585988, 0.554157)
  20 => (0.43467, 0.161223)
  12 => (0.00191717, 0.70299)
  24 => (0.562014, 0.246805)
  8  => (0.446736, 0.717801)
  17 => (0.00141105, 0.639474)
  23 => (0.953932, 0.352213)
  1  => (0.234991, 0.0601079)
  22 => (0.559084, 0.538079)
  19 => (0.650726, 0.439868)
  6  => (0.966998, 0.627038)
  11 => (0.464044, 0.594492)
  9  => (0.512324, 0.0257408)
  14 => (0.559852, 0.55035)
  3  => (0.139261, 0.472539)
  7  => (0.560707, 0.767171)
  25 => (0.0731719, 0.791413)
  4  => (0.389615, 0.253273)
  15 => (0.597794, 0.761001)
  ⋮  => 
```

Here, we need to use a `Dict` to represent the points since the original triangulation may not include all points. With this, we always know that `get_generator(vorn, i) == get_point(tri, i)`. If we want the coordinates of the `i`th generator, we could use 

```julia-repl
julia> get_generator(vorn, 16)
(0.5859881889826812, 0.5541565595408635)

julia> get_generator(vorn, 7, 22, 1)
((0.560707474714347, 0.7671706685373769), (0.5590841028546631, 0.5380792727814472), (0.2349910552466512, 0.060107935323076456))
```

The `polygons` field is a `Dict` that maps the indices of the generators to their polygon:

```julia-repl
julia> DelaunayTriangulation.get_polygons(vorn)
Dict{Int64, Vector{Int64}} with 25 entries:
  5  => [24, 29, 15, 8, 24]
  16 => [24, 37, 26, 20, 16, 29, 24]
  20 => [30, 38, 14, 6, 30]
  12 => [3, -3, -9, 39, 3]
  24 => [33, 6, 14, 9, 18, 2, 33]
  8  => [35, 12, 5, 21, 31, 10, 35]
  17 => [39, -9, -8, 17, 22, 39]
  23 => [11, 18, 9, -2, -5, 11]
  1  => [23, 17, -8, -6, 32, 34, 23]
  22 => [33, 2, 16, 20, 1, 4, 33]
  19 => [29, 16, 2, 18, 11, 15, 29]
  6  => [15, 11, -5, -7, 8, 15]
  11 => [1, 7, 10, 31, 36, 19, 4, 1]
  9  => [38, 32, -6, -2, 9, 14, 38]
  14 => [20, 26, 7, 1, 20]
  3  => [19, 36, 22, 17, 23, 19]
  7  => [13, 25, 12, 35, 13]
  25 => [27, 28, -1, -3, 3, 27]
  4  => [4, 19, 23, 34, 30, 6, 33, 4]
  15 => [8, -7, -4, 25, 13, 37, 24, 8]
  ⋮  => ⋮
```

So, if we wanted the polygon for the fifth generator, we would use:

```julia-repl
julia> get_polygon(vorn, 5)
5-element Vector{Int64}:
 24
 29
 15
  8
 24
```

These polygons are always sorted to be counter-clockwise, with the first element equalling the last also.

The `triangulation` field is simply `tri`. The `boundary_polygons` field is essentially the same as `unbounded_polygons` in most cases, and is only relevant later on when we consider clipped tessellations (for unbounded tessellations, this field is always empty). The `cocircular_circumcenters` field is just for keeping track of cocircular circumcenters to avoid adding duplicate points into the tessellation; it is an implementation detail, you shouldn't need to use it.

The `polygon_points` field stores the coordinates of the points used for the polygons:

```julia-repl
julia> DelaunayTriangulation.get_polygon_points(vorn)
39-element Vector{Tuple{Float64, Float64}}:
 (0.5006477732564246, 0.5478946133538716)
 (0.5008095733914651, 0.3918411792218543)
 (0.10565601252653685, 0.692314380288574)
 (0.4270143529620006, 0.4238420900557963)
 (0.33079352748966667, 0.7624750560470457)
 (0.4753703522300822, 0.23819561748527615)
 (0.5226957825549574, 0.5957485305429776)
 (0.8205472522054272, 0.799163351762504)
 (0.8189082811407147, 0.07294487509120647)
 (0.4026965184932325, 0.6487501595269026)
 ⋮
 (0.330822898041236, 0.16744471621902937)
 (0.3146257041841115, 0.6363882426877013)
 (0.34474111887983666, -0.1904242188538559)
 (0.4811272594567467, 0.391643188260307)
 (0.2982325567607496, 0.16795358951510841)
 (0.5199378327863162, 0.7050503627303082)
 (0.2840418508918076, 0.5804155201166311)
 (0.5862985419779859, 0.657897964228626)
 (0.4650876739497739, 0.08866174194893103)
 (0.10972835538971261, 0.6703707990569139)
```

If we wanted the coordinates for the sixth vertex, we would do:

```julia-repl
julia> get_polygon_point(vorn, 6)
(0.4753703522300822, 0.23819561748527615)
```

The `triangle_to_circumcenter` field is similar to `circumcenter_to_triangle`, except this takes triangles to the index of their associated circumcenter:

```julia-repl
julia> DelaunayTriangulation.get_triangle_to_circumcenter(vorn)
Dict{Tuple{Int64, Int64, Int64}, Int64} with 48 entries:
  (22, 14, 11) => 1
  (8, 13, 7)   => 35
  (22, 19, 16) => 16
  (25, 18, -1) => -1
  (3, 17, 1)   => 17
  (17, 21, 12) => 39
  (14, 13, 11) => 7
  (23, 9, -1)  => -2
  (6, 15, 5)   => 8
  (23, 24, 9)  => 9
  (19, 23, 6)  => 11
  (12, 25, -1) => -3
  (18, 15, -1) => -4
  (16, 19, 5)  => 29
  (13, 15, 7)  => 13
  (4, 3, 1)    => 23
  (19, 6, 5)   => 15
  (21, 11, 8)  => 31
  (6, 23, -1)  => -5
  (18, 25, 2)  => 28
  ⋮            => ⋮
```

Lastly, the `unbounded_polygons` field stores the indices of all generators whose associated polygon is unbounded:

```julia-repl
julia> DelaunayTriangulation.get_unbounded_polygons(vorn)
Set{Int64} with 9 elements:
  6
  15
  25
  18
  9
  12
  17
  23
  1
```

## Relevant Docstrings 

Here are some relevant docstrings for the construction of the tessellation.

```@docs 
initialise_voronoi_tessellation 
prepare_add_voronoi_polygon 
get_next_triangle_for_voronoi_polygon
connect_circumcenters!
add_edge_to_voronoi_polygon!
close_voronoi_polygon!
add_voronoi_polygon!
```