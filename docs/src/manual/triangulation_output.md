```@meta 
CurrentModule = DelaunayTriangulation
```

# Triangulation Output 

In this section, we discuss the output given from [`triangulate`](@ref). Let's take a curve-bounded domain and inspect its output in detail; for information about this domain in particular, see the [curve-bounded domain tutorial](../tutorials/curve_bounded.md). First, here is the triangulation.

```@example curvout
using DelaunayTriangulation
using StableRNGs 
using DelaunayTriangulation: EllipticalArc 
curve = [
    [
        [1, 2, 3], [EllipticalArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0), 2, 1 / 2, 0.0)]
    ],
    [
        [BSpline([(0.0, 0.4), (1.0, 0.2), (0.0, 0.1), (-1.0, 0.2), (0.0, 0.4)])]
    ],
    [
        [4, 5, 6, 7, 4]
    ],
    [
        [BezierCurve([(0.0, -2.0), (0.0, -2.5), (-1.0, -2.5), (-1.0, -3.0)])], [CatmullRomSpline([(-1.0, -3.0), (0.0, -4.0), (1.0, -3.0), (0.0, -2.0)])]
    ],
    [
        [12, 11, 10, 12]
    ],
    [
        [CircularArc((1.1, -3.0), (1.1, -3.0), (0.0, -3.0), positive=false)]
    ]
]
points = [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-2.0, -5.0), (2.0, -5.0), (2.0, -1 / 10), (-2.0, -1 / 10), (-1.0, -3.0), (0.0, -4.0), (0.0, -2.3), (-0.5, -3.5), (0.9, -3.0)]
rng = StableRNG(123)
tri = triangulate(points; boundary_nodes=curve, rng)
refine!(tri; max_area=1e-2, rng)
tri
```

Now let's inspect `tri`. If we look at the fields in `tri`, we see that there is a lot of information stored in `tri`:

```@example curvout
propertynames(tri)
```

Note that each field `X` has an associated accessor `DelaunayTriangulation.get_X`. Let's go through each field. 

## Geometry Fields

First, we list the fields relating to the actual geometry.

### [`get_points(tri)`](@ref get_points)

This field stores the points that were triangulated. This may be a mutated version (not a copy, though, `tri.points === points` above) of the points provided into `tri`, note.

```@example curvout 
get_points(tri)
```

In some cases, the points in this vector will not all appear in `tri`, which is why it is recommend you work with the vertices instead, via [`each_vertex(tri)`](@ref each_vertex) or, for the solid or ghost vertices use [`each_solid_vertex(tri)`](@ref each_solid_vertex) or [`each_ghost_vertex(tri)`](@ref each_ghost_vertex), respectively.

### [`get_triangles(tri)`](@ref get_triangles)

This field stores all the triangles in the triangulation, including both solid and ghost triangles.

```@example curvout 
get_triangles(tri)
```

This is also the same as [`each_triangle(tri)`](@ref each_triangle). If you only wanted the solid triangles, you could use [`each_solid_triangle(tri)`](@ref each_solid_triangle), and similarly for the ghost triangles with [`each_ghost_triangle(tri)`](@ref each_ghost_triangle). All the triangles in these sets are defined to be counter-clockwise:

```@example curvout 
using DelaunayTriangulation: is_positively_oriented, triangle_orientation
all(T -> is_positively_oriented(triangle_orientation(get_point(tri, triangle_vertices(T)...)...)), each_solid_triangle(tri))
```

The above check is not true for all the ghost triangles since the domain is non-convex.

### [`get_boundary_nodes(tri)`](@ref get_boundary_nodes)

This field stores the boundary of the triangulation.

```@example curvout 
get_boundary_nodes(tri)
```

Notice that these boundary nodes are not the same as those provided into `triangulate` above. Instead of having `curve`, `triangulate` constructs a new `boundary_nodes` vector and populates that with the piecewise linear approximation. You could work with these nodes using [`get_boundary_nodes`](@ref) again. For examples, the nodes associated with the first curve are given by 

```@example curvout 
get_boundary_nodes(tri, 1)
```

The second section of this curve is given by 

```@example curvout 
get_boundary_nodes(tri, 1, 2)
```

### [`get_interior_segments(tri)`](@ref get_interior_segments)

This field stores the interior segments of the triangulation. These are the segments that are not part of the boundary.

```@example curvout
get_interior_segments(tri)
```

In our case, there are no such segments, but there would be if we had used `segments` when calling `triangulate`. We could add a segment now, for example

```@example curvout 
r = DelaunayTriangulation.num_points(tri)
add_point!(tri, -3/2, -4.0, concavity_protection = true)
add_point!(tri, -3/2, -1.0, concavity_protection = true)
add_segment!(tri, r + 1, r + 2)
```

Now we see that we have some segments:

```@example curvout 
get_interior_segments(tri)
```

Note that, for any segment `(i, j)`, only one of `(i, j)` or `(j, i)` will appear in the list of interior segments.

### [`get_all_segments(tri)`](@ref get_all_segments)

In contrast to the `interior_segments` field, which only stores the interior segments, this field stores both the interior and boundary segments.

```@example curvout
get_all_segments(tri)
```

Just as for the interior segments, only one of `(i, j)` or `(j, i)` will appear in the list of all segments.

### [`get_weights(tri)`](@ref get_weights)

This field stores the weights associated with each point in the triangulation. Currently, this field is not used properly anywhere and should not be inspected until weighted triangulations are fully implemented. By default, you will see a [`ZeroWeight()`](@ref ZeroWeight):

```@example curvout 
get_weights(tri)
```

## Topology Fields 

Now we list the fields relating to the topology of the triangulation itself.

### [`get_adjacent(tri)`](@ref get_adjacent)

This field stores the adjacent map, mapping each edge `(u, v)` in the triangulation to the vertex `w` such that `(u, v, w)` is a positively oriented triangle in the triangulation. In cases where there is no such triangle, the vertex `0` is returned. For this triangulation, we have:

```@example curvout 
get_adjacent(tri)
```

For example, the mapping `(1648, 2064) => 803` implies that the triangle `(1648, 2064, 803)` is in the triangulation and is positively oriented, which we can verify:

```@example curvout 
DelaunayTriangulation.contains_triangle(tri, 1648, 2064, 803)
```

It is important to note that, for any triangle `(u, v, w)`, the mappings `(u, v) => w`, `(v, w) => u`, and `(w, u) => v` will all appear in the adjacent map. For example:

```@example curvout 
get_adjacent(tri, 1648, 2064), get_adjacent(tri, 2064, 803), get_adjacent(tri, 803, 1648)
```

You can use `get_adjacent(tri, u, v)` to find this vertex `w`, as we demonstrated above. For cases where the edge does not exist, you will get `0`:

```@example curvout 
get_adjacent(tri, 1, 2)
```

One last thing to notice is that some of the vertices `w` will be ghost vertices, and similarly some of the edges `(u, v)` will be ghost edges. For example,

```@example curvout 
get_adjacent(tri, 423, 97)
```

This vertex `-3` means that `(423, 97)` is an edge of the boundary associated with the ghost vertex `-3`.

### [`get_adjacent2vertex(tri)`](@ref get_adjacent2vertex)

The `adjacent2vertex` field is similar to `adjacent`. The difference is that, instead of mapping edges to vertices, `adjacent2vertex` maps vertices `w` to the set of all edges `(u, v)` such that `(u, v, w)` is a positively oriented triangle in the triangulation. For this triangulation, we have:

```@example curvout
get_adjacent2vertex(tri)
```

An example of this mapping is:

```@example curvout 
get_adjacent2vertex(tri, 719)
```

This output means that `(2057, 1625, 719)`, `(1055, 1680, 719)`, `(1625, 1649, 719)`, `(1649, 1162, 719)`, `(1162, 1055, 719)`, and `(1680, 2057, 719)` are all positively oriented triangles in the triangulation, and these are the only triangles that contain the vertex `719`. In contrast to `get_adjacent`, calling `get_adjacent2vertex` on a vertex not in the triangulation will throw a `KeyError`. For ghost vertices, you will get the set of all edges on the boundary associated with that vertex, for example

```@example curvout 
get_adjacent2vertex(tri, -1)
```

gives the set of all edges on the boundary associated with the ghost vertex `-1`. It is important to note that the edges in this set are not returned in any particular order, and so you should not rely on the order of the edges in the set.

### [`get_graph(tri)`](@ref get_graph)

The last field relating to topology is the `graph`, which stores the graph information of the underlying triangulation. In particular, it stores the information about which vertices share an edge. All together, we have three different types of connectivity information being stored in a triangulation:

1. [`Adjacent`](@ref Adjacent): The map from edges to vertices.
2. [`Adjacent2Vertex`](@ref Adjacent2Vertex): The map from vertices to edges.
3. [`Graph`](@ref Graph): The map from vertices to vertices.

For this triangulation, we have:

```@example curvout
get_graph(tri)
```

This output by itself of course isn't too useful. The returned `Graph` does have its own fields,

```@example curvout 
propertynames(get_graph(tri))
```

which you can expect using [`DelaunayTriangulation.get_vertices(tri)`](@ref get_vertices), [`DelaunayTriangulation.get_edges(tri)`](@ref get_edges), and [`get_neighbours(tri)`](@ref get_neighbours). Let us go through each of these fields one at a time.

- For the vertices, you should not be inspecting this field directly. Instead, we provide [`each_vertex(tri)`](@ref each_vertex), [`each_solid_vertex(tri)`](@ref each_solid_vertex), and [`each_ghost_vertex(tri)`](@ref each_ghost_vertex) to iterate over the vertices in the triangulation. These functions return iterators. For example,
    
```@example curvout
each_solid_vertex(tri)
```

```@example curvout 
collect(each_solid_vertex(tri))
```

Since these iterators are derived from `Set`s, you should not rely on the particular ordering of vertices returned.

- For the edges, again we recommend you do not inspect this field directly. Instead, we provide [`each_edge(tri)`](@ref each_edge), [`each_solid_edge(tri)`](@ref each_solid_edge), and [`each_ghost_edge(tri)`](@ref each_ghost_edge) to iterate over the edges in the triangulation. These functions return iterators. For example,

```@example curvout
each_edge(tri)
```

```@example curvout 
collect(each_edge(tri))
```

Since these iterators are derived from `Set`s, you should not rely on the particular ordering of edges returned.

- The neighbours are typically the more useful part of `Graph`. Looking at [`get_neighbours`](@ref) by itself, you obtain a mapping:

```@example curvout 
get_neighbours(tri)
```

This mapping shows the neighbours for each vertex in the triangulation. For example,

```@example curvout 
get_neighbours(tri, 1)
```

shows that the vertices of `1` are `-1`, `18`, `524`, and `-2`. These two ghost vertices, `-1` and `-2`, imply that `1` is on the boundary of the triangulation and is also at the corner of the two sections of the boundary associated with the ghost vertices `-1` and `-2`. Similarly,

```@example curvout
get_neighbours(tri, -2)
```

shows the set of all vertices that are on the boundary associated with the ghost vertex `-2`. Once again, the vertices in this set are not returned in any particular order, and so you should not rely on the order of the vertices in the set. If you did want to traverse the boundary in some order, you would need to use `get_adjacent` to do so, or use the `boundary_nodes`; see, for example, the `get_triangulation_area` function in [this tutorial](../tutorials/constrained_multiply_connected.md). (If you have a node on the boundary you are interested in, [`DelaunayTriangulation.get_left_boundary_node`](@ref) and [`DelaunayTriangulation.get_right_boundary_node`](@ref) are available.)

## Boundary Handling Fields 

The next set of fields relate to how the boundary is handled in the triangulation.

### [`get_boundary_curves(tri)`](@ref get_boundary_curves)

This field is relevant only for curve-bounded domains, like the triangulation in this example. It stores the curves that were used to define the boundary of the triangulation. For our triangulation, we have 

```@example curvout 
get_boundary_curves(tri)
```

This output is not just giving the names of the curves - the actual curves themselves are stored. For example,

```@example curvout 
get_boundary_curves(tri)[3]
```

is the actual `BSpline` that we have provided. The curves are stored in the same order as they were provided, so that the `j`th curve is associated with the `j`th section (i.e., the ghost vertex `-j`) of the boundary. The parts of the boundary that are a sequence of straight lines between vertices are represented by a `PiecewiseLinear` curve. For a triangulation that is not curve-bounded, the output is an empty `Tuple`. For example:

```@example curvout
tri2 = triangulate(rand(2, 50))
get_boundary_curves(tri2)
```

### [`get_boundary_edge_map(tri)`](@ref get_boundary_edge_map)

This field is used to give information about where boundary edges are located in the triangulation. The precise definition is a bit cumbersome to describe: For a given boundary edge `(u, v)`, the `boundary_edge_map` will map `(u, v)` to a tuple of the form `(pos, ℓ)`, so that `pos` is the position of the edge in the boundary, and `ℓ` is the length of the boundary up to the edge. In particular, if `get_boundary_edge_map(tri, u, v) == (pos, ℓ)` and `bn = get_boundary_nodes(tri, pos)`, then `get_boundary_nodes(bn, ℓ) = u` and `get_boundary_nodes(bn, ℓ + 1) = v`. For our triangulation, we have

```@example curvout
get_boundary_edge_map(tri)
```

Since our domain has multiple curves, the `pos` values are `Tuple`s of the form `(m, n)`, where `m` is the curve and `n` is the section of that curve. So, for example, we have

```@example curvout
pos, ℓ = get_boundary_edge_map(tri, 1493, 763)
```

This means that `(1493, 763)` is the `45`th edge of the first section on the third boundary:

```@example curvout
bn = get_boundary_nodes(tri, pos) # notice that we can pass Tuples (m, n) as a single argument
get_boundary_nodes(bn, ℓ), get_boundary_nodes(bn, ℓ + 1)
```

For simpler domains, `pos` may be either a vector of vertices (in the case of a contiguous boundary) or a single vertex (in the case of a single boundary). For the contiguous case, we have:

```@example curvout
tri3 = triangulate_rectangle(0, 1, 0, 1, 10, 10, single_boundary = true) 
get_boundary_edge_map(tri3)
```

This may seem strange, but we note that `get_boundary_nodes(tri, pos)`, when `pos` is a vector of vertices, just returns `pos` once again so that `get_boundary_nodes(get_boundary_nodes(tri, pos), ℓ) == u` is always true, regardless of the boundary type. For example:

```@example curvout 
pos, ℓ = get_boundary_edge_map(tri3, 1, 2)
bn = get_boundary_nodes(tri3, pos)
```

For the sectioned case, we have:

```@example curvout 
tri4 = triangulate_rectangle(0, 1, 0, 1, 10, 10, single_boundary = false) 
get_boundary_edge_map(tri4)
```

When there is no constrained boundary, the boundary edge map will be empty:

```@example curvout
tri5 = triangulate(rand(2, 50))
get_boundary_edge_map(tri5)
```

### [`get_ghost_vertex_map(tri)`](@ref get_ghost_vertex_map)

This field is used to identify to what curve and to what section a ghost vertex is associated. For example, we have

```@example curvout
get_ghost_vertex_map(tri)
```

We see, for example, that the ghost vertex `-7` comes from the first section of the fifth boundary curve. In general, the mappings are of the form `g => pos`, where `g` is the ghost vertex and `pos` is such that `get_boundary_nodes(tri, pos)` is the boundary curve associated with `g`. For example, using [`map_ghost_vertex`](@ref),

```@example curvout
pos = map_ghost_vertex(tri, -7)
bn = get_boundary_nodes(tri, pos)
```

The forms of `pos` for the case of a contiguous boundary, a sectioned boundary, and no boundary are the same as for `get_boundary_edge_map`. For example,

```@example curvout
get_ghost_vertex_map(tri3) # contiguous
```

```@example curvout
get_ghost_vertex_map(tri4) # sectioned
```

```@example curvout
get_ghost_vertex_map(tri5) # no boundary
```

### [`get_ghost_vertex_ranges(tri)`](@ref get_ghost_vertex_ranges)

This field is used to identify all the ghost vertices associated with a curve that has a specific ghost vertex. For example, we have

```@example curvout
get_ghost_vertex_ranges(tri)
```

This output means, for example, that the ghost vertex `-5` is associated with a curve that has both ghost vertices `-6` and `-5` associated with it. You would use [`get_ghost_vertex_range`](@ref) to get the range of ghost vertices associated with a specific ghost vertex. For example,

```@example curvout
get_ghost_vertex_range(tri, -5)
```

Similarly, the output for the contiguous, sectioned, and no boundary cases are

```@example curvout
get_ghost_vertex_ranges(tri3) # contiguous
```

```@example curvout
get_ghost_vertex_ranges(tri4) # sectioned
```

```@example curvout
get_ghost_vertex_ranges(tri5) # no boundary
```

## Other Fields 

There are some other more general fields to inspect.

### [`get_convex_hull(tri)`](@ref get_convex_hull)

For all triangulations, the `convex_hull` field stores the [`ConvexHull`](@ref) of the triangulation. For our triangulation, we have

```@example curvout
get_convex_hull(tri)
```

To inspect the vertices of the convex hull, you use [`get_convex_hull_vertices`](@ref):

```@example curvout 
get_convex_hull_vertices(tri)
```

If you ever need to reconstruct the convex hull, say after some dynamic updates, you would use [`convex_hull!`](@ref).

### [`get_representative_point_list(tri)`](@ref get_representative_point_list)

This field is related to the need for a representative point for each curve, as described in the [ghost vertices section](../manual/ghost_triangles.md). For our triangulation, we have

```@example curvout
DelaunayTriangulation.get_representative_point_list(tri)
```

This shows, for example, that the sixth boundary curve's representative point, i.e. its pole of inaccessibility, is the point `(0, -3)`. You could also review this using [`DelaunayTriangulation.get_representative_point_coordinates`](@ref):

```@example curvout
DelaunayTriangulation.get_representative_point_coordinates(tri, 6)
```

In the case of a triangulation with no constrained boundary, the representative point list is simply the centroid of the domain. If you ever need to recompute the representative points, you need [`DelaunayTriangulation.compute_representative_points!`](@ref).

### [`get_polygon_hierarchy(tri)`](@ref get_polygon_hierarchy)

This field is used to store information about which boundary curves are contained in other boundary curves. For our triangulation, we have

```@example curvout
DelaunayTriangulation.get_polygon_hierarchy(tri)
```

This field is not intended for public use. One useful method that makes direct use of the polygon hierarchy is [`find_polygon`](@ref), but the hierarchy's main use is in boundary enrichment for curve-bounded domains, and for checking arguments via [`check_args`](@ref).

### [`boundary_enricher(tri)`](@ref get_boundary_enricher)

This field, just as for the `polygon_hierarchy`, is not intended for public use. It is the field used for enriching the boundary for initialising the triangulation of a curve-bounded domain. 

```@example curvout 
DelaunayTriangulation.get_boundary_enricher(tri)
```

### [`get_cache(tri)`](@ref get_cache)

This field is not intended for public use. It is used to provide several caches for use during triangulation, such as reusing arrays for constrained triangulations.

```@example curvout 
DelaunayTriangulation.get_cache(tri)
```