```@meta 
CurrentModule = DelaunayTriangulation
```

# Voronoi Tessellation output

In this section, we discuss the output given from [`voronoi`](@ref). We consider a simple clipped example for examining this output.

```@example vorout 
using DelaunayTriangulation
using StableRNGs 
rng = StableRNG(123)
p1 = (0.0, 1.0)
p2 = (3.0, -1.0)
p3 = (2.0, 0.0)
p4 = (-1.0, 2.0)
p5 = (4.0, 2.0)
p6 = (-2.0, -1.0)
p7 = (2.0, 1.0)
p8 = (5.0, 1.0)
points = [p1, p2, p3, p4, p5, p6, p7, p8]
tri = triangulate(points; rng)
vorn = voronoi(tri, clip = true, rng = rng)
vorn
```

Now let's inspect `vorn`. The fields in `vorn` are:

```@example vorout 
propertynames(vorn)
```

Let's examine this fields one at a time.

## [`get_triangulation(vorn)`](@ref get_triangulation)

This field stores the triangulation object used to generate the Voronoi tessellation.

```@example vorout
DelaunayTriangulation.get_triangulation(vorn)
```

## [`get_generators(vorn)`](@ref get_generators)

This field stores the generators of the Voronoi tessellation, i.e. the points associated with the polygons in the tessellation. These are simply the points in the triangulation, but are stored differently in case some points are not present in the triangulation.

```@example vorout 
DelaunayTriangulation.get_generators(vorn)
```

See that the generators are stored as a `Dict`, with the vertices mapping to their associated coordinates. The preferred way to access the generators is through [`get_generator`](@ref), which can be used to obtain the coordinates for a given vertex. For example:

```@example vorout 
get_generator(vorn, 1)
```

```@example vorout 
get_generator(vorn, 2, 6)
```

If you want to iterate over each generator, you should use [`each_generator`](@ref), which simply returns an unordered iterator over the generators.

```@example vorout
each_generator(vorn)
```

You can get the number of generators using [`num_generators`](@ref):

```@example vorout
DelaunayTriangulation.num_generators(vorn)
```

## [`get_polygon_points(vorn)`](@ref get_polygon_points)

This field is used for storing the coordinates of the Voronoi polygons:

```@example vorout
DelaunayTriangulation.get_polygon_points(vorn)
```

These coordinates are used to define the vertices of these polygons. It is important to note that these points are not guaranteed to be unique if a circumcenter appears on the boundary, or if the tessellation is clipped. To access the coordinate of a polygon associated with a given polygon vertex, use [`get_polygon_point`](@ref):

```@example vorout
get_polygon_point(vorn, 1)
```

If you wanted the total number of polygon vertices, possibly counting duplicate vertices, you can use [`num_polygon_vertices`](@ref):

```@example vorout
num_polygon_vertices(vorn)
```

## [`get_polygons(vorn)`](@ref get_polygons)

This field is used for mapping polygon indices to their vertices, with these vertices referring to coordinates from `get_polygon_points(vorn)`:

```@example vorout
DelaunayTriangulation.get_polygons(vorn)
```

For example, using [`get_polygon`](@ref) with the first polygon:

```@example vorout
get_polygon(vorn, 1)
```

means that the first polygon has vertices `[19, 21, 4, 18, 16, 3, 19]`. All polygons are stored as a circular vector of counter-clockwise oriented vertices. Note also that polygon indices are the same as the generator vertex, so that e.g. this polygon `1` is that associated with generator `1`.  You can get the number of polygons using [`num_polygons`](@ref):

```@example vorout
num_polygons(vorn)
```

Of course, this is the same as the number of generators. If you wanted to get the coordinates instead of the vertices associated with a polygon, you use [`get_polygon_coordinates`](@ref):

```@example vorout
get_polygon_coordinates(vorn, 1)
```

This is almost the same result as doing

```@example vorout
get_polygon_point.(Ref(vorn), get_polygon(vorn, 1))
```

but is more efficient.

When the tessellation is unbounded, this field is slightly different. Let's consider an example.

```@example vorout 
rng2 = StableRNG(321)
points2 = rand(rng2, 2, 10)
tri2 = triangulate(points2; rng = rng2)
vorn2 = voronoi(tri2, rng = rng2)
DelaunayTriangulation.get_polygons(vorn2)
```

Notice that some of these polygons have negative vertices. These vertices are not the same as ghost vertices from triangulations. They do still represent points out at infinity, but their numbering is of no importance - only the fact that they are negative is. All polygons that have a negative vertex are unbounded, and the negative vertex is the vertex that is out at infinity. When we try to get the coordinates of a polygon with a negative vertex, we get an error:

```@repl vorout 
get_polygon_coordinates(vorn2, 5)
```

We need to provide a bounding box to obtain these coordinates. As mentioned, we can obtain a reasonable default for a bounding box using [`polygon_bounds`](@ref), which returns the bounding box in the form `(xmin, xmax, ymin, ymax)`:

```@example vorout
bbox = polygon_bounds(vorn2)
```

If we now use this bounding box to get the coordinates of the polygon, we get a polygon that will have been clipped with this box:

```@example vorout
get_polygon_coordinates(vorn2, 4, bbox)
```

For iterating over polygons, there are two main methods. If you only care about the vertices and not the index associated with a polygon, you can use [`each_polygon`](@ref):

```@example vorout
each_polygon(vorn)
```

Alternatively, [`each_polygon_index`](@ref) can be used 
to iterate over the index of each polygon (equivalently, over each generator vertex), followed by `get_polygon` to access the associated vertices of that polygon:

```@example vorout
for v in each_polygon_index(vorn)
    println(v, ": ", get_polygon(vorn, v))
end
```

## [`get_circumcenter_to_triangle(vorn)`](@ref get_circumcenter_to_triangle)

When a Voronoi tessellation is constructed, the circumcenters of the triangles from the underlying Delaunay triangulation are used to define the polygon vertices. This field is used to associate a given polygon vertex with the triangle that it came from:

```@example vorout
DelaunayTriangulation.get_circumcenter_to_triangle(vorn)
```

We see, for example, that the fifth polygon vertex is derived from the circumcenter of the triangle `(3, 6, 2)`. We can verify this:

```@example vorout 
p = DelaunayTriangulation.triangle_circumcenter(tri, (3, 6, 2)) 
q = get_polygon_point(vorn, 5)
p == q 
```

The negative vertices show what ghost triangle is being associated with an unbounded edge of the tessellation. Since we clipped our tessellation, these no longer appear in the previously mentioned fields, but they are still stored in this field; they are used in clipping for determining the direction of the unbounded edges. If you wanted the triangle associated with a specific vertex, you could use 
 
```@example vorout 
DelaunayTriangulation.get_circumcenter_to_triangle(vorn, 5)
```

## [`get_triangle_to_circumcenter(vorn)`](@ref get_triangle_to_circumcenter)

This map is the inverse of `get_circumcenter_to_triangle(vorn)`, and is used to map a triangle to its circumcenter:

```@example vorout
DelaunayTriangulation.get_triangle_to_circumcenter(vorn)
```

The triangles are stored so that the minimum vertex is always last, maintaining the counter-clockwise orientation of each. You can get the circumcenter index associated with a given triangle using:

```@example vorout 
DelaunayTriangulation.get_triangle_to_circumcenter(vorn, (3, 7, 1))
```

## [`get_unbounded_polygons(vorn)`](@ref get_unbounded_polygons)

This field is used to store the polygons that are unbounded. For our clipped tessellation, this field is empty:

```@example vorout
DelaunayTriangulation.get_unbounded_polygons(vorn)
```

 For our unclipped tessellation, we have

```@example vorout
DelaunayTriangulation.get_unbounded_polygons(vorn2)
```

and we can verify that this is consistent with our observation that each unbounded polygon has a negative vertex:

```@example vorout
S1 = DelaunayTriangulation.get_unbounded_polygons(vorn2)
S2 = Set{Int}()
for v in each_polygon_index(vorn2)
    C = get_polygon(vorn2, v)
    if any(<(0), C)
        push!(S2, v)
    end
end
S1 == S2
```

## [`get_cocircular_circumcenters(vorn)`](@ref get_cocircular_circumcenters)

In cases where pairs of triangles are cocircular, meaning lie on the same circle, their circumcenter will be equal and thus the polygon vertex associated with these two triangles will be the same. This field will store the vertices of these circumcenters that came from cocircular triangles. For our case, there is a single entry:

```@example vorout
DelaunayTriangulation.get_cocircular_circumcenters(vorn)
```

This means that the first polygon vertex is associated with at least two cocircular triangles, and is why, for example, the values from `get_triangle_to_circumcenter(vorn)` are not unique:

```@example vorout
t2c = DelaunayTriangulation.get_triangle_to_circumcenter(vorn)
allunique(values(t2c))
```

and also explains why `get_circumcenter_to_triangle(vorn)` and `get_triangle_to_circumcenter(vorn)` are not perfect inverses of each other. The triangle mapped to from `5` inside the `circumcenter_to_triangle` map is simply the first one that the algorithm happened to encounter. To find the two cocircular triangles, we could use the `triangle_to_circumcenter` map:

```@example vorout
triangles = NTuple{3,Int}[]
for (t, c) in t2c
    c == 1 && push!(triangles, t)
end
```

In fact, there are three triangles that all have the same circumcenter!

```@example vorout 
DelaunayTriangulation.triangle_circumcenter.(Ref(tri), triangles)
```

## [`get_adjacent(vorn)`](@ref get_adjacent)

This map is similar to the [`Adjacent`](@ref) map for triangulations:

```@example vorout
get_adjacent(vorn)
```

In this case, oriented edges are being mapped to the polygon it belongs to. For example,

```@example vorout 
get_adjacent(vorn, 3, 19)
```

means that the edge `(3, 19)` belongs to the first polygon; these vertices `3` and `19` refer to the vertices from the polygon, while the `1` from the output refers to a generator.

## [`get_boundary_polygons(vorn)`](@ref get_boundary_polygons)

This field is used to store the polygons that are on the boundary of the tessellation, and is only relevant for clipped tessellations. For our clipped tessellation, we have:

```@example vorout
DelaunayTriangulation.get_boundary_polygons(vorn)
```

For our unclipped tessellation, the field is empty:

```@example vorout
DelaunayTriangulation.get_boundary_polygons(vorn2)
```
