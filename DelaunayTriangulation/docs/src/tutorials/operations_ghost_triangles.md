```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/operations_ghost_triangles.jl"
```

# Triangulation Operations
## Adding or Clearing Ghost Triangles

An important component to be aware of when considering
dynamic updates to triangulations are the ghost triangles.
As we discussed in the [vertex insertion/deletion example](operations_vertex_insertion_deletion.md),
ghost triangles are needed when we are making updates outside of the boundary
of the current triangulation.

````@example operations_ghost_triangles
using DelaunayTriangulation
using CairoMakie
````

Let us take an example triangulation.

````@example operations_ghost_triangles
points = [(-1.0, -1.0), (1.0, -1.0), (0.0, 1.0)]
tri = triangulate(points)
````

````@example operations_ghost_triangles
fig, ax, sc = triplot(tri, show_ghost_edges = true)
fig
````

The ghost triangles are represented by the convex regions
bounded by the blue lines. By default, `triangulate` will keep these
ghost triangles. If you want to remove them, you'd have to use
[`delete_ghost_triangles!`](@ref) or use `delete_ghosts = true` inside `triangulate`.

If you do need to query whether your triangulation already
has ghost triangles, use

````@example operations_ghost_triangles
DelaunayTriangulation.has_ghost_triangles(tri)
````

To clear the ghost triangles, use

````@example operations_ghost_triangles
delete_ghost_triangles!(tri)
DelaunayTriangulation.has_ghost_triangles(tri)
````

An important note for us to make is that ghost triangles
are not just there as a concept, but they are actually
physically stored. Adding them back with [`add_ghost_triangles!`](@ref), we have:

````@example operations_ghost_triangles
add_ghost_triangles!(tri)
get_triangles(tri)
````

See that there is not just the triangle `(1, 2, 3)`, but
also `(3, 2, -1)`, `(1, 3, -1)`, and `(2, 1, -1)` (where
the ghost triangles are also oriented counter-clockwise). For example,

````@example operations_ghost_triangles
get_adjacent(tri, 3, 2)
````

````@example operations_ghost_triangles
get_adjacent(tri, 3, -1)
````

````@example operations_ghost_triangles
get_adjacent(tri, -1, 2)
````

````@example operations_ghost_triangles
get_neighbours(tri, -1)
````

````@example operations_ghost_triangles
get_adjacent2vertex(tri, -1)
````

If we delete them, they are no longer there.

````@example operations_ghost_triangles
delete_ghost_triangles!(tri)
get_triangles(tri)
````

As a last note, we remark that the ghost vertices that define the
vertex of these ghost triangles is still there regardless of whether
the triangulation has ghost triangles. Thus, for example, the
following still work

````@example operations_ghost_triangles
get_neighbours(tri, -1)
````

````@example operations_ghost_triangles
get_adjacent2vertex(tri, -1)
````

````@example operations_ghost_triangles
get_adjacent(tri, 3, 2)
````

You can remove them from the graph, using

````@example operations_ghost_triangles
DelaunayTriangulation.delete_ghost_vertices_from_graph!(tri)
````

so that e.g. `get_neighbours(tri, -1)` is then an error. This
will still not remove them from the `adjacent` and `adjacent2vertex`
maps, but it does mean for example that

````@example operations_ghost_triangles
collect(each_solid_vertex(tri)) == collect(each_vertex(tri))
````

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/operations_ghost_triangles.jl).

```julia
using DelaunayTriangulation
using CairoMakie

points = [(-1.0, -1.0), (1.0, -1.0), (0.0, 1.0)]
tri = triangulate(points)

fig, ax, sc = triplot(tri, show_ghost_edges = true)
fig

DelaunayTriangulation.has_ghost_triangles(tri)

delete_ghost_triangles!(tri)
DelaunayTriangulation.has_ghost_triangles(tri)

add_ghost_triangles!(tri)
get_triangles(tri)

get_adjacent(tri, 3, 2)

get_adjacent(tri, 3, -1)

get_adjacent(tri, -1, 2)

get_neighbours(tri, -1)

get_adjacent2vertex(tri, -1)

delete_ghost_triangles!(tri)
get_triangles(tri)

get_neighbours(tri, -1)

get_adjacent2vertex(tri, -1)

get_adjacent(tri, 3, 2)

DelaunayTriangulation.delete_ghost_vertices_from_graph!(tri)

collect(each_solid_vertex(tri)) == collect(each_vertex(tri))
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

