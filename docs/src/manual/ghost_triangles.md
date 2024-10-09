```@meta 
CurrentModule = DelaunayTriangulation
```

# Ghost (Negative) Vertices

In this section, we will describe how ghost triangles and ghost vertices are used in this package to represent boundary information. This is a crucial part of the package, and especially for understanding output from triangulations, as it is how we represent regions exterior to a boundary and traverse this exterior.

Ghost vertices are negative vertices that are associated with a part of a boundary. Mathematically speaking, they are typically treated as points out at infinity, and each edge on the boundary adjoins a ghost vertex. This way, all edges have two adjoining vertices and we can associate any point in space, including points outside of the triangulation, with a triangle. As an example, if `tri` is a triangulation and `get_adjacent(tri, u, v) == -1`, then this means that `(u, v)` is an edge on the boundary of the triangulation. This number `-1` is defined as `DelaunayTriangulation.𝒢` internally.

In the case of a single contiguous boundary, the only possible ghost vertex is `-1`. When it comes to considering a boundary with mutliple sections or multiple boundaries, then we need to have multiple ghost vertices to refer to eachs ection separately. We accomplish this by simply subtracting 1 from the current ghost vertex for each new section. For example, if the boundary node vector is 

```julia 
boundary_nodes = [
    [section_1, section_2, section_3],
    [section_4, section_5],
    [section_6],
    [section_7, section_8, section_9, section_10]
]
```

then the `i`th section will be associated with the ghost vertex `-i`. We note that this setup does not allow for all nodes to be uniquely associated with a ghost vertex, since nodes connecting between two sections will necessarily have two ghost vertices associated with them. Instead, each edge has only one ghost vertex associated with it. This issue causes problems when using `get_adjacent` to step along a boundary and is why you will notice, if you look at the internals for `get_adjacent(::Triangulation, ::Any)`, the implementation is somewhat verbose. These issues are also why the `Triangulation` data structure has the `ghost_vertex_map` and `ghost_vertex_ranges` fields, with the former allowing `get_boundary_nodes(tri, DelaunayTriangulation.map_ghost_vertex(tri, g))` to be used to get the boundary nodes associated with a ghost vertex `g`. The issue with a curve having multiple ghost vertices is handled with `ghost_vertex_ranges`, allowing us to identify all possible ghost vertices associated with a curve that has a specific ghost vertex via `DelaunayTriangulation.get_ghost_vertex_range(tri, g)`. For example, in the `boundary_nodes` example above, `-2` would map to `-3:-1`.

# Ghost Triangles 

Ghost triangles are special triangles that have a solid edge `(u, v)` and a vertex `g` associated with some ghost vertex. These ghost triangles are needed to make point location work when points are outside of the triangulation, provided that we associate the ghost vertex `g` with a physical point. For the first boundary, the physical point just has to be somewhat in the centre of the domain, which we define using a centroid when building the triangulation and the pole of inaccessibility once the triangulation is built. With this physical point, ghost edges `(u, g)` are then interpreted to be of infinite extent, pointing from `u` out to infinity, but collinear with this central point.

The physical point associated with a ghost vertex `g`, and the interpretation of a ghost edge `(u, g)`, also depends on whether the curve's orientation. For all curves, the pole of inaccessibility defines the physical point for the single curve, except for the first ghost vertex in which case the pole of inaccessibility is defined relative to the entire domain. The ghost edge is considered to be of infinite extent, as mentioned above, for curves that are counter-clockwise. For curves that are oriented clockwise, the ghost edge has a finite interpretation, simply connecting the point `u` with the associated curve's pole of inaccessibility. We give an example of this below.

```@example ghosttrifigex
using DelaunayTriangulation
using CairoMakie
points = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
boundary_nodes = [[[1, 2, 3, 4, 1]], [[CircularArc((0.0, 0.5), (0.0, 0.5), (0.0, 0.0), positive=false)]]]
tri = triangulate(points; boundary_nodes, coarse_n=32)
fig, ax, sc = triplot(tri, show_ghost_edges=true)
c1 = DelaunayTriangulation.get_representative_point_coordinates(tri, 1)
c2 = DelaunayTriangulation.get_representative_point_coordinates(tri, 2)
scatter!(ax, [c1, c2], color = [:red, :magenta])
fig
```

As you can see, the outer boundary has ghost edges (shown in blue) going out to infinity, oriented with the pole of inaccessibility of the entire domain (shown in red). The ghost edges along the circular boundary are finite and simply connect with the pole of inaccessibility of the circle (shown in magenta). 

For more complex domains, in particular non-convex domains, the ghost edges start to overlap and they become less useful, which unfortunately slows down point location (see [`find_triangle`](@ref)'s docstring).
