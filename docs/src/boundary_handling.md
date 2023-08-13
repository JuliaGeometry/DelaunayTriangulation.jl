```@meta
CurrentModule = DelaunayTriangulation
```

# Ghost Triangles, Boundary Handling, and Boundary Specifcations

Here we will give a description about how we represent boundaries, and also how we use ghost triangles.

## Boundary Indices 

We use negative indices to denote vertices belonging to a boundary. For example, if `tri` is a triangulation and `get_adjacent(adj, u, v) == -1`, then this means that `(u, v)` is an edge on the boundary. We call this `-1` a ghost vertex or a boundary index (positive vertices could also be called solid vertices), with `-1` defined from `DelaunayTriangulation.BoundaryIndex`, which we discuss more in the next section.

In the case of a single contiguous outer boundary, having `-1` as the only boundary index is simple and works fine. If we have multiple segments or multiple boundaries, then we need to somehow have multiple boundary indices so that we can refer to each segment separately. We accomplish this by simply subtracting 1 from the current boundary index for each new segment. This is handled by `add_boundary_information!`. For example, if we had

```julia
bn = [
    [segment_1, segment_2, segment_3],
    [segment_4, segment_5],
    [segment_6],
    [segment_7, segment_8, segment_9, segment_10]
]
```

then the `i`th segment will map to `-i`. Note that this always means that `-1` belongs to the outer-most boundary. There is a possible issue we can have with this when, for example, stepping around a boundary, since nodes will occur in two segments and hence nodes may not necessarily have a unique boundary index assigned to them. To handle this case, allowing us to check for all possible boundary indices when stepping around a boundary (or however else we might want to use `Adjacent` or similar), `get_adjacent` can call a safer version `_safe_get_adjacent`, which checks the boundary indices using `boundary_index_ranges`.

These issues are why the `Triangulation` data structure has the `boundary_map` and `boundary_index_ranges` fields. `boundary_map`, constructed with `construct_boundary_map`, is used to map a boundary index to its segment in the set of boundary nodes, so that `get_boundary_nodes(boundary_nodes, map_boundary_index(boundary_map, g))` gives the nodes corresponding to the segment which has boundary index `g`. To handle the issue with a curve having multiple boundary indices, we use `boundary_index_ranges`, constructed with `construct_boundary_index_ranges`, to map a boundary index `g` to all other indices that could be found on the associated curve. For example, in the `bn` example above, `-2` would map to `-3:-1`.

## Ghost Triangles 

Ghost triangles are a special triangle that has a solid edge `(u, v)` and a vertex `g` associated with some boundary index. These ghost triangles are needed to make point location actually work properly when points are outside of the triangle, provided we associate the ghost vertex `g` with a physical point. For the outer-most boundary, this physical point just has to be somewhat in the center of the domain, which we define using a centroid when building the triangulation and the pole of inaccessibility once we have built the entire triangulation. With this physical point, ghost edges `(u, g)` are then interpreted to be of infinite extent, pointing from `u` out to infinity, but collinear with this central point. 

In the case of an inner boundary, these ghost edges are of finite extent, simply connecting with the point `g` which is define via the pole of inaccessibility.

For the Bowyer-Watson algorithm, we need a definition for the circumcircle of a ghost triangle. For ghost triangles that belong to inner boundaries, we simply use the triangle that connects the points, since there is no issue with infinity here. For the outer ghost triangles, we need to be careful. Imagine taking a triangle and pulling away one of its vertices to infinity. The circumcircle would keep growing until it eventually covers the entire space on the side of the fixed edge that the point was on. In particular, the circle becomes the line through the fixed edge, dividing the plane into two half-planes. We then say that a point is in the circumcircle of an outer ghost triangle if it is in the open half-plane on the other side of the edge from the triangulation, or if it is on the edge itself (but is not one of the vertices). The union of the open half-plane and this open edge is called the outer half-plane of the edge.

## Boundary Specification 

When considering constrained triangulations, we allow for a set of boundary nodes to be provided. These nodes can be defined to represent three possible domains:

1. A contiguous boundary with no holes. For example, a circle.
2. A single boundary with no holes, but the boundary is split into multiple segments that can be identified separately via boundary indices. For example, a square with a boundary index for each side.
3. A boundary comprising multiple disjoint curves, i.e. a multiply-connected domain. For example, a square with a circular hole inside. The most support is provided for the case where no interior hole contains another interior, but it is allowed (see the discussion in the constrained triangulation section).

We do also provide support for much more complex geometries, such as disjoint regions, but the support is only for constructing the mesh. Limited support exists, for example, for point location in such regions. See the examples in the constrained triangulation section. There might also be issues with domains that are non-convex.

The way to represent boundaries can be customised as needed (see the Interfaces section), but by default we provide the following specifications; this specification is what we use in the `convert_boundary_points_to_indices` function. Let `bn` be the set of boundary nodes and `I` the integer type.

### Contiguous boundary 

For a contiguous boundary, `bn` should be a single `Vector{I}` such that `bn[begin] == bn[end]`, defining a boundary in counter-clockwise order. 

### Single boundary split into segments 

Now consider a single boundary with no holes, but given by multiple segments, say into `ns` segments. Then `bn` should be a `Vector{Vector{I}}` with `length(bn) = ns`. The vectors should be defined such that the endpoints connect, meaning `bn[i][end] == bn[i+1][begin]` for `i < ns` or `bn[end][end] == bn[begin][begin]`. Additionally, the nodes should be defined in counter-clockwise order. 

### Multiply-connected domain 

Now suppose we have a boundary that is split into multiple curves, say into `nc` curves. In this case, `bn` is a `Vector{Vector{Vector{I}}}` where each individual curve `bn[i]` is assumed to represent a single boundary split into segments, with each `bn[i]` matching the previous specification. In this case, we provide the most support for the case where `bn[1]` is the outermost boundary curve, and all curves `bn[i]` for `i > 1` are inside `bn[1]`. With this, `bn[1]` should be counter-clockwise, while `bn[i]` for `i > 1` should be clockwise. While it is possible to place interior curves inside other interior curves, e.g. `bn[j]` could be inside `bn[i]` as long as `bn[j]` is counter-clockwise, less support for point location and other features is provided for this case. See the constrained triangulation section.