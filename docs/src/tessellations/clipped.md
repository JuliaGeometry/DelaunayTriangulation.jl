```@meta
CurrentModule = DelaunayTriangulation
```

# Clipped Voronoi Tessellations

Often, it is useful to chop off the unbounded polygons in a tessellation, truncating them to some boundary. Usually this is just a box. We currently provide support for chopping to a convex hull, but the aim is to eventually support chopping to any boundary (the same algorithm we use should apply here, but there are some special cases that are really quite difficult). 

The same function `voronoi` is used for this, just setting the second argument to be `true`. The algorithm is described at the end of this section.

## Example

Let us give an example.

```julia
using DelaunayTriangulation, CairoMakie
pts = randn(2, 250)
tri = triangulate(pts)
vorn = voronoi(tri)
vorn_clip = voronoi(tri, true)
cmap = Makie.cgrad(:jet)
colors = get_polygon_colors(vorn, cmap)
fig = Figure()
ax = Axis(fig[1, 1], aspect=1)
voronoiplot!(ax, vorn, strokecolor=:red, strokewidth=0.2, polygon_color=colors, show_generators=false)
triplot!(ax, tri, strokewidth=0.0, strokecolor=(:black, 0.4), show_all_points=false)
xlims!(ax, -5, 5)
ylims!(ax, -5, 5)
ax = Axis(fig[1, 2], aspect=1)
voronoiplot!(ax, vorn_clip, strokecolor=:red, strokewidth=0.2, polygon_color=colors, show_generators=false)
triplot!(ax, tri, strokewidth=0.0, strokecolor=(:black, 0.4), show_all_points=false)
xlims!(ax, -5, 5)
ylims!(ax, -5, 5)
```

```@raw html
<figure>
    <img src='../figs/bounded.png', alt='Clipped Voronoi Tessellation'><br>
</figure>
```

As we can see, all the polygons have now been chopped so that the entire tessellation fits into the original boundary of the dual triangulation. Also, the `unbounded_polygons` and `boundary_polygons` fields have been updated:

```julia-repl
julia> DelaunayTriangulation.get_unbounded_polygons(vorn_clip)
Set{Int64}()

julia> DelaunayTriangulation.get_boundary_polygons(vorn_clip)
Set{Int64} with 31 elements:
  169
  56
  200
  195
  72
  180
  221
  8
  37
  187
  32
  6
  171
  190
  69
  219
  73
  ⋮
```

## Algorithm

Since the code for clipping is quite involved, most likely more involved than it should be, it is probably worthwhile describing the actual implementation here. I want to one day refine this approach, and allow it to work on any domain. The approach is mostly based on the ideas in the paper "Efficient Computation of Clipped Voronoi Diagram for Mesh Generation"
by Yan, Wang, Lévy, and Liu. I imagine if I had access to their code, everything would be a lot nicer (that seems to be true in many computational geometry papers).

The function that performs the clipping is `clip_voronoi_tessellation!`:

```@docs 
clip_voronoi_tessellation!
```

This function starts off by finding all the intersections and then adding them into the point set, then using these new points to clip the polygons. We describe these steps below.

### Finding all intersections 

The function `find_all_intersections` finds the intersections of the polygons with the boundary.

```@docs 
find_all_intersections
initialise_clipping_arrays
```

The idea is to use a first-in first-out (FIFO) queue to prioritise the edges and polygons to consider. In particular, we pick some random edge and then use its midpoint to find what polygon that midpoint is in (via `jump_and_march`). This edge and the polygon are then added into the queue to be processed. Starting with this edge-polygon pair, we keep looking for intersections until we have processed all boundary edges and there are no more pairs in the queue. The function that does this processing of queued pairs is `dequeue_and_process!`:

```@docs 
dequeue_and_process!
enqueue_new_edge!
```

Essentially, this function first checks if the pair has already been processed and, if so, returns. Otherwise, the first step is to go into `process_polygon!` to find the intersections of the polygons with the nearby boundary:

```@docs 
process_polygon!
```

The goal of `process_polygon!` is to iterate over the edges of the current polygon to look for intersections, searching (1) the current edge, (2) the edge to the left of the current edge, and (3) the edge to the right of the current edge. Suppose we have an edge $e_{uv}$ of the polygon being considered. There are four possibilities:

1. First, both $u$ and $v$ could correspond to ghost vertices, meaning no intersection with the boundary is possible.
2. Alternatively, $e_{uv}$ could be a ray going inwards, meaning $u$ is a ghost vertex while $v$ corresponds to an actual circumcenter. In this case, `process_ray_intersection!` is used to look for an intersection of $e_{uv}$ with the edge of the ghost triangle corresponding to $u$. If an intersection is found, we do need to be careful of rays that intersect multiple boundary edges (as is common when we have very few triangles in the underlying triangulation). This check is done with `process_ray_intersection_with_other_edges!`, which simply looks at the left and right edges along with the current edge. 
3. Same as the above, except $e_{uv}$ could be a ray going inwards, meaning $u$ now corresponds to an actual cicrcumcenter while $v$ is a ghost vertex.
4. Lastly, $e_{uv}$ could be a fniite segment, in which case we can use simple intersection formula to test for intersections with the left and right edges along with the current edge. This check is done via `process_segment_intersection!`.

In these steps, when there are no intersections we also check if it was because a segment lies entirely outside of the domain. We keep track of these points as these will need to be deleted later.

Some of the relevant docstrings for working with `process_polygon!` are below.

```@docs 
is_segment_between_two_ghosts
is_ray_going_in
is_ray_going_out
is_finite_segment
process_segment_intersection!
process_ray_intersection!
add_to_intersected_edge_cache!
add_segment_intersection!
segment_intersection_coordinates
intersection_of_edge_and_bisector_ray
classify_and_compute_segment_intersection
process_ray_intersection_with_other_edges!
```

Once we have gone through `process_polygon!`, we need to assign the identified intersections to the current edge, the left edge, or to the right edge. This is done with `classify_intersections!`.

```@docs 
classify_intersections!
```

Once this is done, we need to actually consider the intersection points. The two goals here are: (1) Check if a corner point inside the polygon has to be added, and (2) look for other incident polygons to process. This processing is done via `process_intersection_points!`:

```@docs 
process_intersection_points!
```

Basically, this function first looks at the left and right intersections. If we have intersections on the left and on the current edge, then we will also have the vertex shared by the two edges included as a point we need to chop to, unless the current polygon being considered corresponds to a generator not on a boundary. We check this for each intersection and update the intersections with `add_segment_intersection!` accordingly. In such a case, we add the current edge, together with the indices of the left edge (corresponding to polygons) to the queue, ready for the next iteration. The same is done for the right edge. After handling this case, we check all the intersections without worry for corner points. For each intersection, we take the polygon on the other side of the intersecting segment and add it to the queue together with the boundary ede that was intersected. 

We do all the above steps until we run out of edges to process and the queue is empty, adding all the intersection points into the point list using `add_intersection_points!`:

```@docs 
add_intersection_points!
```

### Clipping the polygons

The next step is to clip the polygons to the boundary using the computed intersections. This is handled via `clip_all_polygons!`:

```@docs 
clip_all_polygons!
```

This function iterates over all the identified boundary polygons and calls `clip_polygon!` on each:

```@docs 
clip_polygon!
```

Consider a single polygon. The `clip_polygon!` function starts by deleting the adjacencies for this polygon, ready for updating later. We then remove all ghost vertices from the vertex list, and all circumcenters identified from the previous step that lie outside of the domain. We then add all the vertices. With this processing, the vertex list is no longer sorted, so `sort_convex_polygon!` is called to get a counter-clockwise representation. Adding in the adjacencies via `add_polygon_adjacent!`, this completes the clipping of this polygon. 

With all polygons processed, we add in all that we processed into the `boundary_polygons` field via `add_all_boundary_polygons!`:

```@docs 
add_all_boundary_polygons!
```

With this last step, the clipping is complete.

