```@meta
CurrentModule = DelaunayTriangulation
```

# Constrained Triangulations

Here we describe the algorithm used for computing a constrained triangulation. We assume that we have some point set $\mathcal P$, a set of edges $\mathcal E$ to be inserted, and some boundary edges $\mathcal B$. The algorithm we implement is given [here](https://doi.org/10.1016/j.comgeo.2015.04.006) and is built upon the basic idea of Chew's algorithm for triangulating convex polygons, namely solving the point location problem by deleting an associated cavity in a clever way.

First, suppose we have computed the unconstrained Delaunay triangulation $\mathcal D\mathcal T(\mathcal P)$ of our point set. The algorithm for then computing the constrained Delaunay triangulation $\mathcal D\mathcal T(\mathcal P, \mathcal E, \mathcal B)$ then works incrementally, inserting edges one at a time. Let us, then, describe the procedure for inserting some edge $e \in \mathcal C$, where $\mathcal C = \mathcal E \cup \mathcal B$. We break the discussion into sections. 

## Segment Location 

Similar to how in adding points into a triangulation the first step is point location, here the first step is _segment location_. Here, our aim is to find all triangles that intersect the edge $e$. This is easy to do if we simply remember that the jump-and-march algorithm walks through all triangles that 
intersect an initial scan line until stopping: this is exactly what we want. So, what we have done is modify our jump-and-march code such that the history of triangles walked through is recorded. This information is recorded into a `PointLocationHistory` struct. With this, we also note that since the jump-and-march algorithm will require rotating around an initial point (one of the indices of $e$), we should want to minimise the number of triangles we may need to rate around. Thus, we rotate $e$ such that its initial vertex is the one with the least degree, hence we have less triangles to rotate around initially, thus reducing the time spent searching. From here on, we let $e$ be this rotated form. This segment location is handled via the following function:

```@docs 
locate_intersecting_triangles
```

This function also returns information about any segments that are collinear with `e` and all vertices of the intersecting triangles to the left and to the right of `e`. For the collinear segments, these are processed via the functions `fix_segments!`, `connect_segments!`, and `extend_segments!`, breaking `e` into a smaller set of segments so that no segments are collinear anymore. We will assume that there are no collinear segments for simplicity. For the vertices to the left and to the right of `e`, these are needed as they define the outline of points to be deleted on each side of `e`, thus giving a polygonal cavity (possibly self-intersecting, but this detail doesn't actually matter) that we can triangulate individually. For example, consider the triangulation:

```julia 
a = (0.0, 0.0)
b = (0.0, 1.0)
c = (0.0, 2.5)
d = (2.0, 0.0)
e = (6.0, 0.0)
f = (8.0, 0.0)
g = (8.0, 0.5)
h = (7.5, 1.0)
i = (4.0, 1.0)
j = (4.0, 2.5)
k = (8.0, 2.5)
pts = [a, b, c, d, e, f, g, h, i, j, k]
tri = triangulate(pts; delete_ghosts=false, randomise=false)
fig, ax, sc = triplot(tri)
lines!(ax, [get_point(tri, 2, 7)...], color=:blue, linewidth=2)
```

```@raw html
<figure>
    <img src='../figs/segment_example.png', alt='An edge through a triangulation'><br>
</figure>
```

When we perform segment location on this example on the highlighted segment `(2, 7)`, we find:

```julia-repl
julia> e = (2, 7)
(2, 7)

julia> intersecting_triangles, collinear_segments, left_vertices, right_vertices = DelaunayTriangulation.locate_intersecting_triangles(tri, e);

julia> intersecting_triangles
8-element Vector{Tuple{Int64, Int64, Int64}}:
 (4, 3, 2)
 (3, 4, 10)
 (10, 4, 9)
 (9, 4, 5)
 (9, 5, 10)
 (10, 5, 8)
 (8, 5, 6)
 (8, 6, 7)

julia> collinear_segments
Tuple{Int64, Int64}[]

julia> left_vertices
7-element Vector{Int64}:
  7
  8
 10
  9
 10
  3
  2

julia> right_vertices
5-element Vector{Int64}:
 2
 4
 5
 6
 7
```

The intersecting triangles gives the triangles intersected by `e` in order of occurrence. For `left_vertices`, these are given in counter-clockwise order, with the left- and right-most elements being the indices of `e`. See that the vertex `10` is repeated in `left_vertices`. This is because the point `9` creates a sort of dangling edge in the polygonal cavity that we have to delete, and we need to somehow know to insert this point back into the triangulation. What we do, then, is to imagine an ant walking around the polygonal cavity. The ant will walk from `10` to `9`, but then it has to come back down to `10`, so we include it twice to represent these two visits. To see the cavity, let us delete these triangles:

```julia
DelaunayTriangulation.delete_intersected_triangles!(tri, intersecting_triangles)
```

```@raw html
<figure>
    <img src='../figs/segment_example_deleted_triangles.png', alt='An edge through a triangulation with excavated cavities'><br>
</figure>
```

The polygonal cavities on each side of the blue segment are what we need to re-triangulate separately. Let us now describe this procedure.

## Triangulating the Polygonal Cavities

Now we need to triangulate each cavity. We will describe this only for the cavity above `(2, 7)` in the figure below, but the procedure is exactly the same for the other cavity. Let $\mathcal V = (v_1, \ldots, v_m)$ be the sequence of vertices in counter-clockwise order around the cavity when we insert the segment $e=(v_1,v_m)$. In this case, $\mathcal V = (7, 8, 10, 9, 10, 3, 2)$. First, just like in Chew's algorithm, we need to build up a linked-list representing the cavity. This is done via the function 

```@docs 
prepare_vertex_linked_list
```

(Note: The algorithm given by Shewchuk and Brown linked above also allocates a `distance` array storing `orient` determinants for this preparation of the linked list, giving values proportional to the distance from $e$. This is not exactly robust for us, since we need to compute the sign of a difference of two determinants. In particular, let $o_1$ and $o_2$ be two robust estimates for an `orient` determinant. To determine if `o_1 < o_2` is the same as defining a predicate for the sign of `o_1 - o_2`, but this is problematic as, while the computation of `o_1` and `o_2` may be reliable, their difference is not. Instead, we recompute this predicate in a robust manner each time, trading performance for robustness. This predicate is defined by `point_closest_to_line`.)

Once we have prepared the linked list, we need to delete vertices from it in a random order, corresponding to deleting vertices from the polygon in a random order. Like in Chew's algorithm, this is done in such a way that we can reverse the process and automatically get point location without ever needing the jump-and-march algorithm. This deletion is handled via `delete_polygon_vertices_in_random_order!` which simply loops over each vertex, calling `select_random_vertex` and `update_vertex_linked_list!` at each iteration:

```@docs 
delete_polygon_vertices_in_random_order!
select_random_vertex 
update_vertex_linked_list! 
```

An important note is that, while this insertion algorithm works even with dangling edges, self-intersections, etc., it does not work when a point has an interior angle exceeding 360 degrees. We can detect this case by finding a point in the polygon that is closer to $e$ than its two neighbours, which is the only time such an interior angle is possible (see the paper for a proof). Thus, `select_random_vertex` actually keeps sampling vertices to delete in a random order until this is not the case, making use of `vertex_is_closer_than_neighbours`:

```@docs 
vertex_is_closer_than_neighbours 
```

Once this is all done, we ready to start triangulating the cavity. After adding an initial triangle from the three remaining vertices, we add points in one at a time, making use of a function `add_point_cavity_cdt!`. This function is defined by:

```julia
function add_point_cavity_cdt!(tri::Triangulation, u, v, w)
    x = get_adjacent(tri, w, v)
    if !edge_exists(x)
        insert_flag = true
    else
        p, q, r, s = get_point(tri, w, v, x, u) 
        incircle_test = point_position_relative_to_circle(p, q, r, s)
        orient_test = triangle_orientation(tri, u, v, w)
        insert_flag = !is_inside(incircle_test) && is_positively_oriented(orient_test)
    end
    if insert_flag
        add_triangle!(tri, u, v, w; protect_boundary=true, update_ghost_edges=false)
    else
        delete_triangle!(tri, w, v, x; protect_boundary=true, update_ghost_edges=false)
        add_point_cavity_cdt!(tri, u, v, x)
        add_point_cavity_cdt!(tri, u, x, w)
    end
    return nothing
end
```

In particular, `x = get_adjacent(tri, w, v)` is used to find the triangle on the other side of the edge $vw$ from $u$, the point being inserted. The only way that $wv$ should not be deleted is if the triangle $wvx$ does not exist, as detected via `!edge_exists(x)`, or if `u` is not inside the circumcircle of $wvx$ and $u$ is on the correct side of the edge $vw$. This is detected by the computation of `insert_flag`. If `insert_flag`, just add the triangle. Otherwise, we need to delete the triangle $wvx$ and $uvw$ as they are no longer constrained Delaunay. This is done by flipping $vw$ onto $ux$. Once we have done this for each point, we have successfuly triangulated the cavity. 

This is all handled via the `triangulate_cavity_cdt` function:

```@docs 
triangulate_cavity_cdt
```

For our example, what we find is (the triangulation `tri` was updated to put the missing triangles back in from the last piece of code, note):

```julia 
julia> add_edge!(tri, 2, 7) # calls triangulate_cavity_cdt on each cavity
```

```@raw html
<figure>
    <img src='../figs/segment_example_completed.png', alt='A constrained edge through a triangulation'><br>
</figure>
```

Here we used `add_edge!(tri, 2, 7)`, which does all this pre-processing for us. Similarly, for adding many edges, the method 

```@docs 
triangulate_constrained!
``` 

is useful (`triangulate` calls this internally).

## Excavating Exterior Faces 

When we define boundary curves, we typically want to delete any points and triangles exterior to them. The logic of the method we use for this is simple. Basically, we "plant" a seed in an exterior face, and use it to infect other points in this exterior face, continuing this spread until all exterior faces are found. (The actual implementation is a lot more involved, obviously). The function that performs this is `delete_holes!`, with relevant docstrings below.

```@docs 
delete_holes!
has_interiors_within_interiors 
find_all_points_to_delete
find_all_triangles_to_delete 
delete_all_exterior_triangles 
clear_deleted_points!
delete_remaining_triangles_connecting_boundary_edges!
```

Let's show how this all works by example. Let us first build a triangulation.

```julia
using DelaunayTriangulation, CairoMakie
curve_1 = [[
    (0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (12.0, 0.0), (12.0, 4.0),
    (12.0, 8.0), (14.0, 10.0), (16.0, 12.0), (16.0, 16.0),
    (14.0, 18.0), (12.0, 20.0), (12.0, 24.0), (12.0, 28.0),
    (8.0, 28.0), (4.0, 28.0), (0.0, 28.0), (-2.0, 26.0), (0.0, 22.0),
    (0.0, 18.0), (0.0, 10.0), (0.0, 8.0), (0.0, 4.0), (-4.0, 4.0),
    (-4.0, 0.0), (0.0, 0.0),
]]
curve_2 = [[
    (4.0, 26.0), (8.0, 26.0), (10.0, 26.0), (10.0, 24.0),
    (10.0, 22.0), (10.0, 20.0), (8.0, 20.0), (6.0, 20.0),
    (4.0, 20.0), (4.0, 22.0), (4.0, 24.0), (4.0, 26.0)
]]
curve_3 = [[(4.0, 16.0), (12.0, 16.0), (12.0, 14.0), (4.0, 14.0), (4.0, 16.0)]]
curve_4 = [[(4.0, 8.0), (10.0, 8.0), (8.0, 6.0), (6.0, 6.0), (4.0, 8.0)]]
curves = [curve_1, curve_2, curve_3, curve_4]
points = [
    (2.0, 26.0), (2.0, 24.0), (6.0, 24.0), (6.0, 22.0), (8.0, 24.0), (8.0, 22.0),
    (2.0, 22.0), (0.0, 26.0), (10.0, 18.0), (8.0, 18.0), (4.0, 18.0), (2.0, 16.0),
    (2.0, 12.0), (6.0, 12.0), (2.0, 8.0), (2.0, 4.0), (4.0, 2.0),
    (-2.0, 2.0), (4.0, 6.0), (10.0, 2.0), (10.0, 6.0), (8.0, 10.0), (4.0, 10.0),
    (10.0, 12.0), (12.0, 12.0), (14.0, 26.0), (16.0, 24.0), (18.0, 28.0),
    (16.0, 20.0), (18.0, 12.0), (16.0, 8.0), (14.0, 4.0), (14.0, -2.0),
    (6.0, -2.0), (2.0, -4.0), (-4.0, -2.0), (-2.0, 8.0), (-2.0, 16.0),
    (-4.0, 22.0), (-4.0, 26.0), (-2.0, 28.0), (6.0, 15.0), (7.0, 15.0),
    (8.0, 15.0), (9.0, 15.0), (10.0, 15.0), (6.2, 7.8),
    (5.6, 7.8), (5.6, 7.6), (5.6, 7.4), (6.2, 7.4), (6.0, 7.6),
    (7.0, 7.8), (7.0, 7.4)]
boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points)
tri = triangulate(points; boundary_nodes, delete_ghosts = false, delete_holes = false)
fig = Figure()
ax = Axis(fig[1, 1], xlabel=L"x",ylabel=L"y",width=600,height=300)
triplot!(ax, tri)
```

```@raw html
<figure>
    <img src='../figs/intermediate_triangulation.png', alt='Triangulation prior to excavation'><br>
</figure>
```

The aim is to delete all triangles outside of the outermost magenta boundary and inside the interior magenta boundaries. First, let us using the boundaries to find points that are next to boundary edges away from the interior faces, accomplished via `find_all_points_to_delete`.

```julia
points_to_delete = DelaunayTriangulation.find_all_points_to_delete(tri)
scatter!(ax, [get_point(tri, i) for i in points_to_delete], color=:blue,markersize=14)
```

```@raw html
<figure>
    <img src='../figs/intermediate_triangulation_highlighted.png', alt='Triangulation prior to excavation with highlighted points for deletion'><br>
</figure>
```

Now let us use these points to find the triangles to delete.

```julia
triangles = DelaunayTriangulation.find_all_triangles_to_delete(tri, points_to_delete)
[poly!(ax, [get_point(tri, u, v, w, u)...], color = (:red, 0.2)) for (u, v, w) in triangles if !DelaunayTriangulation.is_ghost_triangle(u,v,w)]
```

```@raw html
<figure>
    <img src='../figs/intermediate_triangulation_highlighted.png', alt='Triangulation prior to excavation with highlighted points and triangles for deletion'><br>
</figure>
```

We see that all of the triangles are correct in this case, but it's missing some triangles in the top-most interior hole. This is clearer once we delete all the highlighted triangles:

```julia
DelaunayTriangulation.delete_all_exterior_triangles(tri, triangles)
```

```@raw html
<figure>
    <img src='../figs/intermediate_triangulation_highlighted.png', alt='Triangulation partially excavated'><br>
</figure>
```

These remaining triangles to be deleted are handled via `delete_remaining_triangles_connecting_boundary_edges!`.