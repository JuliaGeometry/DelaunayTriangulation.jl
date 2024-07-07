# Clipped Voronoi Tessellations

Now we consider a variant of the Voronoi tessellation called the _clipped Voronoi tessellation_. In the clipped Voronoi tessellation, the Voronoi polygons are clipped to the convex hull of the point set. This is useful when we want to ensure that the Voronoi polygons are bounded and do not extend to infinity. The computation of this tessellation is much more involved than the standard Voronoi tessellation. To be exact, we are interested in computing $\tilde{\mathcal V}(\mathcal P) \equiv \mathcal V(\mathcal P) \cap \mathcal C\mathcal H(\mathcal P)$. An example of a clipped Voronoi tessellation is shown below. We will revisit this example throughout the discussion of the algorithm.

```@setup vornclip11
using DelaunayTriangulation 
using CairoMakie 
using StableRNGs
A, B, C, D, E, F, G, H, I, J, K, L = 
(-6.0, 4.0), (-2.0, 5.0), (0.92, 2.86), (0.0, -2.0), (-5.0, -4.0), (-7.0, -2.0), (-10.0, 2.0), (-4.0, 3.0), (-5.0, 2.0), (-2.0, 1.0), (-3.0, -2.0), (4.0, 1.0)
pts = [A, B, C, D, E, F, G, H, I, J, K, L]
tri = triangulate(pts)
vorn = voronoi(tri)
vorn2 = voronoi(tri, clip=true)
fig = Figure(fontsize = 24)
ax1 = Axis(fig[1, 1], width = 300, height = 300, title = "Unclipped", titlealign = :left)
voronoiplot!(ax1, vorn)
ax2 = Axis(fig[1, 2], width = 300, height = 300, title = "Clipped", titlealign = :left)
voronoiplot!(ax2, vorn2)
xlims!(ax1, -11, 5)
xlims!(ax2, -11, 5)
ylims!(ax1, -5, 6)
ylims!(ax2, -5, 6)
resize_to_layout!(fig)
```

```@example vornclip11
fig # hide
```

At the end of this section, we also discuss the intersection of $\mathcal V(\mathcal P)$ with a rectangle.

## Computing intersections of the Voronoi polygons and the convex hull

The main complication with clipping the Voronoi polygons to the convex hull is in computing all the intersections of the polygons with the convex hull. For this step, we use an iterative approach. The steps are described below. We will start by describing the algorithm in words, and then provide a clearer example.

### Initialising the queue

We maintain a queue of edges to determine where we need to process the intersections. First, we initialise a set $\mathcal E$ containing all boundary edges. We take a random edge $e \in \mathcal E$ and use it to initialise a queue $\mathcal Q$ of polygon edges. To be specific, into $\mathcal Q$ we enqueue $(e, \mathcal V)$, where $\mathcal V_i$ is the polygon incident to the boundary edge $e$; we find this incident polygon by finding the nearest neighbour to the edge's midpoint.

### Processing an edge

If either of $\mathcal Q$ or $\mathcal E$ are not empty, we can consider processing an edge for intersections. If $\mathcal Q = \emptyset$, we take the next edge in $\mathcal E$ and enqueue it and its incident polygon into $\mathcal Q$. Now, dequeue a pair $(e, \mathcal V)$ from $\mathcal Q$. If the pair $(e, \mathcal V)$ has already been processed, we skip this step and dequeue the next pair from $\mathcal Q$.

For the first step of this processing, we need to find any intersections of $e$ with the polygon $\mathcal V$. To do this, we first find the edges $e_\ell$ and $e_r$ left and right of $e$, respectively. We then need to consider each edge of $\mathcal V$. Let $s$ be an edge of $\mathcal V$. How $s$ intersects the edges depends on whether (1) it is an unbounded ray oriented so that it is going out to infinity, (2) it is an unbounded ray oriented so that it is coming in from infinity, or (3) it is a bounded edge. The first two cases are reasonably straightforward - just compute the intersection of the edge with the unbounded rays. For the third case, we compute all the intersections of $e$, $e_\ell$, and $e_r$ with $s$. The reason for needing to check all of $e$, $e_\ell$, and $e_r$ rather than $e$ alone is that there may be issues finding intersections when considering Voronoi polygons near corners of the convex hull.

Once we have processed all of the edges of $\mathcal V$, we store the intersections of the edges with $e$, $e_\ell$, and $e_r$ into the sets $\mathcal E(\mathcal V)$, $\mathcal L(\mathcal V)$, and $\mathcal R(\mathcal V)$, respectively. We use these intersections to determine what edges need to be enqueued into $\mathcal Q$, taking care that we do not miss any corner points. Consider some intersection point $p \in \mathcal L(\mathcal V)$; the cases for the other sets are similar. If $(e_\ell, \mathcal V)$ has not been processed, then we enqueue $(e_\ell, \mathcal V_i)$ and $(e_\ell, \mathcal V_j)$ into $\mathcal Q$, where $v_i$ and $v_j$ are the vertices of $e_\ell$. This ensures that we will find the intersections next to this polygon. After enqueueing these pairs, we still need to protect against corner points. To do this, we note that a corner point only needs to be checked if the incident polygon $\mathcal V$ belongs to the boundary of the convex hull. If it is, then let $v_u$ be the shared vertex between $e_i$ and $e_\ell$. If $v_u$ is exactly the generator of $\mathcal V$, then we have a corner point, and so we need to add this corner point into the set of intersections since it will necessarily be included in the clipped Voronoi tessellation. Once we have processed all the intersections for $e$, $e_\ell$, and $e_r$ as above, we need to then also enqueue all of the edges and their adjacent incident poygons into $\mathcal Q$, provided that these pairs have not already been processed.

Let us give an example that illustrates this part of the procedure clearly.

```@setup vornclip11
tri = triangulate(pts)
lock_convex_hull!(tri)
vorn = voronoi(tri, clip=false)
edges_to_process,
polygon_edge_queue,
boundary_sites,
segment_intersections,
processed_pairs,
intersected_edge_cache,
exterior_circumcenters,
left_edge_intersectors,
right_edge_intersectors,
current_edge_intersectors,
equal_circumcenter_mapping = DelaunayTriangulation.initialise_clipping_arrays(vorn)
e = DelaunayTriangulation.convert_to_edge_adjoining_ghost_vertex(vorn, first(edges_to_process))
DelaunayTriangulation.enqueue_new_edge!(polygon_edge_queue, vorn, e)
e, incident_polygon = popfirst!(polygon_edge_queue)
push!(processed_pairs, (e, incident_polygon))
for cache in (intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors)
    empty!(cache)
end
left_edge, right_edge, e = DelaunayTriangulation.process_polygon!(vorn, e, incident_polygon, boundary_sites, segment_intersections, intersected_edge_cache, exterior_circumcenters, equal_circumcenter_mapping)
e1 = deepcopy(e)
left_edge1 = deepcopy(left_edge)
right_edge1 = deepcopy(right_edge)
segment_intersections1 = copy(segment_intersections)
boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping = DelaunayTriangulation.find_all_intersections(vorn)


fig = Figure(fontsize = 28)
ax1 = Axis(fig[1, 1], width = 300, height = 300)
voronoiplot!(ax1, vorn, show_generators = false)
lines!(ax1, [get_generator(vorn, get_convex_hull_vertices(tri)...)...], color = :magenta, linestyle = :dash, linewidth = 3)
lines!(ax1, [get_generator(vorn, e1...)...], color = :red, linewidth = 3)
lines!(ax1, [get_generator(vorn, left_edge1...)...], color = :green, linewidth = 3)
lines!(ax1, [get_generator(vorn, right_edge1...)...], color = :blue, linewidth = 3)
scatter!(ax1, pts, color = :black, markersize = 14)
text!(ax1, [(-5.7, -3.1)], text = [L"e"], color = :white)
text!(ax1, [(-1.5, -3.3)], text = [L"e_r"], color = :white)
text!(ax1, [(-9.0, -1.8)], text = [L"e_{\ell}"], color = :white)
text!(ax1, [(-10.0, -4.5)], text = [L"\mathcal{V}"], color = :white)
scatter!(ax1, segment_intersections1, color = :orange, markersize = 14)
ax2 = Axis(fig[1, 2], width = 300, height = 300)
voronoiplot!(ax2, vorn, show_generators = false)
lines!(ax2, [get_generator(vorn, get_convex_hull_vertices(tri)...)...], color = :magenta, linestyle = :dash, linewidth = 3)
scatter!(ax2, pts, color = :black, markersize = 14)
scatter!(ax2, segment_intersections, color = :orange, markersize = 14)
resize_to_layout!(fig)
```

```@example vornclip11
fig # hide
```

In the first figure above, we are considering the processing of an edge $e$. The Voronoi polygon $\mathcal V$ we consider incident to $e$ is shown, obtained by finding that the midpoint of $e$ is contained in $\mathcal V$. By just processing the intersections of $\mathcal V$, we find the intersections with $\mathcal C\mathcal H(\mathcal P)$ (shown in magenta) shown in orange. This alone is not enough, though, as we can see that we don't identify that the black dot between $e_\ell$ and $e$ should be included in this intersection. This dot is an example of a corner point we were discussing previously, showing the need for this extra processing step. Note also from the above figure that  unbounded polygons $\mathcal V$ alone are not sufficient for checking all intersections, as we can see that some of the bounded polygons also intersect with the convex hull. Eventually, after processing all edges in this way, we obtain the set of orange points shown in the second figure.

### Clipping the polygons

Once all the intersections have been computed, clipping the polygons is straightforward. For each $\mathcal V$ for which intersections were found, we store with it all intersections found between it and the convex hull. We then remove all vertices from $\mathcal V$ that are outside of the domain and then add in these intersection points, sorting $\mathcal V$ as needed so that it remains a convex polygon. Doing this for all $\mathcal V$ completes the clipping, giving the tessellation shown at the start of this section.

## Clipping a Voronoi polygon to a rectangle

Now let us describe clipping to some rectangle $\mathcal R = [a, b] \times [c, d]$ rather than to the convex hull. The procedure here is very different, and rather than processing the entire tessellation at once we apply the procedure to individual polygons. Before we discuss this procedure in general, we need to discuss two algorithms: The Sutherland-Hodgman algorithm and the Liang-Barsky algorithm.

### The Sutherland-Hodgman algorithm

The Sutherland-Hodgman algorithm is an algorithm for clipping a polygon (called the subject polygon) to a convex polygon (called the clip polygon). The algorithm works by essentially extending each edge of the clip polygon to infinity and then iteratively clipping the subject polygon to each edge of the clip polygon. Let us give an example of this.

```@setup suthodg
using DelaunayTriangulation
using CairoMakie  
using Random
function to_rand_point_verts(___points)
    _points = [Tuple(rand(2)) for _ in 1:500]
    _points = [_points; ___points]
    shuffle!(_points)
    vertices = identity.(indexin(___points, _points)) # identity to convert eltype from Union{Nothing, Int}
    return _points, vertices
end
a = (-4.0, 4.0)
b = (-1.0, 6.0)
c = (3.0, 6.0)
d = (4.0, 4.0)
e = (4.0, -1.0)
f = (2.0, -3.0)
g = (-2.94, -1.32)
points = [g, f, e, d, c, b, a] # ccw
npoints, nvertices = to_rand_point_verts(deepcopy(points))
h = (-2.0, 7.0)
i = (-5.0, 6.0)
j = (-5.0, 2.0)
k = (-4.0, -2.0)
ℓ = (-1.0, -3.0)
m = (2.0, 2.0)
n = (1.0, 5.0)
clip_points = [h, i, j, k, ℓ, m, n]
nclip_points, nclip_vertices = to_rand_point_verts(deepcopy(clip_points))
result = DelaunayTriangulation.clip_polygon(nvertices, npoints, nclip_vertices, nclip_points)
fig = Figure()
ax1 = Axis(fig[1, 1], width = 300, height = 300)
lines!(ax1, [points; g], color = :black)
lines!(ax1, [clip_points; h], color = :red)
poly!(ax1, result, color = (:blue, 0.3))
hidedecorations!(ax1)
resize_to_layout!(fig)
```

```@example suthodg
fig # hide
```

In the figure above, the black polygon shows the subject polygon and the red polygon is the clip polygon. Our aim is to clip the subject polygon to the clip polygon to obtain the blue polygon shown. Letting $\mathcal P_S = \{p_1, \ldots, p_n\}$ and $\mathcal P_C = \{q_1, \ldots, q_m\}$ be the subject and clip polygons listed in counter-clockwise order with $p_1 \neq p_n$ and $q_1 \neq q_m$, respectively, the procedure for this clipping as follows:

1. Define $\mathcal O = \mathcal P_S$ and let $q = q_n$. For each $p \in \mathcal P_C$, we first let $\mathcal I = \mathcal O$ and then reset $\mathcal O = \emptyset$.
2. With this $p$, let $s = \mathcal I_r$, where $|\mathcal I| = r$ and $\mathcal I_j$ denotes the $j$th element of $\mathcal I$. For $i=1,2,\ldots,r$, we consider the edge $e = (s, \mathcal I_i)$. If $\mathcal I_i$ is left of $\overrightarrow{qp}$, then it is outside of $\mathcal P_S$, and so we need to check if there is any intersection, i.e. if $e$ intersects $\overrightarrow{qp}$, which can be easily checked by considering the position of $s$ relative to $\overrightarrow{qp}$. If  $\mathcal I_i$ is not left of $\overrightarrow{qp}$ but $s$ is left, then there must be an intersection of $e$ with $\overrightarrow{qp}$ since $s$ and $\mathcal I_i$ are on opposite sides of $\overrightarrow{qp}$. In either of these cases, the intersection that we find gets pushed into $\mathcal O$, and we then set $s = \mathcal I_i$ and continue onto the next $i$.
3. Once we have processed each $\mathcal I_i$, we set $q = p$ and then continue onto the next $p$ in $\mathcal P_C$.
4. Finally, once each vertex in $\mathcal P_C$ has been processed, the final polygon is defined by the points in $\mathcal O$.

Let's visualise this procedure.

```@setup suthodg
fig = Figure()
ax = Axis(fig[1, 1], width = 300, height = 300)
hidedecorations!(ax)
poly!(ax, [clip_points; h], color = (:red, 0.3))
lines!(ax, [points; g], color = :black, linewidth = 3)
xlims!(ax, -6, 5)
ylims!(ax, -4, 8)
vertices = nvertices 
points = npoints 
clip_vertices = nclip_vertices
clip_points = nclip_points
poly = DelaunayTriangulation.Polygon(vertices, points)
clip_poly = DelaunayTriangulation.Polygon(clip_vertices, clip_points)
orig_clip = deepcopy(clip_poly)
output_list = poly
q = clip_poly[end]
j = 1
t = LinRange(-5, 5, 2500)
for p in clip_poly  
    global j, q, output_list
    j += 1
    ax = Axis(fig[j > 4 ? 2 : 1, mod1(j, 4)], width = 300, height = 300)
    poly!(ax, [orig_clip; orig_clip[1]], color = (:red, 0.3))
    extended_line = [p .+ t .* (p .- q) for t in t]
    lines!(ax, extended_line, color = :blue)
    input_list = output_list 
    isempty(output_list) && break
    output_list = DelaunayTriangulation.clip_polygon_to_edge(input_list, q, p)
    lines!(ax, [output_list; output_list[1]], color = :black, linewidth = 3)
    xlims!(ax, -6, 5)
    ylims!(ax, -4, 8)
    q = p
    hidedecorations!(ax)
end
resize_to_layout!(fig)
```

```@example suthodg
fig # hide
```

In the figure above, we show the individual steps of this algorithm. In the second panel, the blue line shows the extended edge of the clip polygon that we use to slice the subject polygon, clipping it onto the edge. For the next four panels, the blue line never touches the subject polygon, and so nothing happens. The last two panels show the last two clips needed to obtain the final polygon.

### The Liang-Barsky algorithm

Now we describe the Liang-Barsky algorithm, an algorithm for clipping a line segment to a rectangle. Let us take a line $\vb p(t) = \vb p_0 + t\vb d$, where $\vb p_0 = (x_0, y_0)$, $\vb d = (\Delta x, \Delta y)$, and $0 \leq t \leq 1$. A point $(x, y) = \vb p(t)$ is in the rectangle if (1) $a \leq x_0 + t \Delta x \leq b$ and (2) $c \leq y_0 = t \Delta y \leq d$, or equivalently $tp_i \leq q_i$ for $i=1,2,3,4$, where:

1. _Left edge_: $p_1 = -\Delta x$, $q_1 = x_0 - a$.
2. _Right edge_: $p_2 = \Delta x$, $q_2 = b - x_0$.
3. _Bottom edge_: $p_3 = -\Delta y$, $q_3 = y_0 - c$.
4. _Top edge_: $p_4 = \Delta y$, $q_4 = d - y_0$.

Using these inequalities, we can efficiently compute the intersections by processing each side of the rectangle at a time. Starting with $t_1 = 0$ and $t_2 = 1$ defining the current interval for the intersections, for each edge we do the following: Compute the $p_i$ and $q_i$ associated with the edge, and then $r_i = q_i / p_i$. This $r_i$ gives the parameter value for the intersection point of the line and the current edge (possibly extended outside of $0 \leq t \leq 1$). There are three cases to consider:

1. If $p_i = 0$ and $q_i < 0$, then the line is parallel with the edge but outside of the rectangle, and so there are no intersections of the line with the rectangle.
2. If $p_i > 0$ and $r < t_1$, then the the line enters the rectangle earlier than what is proposed by the current interval $[t_1, t_2]$, and so the line is outside of the rectangle and we have no intersections to consider. If $r < t_2$, then we update the intersection interval and let $t_2 = r$.
3. If $p_i < 0$ and $r > t_2$, then just like above we see that the line is outside of the rectangle and so we return no intersections. If $r > t_1$, set $t_1 = r$.

Applying these cases to each edge one at a time, updating $t_1$ and $t_2$ as required, gives the intersection parameters $t_1$ and $t_2$ defining the intersection of the line with the rectangle. If we ever did exit early, then we do not return any intersections.

### Clipping a bounded Voronoi polygon

Now that we have two important algorithms for clipping, let us begin our discussion on clipping the Voronoi polygons in particular, starting with bounded Voronoi polygons. This case is simple - simply apply the Sutherland-Hodgman algorithm to the polygon, using the clip polygon as the rectangle.

### Clipping an unbounded Voronoi polygon

The case of clipping an unbounded Voronoi polygon is more complicated. Rather than trying to come up with an effective way to clip the unbounded rays of the polygon to a rectangle and taking care of all the possible edge cases, we use a more direct approach. Our aim is to convert the unbounded polygon into a finite polygon such that its intersection with the rectangle is the same as the intersection of the unbounded polygon with the rectangle. We do this in steps.

1. First, for the two unbounded edges, we let $u_m$ and $v_m$ be the two midpoints of the edges of the convex hull that the unbounded edges go through.
2. Next, we compute the maximum distance of the box to $u_m$ and $v_m$, letting $m = \max\{\operatorname{dist}(u_m, \mathcal R), \operatorname{dist}(v_m, \mathcal R)\}$.

    Now, starting with `inside = true` and $t = 1$, we do the following until `inside = false`.

3. Replace $t$ by $2t$, and compute $p = u_m + td_1$, where $d_1$ is a unit vector in the direction of the unbounded edge associated with $u_m$, and similarly $q = v_m + td_2$. 
4. Apply the Liang-Barsky algorithm to $\mathcal R$ and the line segment through $p$ and $q$, and check if the line segment is completely outside of $\mathcal R$. If it is, then set `outside = true`.
5. We then need to be careful about the case where the generator associated with the polygon is outside of $\mathcal R$. In this case, the unbounded edge might start outside of the rectangle and eventually find its way inside. To avoid this, we be conservative and check that the length of each ray is greater than the maximum distance from $u_m$ and $v_m$ to the clip rectangle. In particular, compute $\delta_1 = \|p - u_m\|$ and $\delta_2 = \|q - v_m\|$. If $\min\{\delta_1, \delta_2\} < m$, then the edge might possibly be inside $\mathcal R$. If this is the case, we let `inside = true` and continue, and otherwise we let `inside = false`.
6. Once we stop iterating, the final values of $p$ and $q$ are the points that the unbounded edges will now stop at, thus defining a bounded polygon.

The reason we start with $t = 1$ instead of, say, $t = 1$ is so that we start at $t = 2$ in the loop and avoid any possibility of duplicated vertices for polygons completely outside of the box. Note also that an important reason that we need to apply the Liang-Barsky algorithm to the line segment through $p$ and $q$ is so that we can join the two unbounded edges together and be confident that this line is completely outside of this rectangle.  Let's make this polygon growing procedure clearer with an example, where we also pick a complicated example that demonstrates the need for the conservative check in the fifth step.

```@setup suthodg 
points = [(-3,7),(1,6),(-1,3),(-2,4),(3,-2),(5,5),(-4,-3),(3,8)]
to_fl64tup(x) = Float64.(x)
points = to_fl64tup.(points)
bounding_box = (0.0, 5.0, -15.0, 15.0)
i = 7
tri = triangulate(points)
vorn = voronoi(tri)
fig = Figure()
ax = Axis(fig[1, 1], width = 300, height = 300)
voronoiplot!(ax, vorn, show_generators = false, clip = (-20.0, 20.0, -20.0, 20.0))
lines!(ax, [0.0, 5.0, 5.0, 0.0, 0.0], [-15.0, -15.0, 15.0, 15.0, -15.0], color = :red)
scatter!(ax, [points[i]], color = :black, markersize = 14)
xlims!(ax, -20, 8)
ylims!(ax, -20, 16)
hideydecorations!(ax)
resize_to_layout!(fig)
```

```@example suthodg
fig # hide
```

In the figure above, the polygon we are interested is the one corresponding to the black dot, and we want to clip this polygon to the red rectangle. The following figures show how we grow the polygon's unbounded edges to begin.

```@setup suthodg 
a, b, c, d = bounding_box
vertices = get_polygon(vorn, i)
new_vertices, new_points, ghost_vertices = DelaunayTriangulation.get_new_polygon_indices(vorn, vertices)
inside = true
t = 1.0 
u, v = ghost_vertices
u_m, u_r = DelaunayTriangulation._get_ray(vorn, i, u)
v_m, v_r = DelaunayTriangulation._get_ray(vorn, i, v)
u_mx, u_my = DelaunayTriangulation._getxy(u_m)
u_rx, u_ry = DelaunayTriangulation._getxy(u_r)
v_mx, v_my = DelaunayTriangulation._getxy(v_m)
v_rx, v_ry = DelaunayTriangulation._getxy(v_r)
p = (0.0, 0.0)
q = (0.0, 0.0)
dist_to_box = DelaunayTriangulation.maximum_distance_to_box(a, b, c, d, u_m) # this is a squared distance
dist_to_box = max(dist_to_box, DelaunayTriangulation.maximum_distance_to_box(a, b, c, d, v_m))

fig = Figure(fontsize = 24)
for i in 1:5
    global p, q
    t = 2^(i-1)
    p = u_m .+ t .* (u_r .- u_m)
    q = v_m .+ t .* (v_r .- v_m)
    ax = Axis(fig[i > 3 ? 2 : 1, mod1(i, 3)], title = L"t = %$(t)", titlealign = :left, width = 300, height = 300)
    xlims!(ax, -30, 8)
    ylims!(ax, -30, 16)
    unbounded = get_polygon_coordinates(vorn, 7, (-30.0, 30.0, -30.0, 30.0))
    poly!(ax, unbounded, color = (:blue, 0.3))
    lines!(ax, [0.0, 5.0, 5.0, 0.0, 0.0], [-15.0, -15.0, 15.0, 15.0, -15.0], color = :red)
    lines!(ax, [u_m, p], color = :darkgreen, linewidth = 3)
    lines!(ax, [v_m, q], color = :darkgreen, linewidth = 3)
    lines!(ax, [p, q], color = :magenta, linewidth = 3)
    hidedecorations!(ax)
end
ax = Axis(fig[2, 3], width = 300, height = 300, title = "Final polygon", titlealign = :left)
xlims!(ax, -30, 8)
ylims!(ax, -30, 16)
poly!(ax, get_polygon_coordinates(vorn, 7, (-30.0, 30.0, -30.0, 30.0)), color = (:blue, 0.3))
new_points[u] = p
new_points[v] = q
new_vertices[u] = u
new_vertices[v] = v
poly!(ax, [new_points; new_points[1]], color = :darkgreen)
lines!(ax, [0.0, 5.0, 5.0, 0.0, 0.0], [-15.0, -15.0, 15.0, 15.0, -15.0], color = :red)
resize_to_layout!(fig)
```

```@example suthodg
fig # hide
```

In these figures, the blue polygon shows the Voronoi polygon, and the green edges show the approximations to the unbounded rays; the top-most ray starts away from the polygon since the midpoint is further behind the associated circumcenter in this case. The magenta line shows the line segment joining the two approximations - the aim is to grow the green edges long enough such that the magenta line is completely outside of the red rectangle. Let's analyse each figure.

1. At $t = 1$, the line segment is outside of the rectangle, but there is still room to grow the green edges such that they go through the clip rectangle (in particular, the right-most edge). Since the maximum distance of the green edges is not larger than the maximum distance from the green edge's origins to the clip rectangle, our procedure will continue to grow the edges.
2. At $t = 2$, we have the same problem as at $t=1$ in that the magenta edge is still outside of the clip rectangle but the green edges can still grow into the clip rectangle.
3. At $t = 4$, the magenta edge is finally inside of the clip rectangle, but since it intersects it we need to keep growing the green edges until the magenta edge leaves again.
4. The situation at $t=8$ is the same as at $t=4$.
5. Finally, at $t=18$, the magenta edge is completely outside of the clip rectangle, and the length of the two green edges is now greater than the maximum distance from the green edge's origins to the clip rectangle. Thus, we can stop growing the unbounded edges here.
6. The final figure shows the bounded polygon obtained by growing the unbounded edges (with the left side of it out of frame).

Once a bounded polygon representing the unbounded polygon has been obtained, we apply the Sutherland-Hodgman algorithm as before to clip it to the rectangle. For the example above, the polygon we obtain is shown below.

```@setup suthodg
fig = Figure(fontsize = 24)
ax = Axis(fig[1, 1], width = 300, height = 300)
poly!(ax, get_polygon_coordinates(vorn, 7, (a, b, c, d)), color = :darkgreen)
lines!(ax, [0.0, 5.0, 5.0, 0.0, 0.0], [-15.0, -15.0, 15.0, 15.0, -15.0], color = :red)
xlims!(ax, -30, 8)
ylims!(ax, -30, 16)
hidedecorations!(ax)
resize_to_layout!(fig)
```

```@example suthodg
fig # hide
```