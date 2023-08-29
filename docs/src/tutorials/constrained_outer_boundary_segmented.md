```@meta
EditURL = "https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/constrained_outer_boundary_segmented.jl"
```

# Constrained Triangulations
## Segmented Outer Boundary

In this tutorial, we again consider a triangulation with a constrained boundary. In
contrast to the previous tutorial, this outer boundary will be represented as a
chain of multiple paths or segments. This changes nothing geometrically, but it allows
for the identification of separate parts of a boundary. This is useful, for example,
if you want to assign different boundary conditions on different parts of the boundary
for a differential equation problem. To start, let us load in the packages we will need.

````@example constrained_outer_boundary_segmented
using DelaunayTriangulation
using CairoMakie
````

Now, we define some of the points we will be triangulating.

````@example constrained_outer_boundary_segmented
points = [
    (2.0, 8.0), (6.0, 4.0), (2.0, 6.0),
    (2.0, 4.0), (8.0, 2.0)
]
````

We now want to define our boundary. The method for providing a boundary to be
identified into multiple segments is to provided a vector of vectors of indices,
where each vector of indices is a segment. The last index of each segment must match
the first index of the next segment, including the last with the first
segment so that the boundary is closed. Here, we provide three segments.

````@example constrained_outer_boundary_segmented
segment_1 = [(0.0, 0.0), (14.0, 0.0)]
segment_2 = [(14.0, 0.0), (10.0, 4.0), (4.0, 6.0), (2.0, 12.0), (0.0, 14.0)]
segment_3 = [(0.0, 14.0), (0.0, 0.0)]
boundary_points = [segment_1, segment_2, segment_3]
````

We now convert these boundary points to indices using
`convert_boundary_points_to_indices`, and then we triangulate.
We also add a constrained edge.

````@example constrained_outer_boundary_segmented
E = Set(((6, 9),)) # (0, 0) → (4, 6)
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points; existing_points=points)
tri = triangulate(points; boundary_nodes, edges=E)
````

````@example constrained_outer_boundary_segmented
fig, ax, sc = triplot(tri, show_constrained_edges=true, constrained_edge_linewidth=6)
lines!(ax, segment_1, color=:red, linewidth=6)
lines!(ax, segment_2, color=:green, linewidth=6)
lines!(ax, segment_3, color=:blue, linewidth=6)
fig
````

The first segment is in red, the second segment is in green, and the third segment
is in blue. We use boundary indices to identify the segments, where the first segment is identified
by `-1`, the second by `-2`, and the third by `-3`.

Before we go into how the segments can be worked with, let us make a note regarding `get_constrained_edges(tri)`,
now that we have both constrained edges and boundary edges which might technically both be thought of
as being constrained edges. If we look at `get_constrained_edges(tri)`, we get:

````@example constrained_outer_boundary_segmented
get_constrained_edges(tri)
````

This is just the constrained edge we provided, and not the boundary edges. If we instead
want all constrained edges, considering both the boundary edges and the edges provided,
we can instead use `get_all_constrained_edges(tri)`.

````@example constrained_outer_boundary_segmented
get_all_constrained_edges(tri)
````

Let us now explore the several ways available for working with this boundary. First, if we just want to work
with the boundary edges without caring about the order, we can again use the boundary edge map.

````@example constrained_outer_boundary_segmented
get_boundary_edge_map(tri)
````

Remember that the keys are the edges, and the values are `Tuples` that give us (1) the segment index,
and (2) the position of the edge within that segment (more specifically, the position of the first vertex of
the edge). This would be useful if, for example, you don't care about the order of the edges, but you do
care about what segment the edge belongs to so that you can assign a boundary condition for example.
In the previous tutorial, we mentioned that the `Tuples` are of the form `(I, J)` so that the
corresponding edge is identified from `get_boundary_nodes(get_boundary_nodes(tri, I), J)`, but since both
`I` and `J` are integers in this case (since we have segments), we can just use `get_boundary_nodes(tri, (I, J))`.
(The former form is still the most general to support the case of a single boundary.) Here is an example of using
this map to compute the sum
```math
S = \sum_{(i, j) \in \mathcal E} f\left(\frac{x_i+x_j}{2}, \frac{y_i+y_j}{2}\right),
```
where $\mathcal E$ are the set of boundary edges, and
```math
f(x, y) = \begin{cases}
   1 & (x, y) \in \Gamma_1, \\
   \sin(x - y) & (x, y) \in \Gamma_2, \\
   \cos(x + y) & (x, y) \in \Gamma_3,
\end{cases}
```
and $\Gamma_i$ denotes the $i$th segment.

````@example constrained_outer_boundary_segmented
function segment_function(x, y, segment_index)
    f = if abs(segment_index) == 1
        1.0
    elseif abs(segment_index) == 2
        sin(x - y)
    else
        cos(x + y)
    end
    return f
end
function compute_sum(tri)
    bem = get_boundary_edge_map(tri)
    s = 0.0
    for (e, (segment_index, _)) in bem
        u, v = edge_indices(e)
        p, q = get_point(tri, u, v)
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = (px + qx) / 2, (py + qy) / 2
        s += segment_function(mx, my, segment_index)
    end
    return s
end
s = compute_sum(tri)
````

An alternative way to look at each segment is to use `get_adjacent2vertex` with the associated
boundary index.

````@example constrained_outer_boundary_segmented
get_adjacent2vertex(tri, -1)
````

````@example constrained_outer_boundary_segmented
get_adjacent2vertex(tri, -2)
````

````@example constrained_outer_boundary_segmented
get_adjacent2vertex(tri, -3)
````

Note that the provided edges are not in order, but this is helpful for considering specific segments. For
example, if we just wanted to compute the above sum over the second segment, we could do

````@example constrained_outer_boundary_segmented
function compute_sum_2(tri)
    edges = get_adjacent2vertex(tri, -2)
    s = 0.0
    for e in edges
        u, v = edge_indices(e)
        p, q = get_point(tri, u, v)
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = (px + qx) / 2, (py + qy) / 2
        s += segment_function(mx, my, -2)
    end
    return s
end
s = compute_sum_2(tri)
````

If your application instead wanted all the nodes on the segment rather than the edges,
you can look at the neighbours to the boundary index. For example, all the nodes
on the segment segment can be identified using

````@example constrained_outer_boundary_segmented
get_neighbours(tri, -2)
````

Another field is the `boundary_map`, which maps a given boundary index to the associated segment. This is
more so useful for internal methods, but you may sometimes need it.

````@example constrained_outer_boundary_segmented
get_boundary_map(tri)
````

In this case, the `i`th segment just has the boundary index `-i`, but this is typically used
to deal with the case of multiple boundaries so that we know where a boundary index belongs.
## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/new-docs/docs/src/literate_tutorials/constrained_outer_boundary_segmented.jl).

```julia
using DelaunayTriangulation
using CairoMakie

points = [
    (2.0, 8.0), (6.0, 4.0), (2.0, 6.0),
    (2.0, 4.0), (8.0, 2.0)
]

segment_1 = [(0.0, 0.0), (14.0, 0.0)]
segment_2 = [(14.0, 0.0), (10.0, 4.0), (4.0, 6.0), (2.0, 12.0), (0.0, 14.0)]
segment_3 = [(0.0, 14.0), (0.0, 0.0)]
boundary_points = [segment_1, segment_2, segment_3]

E = Set(((6, 9),)) # (0, 0) → (4, 6)
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points; existing_points=points)
tri = triangulate(points; boundary_nodes, edges=E)

fig, ax, sc = triplot(tri, show_constrained_edges=true, constrained_edge_linewidth=6)
lines!(ax, segment_1, color=:red, linewidth=6)
lines!(ax, segment_2, color=:green, linewidth=6)
lines!(ax, segment_3, color=:blue, linewidth=6)
fig

get_constrained_edges(tri)

get_all_constrained_edges(tri)

get_boundary_edge_map(tri)

function segment_function(x, y, segment_index)
    f = if abs(segment_index) == 1
        1.0
    elseif abs(segment_index) == 2
        sin(x - y)
    else
        cos(x + y)
    end
    return f
end
function compute_sum(tri)
    bem = get_boundary_edge_map(tri)
    s = 0.0
    for (e, (segment_index, _)) in bem
        u, v = edge_indices(e)
        p, q = get_point(tri, u, v)
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = (px + qx) / 2, (py + qy) / 2
        s += segment_function(mx, my, segment_index)
    end
    return s
end
s = compute_sum(tri)

get_adjacent2vertex(tri, -1)

get_adjacent2vertex(tri, -2)

get_adjacent2vertex(tri, -3)

function compute_sum_2(tri)
    edges = get_adjacent2vertex(tri, -2)
    s = 0.0
    for e in edges
        u, v = edge_indices(e)
        p, q = get_point(tri, u, v)
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = (px + qx) / 2, (py + qy) / 2
        s += segment_function(mx, my, -2)
    end
    return s
end
s = compute_sum_2(tri)

get_neighbours(tri, -2)

get_boundary_map(tri)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

