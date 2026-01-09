```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/constrained_outer_boundary_segmented.jl"
```

# Constrained Triangulations
## Sectioned Outer Boundary

In this tutorial, we again consider a triangulation with a constrained boundary. In
contrast to the previous tutorial, this outer boundary will be represented as a
chain of multiple paths or sections. This changes nothing geometrically, but it allows
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
    (2.0, 4.0), (8.0, 2.0),
]
````

We now want to define our boundary. The method for providing a boundary to be
identified into multiple sections is to provided a vector of vectors of indices,
where each vector of indices is a section. The last index of each section must match
the first index of the next section, including the last with the first
section so that the boundary is closed. Here, we provide three .

````@example constrained_outer_boundary_segmented
section_1 = [(0.0, 0.0), (14.0, 0.0)]
section_2 = [(14.0, 0.0), (10.0, 4.0), (4.0, 6.0), (2.0, 12.0), (0.0, 14.0)]
section_3 = [(0.0, 14.0), (0.0, 0.0)]
boundary_points = [section_1, section_2, section_3]
````

We now convert these boundary points to indices using
[`convert_boundary_points_to_indices`](@ref), and then we triangulate.
We also add a constrained edge.

````@example constrained_outer_boundary_segmented
E = Set(((6, 9),)) # (0, 0) → (4, 6)
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points; existing_points = points)
tri = triangulate(points; boundary_nodes, segments = E)
````

````@example constrained_outer_boundary_segmented
fig, ax, sc = triplot(tri, show_constrained_edges = true, constrained_edge_linewidth = 6)
lines!(ax, section_1, color = :red, linewidth = 6)
lines!(ax, section_2, color = :green, linewidth = 6)
lines!(ax, section_3, color = :blue, linewidth = 6)
fig
````

The first section is in red, the second section is in green, and the third section
is in blue. We use ghost vertices to identify the sections, where the first section is identified
by `-1`, the second by `-2`, and the third by `-3`.

Before we go into how the sections can be worked with, let us make a note regarding [`get_interior_segments(tri)`](@ref get_interior_segments),
now that we have both constrained segments and boundary segments which might technically both be thought of
as being constrained edges. If we look at `get_interior_segments(tri)`, we get:

````@example constrained_outer_boundary_segmented
get_interior_segments(tri)
````

This is just the constrained segment we provided, and not the boundary segment. If we instead
want all constrained segments, considering both the boundary segments and the segments provided,
we can instead use [`get_all_segments(tri)`](@ref get_all_segments).

````@example constrained_outer_boundary_segmented
get_all_segments(tri)
````

Let us now explore the several ways available for working with this boundary. First, if we just want to work
with the boundary edges without caring about the order, we can again use the boundary edge map.

````@example constrained_outer_boundary_segmented
get_boundary_edge_map(tri)
````

Remember that the keys are the edges, and the values are `Tuples` that give us (1) the section index,
and (2) the position of the edge within that section (more specifically, the position of the first vertex of
the edge). This would be useful if, for example, you don't care about the order of the edges, but you do
care about what section the edge belongs to so that you can assign a boundary condition for example.
In the previous tutorial, we mentioned that the `Tuples` are of the form `(I, J)` so that the
corresponding edge is identified from `get_boundary_nodes(get_boundary_nodes(tri, I), J)`, but since both
`I` and `J` are integers in this case (since we have sections), we can just use `get_boundary_nodes(tri, (I, J))`.
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
and $\Gamma_i$ denotes the $i$th section.

````@example constrained_outer_boundary_segmented
function section_function(x, y, section_index)
    f = if abs(section_index) == 1
        1.0
    elseif abs(section_index) == 2
        sin(x - y)
    else
        cos(x + y)
    end
    return f
end
function compute_sum(tri)
    bem = get_boundary_edge_map(tri)
    s = 0.0
    for (e, (section_index, _)) in bem
        u, v = edge_vertices(e)
        p, q = get_point(tri, u, v)
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = (px + qx) / 2, (py + qy) / 2
        s += section_function(mx, my, section_index)
    end
    return s
end
s = compute_sum(tri)
````

An alternative way to look at each section is to use [`get_adjacent2vertex`](@ref) with the associated
ghost vertex.

````@example constrained_outer_boundary_segmented
get_adjacent2vertex(tri, -1)
````

````@example constrained_outer_boundary_segmented
get_adjacent2vertex(tri, -2)
````

````@example constrained_outer_boundary_segmented
get_adjacent2vertex(tri, -3)
````

Note that the provided edges are not in order, but this is helpful for considering specific sections. For
example, if we just wanted to compute the above sum over the second section, we could do

````@example constrained_outer_boundary_segmented
function compute_sum_2(tri)
    edges = get_adjacent2vertex(tri, -2)
    s = 0.0
    for e in edges
        u, v = edge_vertices(e)
        p, q = get_point(tri, u, v)
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = (px + qx) / 2, (py + qy) / 2
        s += section_function(mx, my, -2)
    end
    return s
end
s = compute_sum_2(tri)
````

If your application instead wanted all the nodes on the section rather than the edges,
you can look at the neighbours to the ghost vertex. For example, all the nodes
on the section section can be identified using

````@example constrained_outer_boundary_segmented
get_neighbours(tri, -2)
````

Another field is the `ghost_vertex_map`, which maps a given ghost vertex to the associated section. This is
more so useful for internal methods, but you may sometimes need it.

````@example constrained_outer_boundary_segmented
get_ghost_vertex_map(tri)
````

In this case, the `i`th section just has the ghost vertex `-i`, but this is typically used
to deal with the case of multiple boundaries so that we know where a ghost vertex belongs.

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/constrained_outer_boundary_segmented.jl).

```julia
using DelaunayTriangulation
using CairoMakie

points = [
    (2.0, 8.0), (6.0, 4.0), (2.0, 6.0),
    (2.0, 4.0), (8.0, 2.0),
]

section_1 = [(0.0, 0.0), (14.0, 0.0)]
section_2 = [(14.0, 0.0), (10.0, 4.0), (4.0, 6.0), (2.0, 12.0), (0.0, 14.0)]
section_3 = [(0.0, 14.0), (0.0, 0.0)]
boundary_points = [section_1, section_2, section_3]

E = Set(((6, 9),)) # (0, 0) → (4, 6)
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points; existing_points = points)
tri = triangulate(points; boundary_nodes, segments = E)

fig, ax, sc = triplot(tri, show_constrained_edges = true, constrained_edge_linewidth = 6)
lines!(ax, section_1, color = :red, linewidth = 6)
lines!(ax, section_2, color = :green, linewidth = 6)
lines!(ax, section_3, color = :blue, linewidth = 6)
fig

get_interior_segments(tri)

get_all_segments(tri)

get_boundary_edge_map(tri)

function section_function(x, y, section_index)
    f = if abs(section_index) == 1
        1.0
    elseif abs(section_index) == 2
        sin(x - y)
    else
        cos(x + y)
    end
    return f
end
function compute_sum(tri)
    bem = get_boundary_edge_map(tri)
    s = 0.0
    for (e, (section_index, _)) in bem
        u, v = edge_vertices(e)
        p, q = get_point(tri, u, v)
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = (px + qx) / 2, (py + qy) / 2
        s += section_function(mx, my, section_index)
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
        u, v = edge_vertices(e)
        p, q = get_point(tri, u, v)
        px, py = getxy(p)
        qx, qy = getxy(q)
        mx, my = (px + qx) / 2, (py + qy) / 2
        s += section_function(mx, my, -2)
    end
    return s
end
s = compute_sum_2(tri)

get_neighbours(tri, -2)

get_ghost_vertex_map(tri)
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

