```@meta
EditURL = "https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/operations_flip_edge.jl"
```

# Triangulation Operations
## Edge Flipping

This tutorial shows we can flip edges in a triangulation. Edge flipping, discussed in more detail
in the [edge flipping application](../applications/edge_flipping.md), is the flipping of some
edge `(i, j)` to the edge `(k, ℓ)`, where `(i, j)` and `(k, ℓ)` are diagonals of the quadrilateral
formed by $p_ip_jp_kp_l$. The edge flip is illustrated below.

````@example operations_flip_edge
using CairoMakie #hide
points = [(0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)] #hide
T1 = points[[1, 2, 3]] #hide
T2 = points[[1, 3, 4]] #hide
T3 = points[[1, 2, 4]] #hide
T4 = points[[2, 3, 4]] #hide
fig = Figure() #hide
ax1 = Axis(fig[1, 1], width = 600, height = 400) #hide
ax2 = Axis(fig[1, 2], width = 600, height = 400) #hide
poly!(ax1, [T1; T2], color=(:white, 0.0), strokewidth=3) #hide
poly!(ax2, [T3; T4], color=(:white, 0.0), strokewidth=3) #hide
for ax in (ax1, ax2) #hide
    hidedecorations!(ax) #hide
    hidespines!(ax) #hide
    text!(ax, [(0.05, -1.1)]; text=[L"p_j"], fontsize=43) #hide
    text!(ax, [(0.9, 0.1)]; text=[L"p_k"], fontsize=43) #hide
    text!(ax, [(0.05, 1.0)]; text=[L"p_i"], fontsize=43) #hide
    text!(ax, [(-1.05, 0.05)]; text=[L"p_\ell"], fontsize=43) #hide
    xlims!(ax, -1.1, 1.1) #hide
    ylims!(ax, -1.3, 1.3) #hide
end #hide
resize_to_layout!(fig) #hide
fig #hide
````

Note that this edge flip only makes sense if the quadrilateral
is convex. If the quadrilateral is not convex, then the edge flip
will not be valid; no checks are made for whether the
quadrilateral is convex inside the `flip_edge!` function.

Let us now showcase how we can flip edges. First, we load in the
packages we need.

````@example operations_flip_edge
using DelaunayTriangulation
using CairoMakie
````

Let us now define our initial triangulation.

````@example operations_flip_edge
points = [(0.0, 0.0), (0.8, 0.0), (1.3, 1.0), (0.0, 1.0)]
tri = triangulate(points, delete_ghosts = false);
nothing #hide
````

Note that we use `delete_ghosts=false` so that we can safely flip the edges
while taking care of boundary effects.

Now, flipping the edge is simple. We simply provide the indices `i` and `j`
for the edge we want to flip. Let us flip the edge `(2, 4)`.

````@example operations_flip_edge
fig, ax, sc = triplot(tri, axis = (title = "Before flipping",))
ax2 = Axis(fig[1, 2], title = "After flipping")
flip_edge!(tri, 2, 4)
triplot!(ax2, tri)
fig
````

As simple as that. Note that no checks are made for whether the edge is actually in the
triangulation, on the boundary, or if the associated quadrilateral is convex. It is
up to you to check this if needed; one way to check would be to use `DelaunayTriangulation.is_legal`,
as is done inside `legalise_edge!` -- see the [next tutorial](operations_legalise_edge.md).
## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/new-docs/docs/src/literate_tutorials/operations_flip_edge.jl).

```julia
using DelaunayTriangulation
using CairoMakie

points = [(0.0, 0.0), (0.8, 0.0), (1.3, 1.0), (0.0, 1.0)]
tri = triangulate(points, delete_ghosts = false);

fig, ax, sc = triplot(tri, axis = (title = "Before flipping",))
ax2 = Axis(fig[1, 2], title = "After flipping")
flip_edge!(tri, 2, 4)
triplot!(ax2, tri)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

