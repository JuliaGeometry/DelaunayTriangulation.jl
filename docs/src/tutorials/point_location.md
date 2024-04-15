```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/point_location.jl"
```

# Point Location

In this tutorial, we demonstrate how triangulations can be
used to perform point location. The problem of interest is: Given
a point `p` and a triangulation `tri`, what triangle `T` in `tri`
contains `p`? We provide a function [`jump_and_march`](@ref) for this task,
implementing the algorithm of [Mücke, Saias, and Zhu (1999)](https://doi.org/10.1016/S0925-7721(98)00035-2).
The algorithm has been slightly modified to allow for regions with holes.
Support is also provided for non-convex and disjoint domains, but the algorithm
is significantly slower in these cases and requires some special case. (The approach for these cases is, basically,
to just keep trying new points to start the algorithm from until it works, but you the user must specify
a keyword argument `concavity_protection` to make an extra check to guarantee even greater safety.)

## Unconstrained example
We start with a simple example, demonstrating point location
on an unconstrained triangulation.

````@example point_location
using DelaunayTriangulation
using CairoMakie
using StableRNGs

points = [
    (-3.0, 6.0), (5.0, 1.0), (-5.0, 3.0), (2.0, -3.0),
    (5.0, 8.0), (0.0, 0.0), (2.0, 5.0), (-3.0, 1.0),
    (-2.0, -1.0), (-1.0, 4.0)
]
tri = triangulate(points)
q = (3.0, 3.0)
fig, ax, sc = triplot(tri)
scatter!(ax, q)
fig
````

The aim is to, from `tri`, find which triangle contains the point `q` shown.
Using the `jump_and_march` function, this is simple.

````@example point_location
V = jump_and_march(tri, q)
````

The result means that the triangle `(2, 7, 6)` contains the point, as we can easily check:

````@example point_location
DelaunayTriangulation.point_position_relative_to_triangle(tri, V, q)
````

When we provide no keyword arguments, the default behaviour of `jump_and_march` is to first
sample some number of points (defaults to $\lceil \sqrt[3]{n}\rceil$, where $n$ is the number of points),
and then start at the point that is closest to `q` out of those sampled, then marching along the triangulation
until `q` is found. This number of samples can be changed using the `m` keyword argument. For example,

````@example point_location
V = jump_and_march(tri, q, m=10)
````

means that we get a sample of size 10, and start at whichever point is the closest.
(For technical reasons, this sample is with replacement, so it is possible that the same point is sampled more than once.)
You could also instead specify the point to start at using the `k` keyword argument, in which case no points are sampled.
For example,

````@example point_location
V = jump_and_march(tri, q, k=6)
````

starts the algorithm at the point `6`.

Note also that the triangles found from `jump_and_march` do not have to be given in the same order as they appear
in the triangulation. For example, if a triangle `(i, j, k)` contains the point `q`, then any of `(i, j, k)`, `(j, k, i)`,
or `(k, i, j)` could be returned.

The point `q` does not have to be in the triangulation. For example, consider the following point.

````@example point_location
q = (-5.0, 8.0)
fig, ax, sc = triplot(tri)
scatter!(ax, q)
fig
````

We obtain:

````@example point_location
V = jump_and_march(tri, q)
````

See that the result is a ghost triangle `(1, 5, -1)`. As discussed in the [manual](../manual/ghost_triangles.md),
this can be interpreted as meaning that `q` is between the two lines through the points `1` and `5` that
start at a central point of the triangulation. (The index `-1` is just the ghost vertex.) This can
be visualised.

````@example point_location
fig, ax, sc = triplot(tri, show_ghost_edges=true)
scatter!(ax, q)
fig
````

## Region with concave boundaries and holes
Now we give an example of point location for a reason with holes. Since the
case where all boundaries are convex is reasonably straight forward, here we consider
concave boundaries and discuss methods for improving the speed of the algorithm in this case.
First, let us give our example triangulation.

````@example point_location
a, b, c = (0.0, 8.0), (0.0, 6.0), (0.0, 4.0)
d, e, f = (0.0, 2.0), (0.0, 0.0), (2.0, 0.0)
g, h, i = (4.0, 0.0), (6.0, 0.0), (8.0, 0.0)
j, k, ℓ = (8.0, 1.0), (7.0, 2.0), (5.0, 2.0)
m, n, o = (3.0, 2.0), (2.0, 3.0), (2.0, 5.0)
p, q, r = (2.0, 7.0), (1.0, 8.0), (1.0, 2.2)
s, t, u = (0.4, 1.4), (1.2, 1.8), (2.8, 0.6)
v, w, z = (3.4, 1.2), (1.6, 1.4), (1.6, 2.2)
outer = [[a, b, c, d, e], [e, f, g, h, i, j, k, ℓ], [ℓ, m, n, o, p, q, a]]
inner = [[r, z, v, u, w, t, s, r]]
boundary_nodes, points = convert_boundary_points_to_indices([outer, inner])
rng = StableRNG(125123)
tri = triangulate(points; rng, boundary_nodes)
refine!(tri; max_area=0.01get_area(tri), rng);
nothing #hide
````

The issue with concavity is that the ghost triangles can no longer be sensibly defined.
To demonstrate this, see the following plot:

````@example point_location
fig, ax, sc = triplot(tri, show_ghost_edges=true)
fig
````

The ghost edges now intersect the boundary, which doesn't make sense, and creates difficulties.
Let us now demonstrate how the function still works here. We try finding the blue points shown below.

````@example point_location
qs = [
    (4.0, 5.0), (1.0, 5.6), (0.2, 5.0),
    (0.0, -1.0), (0.5, 3.5), (2.5, 1.5),
    (1.0, 2.0), (4.5, 1.0), (6.0, 1.5),
    (0.5, 8.5), (1.0, 7.5), (1.2, 1.6)
]
fig, ax, sc = triplot(tri, show_ghost_edges=false)
scatter!(ax, qs, color=:blue, markersize=16)
fig
````

Now let's find the triangles.

````@example point_location
Vs = [jump_and_march(tri, q; rng) for q in qs]
````

While we do find some triangles, they may not all be correct. For example,
the triangle found for `(1.2, 1.6)` is

````@example point_location
Vs[end]
````

but the point `(1.2, 1.6)` is actually inside the triangulation. We can even see
this if we run `jump_and_march` again:

````@example point_location
V = jump_and_march(tri, (1.2, 1.6); rng)
````

To protect against this, you need to use `concavity_protection=true`, which
will enable a check to be made that the point is actually outside the triangulation whenever
a ghost triangle is to be returned. If the check finds this to not be the case, it
restarts. With these results, we now compute:

````@example point_location
Vs = [jump_and_march(tri, q; rng, concavity_protection=true) for q in qs]
````

Here is how we can actually test that these results are now correct. We cannot directly
use [`DelaunayTriangulation.point_position_relative_to_triangle`](@ref) because it does not
know that the ghost triangles are invalid. Instead, we find the distance of each point to the
triangulation's boundary using [`DelaunayTriangulation.dist`](@ref) so that we can classify it as being inside or outside of the triangulation,
and then check the type of the found triangle.

````@example point_location
δs = [DelaunayTriangulation.dist(tri, q) for q in qs]
results = Vector{Bool}(undef, length(qs))
for (j, (q, δ, V)) in (enumerate ∘ zip)(qs, δs, Vs)
    cert = DelaunayTriangulation.point_position_relative_to_triangle(tri, V, q)
    is_ghost = DelaunayTriangulation.is_ghost_triangle(V)
    is_outside = DelaunayTriangulation.is_outside(cert)
    if δ ≥ 0.0
        results[j] = !is_outside && !is_ghost
    else # δ < 0.0 ⟹ outside
        results[j] = !is_outside && is_ghost
    end
end
results
````

As we see, the triangles are now all correct.

## Disjoint domains
Now we continue the previous example by adding in another set of
domains that are disjoint to the current domain, thus allowing us to
demonstrate how `jump_and_march` applies here. The new domain is below,
along with the points we will be searching for.

````@example point_location
m₁, n₁, o₁ = (6.0, 8.0), (8.0, 8.0), (8.0, 4.0)
p₁, q₁, r₁ = (10.0, 4.0), (6.0, 6.0), (8.0, 6.0)
s₁, t₁, u₁ = (9.0, 7.0), (4.0, 4.0), (5.0, 4.0)
v₁, w₁ = (5.0, 3.0), (4.0, 3.0)
new_domain₁ = [[m₁, q₁, o₁, p₁, r₁, s₁, n₁, m₁]]
new_domain₂ = [[t₁, w₁, v₁, u₁, t₁]]
boundary_nodes, points = convert_boundary_points_to_indices(
    [outer, inner, new_domain₁, new_domain₂]
)
rng = StableRNG(125123)
tri = triangulate(points; rng, boundary_nodes)
refine!(tri; max_area=0.001get_area(tri), rng)
qs = [
    (0.6, 6.4), (1.4, 0.8), (3.1, 2.9),
    (6.3, 4.9), (4.6, 3.5), (7.0, 7.0),
    (8.9, 5.1), (5.8, 0.8), (1.0, 1.5),
    (1.5, 2.0), (8.15, 6.0)
]
fig, ax, sc = triplot(tri)
scatter!(ax, qs, color=:blue, markersize=16)
fig
````

Here are the `jump_and_march` results.

````@example point_location
Vs = [jump_and_march(tri, q; rng, concavity_protection=true) for q in qs]
````

Again, we can verify that these are all correct as follows. Without `concavity_protection=true`,
these would not be all correct.

````@example point_location
δs = [DelaunayTriangulation.dist(tri, q) for q in qs]
results = Vector{Bool}(undef, length(qs))
for (j, (q, δ, V)) in (enumerate ∘ zip)(qs, δs, Vs)
    cert = DelaunayTriangulation.point_position_relative_to_triangle(tri, V, q)
    is_ghost = DelaunayTriangulation.is_ghost_triangle(V)
    is_outside = DelaunayTriangulation.is_outside(cert)
    if δ ≥ 0.0
        results[j] = !is_outside && !is_ghost
    else # δ < 0.0 ⟹ outside
        results[j] = !is_outside && is_ghost
    end
end
results
````

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/point_location.jl).

```julia
using DelaunayTriangulation
using CairoMakie
using StableRNGs

points = [
    (-3.0, 6.0), (5.0, 1.0), (-5.0, 3.0), (2.0, -3.0),
    (5.0, 8.0), (0.0, 0.0), (2.0, 5.0), (-3.0, 1.0),
    (-2.0, -1.0), (-1.0, 4.0)
]
tri = triangulate(points)
q = (3.0, 3.0)
fig, ax, sc = triplot(tri)
scatter!(ax, q)
fig

V = jump_and_march(tri, q)

DelaunayTriangulation.point_position_relative_to_triangle(tri, V, q)

V = jump_and_march(tri, q, m=10)

V = jump_and_march(tri, q, k=6)

q = (-5.0, 8.0)
fig, ax, sc = triplot(tri)
scatter!(ax, q)
fig

V = jump_and_march(tri, q)

fig, ax, sc = triplot(tri, show_ghost_edges=true)
scatter!(ax, q)
fig

a, b, c = (0.0, 8.0), (0.0, 6.0), (0.0, 4.0)
d, e, f = (0.0, 2.0), (0.0, 0.0), (2.0, 0.0)
g, h, i = (4.0, 0.0), (6.0, 0.0), (8.0, 0.0)
j, k, ℓ = (8.0, 1.0), (7.0, 2.0), (5.0, 2.0)
m, n, o = (3.0, 2.0), (2.0, 3.0), (2.0, 5.0)
p, q, r = (2.0, 7.0), (1.0, 8.0), (1.0, 2.2)
s, t, u = (0.4, 1.4), (1.2, 1.8), (2.8, 0.6)
v, w, z = (3.4, 1.2), (1.6, 1.4), (1.6, 2.2)
outer = [[a, b, c, d, e], [e, f, g, h, i, j, k, ℓ], [ℓ, m, n, o, p, q, a]]
inner = [[r, z, v, u, w, t, s, r]]
boundary_nodes, points = convert_boundary_points_to_indices([outer, inner])
rng = StableRNG(125123)
tri = triangulate(points; rng, boundary_nodes)
refine!(tri; max_area=0.01get_area(tri), rng);

fig, ax, sc = triplot(tri, show_ghost_edges=true)
fig

qs = [
    (4.0, 5.0), (1.0, 5.6), (0.2, 5.0),
    (0.0, -1.0), (0.5, 3.5), (2.5, 1.5),
    (1.0, 2.0), (4.5, 1.0), (6.0, 1.5),
    (0.5, 8.5), (1.0, 7.5), (1.2, 1.6)
]
fig, ax, sc = triplot(tri, show_ghost_edges=false)
scatter!(ax, qs, color=:blue, markersize=16)
fig

Vs = [jump_and_march(tri, q; rng) for q in qs]

Vs[end]

V = jump_and_march(tri, (1.2, 1.6); rng)

Vs = [jump_and_march(tri, q; rng, concavity_protection=true) for q in qs]

δs = [DelaunayTriangulation.dist(tri, q) for q in qs]
results = Vector{Bool}(undef, length(qs))
for (j, (q, δ, V)) in (enumerate ∘ zip)(qs, δs, Vs)
    cert = DelaunayTriangulation.point_position_relative_to_triangle(tri, V, q)
    is_ghost = DelaunayTriangulation.is_ghost_triangle(V)
    is_outside = DelaunayTriangulation.is_outside(cert)
    if δ ≥ 0.0
        results[j] = !is_outside && !is_ghost
    else # δ < 0.0 ⟹ outside
        results[j] = !is_outside && is_ghost
    end
end
results

m₁, n₁, o₁ = (6.0, 8.0), (8.0, 8.0), (8.0, 4.0)
p₁, q₁, r₁ = (10.0, 4.0), (6.0, 6.0), (8.0, 6.0)
s₁, t₁, u₁ = (9.0, 7.0), (4.0, 4.0), (5.0, 4.0)
v₁, w₁ = (5.0, 3.0), (4.0, 3.0)
new_domain₁ = [[m₁, q₁, o₁, p₁, r₁, s₁, n₁, m₁]]
new_domain₂ = [[t₁, w₁, v₁, u₁, t₁]]
boundary_nodes, points = convert_boundary_points_to_indices(
    [outer, inner, new_domain₁, new_domain₂]
)
rng = StableRNG(125123)
tri = triangulate(points; rng, boundary_nodes)
refine!(tri; max_area=0.001get_area(tri), rng)
qs = [
    (0.6, 6.4), (1.4, 0.8), (3.1, 2.9),
    (6.3, 4.9), (4.6, 3.5), (7.0, 7.0),
    (8.9, 5.1), (5.8, 0.8), (1.0, 1.5),
    (1.5, 2.0), (8.15, 6.0)
]
fig, ax, sc = triplot(tri)
scatter!(ax, qs, color=:blue, markersize=16)
fig

Vs = [jump_and_march(tri, q; rng, concavity_protection=true) for q in qs]

δs = [DelaunayTriangulation.dist(tri, q) for q in qs]
results = Vector{Bool}(undef, length(qs))
for (j, (q, δ, V)) in (enumerate ∘ zip)(qs, δs, Vs)
    cert = DelaunayTriangulation.point_position_relative_to_triangle(tri, V, q)
    is_ghost = DelaunayTriangulation.is_ghost_triangle(V)
    is_outside = DelaunayTriangulation.is_outside(cert)
    if δ ≥ 0.0
        results[j] = !is_outside && !is_ghost
    else # δ < 0.0 ⟹ outside
        results[j] = !is_outside && is_ghost
    end
end
results
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

