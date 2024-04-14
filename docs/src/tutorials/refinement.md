```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/refinement.jl"
```

# Mesh Refinement

In this tutorial, we show how we can use mesh refinement
to improve the quality of a mesh by inserting more points.
In this package, we allow for constant area and angle constraints,
and also for custom constraints based on a user-provided function.
You can also limit the maximum number of points.
Moreover, any type of triangulation can be provided, regardless of whether
the triangulation is unconstrained or triangulations, and regardless of the
number of holes and domains. Curve-bounded domains, not included in this tutorial,
can also be refined as discussed in [this tutorial](../tutorials/curve_bounded.md).

Let us start by loading in the packages we will need.

````@example refinement
using DelaunayTriangulation
using CairoMakie
using StableRNGs
````

## Unconstrained triangulation
Let us start with a simple example, refining an unconstrained triangulation.
We will constrain the triangulation such that the minimum angle is
30 degrees, and the maximum area of a triangulation is 1% of the triangulation's
total area. Note that below we need to make sure `points` is mutable, else
it is not possible to push points into the triangulation. Here we use a tuple, but you
could also use e.g. an `ElasticMatrix` from [ElasticArrays.jl](https://github.com/JuliaArrays/ElasticArrays.jl).

````@example refinement
rng = StableRNG(123)
x = rand(rng, 50)
y = rand(rng, 50)
points = tuple.(x, y)
tri = triangulate(points; rng)
orig_tri = deepcopy(tri)
A = get_area(tri)
refine!(tri; min_angle=30.0, max_area=0.01A, rng)
````

The [`refine!`](@ref) function operates on `tri` in-place. If we wanted to review the
statistics of the refined mesh, we can use [`statistics`](@ref):

````@example refinement
statistics(tri)
````

As we can see, the maximum area of a triangle
is about 0.0064, which is indeed less than 1% of the triangulation's area,
which is about 0.0067. Moreover, the smallest angle is indeed greater than 30.

Let us compare the triangulation pre- and post-refinement.

````@example refinement
fig, ax, sc = triplot(orig_tri, axis=(title="Pre-refinement",))
ax = Axis(fig[1, 2], title="Post-refinement")
triplot!(ax, tri)
fig
````

The triangulation is now much finer. There are still some parts with
many more triangles than other regions, but this is most nearly a boundary
or where was a cluster of random points. If we wanted, we could refine again
to try and improve this.

````@example refinement
refine!(tri; min_angle=30.0, max_area=0.001A, rng) # 0.1% instead of 1%
fig, ax, sc = triplot(tri)
fig
````

The quality has now been improved. We could also try improving the minimum
angle further, but even 30 is a bit closer to the limit of convergence (which is
about 33.9 degrees). For example, if we try a minimum angle of 35 degrees,
the algorithm just doesn't even converge, instead it reaches the maximum
number of points.

````@example refinement
test_tri = deepcopy(tri)
refine!(test_tri; min_angle=35.0, max_area=0.001A, max_points = 5_000, rng) # 20_000 so that it doesn't just keep going
statistics(test_tri)
````

As we can see, the smallest angle is about 29 degrees instead of
35 degrees, and there are now 5000 points in the triangulation. The
resulting triangulation is given below:

````@example refinement
fig, ax, sc = triplot(test_tri)
fig
````

This is certainly not a suitable triangulation.

One useful figure to look at for these plots are histograms
that look at the areas and angles. Looking to `tri`, we can plot
these as follows:

````@example refinement
stats = statistics(tri)
fig = Figure(fontsize=33)
areas = get_all_stat(stats, :area) ./ A
angles = first.(get_all_stat(stats, :angles)) # the first is the smallest
ax = Axis(fig[1, 1], xlabel="A/A(Ω)", ylabel="Count", title="Area histogram", width=400, height=400, titlealign=:left)
hist!(ax, areas, bins=0:0.0001:0.0005)
ax = Axis(fig[1, 2], xlabel="θₘᵢₙ", ylabel="Count", title="Angle histogram", width=400, height=400, titlealign=:left)
hist!(ax, rad2deg.(angles), bins=20:2:60)
vlines!(ax, [30.0], color=:red)
resize_to_layout!(fig)
fig
````

We see that indeed many of the triangle areas are very small, and the angles
are all greater than 30 degrees.

## Constrained triangulation and custom constraints
We now give an example of a constrained triangulation being refined. This is the most common
case where the mesh refinement is needed. For this example, we consider an example
with holes, but note that any triangulation can be refined, regardless of the type.
Here is the triangulation we consider.

````@example refinement
n = 100
θ = LinRange(0, 2π, n + 1)
θ = [θ[1:n]; 0]
rev_θ = reverse(θ) # to go from ccw to cw
r₁ = 10.0
r₂ = 5.0
r₃ = 2.5
outer_x, outer_y = r₁ * cos.(θ), r₁ * sin.(θ)
inner_x, inner_y = r₂ * cos.(rev_θ), r₂ * sin.(rev_θ)
innermost_x, innermost_y = r₃ * cos.(θ), r₃ * sin.(θ)
x = [[outer_x], [inner_x], [innermost_x]]
y = [[outer_y], [inner_y], [innermost_y]]
boundary_nodes, points = convert_boundary_points_to_indices(x, y)
rng = StableRNG(456)
tri = triangulate(points; boundary_nodes, rng)
fig, ax, sc = triplot(tri)
fig
````

Let us now refine this triangulation.

````@example refinement
A = get_area(tri)
refine!(tri; min_angle=27.3, max_area=0.01A, rng)
fig, ax, sc = triplot(tri)
fig
````

We inspect the plot, and we might think that it's perhaps not fine enough.
Let's use finer constraints and see what happens. Since
`refine!` operates on `tri` in-place, refining it again
with the constraints below is going to take roughly
the same amount of time as if we had refined it with these
constraints in the first place.

````@example refinement
refine!(tri; min_angle=33.9, max_area=0.001A, rng)
fig, ax, sc = triplot(tri)
fig
````

This is indeed much better, but notice that the inner hole
is much more fine than the outer. This is because we are applying the same
area constraint inside and outside, when really we should try and take note
of the total contribution to the area that each part of the domain gives.
To refine with this in mind, we need to use custom constraints. The
function that we use for constraining the area takes the form
`f(tri, T)`, where `T` is the triangle's vertices and `tri` is the
triangulation. It should return `true` if the triangle
should be refined, and `false` otherwise. Let us define a function such that,
instead of applying constraints so that the triangles are limited to 1% of the total
triangulation area, we do 0.5% or 0.1% of the area of the inner or outer
domain, respectively.

````@example refinement
outer_area = π * (r₁^2 - r₂^2)
inner_area = π * r₃^2
function in_inner(p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    cx, cy = (px + qx + rx) / 3, (py + qy + ry) / 3
    rad2 = cx^2 + cy^2
    return rad2 ≤ r₃^2
end
function area_constraint(_tri, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(_tri, i, j, k)
    A = DelaunayTriangulation.triangle_area(p, q, r)
    return in_inner(p, q, r) ? (A ≥ 0.005inner_area) : (A ≥ 0.001outer_area)
end
````

Let's now refine. We recompute the triangulation so that we can see
the new results.

````@example refinement
boundary_nodes, points = convert_boundary_points_to_indices(x, y)
rng = StableRNG(456)
tri = triangulate(points; boundary_nodes, rng)
refine!(tri; min_angle=30.0, custom_constraint=area_constraint, rng)
fig, ax, sc = triplot(tri)
fig
````

This is now much better, and the two parts of the domain are
appropriately refined. Let us extend our custom constraint function to also
require that any triangle has minimum angle less than 33 degrees
inside the innermost domain, and less than 20 degrees outside the innermost domain.

````@example refinement
function angle_constraint(_tri, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(_tri, i, j, k)
    θ = rad2deg(minimum(DelaunayTriangulation.triangle_angles(p, q, r)))
    return in_inner(p, q, r) ? (θ ≤ 33.9) : (θ ≤ 20.0)
end
function custom_constraint(_tri, T)
    return area_constraint(_tri, T) || angle_constraint(_tri, T)
end
boundary_nodes, points = convert_boundary_points_to_indices(x, y)
rng = StableRNG(456)
tri = triangulate(points; boundary_nodes, rng)
refine!(tri; custom_constraint, rng)
fig, ax, sc = triplot(tri)
fig
````

Indeed, the inner domain is much finer. These examples could be extended
to more complicated cases, for example using adaptive mesh refinement for a numerical
PDE solution so that triangles are refined based on some *a posteriori* error estimate, implemented
using a custom area constraint like above, or even some refinement based on the triangle's
location in space in case of some geospatial application.

## Domains with small angles
In the examples considered, none of the boundaries had small angles. When domains have
small angles, it is not always possible to satisfy the minimum angle constraints, but the
algorithm will still try its best to refine in these locations. Let's consider a complicated
example with many small angles. We consider the boundary of Switzerland, as obtained in this
[NaturalNeighbours.jl example](https://danielvandh.github.io/NaturalNeighbours.jl/stable/swiss/).

````@example refinement
using Downloads
using DelimitedFiles
boundary_url = "https://gist.githubusercontent.com/DanielVandH/13687b0918e45a416a5c93cd52c91449/raw/a8da6cdc94859fd66bcff85a2307f0f9cd57a18c/boundary.txt"
boundary_dir = Downloads.download(boundary_url)
boundary = readdlm(boundary_dir, skipstart=6)
boundary_points = [(boundary[i, 1], boundary[i, 2]) for i in axes(boundary, 1)]
reverse!(boundary_points)
````

Here is the boundary.

````@example refinement
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
rng = StableRNG(789)
tri = triangulate(points; boundary_nodes, rng)
fig, ax, sc = triplot(tri)
fig
````

Now let's refine.

````@example refinement
A = get_area(tri)
refine!(tri; min_angle=30.0, max_area=0.001A, rng)
````

````@example refinement
fig, ax, sc = triplot(tri)
fig
````

We see that the triangulation is now adequately refined. There are
still triangles near the boundaries whose minimum angle is less
than 30 degrees, though, because of the angles that boundary edges
meet at in some places. Most of the triangles will satisfy the
constraint, though, as we show below.

````@example refinement
stats = statistics(tri)
angles = first.(get_all_stat(stats, :angles)) # the first is the smallest
fig, ax, sc = scatter(rad2deg.(angles))
hlines!(ax, [30.0], color = :red, linewidth = 4)
fig
````

As we can see, the vast majority of the triangles satisfy the constraint,
but there are still some that do not. Here is another set of results with a lower minimum angle constraint.

````@example refinement
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
rng = StableRNG(789)
tri = triangulate(points; boundary_nodes, rng)
refine!(tri; min_angle=18.73, max_area=0.001A, rng)
fig = Figure(fontsize = 43)
ax = Axis(fig[1, 1], width = 600, height = 400)
triplot!(tri)
ax = Axis(fig[1, 2], width = 600, height = 400)
stats = statistics(tri)
angles = first.(get_all_stat(stats, :angles)) # the first is the smallest
scatter!(ax, rad2deg.(angles))
hlines!(ax, [18.73], color = :red, linewidth = 4)
resize_to_layout!(fig)
fig
````

In this case, all the triangles satisfy the constraint, of course
at the expense of some other triangles having lesser quality.
## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/refinement.jl).

```julia
using DelaunayTriangulation
using CairoMakie
using StableRNGs

rng = StableRNG(123)
x = rand(rng, 50)
y = rand(rng, 50)
points = tuple.(x, y)
tri = triangulate(points; rng)
orig_tri = deepcopy(tri)
A = get_area(tri)
refine!(tri; min_angle=30.0, max_area=0.01A, rng)

statistics(tri)

fig, ax, sc = triplot(orig_tri, axis=(title="Pre-refinement",))
ax = Axis(fig[1, 2], title="Post-refinement")
triplot!(ax, tri)
fig

refine!(tri; min_angle=30.0, max_area=0.001A, rng) # 0.1% instead of 1%
fig, ax, sc = triplot(tri)
fig

test_tri = deepcopy(tri)
refine!(test_tri; min_angle=35.0, max_area=0.001A, max_points = 5_000, rng) # 20_000 so that it doesn't just keep going
statistics(test_tri)

fig, ax, sc = triplot(test_tri)
fig

stats = statistics(tri)
fig = Figure(fontsize=33)
areas = get_all_stat(stats, :area) ./ A
angles = first.(get_all_stat(stats, :angles)) # the first is the smallest
ax = Axis(fig[1, 1], xlabel="A/A(Ω)", ylabel="Count", title="Area histogram", width=400, height=400, titlealign=:left)
hist!(ax, areas, bins=0:0.0001:0.0005)
ax = Axis(fig[1, 2], xlabel="θₘᵢₙ", ylabel="Count", title="Angle histogram", width=400, height=400, titlealign=:left)
hist!(ax, rad2deg.(angles), bins=20:2:60)
vlines!(ax, [30.0], color=:red)
resize_to_layout!(fig)
fig

n = 100
θ = LinRange(0, 2π, n + 1)
θ = [θ[1:n]; 0]
rev_θ = reverse(θ) # to go from ccw to cw
r₁ = 10.0
r₂ = 5.0
r₃ = 2.5
outer_x, outer_y = r₁ * cos.(θ), r₁ * sin.(θ)
inner_x, inner_y = r₂ * cos.(rev_θ), r₂ * sin.(rev_θ)
innermost_x, innermost_y = r₃ * cos.(θ), r₃ * sin.(θ)
x = [[outer_x], [inner_x], [innermost_x]]
y = [[outer_y], [inner_y], [innermost_y]]
boundary_nodes, points = convert_boundary_points_to_indices(x, y)
rng = StableRNG(456)
tri = triangulate(points; boundary_nodes, rng)
fig, ax, sc = triplot(tri)
fig

A = get_area(tri)
refine!(tri; min_angle=27.3, max_area=0.01A, rng)
fig, ax, sc = triplot(tri)
fig

refine!(tri; min_angle=33.9, max_area=0.001A, rng)
fig, ax, sc = triplot(tri)
fig

outer_area = π * (r₁^2 - r₂^2)
inner_area = π * r₃^2
function in_inner(p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    cx, cy = (px + qx + rx) / 3, (py + qy + ry) / 3
    rad2 = cx^2 + cy^2
    return rad2 ≤ r₃^2
end
function area_constraint(_tri, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(_tri, i, j, k)
    A = DelaunayTriangulation.triangle_area(p, q, r)
    return in_inner(p, q, r) ? (A ≥ 0.005inner_area) : (A ≥ 0.001outer_area)
end

boundary_nodes, points = convert_boundary_points_to_indices(x, y)
rng = StableRNG(456)
tri = triangulate(points; boundary_nodes, rng)
refine!(tri; min_angle=30.0, custom_constraint=area_constraint, rng)
fig, ax, sc = triplot(tri)
fig

function angle_constraint(_tri, T)
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(_tri, i, j, k)
    θ = rad2deg(minimum(DelaunayTriangulation.triangle_angles(p, q, r)))
    return in_inner(p, q, r) ? (θ ≤ 33.9) : (θ ≤ 20.0)
end
function custom_constraint(_tri, T)
    return area_constraint(_tri, T) || angle_constraint(_tri, T)
end
boundary_nodes, points = convert_boundary_points_to_indices(x, y)
rng = StableRNG(456)
tri = triangulate(points; boundary_nodes, rng)
refine!(tri; custom_constraint, rng)
fig, ax, sc = triplot(tri)
fig

using Downloads
using DelimitedFiles
boundary_url = "https://gist.githubusercontent.com/DanielVandH/13687b0918e45a416a5c93cd52c91449/raw/a8da6cdc94859fd66bcff85a2307f0f9cd57a18c/boundary.txt"
boundary_dir = Downloads.download(boundary_url)
boundary = readdlm(boundary_dir, skipstart=6)
boundary_points = [(boundary[i, 1], boundary[i, 2]) for i in axes(boundary, 1)]
reverse!(boundary_points)

boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
rng = StableRNG(789)
tri = triangulate(points; boundary_nodes, rng)
fig, ax, sc = triplot(tri)
fig

A = get_area(tri)
refine!(tri; min_angle=30.0, max_area=0.001A, rng)

fig, ax, sc = triplot(tri)
fig

stats = statistics(tri)
angles = first.(get_all_stat(stats, :angles)) # the first is the smallest
fig, ax, sc = scatter(rad2deg.(angles))
hlines!(ax, [30.0], color = :red, linewidth = 4)
fig

boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
rng = StableRNG(789)
tri = triangulate(points; boundary_nodes, rng)
refine!(tri; min_angle=18.73, max_area=0.001A, rng)
fig = Figure(fontsize = 43)
ax = Axis(fig[1, 1], width = 600, height = 400)
triplot!(tri)
ax = Axis(fig[1, 2], width = 600, height = 400)
stats = statistics(tri)
angles = first.(get_all_stat(stats, :angles)) # the first is the smallest
scatter!(ax, rad2deg.(angles))
hlines!(ax, [18.73], color = :red, linewidth = 4)
resize_to_layout!(fig)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

