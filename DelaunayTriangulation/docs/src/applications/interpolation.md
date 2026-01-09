```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_applications/interpolation.jl"
```

# Interpolation

For our first application, we consider applying Voronoi tessellations to interpolation.
The discussion here is based on the implementation of what is known as _natural neighbour interpolation_ in
[NaturalNeighbours.jl](https://github.com/DanielVandH/NaturalNeighbours.jl). NaturalNeighbours.jl considers a diverse set of interpolants for this purpose,
but here we will focus only on a single interpolant (called the `Sibson()` interpolant in NaturalNeighbours.jl).

The problem we are considering is as follows: Given some data $\mathcal D = \{(\vb x_i, z_i)\}_{i=1}^m \subseteq \mathbb R^2 \times \mathbb R$, we want to
construct a smooth interpolant $f \colon \mathcal C\mathcal H(X) \to \mathbb R$ of the data, where $X = \{\vb x_1, \ldots, \vb x_m\}$, such
that $f(\vb x_i) = z_i$ for all $i = 1, \ldots, m$.

One idea to interpolate this data would be to use a _piecewise linear interpolant_, obtained by defining, inside each triangle $T \in \mathcal D\mathcal T(X)$, a
piecewise linear interpolant using the data points at each of the three vertices. One problem with this is that $\mathcal D\mathcal T(X)$ does not
depend continuously on $X$, meaning that small perturbations of a data point can lead to topological changes in the triangulation.
One way around this is to instead use the Voronoi tessellation to guide the tessellation, since $\mathcal V(X)$ _does_ depend continuously on $X$.
We give an example of this below, where we show how $\mathcal D\mathcal T(X)$ may change significant after a small perturbation while $\mathcal V(X)$
does not at the same time.

````@example interpolation
using DelaunayTriangulation #hide
using CairoMakie #hide
A, B, C, D, E, F, G, H = (0.3, 1.1), (-0.1, 0.8), (0.2, 0.3), (0.6, 0.2), (0.8, 0.8), (0.3, 0.9), (0.5503600264347, 0.6814266789918), (1.1, 0.5) #hide
G2 = (0.5496217775447, 0.7146478790414) #hide
fig = Figure() #hide
ax = Axis(fig[1, 1], width = 400, height = 400, title = "Original") #hide
xlims!(ax, -0.2, 1.3) #hide
ylims!(ax, 0.1, 1.2) #hide
points = [A, B, C, D, E, F, G, H] #hide
triplot!(ax, points, show_points = true) #hide
scatter!(ax, [G], color = :red, markersize = 14) #hide
hidedecorations!(ax) #hide
ax2 = Axis(fig[2, 1], width = 400, height = 400) #hide
voronoiplot!(ax2, points, clip = (-0.2, 1.3, 0.1, 1.2)) #hide
xlims!(ax2, -0.2, 1.3) #hide
ylims!(ax2, 0.1, 1.2) #hide
scatter!(ax2, [G], color = :red, markersize = 14) #hide
hidedecorations!(ax2) #hide
ax3 = Axis(fig[1, 2], width = 400, height = 400, title = "Perturbed") #hide
points = [A, B, C, D, E, F, G2, H] #hide
triplot!(ax3, points, show_points = true) #hide
hidedecorations!(ax3) #hide
xlims!(ax3, -0.2, 1.3) #hide
ylims!(ax3, 0.1, 1.2) #hide
scatter!(ax3, [G2], color = :red, markersize = 14) #hide
ax4 = Axis(fig[2, 2], width = 400, height = 400) #hide
voronoiplot!(ax4, points, clip = (-0.2, 1.3, 0.1, 1.2)) #hide
scatter!(ax4, [G2], color = :red, markersize = 14) #hide
xlims!(ax4, -0.2, 1.3) #hide
ylims!(ax4, 0.1, 1.2) #hide
hidedecorations!(ax4) #hide
resize_to_layout!(fig) #hide
fig #hide
````

This observation motivates the use of the Voronoi tessellation to guide the interpolation.

## Natural Neighbours

Let us start by defining _natural neighbours_. Consider two Voronoi polygons $\mathcal V_i$ and $\mathcal V_j$,
and let $\mathcal F_{ij} = \mathcal V_i \cap \mathcal V_j$. If $\mathcal F_{ij} \neq \emptyset$, we say that $\vb x_i$ and $\vb x_j$ are natural neighbours in
$X$. For a point $\vb  \in X$, we denote its set of natural neighbours by $N(\vb x) \subseteq X$, and the corresponding
indices by $N_i = \{j : \vb x_j \in N(\vb x_i)\}$.

We use natural neighbours to guide our representation of points. In particular, we use _natural neighbour coordinates_,
which are defined as follows: For a point $\vb x \in X$, a vector $\boldsymbol\lambda$ are called a set of _natural neighbour coordinates_ if
$\lambda_{i} \geq 0$ for each $i$, $\boldsymbol\lambda$ is continuous with respect to $\vb x$, and $\lambda_i > 0 \iff \vb x_i \in N(\vb x)$.

## Sibsonian Interpolation
Here, for $\boldsymbol\lambda$ we will discuss _Sibsonian coordinates_. To define these coordinates, take $\mathcal V(X)$ and consider what happens when a
new point $\vb x_0$ is added into it. This new point creates a new tessellation $\mathcal V(X \cup \{\vb x_0\})$, where all the tiles
that change in the new tessellation compared to $\mathcal V(\vb X)$ are those in $N(\vb x_0)$. We let $A(\vb x_0)$ be the area of the
new tile created by $\vb x_0$, and let $A(\vb x_i)$ be the area of the intersection between the original tile for $\vb x_i$ and the new tile from
$\vb x_0$. The Sibson coordinates are then
```math
\lambda_i(\vb x_0) = \frac{A(\vb x_i)}{A(\vb x_0)}.
```

Using this definition of Sibsonian coordinates, the Sibson interpolant is defined by
```math
f(\vb x_0) = \sum_{i \in N_0} \lambda_i(\vb x_0)z_i.
```
This interpolant is $C^1$ continuous in $\mathcal C\mathcal H(X) \setminus X$, with derivative discontinuities at the data sites.

## Implementation
Let us now implement the Sibson interpolant. We note that, in this implementation, we ignore some edge cases; these are handled properly in the
implementation in NaturalNeighbours.jl. Moreover, our implementation will not be the most efficient, but will be enough
for the purposes of this demonstration. We assume that all points are contained in the interior of $\mathcal C\mathcal H(X)$.

### Bowyer-Watson envelope
The first issue to deal with in our implementation is the computation of the _Bowyer-Watson envelope_. For a point $\vb x$, its
Bowyer-Watson envelope is the boundary of the set of all triangles in $\mathcal D\mathcal T(X)$ whose circumcircles contain $\vb x$.
This envelope tells us the region in which any changes to the triangulation and to the tessellation can occur.
A reasonably straight forward way to implement this is to simply add $\vb x$ into $\mathcal D\mathcal T(X)$ and take all the triangles
containing $\vb x$ as a vertex in $\mathcal D\mathcal T(X \cup \{\vb x\})$. We then remove $\vb x$ from $\mathcal D\mathcal T(X \cup \{\vb x\})$.
We implement this below.[^1]

[^1]: This is more expensive than we need. In NaturalNeighbours.jl, we use the `peek` keyword in `triangulate` to avoid making any changes to the triangulation itself, and use the `InsertionEventHistory` to track all changes made.

````@example interpolation
using DelaunayTriangulation
function compute_envelope(tri::Triangulation, point)
    r = DelaunayTriangulation.num_points(tri)
    x, y = getxy(point)
    add_point!(tri, x, y)
    envelope_vertices = DelaunayTriangulation.get_surrounding_polygon(tri, r + 1)
    push!(envelope_vertices, envelope_vertices[begin])
    envelope_points = [get_point(tri, i) for i in envelope_vertices]
    delete_point!(tri, r + 1)
    DelaunayTriangulation.pop_point!(tri)
    return envelope_vertices, envelope_points
end
````

Let's now check that this function works.

````@example interpolation
using CairoMakie
using StableRNGs
using ElasticArrays
rng = StableRNG(999)
points = ElasticMatrix(randn(rng, 2, 50)) # so that the points are mutable
tri = triangulate(points; rng)
envelope_vertices, envelope_points = compute_envelope(tri, (0.5, 0.5))

fig = Figure(fontsize = 24)
ax = Axis(fig[1, 1], width = 400, height = 400)
triplot!(ax, tri, show_points = true)
ax2 = Axis(fig[1, 2], width = 400, height = 400)
add_point!(tri, 0.5, 0.5)
triplot!(ax2, tri, show_points = true)
poly!(ax2, envelope_points, color = (:red, 0.2))
resize_to_layout!(fig)
fig
````

As we can see, the red region we have computed from our envelope is indeed the envelope we need.

### Computing the Sibsonian coordinates
When we are interpolating at a point $\vb x$, remember that we need to know the area of the
Voronoi polygon that would be produced when $\vb x$ is inserted into $\mathcal V(X)$. To compute this area,
we need to know how we can compute it using the Bowyer-Watson envelope. Remember that the
Voronoi polygons are obtained by drawing lines between the circumcenters of neighbouring triangles.
This can be done using the Bowyer-Watson envelope: The edges of the envelope, together with the new point,
define the triangles that would be produced if it were to be added into the triangulation, and so we can join the circumcenters
of these triangles to compute the new Voronoi polygon. We can then use this polygon, together with the
original tessellation, to compute the Sibsonian coordinates.

Let us first discuss the area of the Voronoi polygons in $N(\vb x)$ ($\vb x$'s natural neighbours, i.e. the vertices of the envelope) before $\vb x$ is inserted. We only need
to compute the part of the area that is contained within the envelope, since everything outside of that
envelope is unchanged. If we included the entire area, then the area that we subtract off for the intersection we compute later
would just cancel it out anyway. Let's zoom in on the envelope and consider a specific example of how we can do this computation.

````@example interpolation
fig = Figure(fontsize = 24)
ax = Axis(fig[1, 1], width = 400, height = 400)
triplot!(ax, tri, show_points = true)
lines!(ax, envelope_points, color = :red)
j = 7 # example vertex
v = envelope_vertices[j]
scatter!(ax, [get_point(tri, v)], color = :blue)
first_neighbour = envelope_vertices[j - 1]
next_triangle = get_adjacent(tri, first_neighbour, v)
next_triangle_2 = get_adjacent(tri, v, envelope_vertices[j + 1])
last_neighbour = envelope_vertices[j + 1]
polygon_points = [
    get_point(tri, v),
    (get_point(tri, v) .+ get_point(tri, last_neighbour)) ./ 2,
    DelaunayTriangulation.triangle_circumcenter(tri, (v, envelope_vertices[j + 1], next_triangle_2)),
    DelaunayTriangulation.triangle_circumcenter(tri, (first_neighbour, v, next_triangle)),
    (get_point(tri, v) .+ get_point(tri, first_neighbour)) ./ 2,
]
poly!(ax, polygon_points, color = (:blue, 0.5), strokecolor = :blue, strokewidth = 2)
xlims!(ax, -0.5, 1.4)
ylims!(ax, -0.15, 1.4)
resize_to_layout!(fig)
fig
````

The relevant polygon is shown above in blue, associated with the generator shown by the blue point. We need to compute the area of this polygon.
This is simple using the [shoelace formula](https://en.wikipedia.org/wiki/Shoelace_formula). Our
implementation of this is given below.

````@example interpolation
function polygon_area(points) # this is the first formula in the "Other formulae" section of the above Wikipedia article
    n = DelaunayTriangulation.num_points(points)
    p, q, r, s = get_point(points, 1, 2, n, n - 1)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    sx, sy = getxy(s)
    area = px * (qy - ry) + rx * (py - sy)
    for i in 2:(n - 1)
        p, q, r = get_point(points, i, i + 1, i - 1)
        px, py = getxy(p)
        qx, qy = getxy(q)
        rx, ry = getxy(r)
        area += px * (qy - ry)
    end
    return area / 2
end
function pre_insertion_area(tri::Triangulation, i, envelope_vertices) # area from the envelope[i]th generator
    poly_points = NTuple{2, Float64}[]
    u = envelope_vertices[i]
    prev_index = i == 1 ? length(envelope_vertices) - 1 : i - 1
    next_index = i == length(envelope_vertices) ? 1 : i + 1
    first_neighbour = envelope_vertices[prev_index]
    last_neighbour = envelope_vertices[next_index]
    v = last_neighbour
    (ux, uy), (vx, vy) = get_point(tri, u, v)
    mx1, my1 = (ux + vx) / 2, (uy + vy) / 2
    push!(poly_points, (mx1, my1))
    while v ≠ first_neighbour
        w = get_adjacent(tri, u, v)
        cx, cy = DelaunayTriangulation.triangle_circumcenter(tri, (u, v, w))
        push!(poly_points, (cx, cy))
        v = w
    end
    vx, vy = get_point(tri, v)
    mx, my = (ux + vx) / 2, (uy + vy) / 2
    push!(poly_points, (mx, my), (mx1, my1))
    return polygon_area(poly_points)
end
````

The details for the post-insertion area are similar, but now the triangles that we take the circumcenters of
are those where the edges instead join with the inserted vertex. The function we use is below.

````@example interpolation
function post_insertion_area(tri::Triangulation, i, envelope_vertices, point)
    u = envelope_vertices[i]
    prev_index = i == 1 ? length(envelope_vertices) - 1 : i - 1
    next_index = i == length(envelope_vertices) ? 1 : i + 1
    first_neighbour = envelope_vertices[prev_index]
    last_neighbour = envelope_vertices[next_index]
    p, q, r = get_point(tri, u, first_neighbour, last_neighbour)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    mpq = (px + qx) / 2, (py + qy) / 2
    mpr = (px + rx) / 2, (py + ry) / 2
    g1 = DelaunayTriangulation.triangle_circumcenter(p, r, point)
    !all(isfinite, g1) && return NaN # point is one of p and r, i.e. we are interpolating at a data site
    g2 = DelaunayTriangulation.triangle_circumcenter(q, p, point)
    !all(isfinite, g2) && return NaN
    points = (mpq, mpr, g1, g2, mpq)
    return polygon_area(points)
end
````

Now that we can compute the pre- and post-insertion areas, we can start computing the Sibsonian coordinates.

````@example interpolation
function compute_sibson_coordinates(tri::Triangulation, envelope_vertices, point)
    coordinates = zeros(length(envelope_vertices) - 1)
    w = 0.0
    for i in firstindex(envelope_vertices):(lastindex(envelope_vertices) - 1)
        pre = max(0.0, pre_insertion_area(tri, i, envelope_vertices))
        post = max(0.0, post_insertion_area(tri, i, envelope_vertices, point))
        if isnan(post) # need to return the the vector λ = [1] since we are exactly at a data site
            return [1.0]
        end
        coordinates[i] = max(pre - post, 0.0) # take care of any precision issues
        w += coordinates[i]
    end
    coordinates ./= w
    return coordinates
end
````

This function gives our $\boldsymbol\lambda$ vector. Notice that, in the computation of these coordinates,
we need needed to have $\mathcal V(X)$ directly or make use of the data $z_i$.

## Evaluating the Sibsonian interpolant
Now we can evaluate our Sibson interpolant. The following function does this for us.

````@example interpolation
function evaluate_sibson_interpolant(tri::Triangulation, z, point)
    envelope_vertices, _ = compute_envelope(tri, point)
    λ = compute_sibson_coordinates(tri, envelope_vertices, point)
    if length(λ) == 1
        for i in each_solid_vertex(tri)
            get_point(tri, i) == point && return z[i]
        end
    else
        itp = 0.0
        for (λ, k) in zip(λ, envelope_vertices)
            itp += λ * z[k]
        end
        return itp
    end
end
````

Let's now use this function to interpolate some data.

````@example interpolation
f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
trit = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 30, 30)
zz = [f(x, y) for (x, y) in DelaunayTriangulation.each_point(trit)]
xx = LinRange(0.001, 0.999, 20) # handling points on the boundary requires more care than we have discussed here
yy = LinRange(0.001, 0.999, 20)
fig = Figure(fontsize = 24)
ax = Axis(fig[1, 1], xlabel = L"x", ylabel = L"y", title = "True function", titlealign = :left, width = 400, height = 400)
contourf!(ax, xx, yy, f.(xx, yy'))
ax2 = Axis(fig[1, 2], xlabel = L"x", ylabel = L"y", title = "Interpolant", titlealign = :left, width = 400, height = 400)
zi = [evaluate_sibson_interpolant(trit, zz, (xᵢ, yᵢ)) for xᵢ in xx, yᵢ in yy]
contourf!(ax2, xx, yy, zi)
resize_to_layout!(fig)
fig
````

Works perfectly!

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_applications/interpolation.jl).

```julia
using DelaunayTriangulation
function compute_envelope(tri::Triangulation, point)
    r = DelaunayTriangulation.num_points(tri)
    x, y = getxy(point)
    add_point!(tri, x, y)
    envelope_vertices = DelaunayTriangulation.get_surrounding_polygon(tri, r + 1)
    push!(envelope_vertices, envelope_vertices[begin])
    envelope_points = [get_point(tri, i) for i in envelope_vertices]
    delete_point!(tri, r + 1)
    DelaunayTriangulation.pop_point!(tri)
    return envelope_vertices, envelope_points
end

using CairoMakie
using StableRNGs
using ElasticArrays
rng = StableRNG(999)
points = ElasticMatrix(randn(rng, 2, 50)) # so that the points are mutable
tri = triangulate(points; rng)
envelope_vertices, envelope_points = compute_envelope(tri, (0.5, 0.5))

fig = Figure(fontsize = 24)
ax = Axis(fig[1, 1], width = 400, height = 400)
triplot!(ax, tri, show_points = true)
ax2 = Axis(fig[1, 2], width = 400, height = 400)
add_point!(tri, 0.5, 0.5)
triplot!(ax2, tri, show_points = true)
poly!(ax2, envelope_points, color = (:red, 0.2))
resize_to_layout!(fig)
fig

fig = Figure(fontsize = 24)
ax = Axis(fig[1, 1], width = 400, height = 400)
triplot!(ax, tri, show_points = true)
lines!(ax, envelope_points, color = :red)
j = 7 # example vertex
v = envelope_vertices[j]
scatter!(ax, [get_point(tri, v)], color = :blue)
first_neighbour = envelope_vertices[j - 1]
next_triangle = get_adjacent(tri, first_neighbour, v)
next_triangle_2 = get_adjacent(tri, v, envelope_vertices[j + 1])
last_neighbour = envelope_vertices[j + 1]
polygon_points = [
    get_point(tri, v),
    (get_point(tri, v) .+ get_point(tri, last_neighbour)) ./ 2,
    DelaunayTriangulation.triangle_circumcenter(tri, (v, envelope_vertices[j + 1], next_triangle_2)),
    DelaunayTriangulation.triangle_circumcenter(tri, (first_neighbour, v, next_triangle)),
    (get_point(tri, v) .+ get_point(tri, first_neighbour)) ./ 2,
]
poly!(ax, polygon_points, color = (:blue, 0.5), strokecolor = :blue, strokewidth = 2)
xlims!(ax, -0.5, 1.4)
ylims!(ax, -0.15, 1.4)
resize_to_layout!(fig)
fig

function polygon_area(points) # this is the first formula in the "Other formulae" section of the above Wikipedia article
    n = DelaunayTriangulation.num_points(points)
    p, q, r, s = get_point(points, 1, 2, n, n - 1)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    sx, sy = getxy(s)
    area = px * (qy - ry) + rx * (py - sy)
    for i in 2:(n - 1)
        p, q, r = get_point(points, i, i + 1, i - 1)
        px, py = getxy(p)
        qx, qy = getxy(q)
        rx, ry = getxy(r)
        area += px * (qy - ry)
    end
    return area / 2
end
function pre_insertion_area(tri::Triangulation, i, envelope_vertices) # area from the envelope[i]th generator
    poly_points = NTuple{2, Float64}[]
    u = envelope_vertices[i]
    prev_index = i == 1 ? length(envelope_vertices) - 1 : i - 1
    next_index = i == length(envelope_vertices) ? 1 : i + 1
    first_neighbour = envelope_vertices[prev_index]
    last_neighbour = envelope_vertices[next_index]
    v = last_neighbour
    (ux, uy), (vx, vy) = get_point(tri, u, v)
    mx1, my1 = (ux + vx) / 2, (uy + vy) / 2
    push!(poly_points, (mx1, my1))
    while v ≠ first_neighbour
        w = get_adjacent(tri, u, v)
        cx, cy = DelaunayTriangulation.triangle_circumcenter(tri, (u, v, w))
        push!(poly_points, (cx, cy))
        v = w
    end
    vx, vy = get_point(tri, v)
    mx, my = (ux + vx) / 2, (uy + vy) / 2
    push!(poly_points, (mx, my), (mx1, my1))
    return polygon_area(poly_points)
end

function post_insertion_area(tri::Triangulation, i, envelope_vertices, point)
    u = envelope_vertices[i]
    prev_index = i == 1 ? length(envelope_vertices) - 1 : i - 1
    next_index = i == length(envelope_vertices) ? 1 : i + 1
    first_neighbour = envelope_vertices[prev_index]
    last_neighbour = envelope_vertices[next_index]
    p, q, r = get_point(tri, u, first_neighbour, last_neighbour)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    mpq = (px + qx) / 2, (py + qy) / 2
    mpr = (px + rx) / 2, (py + ry) / 2
    g1 = DelaunayTriangulation.triangle_circumcenter(p, r, point)
    !all(isfinite, g1) && return NaN # point is one of p and r, i.e. we are interpolating at a data site
    g2 = DelaunayTriangulation.triangle_circumcenter(q, p, point)
    !all(isfinite, g2) && return NaN
    points = (mpq, mpr, g1, g2, mpq)
    return polygon_area(points)
end

function compute_sibson_coordinates(tri::Triangulation, envelope_vertices, point)
    coordinates = zeros(length(envelope_vertices) - 1)
    w = 0.0
    for i in firstindex(envelope_vertices):(lastindex(envelope_vertices) - 1)
        pre = max(0.0, pre_insertion_area(tri, i, envelope_vertices))
        post = max(0.0, post_insertion_area(tri, i, envelope_vertices, point))
        if isnan(post) # need to return the the vector λ = [1] since we are exactly at a data site
            return [1.0]
        end
        coordinates[i] = max(pre - post, 0.0) # take care of any precision issues
        w += coordinates[i]
    end
    coordinates ./= w
    return coordinates
end

function evaluate_sibson_interpolant(tri::Triangulation, z, point)
    envelope_vertices, _ = compute_envelope(tri, point)
    λ = compute_sibson_coordinates(tri, envelope_vertices, point)
    if length(λ) == 1
        for i in each_solid_vertex(tri)
            get_point(tri, i) == point && return z[i]
        end
    else
        itp = 0.0
        for (λ, k) in zip(λ, envelope_vertices)
            itp += λ * z[k]
        end
        return itp
    end
end

f = (x, y) -> sin(x * y) - cos(x - y) * exp(-(x - y)^2)
trit = triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 30, 30)
zz = [f(x, y) for (x, y) in DelaunayTriangulation.each_point(trit)]
xx = LinRange(0.001, 0.999, 20) # handling points on the boundary requires more care than we have discussed here
yy = LinRange(0.001, 0.999, 20)
fig = Figure(fontsize = 24)
ax = Axis(fig[1, 1], xlabel = L"x", ylabel = L"y", title = "True function", titlealign = :left, width = 400, height = 400)
contourf!(ax, xx, yy, f.(xx, yy'))
ax2 = Axis(fig[1, 2], xlabel = L"x", ylabel = L"y", title = "Interpolant", titlealign = :left, width = 400, height = 400)
zi = [evaluate_sibson_interpolant(trit, zz, (xᵢ, yᵢ)) for xᵢ in xx, yᵢ in yy]
contourf!(ax2, xx, yy, zi)
resize_to_layout!(fig)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

