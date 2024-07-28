```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/curve_bounded.jl"
```

# Triangulating Curve-Bounded Domains

In this tutorial, we show how we can triangulate a domain that
is defined by curves rather than by straight lines. This tutorial makes use
of functions introduced in the [refinement tutorial](../tutorials/refinement.md).

Let us start by loading in the packages we will need.

````@example curve_bounded
using DelaunayTriangulation
using DelaunayTriangulation: EllipticalArc # CairoMakie also exports this
using CairoMakie
using StableRNGs
using LinearAlgebra
````

## Curves
In this package, only a small subset of curves are provided (although
arbitrary curves could be provided by the user, as discussed in
[here](../api/curves.md). The curves, defined as parametric curves over
a parameter $t ∈ [0, 1]$, are:

- [`LineSegment`](@ref): Just a simple line segment between two points.
- [`CircularArc`](@ref): A circular arc defined between two points and a center.
- [`EllipticalArc`](@ref): An elliptical arc defined by two points, a center, a major radius, a minor radius, and an angle of rotation.
- [`BezierCurve`](@ref): A Bézier curve defined by a set of control points.
- [`BSpline`](@ref): A B-spline curve defined by a set of control points. Defaults to a cubic B-spline.
- [`CatmullRomSpline`](@ref): A Catmull-Rom spline defined by a set of control points.

This set of curves is sufficient to represent a large number of curves.

## A Circular Domain
We start with a simple domain: a circle. In particular, we consider a circle of
radius $r = 2$ centered at $(x, y) = (1/2, 2)$. In past tutorials, we would have defined a
domain like this by a set of points around a boundary, e.g.:

````@example curve_bounded
n = 50
r = 2.0
xc, yc = 1 / 2, 2.0
θ = range(0, 2π, length = n + 1) |> collect;
θ[end] = θ[begin];
x = xc .+ r * cos.(θ)
y = yc .+ r * sin.(θ);
nothing #hide
````

One problem with this approach is that we had to decide what value of $n$ to use for discretising this
boundary. Instead, we can use `CircularArc`. For this, we use:

````@example curve_bounded
p = (xc + r, yc)
c = (xc, yc)
arc = CircularArc(p, p, c)
````

Here, the syntax is `CircularArc(first_point, last_point, centre)`. Since the circle is closed, we
use `p` for both `first_point` and `last_point`. Notice that the `arc` is a function. In particular,

````@example curve_bounded
typeof(arc) |> supertype |> supertype
````

If we wanted to look at this circle, we would need to evaluate it at a set of $t \in [0, 1]$.

````@example curve_bounded
t = LinRange(0, 1, 2500)
points = arc.(t)
````

````@example curve_bounded
fig, ax, sc = lines(points)
````

Let's now triangulate this domain. We need to put the arc into its own vector, and we still need to pass a set of
points into `triangulate`:

````@example curve_bounded
points = NTuple{2, Float64}[]
rng = StableRNG(123)
tri = triangulate(points; boundary_nodes = [arc], rng)
````

````@example curve_bounded
fig, ax, sc = triplot(tri)
fig
````

Notice that the domain doesn't look like a circle yet. This is because using `triangulate` on the
curve just by itself isn't enough. In fact, the triangulation returned in this case is simply one where:

- No point is contained in the interior of any boundary edge's diametral circle.
- The total variation between any two neighbouring points is less than $\pi/2$, meaning that the absolute change in the angle along the parametric curve between these points is less than $\pi/2$.

This is probably not what we actually want, though. Instead, we need to refine the domain using mesh refinement.
The syntax for this is the same as in the [refinement tutorial](../tutorials/refinement.md):

````@example curve_bounded
refine!(tri; max_area = 1.0e-1, rng)
fig, ax, sc = triplot(tri)
fig
````

Much better! We have now triangulated our first curve-bounded domain.

## A Boundary Defined by Multiple Parametric Curves
We now take a domain which is defined by three separate curves: a circular arc, a B-spline, and a piecewise linear curve.
For piecewise linear curves, we use the same method as we would in previous tutorials, where instead of using coordinates to define
the boundary we use numbers that refer to points in the point set. To start, the point set we will be using is

````@example curve_bounded
points = [(0.1, 0.1), (0.15, 0.15), (0.23, 0.23), (0.009, 0.11), (0.0, -2.0), (0.2, -1.7), (0.000591, 0.00019), (0.111, -0.005), (-0.0001, -0.00991), (1.0, 0.0)]
````

Now, for the boundary, we will take:

- A circular arc defined between $(1, 0)$ and $(0, 1)$ centred at $(0, 0)$.
- A B-spline with control points $(0,1)$, $(-1, 2)$, $(-2,0)$, $(-2,-1)$, and $(0,-2)$.
- A piecewise linear curve defined between the fixth, sixth, and tenth points in `points`.

We can define these curves as follows:

````@example curve_bounded
arc = CircularArc((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))
bspl = BSpline([(0.0, 1.0), (-1.0, 2.0), (-2.0, 0.0), (-2.0, -1.0), (0.0, -2.0)])
pce = [5, 6, 10]
````

Notice that we still must make sure that the curves connect, and that together the curves
define a positively-oriented boundary. The domain we get from this looks like:

````@example curve_bounded
t = LinRange(0, 1, 1000)
pts = vcat(arc.(t), bspl.(t), points[pce])
fig, ax, sc = lines(pts)
fig
````

Let's now get a triangulation of this domain. We will use a custom constraint to force triangles closer
to the origin to be smaller than those outside of it.

````@example curve_bounded
curve = [[arc], [bspl], pce]
rng = StableRNG(123)
tri = triangulate(points; boundary_nodes = curve, rng)
refine!(
    tri; max_area = 1.0e-2, rng, custom_constraint = (_tri, T) -> begin
        i, j, k = triangle_vertices(T)
        p, q, r = get_point(_tri, i, j, k)
        c = (p .+ q .+ r) ./ 3
        return norm(c) < 1 / 2 && DelaunayTriangulation.triangle_area(p, q, r) > 1.0e-3
    end,
)
fig, ax, sc = triplot(tri)
fig
````

## A Complicated Multiply-Connected Disjoint Domain
For our last example, we take a complicated case with a domain that is disjoint, and where the individual
domains are multiply-connected. Let us give the domain followed by an explanation of how it is defined:

````@example curve_bounded
curve = [
    [
        [1, 2, 3], [EllipticalArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0), 2, 1 / 2, 0.0)],
    ],
    [
        [BSpline([(0.0, 0.4), (1.0, 0.2), (0.0, 0.1), (-1.0, 0.2), (0.0, 0.4)])],
    ],
    [
        [4, 5, 6, 7, 4],
    ],
    [
        [BezierCurve([(0.0, -2.0), (0.0, -2.5), (-1.0, -2.5), (-1.0, -3.0)])], [CatmullRomSpline([(-1.0, -3.0), (0.0, -4.0), (1.0, -3.0), (0.0, -2.0)])],
    ],
    [
        [12, 11, 10, 12],
    ],
    [
        [CircularArc((1.1, -3.0), (1.1, -3.0), (0.0, -3.0), positive = false)],
    ],
]
points = [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-2.0, -5.0), (2.0, -5.0), (2.0, -1 / 10), (-2.0, -1 / 10), (-1.0, -3.0), (0.0, -4.0), (0.0, -2.3), (-0.5, -3.5), (0.9, -3.0)]
t = LinRange(0, 1, 1000)
fig
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, [get_point(points, curve[1][1]...)...], color = :red, label = "(1, 2, 3)")
lines!(ax, curve[1][2][1].(t), color = :red, linestyle = :dashdot, label = "EllipticalArc")
lines!(ax, curve[2][1][1].(t), color = :green, label = "BSpline")
lines!(ax, [get_point(points, curve[3][1]...)...], color = :blue, label = "(4, 5, 6, 7, 4)")
lines!(ax, curve[4][1][1].(t), color = :purple, label = "BezierCurve")
lines!(ax, curve[4][2][1].(t), color = :purple, linestyle = :dashdot, label = "CatmullRomSpline")
lines!(ax, [get_point(points, curve[5][1]...)...], color = :orange, label = "(12, 11, 10, 12)")
lines!(ax, curve[6][1][1].(t), color = :black, label = "CircularArc")
fig[1, 2] = Legend(fig, ax, "Curve")
fig
````

Let's walk through the definition of `curve`.

- The first domain that is defined is the red curve in the above figure, defined in terms of a piecewise linear portion and an elliptical arc.
  The elliptical arc is defined by the points $(2, 0)$, $(-2, 0)$, and $(0, 0)$, with a major radius of $2$, a minor radius of $1/2$, and an angle of rotation of $0$.
- The second domain that is defined is the green curve, defined as a cubic B-spline with control points
  $(0, 0.4)$, $(1, 0.2)$, $(0, 0.1)$, $(-1, 0.2)$, $(0, 0.4)$.
- The third domain that we define is the square blue curve.
- The fourth domain that we define is the purple curve, defined by a Bézier curve and a Catmull-Rom spline.
  The Bézier curve is defined by the control points $(-1, -3)$, $(-1, -2.5)$, $(0, -2.5)$, and $(0, -2)$. The Catmull-Rom spline is defined by the control points
  $(0, -2)$, $(1, -3)$, $(0, -4)$, and $(-1, -3)$.
- The fifth domain that we define is the orange curve, defined as a piecewise linear curve to represent a triangle.
- The last domain that we define is a circular arc. To get the orientation correct, i.e. to make sure the circle is defined clockwise so that the domain
  remains positively oriented, we must use `positive=false`.

In addition to ensuring that the curves are all oriented correctly, care has also been taken in the definition of these curves to make sure that
the curves connect at the correct points.

Let's now triangulate.

````@example curve_bounded
rng = StableRNG(123)
tri = triangulate(copy(points); boundary_nodes = curve, rng) # copying so that we don't mutate for the next section
refine!(tri; max_area = 1.0e-2)
fig, ax, sc = triplot(tri)
fig
````

### Using Custom Constraints to Control Refinement
Let's give another example of using custom constraints to better control the refinement within different domains. Referencing the
previous figure where we showed each domain by colour, let us try and use a coarse mesh in the region bounded between the red and green
curves, a slightly denser mesh bounded between the blue and black curves, and finally a dense mesh in the region bounded between the
purple and orange curves. To do this, we must have a method for deciding which region a given point resides in. This is what
`find_polygon` for is, as we also use in the [point-in-polygon tutorial](../tutorials/point_in_polygon.md).
To write this function, we note that the indices of these polygons are 1, 3, and 4 for the red, blue, and purple regions,
respectively.

````@example curve_bounded
poly_constraint = (_tri, T) -> begin
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(_tri, i, j, k)
    c = (p .+ q .+ r) ./ 3
    idx = find_polygon(_tri, c)
    if idx ∉ (1, 3, 4)
        return true
    end
    max_area = if idx == 1 # coarse
        1.0e-1
    elseif idx == 3 # medium
        1.0e-2
    else # dense
        1.0e-3
    end
    area = DelaunayTriangulation.triangle_area(p, q, r)
    return area > max_area
end
rng = StableRNG(123)
tri = triangulate(points; boundary_nodes = curve, rng)
refine!(tri; custom_constraint = poly_constraint, rng)
fig, ax, sc = triplot(tri)
fig
````

## Defining a New Parametric Curve
Let us now give an example where we define a domain by a parametric curve that is not provided natively by this package. For this
example, we consider the astroid, where
```math
\begin{aligned}
x(t) = \cos^3(2\pi t), \\
y(t) = \sin^3(2\pi t).
\end{aligned}
```
Following the docstring of [`DelaunayTriangulation.AbstractParametricCurve`](@ref), we know that to define this curve to be compatible
with this package we need:

- To represent the curve as a callable struct that subtypes `AbstractParametricCurve` that maps `Float64 -> NTuple{2,Float64}`.
- To define [`DelaunayTriangulation.differentiate`](@ref), [`DelaunayTriangulation.twice_differentiate`](@ref), and [`DelaunayTriangulation.thrice_differentiate`](@ref).
- Have a `lookup_table` field that maps `lookup_table[i]` to `(i - 1) / (length(lookup_table) - 1)`, where `lookup_table` is a `Dict`.
- Have defined the parametric curve according to $0 ≤ t ≤ 1$ (already done).

Let's now meet these requirements.

````@example curve_bounded
struct Astroid <: DelaunayTriangulation.AbstractParametricCurve
    lookup_table::Vector{NTuple{2, Float64}}
end
function (c::Astroid)(t)
    if t == 0.0 || t == 1.0
        return (1.0, 0.0)
    end
    x = cos(2π * t)^3
    y = sin(2π * t)^3
    return (x, y)
end
function DelaunayTriangulation.differentiate(c::Astroid, t)
    x = -6π * sin(2π * t) * cos(2π * t)^2
    y = 6π * sin(2π * t)^2 * cos(2π * t)
    return (x, y)
end
function DelaunayTriangulation.twice_differentiate(c::Astroid, t)
    x = -12π^2 * cos(2π * t) * (cos(2π * t)^2 - 2sin(2π * t)^2)
    y = -12π^2 * sin(2π * t) * (sin(2π * t)^2 - 2cos(2π * t)^2)
    return (x, y)
end
function DelaunayTriangulation.thrice_differentiate(c::Astroid, t)
    x = 24π^3 * sin(2π * t) * (7cos(2π * t)^2 - 2sin(2π * t)^2)
    y = 24π^3 * cos(2π * t) * (2cos(2π * t)^2 - 7sin(2π * t)^2)
    return (x, y)
end
````

Let's now define an astroid curve and triangulate it.

````@example curve_bounded
function Astroid(n::Int)
    lookup_table = Vector{NTuple{2, Float64}}(undef, n)
    c = Astroid(lookup_table)
    for i in 1:n
        lookup_table[i] = c((i - 1) / (n - 1))
    end
    return Astroid(lookup_table)
end
rng = StableRNG(123)
curve = Astroid(1000)
tri = triangulate(NTuple{2, Float64}[]; boundary_nodes = [curve], rng)
refine!(tri; max_area = 1.0e-2)
fig, ax, sc = triplot(tri)
fig
````

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/curve_bounded.jl).

```julia
using DelaunayTriangulation
using DelaunayTriangulation: EllipticalArc # CairoMakie also exports this
using CairoMakie
using StableRNGs
using LinearAlgebra

n = 50
r = 2.0
xc, yc = 1 / 2, 2.0
θ = range(0, 2π, length = n + 1) |> collect;
θ[end] = θ[begin];
x = xc .+ r * cos.(θ)
y = yc .+ r * sin.(θ);

p = (xc + r, yc)
c = (xc, yc)
arc = CircularArc(p, p, c)

typeof(arc) |> supertype |> supertype

t = LinRange(0, 1, 2500)
points = arc.(t)

fig, ax, sc = lines(points)

points = NTuple{2, Float64}[]
rng = StableRNG(123)
tri = triangulate(points; boundary_nodes = [arc], rng)

fig, ax, sc = triplot(tri)
fig

refine!(tri; max_area = 1.0e-1, rng)
fig, ax, sc = triplot(tri)
fig

points = [(0.1, 0.1), (0.15, 0.15), (0.23, 0.23), (0.009, 0.11), (0.0, -2.0), (0.2, -1.7), (0.000591, 0.00019), (0.111, -0.005), (-0.0001, -0.00991), (1.0, 0.0)]

arc = CircularArc((1.0, 0.0), (0.0, 1.0), (0.0, 0.0))
bspl = BSpline([(0.0, 1.0), (-1.0, 2.0), (-2.0, 0.0), (-2.0, -1.0), (0.0, -2.0)])
pce = [5, 6, 10]

t = LinRange(0, 1, 1000)
pts = vcat(arc.(t), bspl.(t), points[pce])
fig, ax, sc = lines(pts)
fig

curve = [[arc], [bspl], pce]
rng = StableRNG(123)
tri = triangulate(points; boundary_nodes = curve, rng)
refine!(
    tri; max_area = 1.0e-2, rng, custom_constraint = (_tri, T) -> begin
        i, j, k = triangle_vertices(T)
        p, q, r = get_point(_tri, i, j, k)
        c = (p .+ q .+ r) ./ 3
        return norm(c) < 1 / 2 && DelaunayTriangulation.triangle_area(p, q, r) > 1.0e-3
    end,
)
fig, ax, sc = triplot(tri)
fig

curve = [
    [
        [1, 2, 3], [EllipticalArc((2.0, 0.0), (-2.0, 0.0), (0.0, 0.0), 2, 1 / 2, 0.0)],
    ],
    [
        [BSpline([(0.0, 0.4), (1.0, 0.2), (0.0, 0.1), (-1.0, 0.2), (0.0, 0.4)])],
    ],
    [
        [4, 5, 6, 7, 4],
    ],
    [
        [BezierCurve([(0.0, -2.0), (0.0, -2.5), (-1.0, -2.5), (-1.0, -3.0)])], [CatmullRomSpline([(-1.0, -3.0), (0.0, -4.0), (1.0, -3.0), (0.0, -2.0)])],
    ],
    [
        [12, 11, 10, 12],
    ],
    [
        [CircularArc((1.1, -3.0), (1.1, -3.0), (0.0, -3.0), positive = false)],
    ],
]
points = [(-2.0, 0.0), (0.0, 0.0), (2.0, 0.0), (-2.0, -5.0), (2.0, -5.0), (2.0, -1 / 10), (-2.0, -1 / 10), (-1.0, -3.0), (0.0, -4.0), (0.0, -2.3), (-0.5, -3.5), (0.9, -3.0)]
t = LinRange(0, 1, 1000)
fig
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, [get_point(points, curve[1][1]...)...], color = :red, label = "(1, 2, 3)")
lines!(ax, curve[1][2][1].(t), color = :red, linestyle = :dashdot, label = "EllipticalArc")
lines!(ax, curve[2][1][1].(t), color = :green, label = "BSpline")
lines!(ax, [get_point(points, curve[3][1]...)...], color = :blue, label = "(4, 5, 6, 7, 4)")
lines!(ax, curve[4][1][1].(t), color = :purple, label = "BezierCurve")
lines!(ax, curve[4][2][1].(t), color = :purple, linestyle = :dashdot, label = "CatmullRomSpline")
lines!(ax, [get_point(points, curve[5][1]...)...], color = :orange, label = "(12, 11, 10, 12)")
lines!(ax, curve[6][1][1].(t), color = :black, label = "CircularArc")
fig[1, 2] = Legend(fig, ax, "Curve")
fig

rng = StableRNG(123)
tri = triangulate(copy(points); boundary_nodes = curve, rng) # copying so that we don't mutate for the next section
refine!(tri; max_area = 1.0e-2)
fig, ax, sc = triplot(tri)
fig

poly_constraint = (_tri, T) -> begin
    i, j, k = triangle_vertices(T)
    p, q, r = get_point(_tri, i, j, k)
    c = (p .+ q .+ r) ./ 3
    idx = find_polygon(_tri, c)
    if idx ∉ (1, 3, 4)
        return true
    end
    max_area = if idx == 1 # coarse
        1.0e-1
    elseif idx == 3 # medium
        1.0e-2
    else # dense
        1.0e-3
    end
    area = DelaunayTriangulation.triangle_area(p, q, r)
    return area > max_area
end
rng = StableRNG(123)
tri = triangulate(points; boundary_nodes = curve, rng)
refine!(tri; custom_constraint = poly_constraint, rng)
fig, ax, sc = triplot(tri)
fig

struct Astroid <: DelaunayTriangulation.AbstractParametricCurve
    lookup_table::Vector{NTuple{2, Float64}}
end
function (c::Astroid)(t)
    if t == 0.0 || t == 1.0
        return (1.0, 0.0)
    end
    x = cos(2π * t)^3
    y = sin(2π * t)^3
    return (x, y)
end
function DelaunayTriangulation.differentiate(c::Astroid, t)
    x = -6π * sin(2π * t) * cos(2π * t)^2
    y = 6π * sin(2π * t)^2 * cos(2π * t)
    return (x, y)
end
function DelaunayTriangulation.twice_differentiate(c::Astroid, t)
    x = -12π^2 * cos(2π * t) * (cos(2π * t)^2 - 2sin(2π * t)^2)
    y = -12π^2 * sin(2π * t) * (sin(2π * t)^2 - 2cos(2π * t)^2)
    return (x, y)
end
function DelaunayTriangulation.thrice_differentiate(c::Astroid, t)
    x = 24π^3 * sin(2π * t) * (7cos(2π * t)^2 - 2sin(2π * t)^2)
    y = 24π^3 * cos(2π * t) * (2cos(2π * t)^2 - 7sin(2π * t)^2)
    return (x, y)
end

function Astroid(n::Int)
    lookup_table = Vector{NTuple{2, Float64}}(undef, n)
    c = Astroid(lookup_table)
    for i in 1:n
        lookup_table[i] = c((i - 1) / (n - 1))
    end
    return Astroid(lookup_table)
end
rng = StableRNG(123)
curve = Astroid(1000)
tri = triangulate(NTuple{2, Float64}[]; boundary_nodes = [curve], rng)
refine!(tri; max_area = 1.0e-2)
fig, ax, sc = triplot(tri)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

