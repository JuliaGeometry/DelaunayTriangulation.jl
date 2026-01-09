```@meta
EditURL = "https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/custom_primitive.jl"
```

# Using Custom Structs for Primitives and Boundaries

In this package, we allow for custom structs to be used for defining
points, edges, triangles, and boundaries (and parametric curves, as done for
example in the final example of [this tutorial](../tutorials/curve_bounded.md)).
This tutorial demonstrates how this can be done, showing what methods need to be
defined to use certain parts of the package. We note that there are some good defaults
already defined for most methods that can be overloaded for these primitives, but for
completely new structs you are required to define many new methods.
The packages we will be using are loaded below.

````@example custom_primitive
using DelaunayTriangulation
using CairoMakie
using Random
using StableRNGs
const DT = DelaunayTriangulation;
nothing #hide
````

Let us now define our custom structs.

````@example custom_primitive
struct CustomPoint
    x::Float64
    y::Float64
end
struct CustomPoints
    points::Vector{CustomPoint}
end
struct CustomSegment
    i::Int32
    j::Int32
end
struct CustomSegments
    segments::Set{CustomSegment}
end
struct CustomTriangle
    i::Int32
    j::Int32
    k::Int32
end
struct CustomTriangles
    triangles::Vector{CustomTriangle}
end
struct CustomPolygonSegment
    segments::Vector{CustomSegment}
end
struct CustomPolygon
    segments::Vector{CustomPolygonSegment}
end
struct CustomPolygons{N}
    polygons::NTuple{N, CustomPolygon}
end
````

Now, depending on your application you might not need to define all possible methods. For example, if you just want an unconstrained triangulation, all you need are `CustomPoint`
and `CustomPoints`. So, let's define our methods for increasing complexity. First, for unconstrained triangulations, all we need to define are:

````@example custom_primitive
DT.getx(p::CustomPoint) = p.x
DT.gety(p::CustomPoint) = p.y
DT.number_type(::Type{CustomPoint}) = Float64

Base.eachindex(points::CustomPoints) = Base.eachindex(points.points)
Base.iterate(points::CustomPoints, state...) = Base.iterate(points.points, state...)
Base.length(points::CustomPoints) = length(points.points)
Base.getindex(points::CustomPoints, i) = points.points[i]
DT.number_type(::Type{CustomPoints}) = Float64
DT.is_planar(::CustomPoints) = true

DT.geti(T::CustomTriangle) = T.i
DT.getj(T::CustomTriangle) = T.j
DT.getk(T::CustomTriangle) = T.k
DT.number_type(::Type{CustomTriangle}) = Int32
DT.construct_triangle(::Type{CustomTriangle}, i, j, k) = CustomTriangle(i, j, k)

Base.iterate(triangles::CustomTriangles, state...) = Base.iterate(triangles.triangles, state...)
Base.sizehint!(triangles::CustomTriangles, n) = sizehint!(triangles.triangles, n)
Base.eltype(::Type{CustomTriangles}) = CustomTriangle
Base.push!(triangles::CustomTriangles, T::CustomTriangle) = push!(triangles.triangles, T)
Base.delete!(triangles::CustomTriangles, T::CustomTriangle) = deleteat!(triangles.triangles, findfirst(==(T), triangles.triangles))
Base.length(triangles::CustomTriangles) = length(triangles.triangles)
CustomTriangles() = CustomTriangles(Vector{CustomTriangle}())
````

Now let's suppose we want to add in some segments. For this, we also need the following methods.

````@example custom_primitive
DT.construct_edge(::Type{CustomSegment}, i, j) = CustomSegment(i, j)
DT.initial(e::CustomSegment) = e.i
DT.terminal(e::CustomSegment) = e.j

Base.iterate(segments::CustomSegments, state...) = Base.iterate(segments.segments, state...)
Base.eltype(::Type{CustomSegments}) = CustomSegment
Base.push!(segments::CustomSegments, e::CustomSegment) = push!(segments.segments, e)
Base.delete!(segments::CustomSegments, e::CustomSegment) = delete!(segments.segments, e)
Base.rand(rng::AbstractRNG, segments::CustomSegments) = rand(rng, segments.segments)
Base.length(segments::CustomSegments) = length(segments.segments)
CustomSegments() = CustomSegments(Set{CustomSegment}())
````

Next, we want to allow for defining a boundary. For this, we need the following methods.

````@example custom_primitive
DT.has_multiple_curves(::CustomPolygons{N}) where {N} = N > 1
DT.has_multiple_curves(::CustomPolygon) = false
DT.has_multiple_curves(::CustomPolygonSegment) = false
DT.has_multiple_sections(::CustomPolygons) = true
DT.has_multiple_sections(::CustomPolygon) = true
DT.has_multiple_sections(::CustomPolygonSegment) = false
DT.num_curves(::CustomPolygons{N}) where {N} = N
DT.num_sections(poly::CustomPolygon) = length(poly.segments)
DT.num_boundary_edges(seg::CustomPolygonSegment) = length(seg.segments)
DT.get_boundary_nodes(poly::CustomPolygons, m::Integer) = poly.polygons[m] # go down to the mth polygon
DT.get_boundary_nodes(poly::CustomPolygon, m::Integer) = poly.segments[m] # go down to the mth segment
DT.get_boundary_nodes(seg::CustomPolygonSegment, m::Integer) = m > length(seg.segments) ? DT.terminal(seg.segments[m - 1]) : DT.initial(seg.segments[m]) # go down to the mth edge and extract the left node
DT.get_boundary_nodes(poly::CustomPolygons, (m, n)::NTuple{2, Int32}) = DT.get_boundary_nodes(DT.get_boundary_nodes(poly, m), n)
````

We now have all that we need for defining our custom primitives for constrained triangulations. We can go further and define methods
for working with Voronoi tessellations and centroidal Voronoi tessellations. For these methods, note that these only apply to unconstrained
triangulations, so you of course would not have to define the methods we have just defined for segments and boundaries. For Voronoi tessellations,
no extra methods are needed, except for

````@example custom_primitive
Base.empty!(segments::CustomSegments) = empty!(segments.segments)
````

if we specify `EdgesType` inside `triangulate` together with `clip=true` (and not otherwise). For centroidal Voronoi tessellations, we need

````@example custom_primitive
DT.set_point!(points::CustomPoints, i, x, y) = points.points[i] = CustomPoint(x, y)
````

We also need to consider methods needed for mesh refinement. We could also consider how we can define a method for a new parametric curve,
but an example of this is already given in the [curve bounded tutorial](../tutorials/curve_bounded.md), so we will only consider methods
with piecewise linear boundaries. For refinement, we need the following methods:

````@example custom_primitive
Base.pop!(points::CustomPoints) = pop!(points.points)
DT.push_point!(points::CustomPoints, x, y) = push!(points.points, CustomPoint(x, y))

DT.contains_edge(e::CustomSegment, Es::CustomSegments) = e ∈ Es.segments

Base.empty!(triangles::CustomTriangles) = empty!(triangles.triangles)

function Base.insert!(seg::CustomPolygonSegment, index, node)
    cur_segment = seg.segments[index - 1]
    u, v = edge_vertices(cur_segment)
    seg.segments[index - 1] = CustomSegment(u, node)
    insert!(seg.segments, index, CustomSegment(node, v))
    return seg
end
````

Now we have all that we need. Let's now demonstrate that this works. First, let's define the points, segments, and the boundary.

````@example custom_primitive
p1 = CustomPoint(0.0, 0.0)
p2 = CustomPoint(1.0, 0.0)
p3 = CustomPoint(1.0, 1.0)
p4 = CustomPoint(0.0, 1.0)
p5 = CustomPoint(0.5, 0.5)
p6 = CustomPoint(0.25, 0.25)
p7 = CustomPoint(0.75, 0.25)
p8 = CustomPoint(0.75, 0.75)
p9 = CustomPoint(0.25, 0.75)
points = CustomPoints([p1, p2, p3, p4, p5, p6, p7, p8, p9])
segments = CustomSegments(Set{CustomSegment}((CustomSegment(2, 7), CustomSegment(8, 3))))
outer_polygon = CustomPolygon(
    [
        CustomPolygonSegment([CustomSegment(1, 2), CustomSegment(2, 3)]),
        CustomPolygonSegment([CustomSegment(3, 4), CustomSegment(4, 1)]),
    ],
)
inner_polygon = CustomPolygon(
    [
        CustomPolygonSegment([CustomSegment(6, 9), CustomSegment(9, 8), CustomSegment(8, 7), CustomSegment(7, 6)]),
    ],
)
polygons = CustomPolygons((outer_polygon, inner_polygon));
nothing #hide
````

Now we triangulate and refine.

````@example custom_primitive
rng = StableRNG(123)
tri = triangulate(
    points; boundary_nodes = polygons, segments,
    IntegerType = Int32,
    EdgeType = CustomSegment,
    TriangleType = CustomTriangle,
    EdgesType = CustomSegments,
    TrianglesType = CustomTriangles,
    rng,
)
refine!(tri; max_area = 1.0e-3, rng)
fig, ax, sc = triplot(tri)
fig
````

Now let's give an example of a centroidal Voronoi tessellation to show that
this all works.

````@example custom_primitive
rng = StableRNG(123)
points = CustomPoints([p1, p2, p3, p4, p5, p6, p7, p8, p9])
tri = triangulate(
    points;
    IntegerType = Int32,
    EdgeType = CustomSegment,
    TriangleType = CustomTriangle,
    EdgesType = CustomSegments,
    TrianglesType = CustomTriangles,
    rng,
)
vorn = voronoi(tri; clip = true, smooth = true, rng)
fig, ax, sc = voronoiplot(vorn)
fig
````

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/custom_primitive.jl).

```julia
using DelaunayTriangulation
using CairoMakie
using Random
using StableRNGs
const DT = DelaunayTriangulation;

struct CustomPoint
    x::Float64
    y::Float64
end
struct CustomPoints
    points::Vector{CustomPoint}
end
struct CustomSegment
    i::Int32
    j::Int32
end
struct CustomSegments
    segments::Set{CustomSegment}
end
struct CustomTriangle
    i::Int32
    j::Int32
    k::Int32
end
struct CustomTriangles
    triangles::Vector{CustomTriangle}
end
struct CustomPolygonSegment
    segments::Vector{CustomSegment}
end
struct CustomPolygon
    segments::Vector{CustomPolygonSegment}
end
struct CustomPolygons{N}
    polygons::NTuple{N, CustomPolygon}
end

DT.getx(p::CustomPoint) = p.x
DT.gety(p::CustomPoint) = p.y
DT.number_type(::Type{CustomPoint}) = Float64

Base.eachindex(points::CustomPoints) = Base.eachindex(points.points)
Base.iterate(points::CustomPoints, state...) = Base.iterate(points.points, state...)
Base.length(points::CustomPoints) = length(points.points)
Base.getindex(points::CustomPoints, i) = points.points[i]
DT.number_type(::Type{CustomPoints}) = Float64
DT.is_planar(::CustomPoints) = true

DT.geti(T::CustomTriangle) = T.i
DT.getj(T::CustomTriangle) = T.j
DT.getk(T::CustomTriangle) = T.k
DT.number_type(::Type{CustomTriangle}) = Int32
DT.construct_triangle(::Type{CustomTriangle}, i, j, k) = CustomTriangle(i, j, k)

Base.iterate(triangles::CustomTriangles, state...) = Base.iterate(triangles.triangles, state...)
Base.sizehint!(triangles::CustomTriangles, n) = sizehint!(triangles.triangles, n)
Base.eltype(::Type{CustomTriangles}) = CustomTriangle
Base.push!(triangles::CustomTriangles, T::CustomTriangle) = push!(triangles.triangles, T)
Base.delete!(triangles::CustomTriangles, T::CustomTriangle) = deleteat!(triangles.triangles, findfirst(==(T), triangles.triangles))
Base.length(triangles::CustomTriangles) = length(triangles.triangles)
CustomTriangles() = CustomTriangles(Vector{CustomTriangle}())

DT.construct_edge(::Type{CustomSegment}, i, j) = CustomSegment(i, j)
DT.initial(e::CustomSegment) = e.i
DT.terminal(e::CustomSegment) = e.j

Base.iterate(segments::CustomSegments, state...) = Base.iterate(segments.segments, state...)
Base.eltype(::Type{CustomSegments}) = CustomSegment
Base.push!(segments::CustomSegments, e::CustomSegment) = push!(segments.segments, e)
Base.delete!(segments::CustomSegments, e::CustomSegment) = delete!(segments.segments, e)
Base.rand(rng::AbstractRNG, segments::CustomSegments) = rand(rng, segments.segments)
Base.length(segments::CustomSegments) = length(segments.segments)
CustomSegments() = CustomSegments(Set{CustomSegment}())

DT.has_multiple_curves(::CustomPolygons{N}) where {N} = N > 1
DT.has_multiple_curves(::CustomPolygon) = false
DT.has_multiple_curves(::CustomPolygonSegment) = false
DT.has_multiple_sections(::CustomPolygons) = true
DT.has_multiple_sections(::CustomPolygon) = true
DT.has_multiple_sections(::CustomPolygonSegment) = false
DT.num_curves(::CustomPolygons{N}) where {N} = N
DT.num_sections(poly::CustomPolygon) = length(poly.segments)
DT.num_boundary_edges(seg::CustomPolygonSegment) = length(seg.segments)
DT.get_boundary_nodes(poly::CustomPolygons, m::Integer) = poly.polygons[m] # go down to the mth polygon
DT.get_boundary_nodes(poly::CustomPolygon, m::Integer) = poly.segments[m] # go down to the mth segment
DT.get_boundary_nodes(seg::CustomPolygonSegment, m::Integer) = m > length(seg.segments) ? DT.terminal(seg.segments[m - 1]) : DT.initial(seg.segments[m]) # go down to the mth edge and extract the left node
DT.get_boundary_nodes(poly::CustomPolygons, (m, n)::NTuple{2, Int32}) = DT.get_boundary_nodes(DT.get_boundary_nodes(poly, m), n)

Base.empty!(segments::CustomSegments) = empty!(segments.segments)

DT.set_point!(points::CustomPoints, i, x, y) = points.points[i] = CustomPoint(x, y)

Base.pop!(points::CustomPoints) = pop!(points.points)
DT.push_point!(points::CustomPoints, x, y) = push!(points.points, CustomPoint(x, y))

DT.contains_edge(e::CustomSegment, Es::CustomSegments) = e ∈ Es.segments

Base.empty!(triangles::CustomTriangles) = empty!(triangles.triangles)

function Base.insert!(seg::CustomPolygonSegment, index, node)
    cur_segment = seg.segments[index - 1]
    u, v = edge_vertices(cur_segment)
    seg.segments[index - 1] = CustomSegment(u, node)
    insert!(seg.segments, index, CustomSegment(node, v))
    return seg
end

p1 = CustomPoint(0.0, 0.0)
p2 = CustomPoint(1.0, 0.0)
p3 = CustomPoint(1.0, 1.0)
p4 = CustomPoint(0.0, 1.0)
p5 = CustomPoint(0.5, 0.5)
p6 = CustomPoint(0.25, 0.25)
p7 = CustomPoint(0.75, 0.25)
p8 = CustomPoint(0.75, 0.75)
p9 = CustomPoint(0.25, 0.75)
points = CustomPoints([p1, p2, p3, p4, p5, p6, p7, p8, p9])
segments = CustomSegments(Set{CustomSegment}((CustomSegment(2, 7), CustomSegment(8, 3))))
outer_polygon = CustomPolygon(
    [
        CustomPolygonSegment([CustomSegment(1, 2), CustomSegment(2, 3)]),
        CustomPolygonSegment([CustomSegment(3, 4), CustomSegment(4, 1)]),
    ],
)
inner_polygon = CustomPolygon(
    [
        CustomPolygonSegment([CustomSegment(6, 9), CustomSegment(9, 8), CustomSegment(8, 7), CustomSegment(7, 6)]),
    ],
)
polygons = CustomPolygons((outer_polygon, inner_polygon));

rng = StableRNG(123)
tri = triangulate(
    points; boundary_nodes = polygons, segments,
    IntegerType = Int32,
    EdgeType = CustomSegment,
    TriangleType = CustomTriangle,
    EdgesType = CustomSegments,
    TrianglesType = CustomTriangles,
    rng,
)
refine!(tri; max_area = 1.0e-3, rng)
fig, ax, sc = triplot(tri)
fig

rng = StableRNG(123)
points = CustomPoints([p1, p2, p3, p4, p5, p6, p7, p8, p9])
tri = triangulate(
    points;
    IntegerType = Int32,
    EdgeType = CustomSegment,
    TriangleType = CustomTriangle,
    EdgesType = CustomSegments,
    TrianglesType = CustomTriangles,
    rng,
)
vorn = voronoi(tri; clip = true, smooth = true, rng)
fig, ax, sc = voronoiplot(vorn)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

