```@meta
EditURL = "https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/main/docs/src/literate_tutorials/custom_primitive.jl"
```

# Using Custom Structs for Primitives and Boundaries

In this package, we allow for custom structs to be used for defining
points, edges, triangles, and boundaries. This tutorial
demonstrates how this can be done, showing what methods need
to be defined to use certain parts of the package. The packages we will
be using are loaded below.

````@example custom_primitive
using DelaunayTriangulation
using CairoMakie
using StableRNGs
const DT = DelaunayTriangulation;
nothing #hide
````

## Structs
We start by defining the structs we will be using.

````@example custom_primitive
#### Structs
struct CustomPoint
    x::Float64
    y::Float64
end
struct CustomPoints
    points::Vector{CustomPoint}
end
struct CustomEdge
    i::Int32
    j::Int32
end
struct CustomEdges
    edges::Set{CustomEdge}
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
    edges::Vector{CustomEdge}
end
struct CustomPolygon
    segments::Vector{CustomPolygonSegment}
end
struct CustomPolygons{N}
    polygons::NTuple{N,CustomPolygon}
end
````

## Methods
Now we define the methods that are needed, breaking them
into what methods are needed for each unconstrained, constrained,
and refined triangulations.

````@example custom_primitive
#### Methods
### Definitions needed for unconstrained triangulations
# Point
DT.getx(p::CustomPoint) = p.x
DT.gety(p::CustomPoint) = p.y
DT.number_type(::Type{CustomPoint}) = Float64

# Points
DT.each_point_index(pts::CustomPoints) = eachindex(pts.points)
DT.num_points(pts::CustomPoints) = length(pts.points)
DT.each_point(pts::CustomPoints) = pts.points
DT.number_type(::Type{CustomPoints}) = DT.number_type(CustomPoint)
DT.getpoint(pts::CustomPoints, i::Integer) = pts.points[i]

# Edge
DT.construct_edge(::Type{CustomEdge}, i, j) = CustomEdge(i, j)
DT.initial(e::CustomEdge) = e.i
DT.terminal(e::CustomEdge) = e.j

# Edges
DT.initialise_edges(::Type{CustomEdges}) = CustomEdges(Set{CustomEdge}())
DT.add_to_edges!(edges::CustomEdges, edge::CustomEdge) = push!(edges.edges, edge)
Base.iterate(edges::CustomEdges, state...) = Base.iterate(edges.edges, state...)
DT.each_edge(edges::CustomEdges) = edges.edges
DT.delete_from_edges!(edges::CustomEdges, edge::CustomEdge) = delete!(edges.edges, edge)
DT.contains_edge(edge::CustomEdge, edges::CustomEdges) = edge ∈ edges.edges
DT.num_edges(edges::CustomEdges) = length(edges.edges)

# Triangle
DT.construct_triangle(::Type{CustomTriangle}, i, j, k) = CustomTriangle(i, j, k)
DT.geti(tri::CustomTriangle) = tri.i
DT.getj(tri::CustomTriangle) = tri.j
DT.getk(tri::CustomTriangle) = tri.k
DT.integer_type(::Type{CustomTriangle}) = Int32

# Triangles
function DT.delete_from_triangles!(tri::CustomTriangles, triangle::CustomTriangle)
    i = findfirst(==(triangle), tri.triangles)
    deleteat!(tri.triangles, i)
    return nothing
end
DT.initialise_triangles(::Type{CustomTriangles}) = CustomTriangles(Vector{CustomTriangle}())
DT.sizehint!(tri::CustomTriangles, n) = sizehint!(tri.triangles, n)
DT.triangle_type(::Type{CustomTriangles}) = CustomTriangle
DT.add_to_triangles!(tri::CustomTriangles, triangle::CustomTriangle) = push!(tri.triangles, triangle)
DT.num_triangles(tri::CustomTriangles) = length(tri.triangles)
Base.iterate(tri::CustomTriangles, state...) = Base.iterate(tri.triangles, state...)
Base.in(T, V::CustomTriangles) = T ∈ V.triangles

### Definitions needed for constrained segments
# Edge
DT.integer_type(::Type{CustomEdge}) = Int32

# Edges
DT.edge_type(::Type{CustomEdges}) = CustomEdge

### Definitions needed for boundary nodes
# Triangles
DT.each_triangle(tri::CustomTriangles) = tri.triangles

# Boundary nodes
DT.has_multiple_curves(::CustomPolygons{N}) where {N} = N > 1
DT.has_multiple_curves(::CustomPolygon) = false
DT.has_multiple_curves(::CustomPolygonSegment) = false
DT.has_multiple_segments(::CustomPolygons) = true
DT.has_multiple_segments(poly::CustomPolygon) = true
DT.has_multiple_segments(seg::CustomPolygonSegment) = false
DT.num_curves(::CustomPolygons{N}) where {N} = N
DT.getboundarynodes(poly::CustomPolygons, m) = poly.polygons[m] # go down to the mth polygon
DT.getboundarynodes(poly::CustomPolygon, m) = poly.segments[m] # go down to the mth segment
DT.getboundarynodes(seg::CustomPolygonSegment, m) = m > length(seg.edges) ? DT.terminal(seg.edges[m-1]) : DT.initial(seg.edges[m]) # go down to the mth edge and extract the left node
DT.getboundarynodes(poly::CustomPolygons, (m, n)::NTuple{2,Int32}) = DT.getboundarynodes(DT.getboundarynodes(poly, m), n)
DT.num_segments(poly::CustomPolygon) = length(poly.segments)
DT.num_boundary_edges(seg::CustomPolygonSegment) = length(seg.edges)

### Definitions needed for refinement
# Points
DT.push_point!(pts::CustomPoints, x, y) = push!(pts.points, CustomPoint(x, y))
DT.pop_point!(pts::CustomPoints) = pop!(pts.points)

# Edges
DT.contains_edge(e::CustomEdge, Es::CustomEdges) = e ∈ Es.edges

# Boundary nodes
function Base.insert!(seg::CustomPolygonSegment, index, node)
    cur_edge = seg.edges[index-1]
    u, v = DT.initial(cur_edge), DT.terminal(cur_edge)
    seg.edges[index-1] = CustomEdge(u, node)
    insert!(seg.edges, index, CustomEdge(node, v))
    return nothing
end

### Definitions needed for centroidal tessellations
# Points
DT.set_point!(points::CustomPoints, i, x, y) = points.points[i] = CustomPoint(x, y)
````

## Example
Let us now give an example just to demonstrate that this all actually works.
First, we build up the points, edges, and boundary nodes.

````@example custom_primitive
# Example
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
````

````@example custom_primitive
edges = CustomEdges(Set{CustomEdge}((CustomEdge(2, 7), CustomEdge(8, 3))))
````

````@example custom_primitive
outer_polygon = CustomPolygon([
    CustomPolygonSegment([CustomEdge(1, 2), CustomEdge(2, 3)]),
    CustomPolygonSegment([CustomEdge(3, 4), CustomEdge(4, 1)]),
])
inner_polygon = CustomPolygon([
    CustomPolygonSegment([CustomEdge(6, 9), CustomEdge(9, 8), CustomEdge(8, 7), CustomEdge(7, 6)]),
])
polygons = CustomPolygons((outer_polygon, inner_polygon))
````

Now we triangulate. To use the custom definitions, we need to pass the
types as keyword arguments.

````@example custom_primitive
rng = StableRNG(123)
tri = triangulate(points; edges=edges, boundary_nodes=polygons,
    IntegerType=Int32,
    EdgeType=CustomEdge,
    TriangleType=CustomTriangle,
    EdgesType=CustomEdges,
    TrianglesType=CustomTriangles,
    delete_ghosts=false,
    rng
)
````

We can then refine the triangulation as usual.

````@example custom_primitive
A = get_total_area(tri)
stats = refine!(tri; max_area=1e-3A, min_angle=27.3, rng);
nothing #hide
````

We can then plot the triangulation.

````@example custom_primitive
fig, ax, sc = triplot(tri)
````

Similarly, we can work with Voronoi tessellations. We use an unconstrained triangulation
so that we can also demonstrate clipped and centroidal tessellations.

````@example custom_primitive
tri = triangulate(points;
    IntegerType=Int32,
    EdgeType=CustomEdge,
    TriangleType=CustomTriangle,
    EdgesType=CustomEdges,
    TrianglesType=CustomTriangles,
    delete_ghosts=false,
    rng
)
vorn = voronoi(tri)
````

````@example custom_primitive
clipped_vorn = voronoi(tri, true)
````

````@example custom_primitive
centroidal_vorn = centroidal_smooth(clipped_vorn; rng)
````

Here are the plots of these tessellations.

````@example custom_primitive
fig = Figure(fontsize=43)
ax1 = Axis(fig[1, 1], title="Tessellation", width=600, height=400)
ax2 = Axis(fig[1, 2], title="Clipped Tessellation", width=600, height=400)
ax3 = Axis(fig[1, 3], title="Centroidal Tessellation", width=600, height=400)
voronoiplot!(ax1, vorn, colormap=:matter, strokewidth=2, markersize=7)
voronoiplot!(ax2, clipped_vorn, colormap=:matter, strokewidth=2, markersize=7)
voronoiplot!(ax3, centroidal_vorn, colormap=:matter, strokewidth=2, markersize=7)
resize_to_layout!(fig)
fig
````

## Just the code
An uncommented version of this example is given below.
You can view the source code for this file [here](https://github.com/DanielVandH/DelaunayTriangulation.jl/tree/new-docs/docs/src/literate_tutorials/custom_primitive.jl).

```julia
using DelaunayTriangulation
using CairoMakie
using StableRNGs
const DT = DelaunayTriangulation;

#### Structs
struct CustomPoint
    x::Float64
    y::Float64
end
struct CustomPoints
    points::Vector{CustomPoint}
end
struct CustomEdge
    i::Int32
    j::Int32
end
struct CustomEdges
    edges::Set{CustomEdge}
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
    edges::Vector{CustomEdge}
end
struct CustomPolygon
    segments::Vector{CustomPolygonSegment}
end
struct CustomPolygons{N}
    polygons::NTuple{N,CustomPolygon}
end

#### Methods
### Definitions needed for unconstrained triangulations
# Point
DT.getx(p::CustomPoint) = p.x
DT.gety(p::CustomPoint) = p.y
DT.number_type(::Type{CustomPoint}) = Float64

# Points
DT.each_point_index(pts::CustomPoints) = eachindex(pts.points)
DT.num_points(pts::CustomPoints) = length(pts.points)
DT.each_point(pts::CustomPoints) = pts.points
DT.number_type(::Type{CustomPoints}) = DT.number_type(CustomPoint)
DT.getpoint(pts::CustomPoints, i::Integer) = pts.points[i]

# Edge
DT.construct_edge(::Type{CustomEdge}, i, j) = CustomEdge(i, j)
DT.initial(e::CustomEdge) = e.i
DT.terminal(e::CustomEdge) = e.j

# Edges
DT.initialise_edges(::Type{CustomEdges}) = CustomEdges(Set{CustomEdge}())
DT.add_to_edges!(edges::CustomEdges, edge::CustomEdge) = push!(edges.edges, edge)
Base.iterate(edges::CustomEdges, state...) = Base.iterate(edges.edges, state...)
DT.each_edge(edges::CustomEdges) = edges.edges
DT.delete_from_edges!(edges::CustomEdges, edge::CustomEdge) = delete!(edges.edges, edge)
DT.contains_edge(edge::CustomEdge, edges::CustomEdges) = edge ∈ edges.edges
DT.num_edges(edges::CustomEdges) = length(edges.edges)

# Triangle
DT.construct_triangle(::Type{CustomTriangle}, i, j, k) = CustomTriangle(i, j, k)
DT.geti(tri::CustomTriangle) = tri.i
DT.getj(tri::CustomTriangle) = tri.j
DT.getk(tri::CustomTriangle) = tri.k
DT.integer_type(::Type{CustomTriangle}) = Int32

# Triangles
function DT.delete_from_triangles!(tri::CustomTriangles, triangle::CustomTriangle)
    i = findfirst(==(triangle), tri.triangles)
    deleteat!(tri.triangles, i)
    return nothing
end
DT.initialise_triangles(::Type{CustomTriangles}) = CustomTriangles(Vector{CustomTriangle}())
DT.sizehint!(tri::CustomTriangles, n) = sizehint!(tri.triangles, n)
DT.triangle_type(::Type{CustomTriangles}) = CustomTriangle
DT.add_to_triangles!(tri::CustomTriangles, triangle::CustomTriangle) = push!(tri.triangles, triangle)
DT.num_triangles(tri::CustomTriangles) = length(tri.triangles)
Base.iterate(tri::CustomTriangles, state...) = Base.iterate(tri.triangles, state...)
Base.in(T, V::CustomTriangles) = T ∈ V.triangles

### Definitions needed for constrained segments
# Edge
DT.integer_type(::Type{CustomEdge}) = Int32

# Edges
DT.edge_type(::Type{CustomEdges}) = CustomEdge

### Definitions needed for boundary nodes
# Triangles
DT.each_triangle(tri::CustomTriangles) = tri.triangles

# Boundary nodes
DT.has_multiple_curves(::CustomPolygons{N}) where {N} = N > 1
DT.has_multiple_curves(::CustomPolygon) = false
DT.has_multiple_curves(::CustomPolygonSegment) = false
DT.has_multiple_segments(::CustomPolygons) = true
DT.has_multiple_segments(poly::CustomPolygon) = true
DT.has_multiple_segments(seg::CustomPolygonSegment) = false
DT.num_curves(::CustomPolygons{N}) where {N} = N
DT.getboundarynodes(poly::CustomPolygons, m) = poly.polygons[m] # go down to the mth polygon
DT.getboundarynodes(poly::CustomPolygon, m) = poly.segments[m] # go down to the mth segment
DT.getboundarynodes(seg::CustomPolygonSegment, m) = m > length(seg.edges) ? DT.terminal(seg.edges[m-1]) : DT.initial(seg.edges[m]) # go down to the mth edge and extract the left node
DT.getboundarynodes(poly::CustomPolygons, (m, n)::NTuple{2,Int32}) = DT.getboundarynodes(DT.getboundarynodes(poly, m), n)
DT.num_segments(poly::CustomPolygon) = length(poly.segments)
DT.num_boundary_edges(seg::CustomPolygonSegment) = length(seg.edges)

### Definitions needed for refinement
# Points
DT.push_point!(pts::CustomPoints, x, y) = push!(pts.points, CustomPoint(x, y))
DT.pop_point!(pts::CustomPoints) = pop!(pts.points)

# Edges
DT.contains_edge(e::CustomEdge, Es::CustomEdges) = e ∈ Es.edges

# Boundary nodes
function Base.insert!(seg::CustomPolygonSegment, index, node)
    cur_edge = seg.edges[index-1]
    u, v = DT.initial(cur_edge), DT.terminal(cur_edge)
    seg.edges[index-1] = CustomEdge(u, node)
    insert!(seg.edges, index, CustomEdge(node, v))
    return nothing
end

### Definitions needed for centroidal tessellations
# Points
DT.set_point!(points::CustomPoints, i, x, y) = points.points[i] = CustomPoint(x, y)

# Example
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

edges = CustomEdges(Set{CustomEdge}((CustomEdge(2, 7), CustomEdge(8, 3))))

outer_polygon = CustomPolygon([
    CustomPolygonSegment([CustomEdge(1, 2), CustomEdge(2, 3)]),
    CustomPolygonSegment([CustomEdge(3, 4), CustomEdge(4, 1)]),
])
inner_polygon = CustomPolygon([
    CustomPolygonSegment([CustomEdge(6, 9), CustomEdge(9, 8), CustomEdge(8, 7), CustomEdge(7, 6)]),
])
polygons = CustomPolygons((outer_polygon, inner_polygon))

rng = StableRNG(123)
tri = triangulate(points; edges=edges, boundary_nodes=polygons,
    IntegerType=Int32,
    EdgeType=CustomEdge,
    TriangleType=CustomTriangle,
    EdgesType=CustomEdges,
    TrianglesType=CustomTriangles,
    delete_ghosts=false,
    rng
)

A = get_total_area(tri)
stats = refine!(tri; max_area=1e-3A, min_angle=27.3, rng);

fig, ax, sc = triplot(tri)

tri = triangulate(points;
    IntegerType=Int32,
    EdgeType=CustomEdge,
    TriangleType=CustomTriangle,
    EdgesType=CustomEdges,
    TrianglesType=CustomTriangles,
    delete_ghosts=false,
    rng
)
vorn = voronoi(tri)

clipped_vorn = voronoi(tri, true)

centroidal_vorn = centroidal_smooth(clipped_vorn; rng)

fig = Figure(fontsize=43)
ax1 = Axis(fig[1, 1], title="Tessellation", width=600, height=400)
ax2 = Axis(fig[1, 2], title="Clipped Tessellation", width=600, height=400)
ax3 = Axis(fig[1, 3], title="Centroidal Tessellation", width=600, height=400)
voronoiplot!(ax1, vorn, colormap=:matter, strokewidth=2, markersize=7)
voronoiplot!(ax2, clipped_vorn, colormap=:matter, strokewidth=2, markersize=7)
voronoiplot!(ax3, centroidal_vorn, colormap=:matter, strokewidth=2, markersize=7)
resize_to_layout!(fig)
fig
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

