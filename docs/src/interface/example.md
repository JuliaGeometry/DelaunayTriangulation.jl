```@meta
CurrentModule = DelaunayTriangulation
```

# Defining a Custom Interface

In the code below, we given an example of how we can provide a fully customised interface for constructing an unconstrained triangulation. One important feature to note in this custom definition is that, rather than writing boundary nodes a sequence of integers defining the polygon vertices, we write them as a sequence of edges built up by individual polygons themselves built up of segments. 

```julia 
using DelaunayTriangulation, CairoMakie, StableRNGs 
const DT = DelaunayTriangulation

## Struct definitions
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

### Defining the methods
## Definitions needed for unconstrained triangulations
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

## Definitions needed for constrained segments 
# Edge 
DT.integer_type(::Type{CustomEdge}) = Int32

# Edges 
DT.edge_type(::Type{CustomEdges}) = CustomEdge

## Definitions needed for boundary nodes 
# Triangles 
DT.each_triangle(tri::CustomTriangles) = tri.triangles

# Boundary Nodes 
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

## Definitions needed for refinement 
# Points 
DT.push_point!(pts::CustomPoints, x, y) = push!(pts.points, CustomPoint(x, y))
DT.pop_point!(pts::CustomPoints) = pop!(pts.points)

# Edges 
DT.contains_edge(e::CustomEdge, Es::CustomEdges) = e ∈ Es.edges

# Boundary Nodes 
function Base.insert!(seg::CustomPolygonSegment, index, node)
    cur_edge = seg.edges[index-1]
    u, v = DT.initial(cur_edge), DT.terminal(cur_edge)
    seg.edges[index-1] = CustomEdge(u, node)
    insert!(seg.edges, index, CustomEdge(node, v))
    return nothing
end

### Example
## Build the points 
p1 = CustomPoint(0.0, 0.0)
p2 = CustomPoint(1.0, 0.0)
p3 = CustomPoint(1.0, 1.0)
p4 = CustomPoint(0.0, 1.0)
p5 = CustomPoint(0.5, 0.5)
p6 = CustomPoint(0.25, 0.25)
p7 = CustomPoint(0.75, 0.25)
p8 = CustomPoint(0.75, 0.75)
p9 = CustomPoint(0.25, 0.75)
p10 = CustomPoint(2.3, 5.5)
p11 = CustomPoint(-0.5, 2.3)
points = CustomPoints([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11])

## The edges
edges = CustomEdges(Set{CustomEdge}((CustomEdge(2, 7), CustomEdge(8, 3))))

## The boundary nodes
outer_polygon = CustomPolygon([
    CustomPolygonSegment([CustomEdge(1, 2), CustomEdge(2, 3)]),
    CustomPolygonSegment([CustomEdge(3, 4), CustomEdge(4, 1)]),
])
inner_polygon = CustomPolygon([
    CustomPolygonSegment([CustomEdge(6, 9), CustomEdge(9, 8), CustomEdge(8, 7), CustomEdge(7, 6)]),
])
polygons = CustomPolygons((outer_polygon, inner_polygon))

## Triangulate
rng = StableRNG(982381238)
tri = triangulate(points; edges=edges, boundary_nodes=polygons,
    IntegerType=Int32,
    EdgeType=CustomEdge,
    TriangleType=CustomTriangle,
    EdgesType=CustomEdges,
    TrianglesType=CustomTriangles,
    randomise=true,
    delete_ghosts=false,
    rng
)

## Refine 
A = get_total_area(tri)
stats = refine!(tri; max_area=1e-3A, min_angle=27.3, rng)
fig, ax, sc = triplot(tri; show_convex_hull=false)
xlims!(ax, 0, 1)
ylims!(ax, 0, 1)
```

```@raw html
<figure>
    <img src='../figs/interface_example.png', alt='Triangulation using a custom interface'><br>
</figure>
```