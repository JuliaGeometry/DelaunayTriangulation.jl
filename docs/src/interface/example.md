```@meta
CurrentModule = DelaunayTriangulation
```

In the code below, we given an example of how we can provide a fully customised interface for cosntructing an unconstrained triangulation. 

```julia 
using DelaunayTriangulation
using CairoMakie 

## Define the new types
DT = DelaunayTriangulation
struct CustomPoint
    x::Float64
    y::Float64
end
DT.getx(p::CustomPoint) = p.x
DT.gety(p::CustomPoint) = p.y

struct CustomPoints
    pts::Vector{CustomPoint}
end
DT.each_point_index(pts::CustomPoints) = eachindex(pts.pts)
DT.num_points(pts::CustomPoints) = length(pts.pts)
DT.each_point(pts::CustomPoints) = pts.pts
DT.getpoint(pts::CustomPoints, n::Integer) = pts.pts[n]
DT.number_type(::Type{CustomPoints}) = Float64

struct CustomEdge
    i::Int32
    j::Int32
end
DT.construct_edge(::Type{CustomEdge}, i, j) = CustomEdge(i, j)
DT.initial(edge::CustomEdge) = edge.i
DT.terminal(edge::CustomEdge) = edge.j

struct CustomEdges
    edges::Vector{CustomEdge}
end
DT.initialise_edges(::Type{CustomEdges}) = CustomEdges(CustomEdge[])
DT.add_to_edges!(edges::CustomEdges, e) = push!(edges.edges, e)
DT.each_edge(edges::CustomEdges) = edges.edges
DT.delete_from_edges!(edges::CustomEdges, e) = deleteat!(edges.edges, findfirst(==(e), edges.edges))
DT.num_edges(edges::CustomEdges) = length(edges.edges)
DT.is_empty(edges::CustomEdges) = isempty(edges.edges)

struct CustomTriangle
    i::Int32
    j::Int32
    k::Int32
end
DT.geti(T::CustomTriangle) = T.i
DT.getj(T::CustomTriangle) = T.j
DT.getk(T::CustomTriangle) = T.k
DT.construct_triangle(::Type{CustomTriangle}, i, j, k) = CustomTriangle(i, j, k)
DT.integer_type(::Type{CustomTriangle}) = Int32

struct CustomTriangles
    tris::Set{CustomTriangle}
end
DT.initialise_triangles(::Type{CustomTriangles}) = CustomTriangles(Set{CustomTriangle}())
Base.sizehint!(tris::CustomTriangles, n) = sizehint!(tris.tris, n)
DT.triangle_type(::Type{CustomTriangles}) = CustomTriangle
DT.add_to_triangles!(tris::CustomTriangles, T) = push!(tris.tris, T)
DT.num_triangles(tris::CustomTriangles) = length(tris.tris)
DT.delete_from_triangles!(tris::CustomTriangles, T) = delete!(tris.tris, T)
Base.in(T, V::CustomTriangles) = T ∈ V.tris
DT.each_triangle(tris::CustomTriangles) = tris.tris

## Triangulate
x = rand(100)
y = rand(100)
pts = CustomPoints(CustomPoint.(x, y))
tri1 = triangulate(pts;
    IntegerType=Int32,
    EdgeType=CustomEdge,
    EdgesType=CustomEdges,
    TriangleType=CustomTriangle,
    TrianglesType=CustomTriangles)
tri2 = triangulate(pts)

fig = Figure(fontsize=24)
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y",aspect=1,width=300,height=300)
triplot!(ax, tri1; show_convex_hull=true, convex_hull_linewidth=4)
ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y",aspect=1,width=300,height=300)
triplot!(ax, tri2; show_convex_hull=true, convex_hull_linewidth=4)
resize_to_layout!(fig)
```

![A triangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/example_image.png?raw=true)