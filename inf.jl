using DelaunayTriangulation, Random, DispatchDoctor
const DT = DelaunayTriangulation;
struct CustomPoint
    x::Float32
    y::Float32
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
    polygons::NTuple{N,CustomPolygon}
end
@stable begin
    DT.getx(p::CustomPoint) = p.x
    DT.gety(p::CustomPoint) = p.y
    DT.number_type(::Type{CustomPoint}) = Float32
    Base.eachindex(points::CustomPoints) = Base.eachindex(points.points)
    Base.iterate(points::CustomPoints, state...) = Base.iterate(points.points, state...)
    Base.length(points::CustomPoints) = length(points.points)
    Base.getindex(points::CustomPoints, i) = points.points[i]
    DT.number_type(::Type{CustomPoints}) = Float32
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
    DT.get_boundary_nodes(seg::CustomPolygonSegment, m::Integer) = m > length(seg.segments) ? DT.terminal(seg.segments[m-1]) : DT.initial(seg.segments[m]) # go down to the mth edge and extract the left node 
    DT.get_boundary_nodes(poly::CustomPolygons, (m, n)::NTuple{2,Int32}) = DT.get_boundary_nodes(DT.get_boundary_nodes(poly, m), n)
    Base.empty!(segments::CustomSegments) = empty!(segments.segments)
    DT.set_point!(points::CustomPoints, i, x, y) = points.points[i] = CustomPoint(x, y)
    Base.pop!(points::CustomPoints) = pop!(points.points)
    DT.push_point!(points::CustomPoints, x, y) = push!(points.points, CustomPoint(x, y))
    DT.contains_edge(e::CustomSegment, Es::CustomSegments) = e ∈ Es.segments
    Base.empty!(triangles::CustomTriangles) = empty!(triangles.triangles)
    function Base.insert!(seg::CustomPolygonSegment, index, node)
        cur_segment = seg.segments[index-1]
        u, v = edge_vertices(cur_segment)
        seg.segments[index-1] = CustomSegment(u, node)
        insert!(seg.segments, index, CustomSegment(node, v))
        return seg
    end
end

p1 = CustomPoint(0.0f0, 0.0f0)
p2 = CustomPoint(1.0f0, 0.0f0)
p3 = CustomPoint(1.0f0, 1.0f0)
p4 = CustomPoint(0.0f0, 1.0f0)
p5 = CustomPoint(0.5f0, 0.5f0)
p6 = CustomPoint(0.25f0, 0.25f0)
p7 = CustomPoint(0.75f0, 0.25f0)
p8 = CustomPoint(0.75f0, 0.75f0)
p9 = CustomPoint(0.25f0, 0.75f0)
points = CustomPoints([p1, p2, p3, p4, p5, p6, p7, p8, p9])
segments = CustomSegments(Set{CustomSegment}((CustomSegment(2, 7), CustomSegment(8, 3))))
outer_polygon = CustomPolygon([
    CustomPolygonSegment([CustomSegment(1, 2), CustomSegment(2, 3)]),
    CustomPolygonSegment([CustomSegment(3, 4), CustomSegment(4, 1)]),
])
inner_polygon = CustomPolygon([
    CustomPolygonSegment([CustomSegment(6, 9), CustomSegment(9, 8), CustomSegment(8, 7), CustomSegment(7, 6)]),
])
polygons = CustomPolygons((outer_polygon, inner_polygon));

tri = triangulate(points; boundary_nodes=polygons, segments, IntegerType=Int32,
    EdgeType=CustomSegment,
    TriangleType=CustomTriangle,
    EdgesType=CustomSegments,
    TrianglesType=CustomTriangles,
)
vorn = voronoi(tri, clip=true)
cvor = centroidal_smooth(vorn)

refine!(tri; max_area = 1f-3)

statistics(tri)