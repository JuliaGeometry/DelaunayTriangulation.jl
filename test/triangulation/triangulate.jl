using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using Test
using DataStructures
using CairoMakie

include("../helper_functions.jl")

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

@testset "Random tests" begin
    for _ in 1:100
        pts = rand(2, 382)
        tri = triangulate(pts)
        validate_triangulation(tri)
        _tri = DT.triangulate_bowyer_watson(pts)
        @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri)) &&
              get_adjacent(tri) == get_adjacent(_tri) &&
              get_adjacent2vertex(tri) == get_adjacent2vertex(_tri) &&
              get_graph(tri) == get_graph(_tri) &&
              get_convex_hull(tri) == get_convex_hull(_tri)
    end
end

@testset "Lots of collinearity" begin
    _tri = triangulate_rectangle(-3.0, 2.0, 5.0, 17.3, 23, 57; single_boundary=true)
    validate_triangulation(_tri)
    for _ in 1:100
        @time tri = triangulate(_tri.points)
        validate_triangulation(tri)
    end
end

@testset "A detailed example" begin
    @time for _ in 1:500
        T = Set{NTuple{3,Int64}}(((2, 9, 8),
            (2, 8, 6),
            (2, 6, 1),
            (1, 6, 10),
            (6, 8, 10),
            (10, 8, 3),
            (10, 3, 4),
            (10, 4, 5),
            (1, 10, 5),
            (8, 9, 3),
            (9, 7, 3),
            (3, 7, 4),
            (1, 5, DT.BoundaryIndex),
            (5, 4, DT.BoundaryIndex),
            (4, 7, DT.BoundaryIndex),
            (7, 9, DT.BoundaryIndex),
            (9, 2, DT.BoundaryIndex),
            (2, 1, DT.BoundaryIndex)))
        T2 = Set{NTuple{3,Int64}}(((2, 9, 8),
            (2, 8, 6),
            (2, 6, 1),
            (1, 6, 10),
            (6, 8, 10),
            (10, 8, 3),
            (10, 3, 4),
            (10, 4, 5),
            (1, 10, 5),
            (8, 9, 3),
            (9, 7, 3),
            (3, 7, 4)))
        adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
            Dict(@_adj(2, 9, 8)...,
                @_adj(2, 8, 6)...,
                @_adj(2, 6, 1)...,
                @_adj(1, 6, 10)...,
                @_adj(6, 8, 10)...,
                @_adj(10, 8, 3)...,
                @_adj(10, 3, 4)...,
                @_adj(10, 4, 5)...,
                @_adj(1, 10, 5)...,
                @_adj(8, 9, 3)...,
                @_adj(9, 7, 3)...,
                @_adj(3, 7, 4)...,
                @_adj(1, 5, DT.BoundaryIndex)...,
                @_adj(5, 4, DT.BoundaryIndex)...,
                @_adj(4, 7, DT.BoundaryIndex)...,
                @_adj(7, 9, DT.BoundaryIndex)...,
                @_adj(9, 2, DT.BoundaryIndex)...,
                @_adj(2, 1, DT.BoundaryIndex)...)))
        adj2 = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
            Dict(@_adj(2, 9, 8)...,
                @_adj(2, 8, 6)...,
                @_adj(2, 6, 1)...,
                @_adj(1, 6, 10)...,
                @_adj(6, 8, 10)...,
                @_adj(10, 8, 3)...,
                @_adj(10, 3, 4)...,
                @_adj(10, 4, 5)...,
                @_adj(1, 10, 5)...,
                @_adj(8, 9, 3)...,
                @_adj(9, 7, 3)...,
                @_adj(3, 7, 4)...,
                (1, 5) => DT.BoundaryIndex,
                (5, 4) => DT.BoundaryIndex,
                (4, 7) => DT.BoundaryIndex,
                (7, 9) => DT.BoundaryIndex,
                (9, 2) => DT.BoundaryIndex,
                (2, 1) => DT.BoundaryIndex)))
        adj2v = DT.Adjacent2Vertex(Dict(DT.BoundaryIndex => Set{NTuple{2,Int64}}(((1, 5),
                (5, 4),
                (4, 7),
                (7, 9),
                (9, 2),
                (2, 1))),
            1 => Set{NTuple{2,Int64}}(((2, 6), (6, 10),
                (10, 5),
                (5, DT.BoundaryIndex),
                (DT.BoundaryIndex, 2))),
            2 => Set{NTuple{2,Int64}}(((9, 8), (8, 6), (6, 1),
                (1, DT.BoundaryIndex),
                (DT.BoundaryIndex, 9))),
            3 => Set{NTuple{2,Int64}}(((9, 7), (7, 4), (4, 10),
                (10, 8), (8, 9))),
            4 => Set{NTuple{2,Int64}}(((5, 10), (10, 3),
                (3, 7),
                (7, DT.BoundaryIndex),
                (DT.BoundaryIndex, 5))),
            5 => Set{NTuple{2,Int64}}(((1, 10), (10, 4),
                (4, DT.BoundaryIndex),
                (DT.BoundaryIndex, 1))),
            6 => Set{NTuple{2,Int64}}(((2, 8), (8, 10),
                (10, 1), (1, 2))),
            7 => Set{NTuple{2,Int64}}(((4, 3), (3, 9),
                (9, DT.BoundaryIndex),
                (DT.BoundaryIndex, 4))),
            8 => Set{NTuple{2,Int64}}(((3, 10), (10, 6),
                (6, 2), (2, 9), (9, 3))),
            9 => Set{NTuple{2,Int64}}(((7, 3), (3, 8), (8, 2),
                (2, DT.BoundaryIndex),
                (DT.BoundaryIndex, 7))),
            10 => Set{NTuple{2,Int64}}(((8, 3), (3, 4), (4, 5),
                (5, 1), (1, 6),
                (6, 8)))))
        adj2v2 = DT.Adjacent2Vertex(Dict(DT.BoundaryIndex => Set{NTuple{2,Int64}}(((1, 5),
                (5, 4),
                (4, 7),
                (7, 9),
                (9, 2),
                (2, 1))),
            1 => Set{NTuple{2,Int64}}(((2, 6), (6, 10),
                (10, 5))),
            2 => Set{NTuple{2,Int64}}(((9, 8), (8, 6),
                (6, 1))),
            3 => Set{NTuple{2,Int64}}(((9, 7), (7, 4),
                (4, 10), (10, 8),
                (8, 9))),
            4 => Set{NTuple{2,Int64}}(((5, 10), (10, 3),
                (3, 7))),
            5 => Set{NTuple{2,Int64}}(((1, 10), (10, 4))),
            6 => Set{NTuple{2,Int64}}(((2, 8), (8, 10),
                (10, 1), (1, 2))),
            7 => Set{NTuple{2,Int64}}(((4, 3), (3, 9))),
            8 => Set{NTuple{2,Int64}}(((3, 10), (10, 6),
                (6, 2), (2, 9),
                (9, 3))),
            9 => Set{NTuple{2,Int64}}(((7, 3), (3, 8),
                (8, 2))),
            10 => Set{NTuple{2,Int64}}(((8, 3), (3, 4),
                (4, 5), (5, 1),
                (1, 6), (6, 8)))))
        A = zeros(Int64, 10, 10)
        A[1, [2, 6, 10, 5]] .= 1
        A[2, [1, 6, 8, 9]] .= 1
        A[3, [8, 9, 7, 4, 10]] .= 1
        A[4, [5, 10, 3, 7]] .= 1
        A[5, [1, 10, 4]] .= 1
        A[6, [2, 8, 10, 1]] .= 1
        A[7, [4, 3, 9]] .= 1
        A[8, [6, 2, 9, 3, 10]] .= 1
        A[9, [2, 8, 3, 7]] .= 1
        A[10, [1, 6, 8, 3, 4, 5]] .= 1
        B = zeros(Int64, 11, 11)
        B[2:end, 2:end] .= A
        B[1, [1, 5, 4, 7, 9, 2].+1] .= 1
        B[[1, 5, 4, 7, 9, 2].+1, 1] .= 1
        graph = DT.Graph((x -> DT.SimpleGraphs.relabel(x,
            Dict(1 => DT.BoundaryIndex,
                (2:11 .=> 1:10)...)))(DT.SimpleGraphs.UndirectedGraph(B)))
        a = [1.5, 4.0]
        b = [0.0, 3.5]
        c = [2.0, 1.5]
        d = [3.0, 2.5]
        e = [2.5, 3.5]
        f = [0.5, 3.0]
        g = [2.5, -2.0]
        h = [0.5, 1.5]
        i = [0.0, 0.5]
        j = [1.5, 3.0]
        pts = [a, b, c, d, e, f, g, h, i, j]
        ch = DT.ConvexHull(pts, [1, 2, 9, 7, 4, 5, 1])
        tri = triangulate(pts; delete_ghosts=false)
        @test DT.compare_triangle_collections(T, get_triangles(tri)) &&
              get_adjacent(tri) == adj &&
              get_adjacent2vertex(tri) == adj2v &&
              get_graph(tri) == graph &&
              get_convex_hull(tri) == ch
        validate_triangulation(tri)
        delete_ghost_triangles!(tri)
        @test DT.compare_triangle_collections(T2, get_triangles(tri)) &&
              get_adjacent(tri) == adj2 &&
              get_adjacent2vertex(tri) == adj2v2 &&
              get_graph(tri) == graph &&
              get_convex_hull(tri) == ch
        validate_triangulation(tri)
    end
end

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
save("$save_path/custom_interface_testing.png", fig)