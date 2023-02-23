using ..DelaunayTriangulation
const DT = DelaunayTriangulation

include("../helper_functions.jl")

x, y = complicated_geometry()
tri = generate_mesh(x, y, 0.1; convert_result=true, add_ghost_triangles=true)
tri_2 = generate_mesh(x[1], y[1], 0.1; convert_result=true, add_ghost_triangles=true)

i, j, k = 58172, 188110, 92887
u, v, w = 101011, 250291, 82277
T = (u, v, w)
DT.add_triangle!(tri, i, j, k)
DT.add_triangle!(tri, T)
for ((a, b, c)) in ((i, j, k), (u, v, w))
    @test DT.get_adjacent(tri, a, b) == c
    @test DT.get_adjacent(tri, b, c) == a
    @test DT.get_adjacent(tri, c, a) == b
    @test (a, b) ∈ DT.get_adjacent2vertex(tri, c)
    @test (b, c) ∈ DT.get_adjacent2vertex(tri, a)
    @test (c, a) ∈ DT.get_adjacent2vertex(tri, b)
    @test a ∈ DT.get_neighbours(tri, b) && a ∈ DT.get_neighbours(tri, c)
    @test b ∈ DT.get_neighbours(tri, a) && b ∈ DT.get_neighbours(tri, c)
    @test c ∈ DT.get_neighbours(tri, a) && c ∈ DT.get_neighbours(tri, b)
    @test (a, b, c) ∈ tri.triangles
end

pts = rand(2, 500)
u, v, w = 37, 38, 73
pts[:, u] .= 0.0
pts[:, v] .= [1.0, 0.0]
pts[:, w] .= [0.0, 1.0]
tri = Triangulation(pts)
DT.add_triangle!(tri, (u, v, w); update_ghost_edges=true)
BI = DT.BoundaryIndex
@test get_triangles(tri) == Set(((u, v, w), (w, v, BI), (v, u, BI), (u, w, BI)))
@test get_adjacent(get_adjacent(tri)).d.d == Dict((u, v) => w,
                                                  (v, w) => u,
                                                  (w, u) => v,
                                                  (w, v) => BI,
                                                  (v, BI) => w,
                                                  (BI, w) => v,
                                                  (v, u) => BI,
                                                  (u, BI) => v,
                                                  (BI, v) => u,
                                                  (u, w) => BI,
                                                  (w, BI) => u,
                                                  (BI, u) => w)
@test get_adjacent2vertex(get_adjacent2vertex(tri)) ==
      Dict(BI => Set{NTuple{2,Int64}}([(w, v), (v, u), (u, w)]),
           u => Set{NTuple{2,Int64}}([(v, w), (BI, v), (w, BI)]),
           v => Set{NTuple{2,Int64}}([(w, u), (BI, w), (u, BI)]),
           w => Set{NTuple{2,Int64}}([(u, v), (v, BI), (BI, u)]))
@test get_graph(get_graph(tri)).E ==
      Set([(BI, v), (u, w), (u, v), (BI, w), (BI, u), (v, w)])
@test get_graph(get_graph(tri)).V == Set([u, v, w, BI])
