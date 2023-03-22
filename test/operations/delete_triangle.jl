using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
import SimpleGraphs: relabel, UndirectedGraph
using DataStructures

include("../test_setup.jl")
save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

include("../helper_functions.jl")

tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)

fig = Figure()
ax = Axis(fig[1, 1])
triplot!(ax, tri; recompute_centers=true, show_ghost_edges=true, ghost_edge_extension_factor=2.0)
xlims!(ax, -2, 22)
ylims!(ax, -2, 22)
all_i = Set{NTuple{3,Int64}}()
for (a, b, c) in [("a1", "f", "g"),
    ("a1", "g", "i"),
    ("a1", "i", "f"),
    ("f", "i", "ℓ"),
    ("f", "v", "e"),
    ("f", "w", "v"),
    ("e", "v", "z")]
    i, j, k = index_map[a], index_map[b], index_map[c]
    push!(all_i, (i, j, k))
    DT.delete_triangle!(tri, i, j, k; protect_boundary=true, update_ghost_edges=false)
    @test !DT.contains_triangle(tri, i, j, k)[2]
    @test !DT.contains_triangle(tri, j, k, i)[2]
    @test !DT.contains_triangle(tri, k, i, j)[2]
    @test DT.get_adjacent(tri, i, j) == DT.DefaultAdjacentValue
    @test DT.get_adjacent(tri, j, k) == DT.DefaultAdjacentValue
    @test DT.get_adjacent(tri, k, i) == DT.DefaultAdjacentValue
    @test !DT.contains_edge(i, j, DT.get_adjacent2vertex(tri, k))
    @test !DT.contains_edge(j, k, DT.get_adjacent2vertex(tri, i))
    @test !DT.contains_edge(k, i, DT.get_adjacent2vertex(tri, j))
end
for T in each_triangle(tri)
    i, j, k = T
    if !DT.contains_triangle(i, j, k, all_i)[2]
        @test DT.contains_triangle(tri, i, j, k)[2]
        @test DT.contains_triangle(tri, j, k, i)[2]
        @test DT.contains_triangle(tri, k, i, j)[2]
        @test DT.get_adjacent(tri, i, j) == k
        @test DT.get_adjacent(tri, j, k) == i
        @test DT.get_adjacent(tri, k, i) == j
        @test DT.contains_edge(i, j, DT.get_adjacent2vertex(tri, k))
        @test DT.contains_edge(j, k, DT.get_adjacent2vertex(tri, i))
        @test DT.contains_edge(k, i, DT.get_adjacent2vertex(tri, j))
    end
end

ax = Axis(fig[1, 2])
triplot!(ax, tri; show_ghost_edges=true, ghost_edge_extension_factor=2.0)
xlims!(ax, -2, 22)
ylims!(ax, -2, 22)

SAVE_FIGURE && save("$save_path/test_delete_triangles.png", fig)

### Simpler test 
tri = example_triangulation()
DT.add_triangle!(tri, 1, 3, 7)
DT.add_triangle!(tri, 5, 2, 8)
DT.add_triangle!(tri, 6, 2, 3)

## Deleting an interior triangle 
i, j, k = 1, 3, 7
DT.delete_triangle!(tri, i, j, k)
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (4, 1, 5),
    (6, 3, 1),
    (4, 6, 1),
    (5, 1, 3),
    (5, 2, 8),
    (6, 2, 3)
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (5, 1) => 3, (3, 5) => 1,
        (1, 7) => DT.DefaultAdjacentValue,
        (7, 3) => DT.DefaultAdjacentValue,
        (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
        (2, 3) => 6, (3, 6) => 2, (6, 2) => 3,
        (4, 5) => DT.BoundaryIndex, (5, 8) => DT.BoundaryIndex,
        (8, 2) => DT.BoundaryIndex, (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex
    )
)
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(3, 5), (6, 3), (5, 4), (4, 6)]),
    2 => Set{NTuple{2,Int64}}([(5, 3), (8, 5), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(2, 5), (5, 1), (1, 6), (6, 2)]),
    4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
    5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2), (2, 8)]),
    6 => Set{NTuple{2,Int64}}([(3, 1), (1, 4), (2, 3)]),
    7 => Set{NTuple{2,Int64}}([]),
    8 => Set{NTuple{2,Int64}}([(5, 2)])
)
true_DG = relabel(UndirectedGraph(
        [
            0 0 1 0 1 1 1 0 1
            0 0 0 1 1 1 1 0 0
            1 0 0 1 0 1 1 0 1
            0 1 1 0 0 1 1 0 0
            1 1 0 0 0 1 1 0 0
            1 1 1 1 1 0 0 0 1
            1 1 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 0 0
            1 0 1 0 0 1 0 0 0
        ]
    ), Dict(1:9 .=> [-1, (1:8)...]))
@test get_triangles(tri) == true_T
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
@test (get_graph ∘ get_graph)(tri) == true_DG

## Deleting a triangle with two boundary edges 
for (i, j, k) in ((5, 2, 8), (2, 8, 5), (8, 5, 2))
    local true_T, true_adj, true_adj2v, true_DG
    _tri = deepcopy(tri)
    DT.delete_triangle!(_tri, i, j, k)
    true_T = Set{NTuple{3,Int64}}([
        (3, 2, 5),
        (4, 1, 5),
        (6, 3, 1),
        (4, 6, 1),
        (5, 1, 3),
        (6, 2, 3)
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (3, 5) => 1,
            (2, 3) => 6, (3, 6) => 2, (6, 2) => 3,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex, (2, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(3, 5), (6, 3), (5, 4), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3), (3, 6)]),
        3 => Set{NTuple{2,Int64}}([(2, 5), (5, 1), (1, 6), (6, 2)]),
        4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(3, 1), (1, 4), (2, 3)]),
    )
    true_DG = relabel(UndirectedGraph(
            [
                0 0 1 0 1 1 1
                0 0 0 1 1 1 1
                1 0 0 1 0 1 1
                0 1 1 0 0 1 1
                1 1 0 0 0 1 1
                1 1 1 1 1 0 0
                1 1 1 1 1 0 0
            ]
        ), Dict(1:7 .=> [-1, (1:6)...]))
    DT.clear_empty_features!(_tri)
    @test get_triangles(_tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(_tri) == true_adj
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(_tri) == true_adj2v
    @test (get_graph ∘ get_graph)(_tri) == true_DG
end
i, j, k = 5, 2, 8
DT.delete_triangle!(tri, i, j, k)

## Deleting a triangle with a single boundary edge 
for (i, j, k) in ((6, 2, 3), (2, 3, 6), (3, 6, 2))
    local true_T, true_adj, true_adj2v, true_DG
    _tri = deepcopy(tri)
    DT.delete_triangle!(_tri, i, j, k)
    true_T = Set{NTuple{3,Int64}}([
        (3, 2, 5),
        (4, 1, 5),
        (6, 3, 1),
        (4, 6, 1),
        (5, 1, 3),
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (3, 5) => 1,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex,
            (2, 3) => DT.BoundaryIndex, (3, 6) => DT.BoundaryIndex,
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(3, 5), (6, 3), (5, 4), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3)]),
        3 => Set{NTuple{2,Int64}}([(2, 5), (5, 1), (1, 6)]),
        4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(3, 1), (1, 4)]),
    )
    true_DG = relabel(UndirectedGraph(
            [
                0 0 1 1 1 1 1
                0 0 0 1 1 1 1
                1 0 0 1 0 1 0
                1 1 1 0 0 1 1
                1 1 0 0 0 1 1
                1 1 1 1 1 0 0
                1 1 0 1 1 0 0
            ]
        ), Dict(1:7 .=> [-1, (1:6)...]))
    DT.clear_empty_features!(_tri)
    @test get_triangles(_tri) == true_T
    @test (get_adjacent ∘ get_adjacent)(_tri) == true_adj
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(_tri) == true_adj2v
    @test (get_graph ∘ get_graph)(_tri) == true_DG
end

### Deleting the only triangle of a triangulation
tri = example_empty_triangulation()
DT.add_triangle!(tri, 1, 2, 3)
true_T = Set{NTuple{3,Int64}}([])
true_adj = DefaultDict(DT.DefaultAdjacentValue, Dict())
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}(),
    1 => Set{NTuple{2,Int64}}(),
    2 => Set{NTuple{2,Int64}}(),
    3 => Set{NTuple{2,Int64}}()
)
true_DG = relabel(UndirectedGraph([0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]), Dict(1:4 .=> [-1, 1, 2, 3]))
DT.delete_triangle!(tri, 1, 2, 3)
@test get_triangles(tri) == true_T
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
@test (get_graph ∘ get_graph)(tri) == true_DG

### Testing all the boundary deletion cases
p1 = @SVector[0.0, 0.0]
p2 = @SVector[1.0, 0.0]
p3 = @SVector[0.0, 1.0]
pts = [p1, p2, p3]
tri = Triangulation(pts)
DT.add_triangle!(tri, 1, 2, 3; update_ghost_edges=true)
p4 = @SVector[1.7, 1.7]
push!(pts, p4)
DT.add_triangle!(tri, 3, 2, 4; update_ghost_edges=true)
p5 = @SVector[1.0, 3.0]
p6 = @SVector[3.0, 1.0]
push!(pts, p5, p6)
DT.add_triangle!(tri, 3, 4, 5; update_ghost_edges=true)
DT.add_triangle!(tri, 4, 2, 6; update_ghost_edges=true)
DT.add_triangle!(tri, 5, 4, 6; update_ghost_edges=true)
DT.delete_triangle!(tri, 5, 4, 6; update_ghost_edges=true)
true_T = Set{NTuple{3,Int64}}([
    (1, 2, 3),
    (3, 2, 4),
    (3, 4, 5),
    (4, 2, 6),
    (2, 1, DT.BoundaryIndex),
    (1, 3, DT.BoundaryIndex),
    (3, 5, DT.BoundaryIndex),
    (5, 4, DT.BoundaryIndex),
    (4, 6, DT.BoundaryIndex),
    (6, 2, DT.BoundaryIndex)
])
true_adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
        (3, 2) => 4, (2, 4) => 3, (4, 3) => 2,
        (3, 4) => 5, (4, 5) => 3, (5, 3) => 4,
        (4, 2) => 6, (2, 6) => 4, (6, 4) => 2,
        (2, 1) => DT.BoundaryIndex, (1, DT.BoundaryIndex) => 2, (DT.BoundaryIndex, 2) => 1,
        (1, 3) => DT.BoundaryIndex, (3, DT.BoundaryIndex) => 1, (DT.BoundaryIndex, 1) => 3,
        (3, 5) => DT.BoundaryIndex, (5, DT.BoundaryIndex) => 3, (DT.BoundaryIndex, 3) => 5,
        (5, 4) => DT.BoundaryIndex, (4, DT.BoundaryIndex) => 5, (DT.BoundaryIndex, 5) => 4,
        (4, 6) => DT.BoundaryIndex, (6, DT.BoundaryIndex) => 4, (DT.BoundaryIndex, 4) => 6,
        (6, 2) => DT.BoundaryIndex, (2, DT.BoundaryIndex) => 6, (DT.BoundaryIndex, 6) => 2
    )))
true_adj2v = DT.Adjacent2Vertex(
    Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 4), (4, 6), (6, 2), (2, 1)]),
        1 => Set{NTuple{2,Int64}}([(2, 3), (3, DT.BoundaryIndex), (DT.BoundaryIndex, 2)]),
        2 => Set{NTuple{2,Int64}}([(3, 1), (4, 3), (6, 4), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 6)]),
        3 => Set{NTuple{2,Int64}}([(1, 2), (2, 4), (4, 5), (DT.BoundaryIndex, 1), (5, DT.BoundaryIndex)]),
        4 => Set{NTuple{2,Int64}}([(3, 2), (5, 3), (2, 6), (DT.BoundaryIndex, 5), (6, DT.BoundaryIndex)]),
        5 => Set{NTuple{2,Int64}}([(3, 4), (DT.BoundaryIndex, 3), (4, DT.BoundaryIndex)]),
        6 => Set{NTuple{2,Int64}}([(4, 2), (DT.BoundaryIndex, 4), (2, DT.BoundaryIndex)])
    )
)
true_DG = DT.Graph{Int64}()
DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 3, 5, 4, 6, 2]...)
DT.add_neighbour!(true_DG, 1, [2, 3]...)
DT.add_neighbour!(true_DG, 2, [1, 3, 4, 6]...)
DT.add_neighbour!(true_DG, 3, [1, 2, 4, 5]...)
DT.add_neighbour!(true_DG, 4, [2, 3, 5, 6]...)
DT.add_neighbour!(true_DG, 5, [3, 4]...)
DT.add_neighbour!(true_DG, 6, [2, 4]...)
@test DT.compare_triangle_collections(get_triangles(tri), true_T)
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
@test (get_graph ∘ get_graph)(tri) == true_DG.graph
DT.delete_triangle!(tri, 4, 2, 6; update_ghost_edges=true)
DT.delete_triangle!(tri, 3, 4, 5; update_ghost_edges=true)
true_T = Set{NTuple{3,Int64}}([(1, 2, 3),
    (2, 1, DT.BoundaryIndex),
    (1, 3, DT.BoundaryIndex),
    (3, 4, DT.BoundaryIndex),
    (4, 2, DT.BoundaryIndex),
    (3, 2, 4)])
true_adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
        (2, 1) => DT.BoundaryIndex, (1, DT.BoundaryIndex) => 2, (DT.BoundaryIndex, 2) => 1,
        (1, 3) => DT.BoundaryIndex, (3, DT.BoundaryIndex) => 1, (DT.BoundaryIndex, 1) => 3,
        (3, 4) => DT.BoundaryIndex, (4, DT.BoundaryIndex) => 3, (DT.BoundaryIndex, 3) => 4,
        (4, 2) => DT.BoundaryIndex, (2, DT.BoundaryIndex) => 4, (DT.BoundaryIndex, 4) => 2,
        (3, 2) => 4, (2, 4) => 3, (4, 3) => 2
    )))
true_adj2v = DT.Adjacent2Vertex(
    Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(2, 1), (1, 3), (3, 4), (4, 2)]),
        1 => Set{NTuple{2,Int64}}([(2, 3), (DT.BoundaryIndex, 2), (3, DT.BoundaryIndex)]),
        2 => Set{NTuple{2,Int64}}([(3, 1), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 4), (4, 3)]),
        3 => Set{NTuple{2,Int64}}([(1, 2), (DT.BoundaryIndex, 1), (4, DT.BoundaryIndex), (2, 4)]),
        4 => Set{NTuple{2,Int64}}([(DT.BoundaryIndex, 3), (2, DT.BoundaryIndex), (3, 2)])
    )
)
true_DG = DT.Graph{Int64}()
DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 3, 4, 2]...)
DT.add_neighbour!(true_DG, 1, [2, 3]...)
DT.add_neighbour!(true_DG, 2, [1, 3, 4]...)
DT.add_neighbour!(true_DG, 3, [1, 2, 4]...)
DT.add_neighbour!(true_DG, 4, [2, 3]...)
DT.clear_empty_features!(tri)
@test DT.compare_triangle_collections(get_triangles(tri), true_T)
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
@test (get_graph ∘ get_graph)(tri) == true_DG.graph
DT.delete_triangle!(tri, 3, 2, 4; update_ghost_edges=true)
true_T = Set{NTuple{3,Int64}}([(1, 2, 3),
    (2, 1, DT.BoundaryIndex),
    (1, 3, DT.BoundaryIndex),
    (3, 2, DT.BoundaryIndex)])
true_adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
        (2, 1) => DT.BoundaryIndex, (1, DT.BoundaryIndex) => 2, (DT.BoundaryIndex, 2) => 1,
        (1, 3) => DT.BoundaryIndex, (3, DT.BoundaryIndex) => 1, (DT.BoundaryIndex, 1) => 3,
        (3, 2) => DT.BoundaryIndex, (2, DT.BoundaryIndex) => 3, (DT.BoundaryIndex, 3) => 2
    )))
true_adj2v = DT.Adjacent2Vertex(
    Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(2, 1), (1, 3), (3, 2)]),
        1 => Set{NTuple{2,Int64}}([(2, 3), (DT.BoundaryIndex, 2), (3, DT.BoundaryIndex)]),
        2 => Set{NTuple{2,Int64}}([(3, 1), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 3)]),
        3 => Set{NTuple{2,Int64}}([(1, 2), (DT.BoundaryIndex, 1), (2, DT.BoundaryIndex)])
    )
)
true_DG = DT.Graph{Int64}()
DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 2, 3]...)
DT.add_neighbour!(true_DG, 1, [2, 3]...)
DT.add_neighbour!(true_DG, 2, [1, 3]...)
DT.add_neighbour!(true_DG, 3, [1, 2]...)
DT.clear_empty_features!(tri)
@test DT.compare_triangle_collections(get_triangles(tri), true_T)
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
@test (get_graph ∘ get_graph)(tri) == true_DG.graph
DT.delete_triangle!(tri, 1, 2, 3; update_ghost_edges=true)
true_T = Set{NTuple{3,Int64}}()
true_adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
true_adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
true_DG = DT.Graph{Int64}()
DT.clear_empty_features!(tri)
@test DT.compare_triangle_collections(get_triangles(tri), true_T)
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
@test (get_graph ∘ get_graph)(tri) == true_DG.graph

### A larger example
p1 = @SVector[-3.32, 3.53]
p2 = @SVector[-5.98, 2.17]
p3 = @SVector[-6.36, -1.55]
p4 = @SVector[-2.26, -4.31]
p5 = @SVector[6.34, -3.23]
p6 = @SVector[-3.24, 1.01]
p7 = @SVector[0.14, -1.51]
p8 = @SVector[0.2, 1.25]
p9 = @SVector[1.0, 4.0]
p10 = @SVector[4.74, 2.21]
p11 = @SVector[2.32, -0.27]
pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
tri = Triangulation(pts)
for (i, j, k) in (
    (1, 2, 6),
    (1, 6, 8),
    (9, 1, 8),
    (9, 8, 10),
    (10, 8, 11),
    (8, 7, 11),
    (8, 6, 7),
    (6, 2, 3),
    (6, 3, 4),
    (6, 4, 7),
    (7, 4, 5),
    (11, 7, 5),
    (10, 11, 5)
)
    DT.add_triangle!(tri, i, j, k; update_ghost_edges=true)
end
[DT.delete_triangle!(tri, i, j, k; update_ghost_edges=true) for (i, j, k) in (
    (1, 8, 9), (9, 8, 10), (10, 11, 5), (5, 7, 4), (4, 6, 3), (3, 6, 2), (2, 6, 1)
)]
true_T = Set{NTuple{3,Int64}}([
    (1, 6, 8),
    (6, 7, 8),
    (6, 4, 7),
    (8, 7, 11),
    (11, 7, 5),
    (10, 8, 11),
    (1, 8, DT.BoundaryIndex),
    (8, 10, DT.BoundaryIndex),
    (10, 11, DT.BoundaryIndex),
    (11, 5, DT.BoundaryIndex),
    (5, 7, DT.BoundaryIndex),
    (7, 4, DT.BoundaryIndex),
    (4, 6, DT.BoundaryIndex),
    (6, 1, DT.BoundaryIndex)
])
true_adj = DT.Adjacent(
    DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (1, 6) => 8, (6, 8) => 1, (8, 1) => 6,
            (6, 7) => 8, (7, 8) => 6, (8, 6) => 7,
            (6, 4) => 7, (4, 7) => 6, (7, 6) => 4,
            (8, 7) => 11, (7, 11) => 8, (11, 8) => 7,
            (11, 7) => 5, (7, 5) => 11, (5, 11) => 7,
            (10, 8) => 11, (8, 11) => 10, (11, 10) => 8,
            (1, 8) => DT.BoundaryIndex, (8, DT.BoundaryIndex) => 1, (DT.BoundaryIndex, 1) => 8,
            (8, 10) => DT.BoundaryIndex, (10, DT.BoundaryIndex) => 8, (DT.BoundaryIndex, 8) => 10,
            (10, 11) => DT.BoundaryIndex, (11, DT.BoundaryIndex) => 10, (DT.BoundaryIndex, 10) => 11,
            (11, 5) => DT.BoundaryIndex, (5, DT.BoundaryIndex) => 11, (DT.BoundaryIndex, 11) => 5,
            (5, 7) => DT.BoundaryIndex, (7, DT.BoundaryIndex) => 5, (DT.BoundaryIndex, 5) => 7,
            (7, 4) => DT.BoundaryIndex, (4, DT.BoundaryIndex) => 7, (DT.BoundaryIndex, 7) => 4,
            (4, 6) => DT.BoundaryIndex, (6, DT.BoundaryIndex) => 4, (DT.BoundaryIndex, 4) => 6,
            (6, 1) => DT.BoundaryIndex, (1, DT.BoundaryIndex) => 6, (DT.BoundaryIndex, 6) => 1
        ))
)
true_adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
for (ij, k) in true_adj
    DT.add_adjacent2vertex!(true_adj2v, k, ij)
end
true_DG = DT.Graph{Int64}()
DT.add_neighbour!(true_DG, DT.BoundaryIndex, 1, 8, 10, 11, 5, 7, 4, 6, 1)
DT.add_neighbour!(true_DG, 1, 6, 8)
DT.add_neighbour!(true_DG, 4, 6, 7)
DT.add_neighbour!(true_DG, 5, 7, 11)
DT.add_neighbour!(true_DG, 6, 1, 8, 7, 4)
DT.add_neighbour!(true_DG, 7, 4, 6, 8, 11, 5)
DT.add_neighbour!(true_DG, 8, 1, 6, 7, 11, 10)
DT.add_neighbour!(true_DG, 10, 8, 11)
DT.add_neighbour!(true_DG, 11, 10, 8, 7, 5)
DT.clear_empty_features!(tri)
@test DT.compare_triangle_collections(get_triangles(tri), true_T)
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj.adjacent
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v.adjacent2vertex
@test (get_graph ∘ get_graph)(tri) == true_DG.graph
