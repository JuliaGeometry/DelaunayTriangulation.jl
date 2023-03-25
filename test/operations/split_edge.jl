using ..DelaunayTriangulation
const DT = DelaunayTriangulation
import SimpleGraphs: relabel, UndirectedGraph
using DataStructures

include("../helper_functions.jl")

tri = example_triangulation()
DT.add_triangle!(tri, 6, 2, 3)
DT.split_triangle!(tri, 1, 3, 5, 7)
DT.flip_edge!(tri, 1, 5)
p9 = @SVector[2.5, 0.5]
p10 = @SVector[1.5, 2.0]
push!(get_points(tri), p9, p10)

## Interior edge 
i, j, r = 1, 7, 8
DT.split_edge!(tri, i, j, r)
DT.split_edge!(tri, j, i, r)
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (1, 8, 4),
    (7, 8, 3),
    (3, 5, 7),
    (6, 3, 1),
    (4, 6, 1),
    (4, 7, 5),
    (6, 2, 3),
    (8, 7, 4),
    (8, 1, 3),
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (1, 8) => 4, (8, 4) => 1, (4, 1) => 8,
        (7, 8) => 3, (8, 3) => 7, (3, 7) => 8,
        (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (4, 7) => 5, (7, 5) => 4, (5, 4) => 7,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
        (8, 1) => 3, (1, 3) => 8, (3, 8) => 1,
        (4, 5) => DT.BoundaryIndex,
        (5, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex,
    ))
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(4, 6), (6, 3), (3, 8), (8, 4)]),
    2 => Set{NTuple{2,Int64}}([(5, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(6, 2), (2, 5), (5, 7), (7, 8), (8, 1), (1, 6)]),
    4 => Set{NTuple{2,Int64}}([(6, 1), (1, 8), (8, 7), (7, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 7), (7, 3), (3, 2)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
    7 => Set{NTuple{2,Int64}}([(3, 5), (5, 4), (4, 8), (8, 3)]),
    8 => Set{NTuple{2,Int64}}([(1, 3), (3, 7), (7, 4), (4, 1)])
)
true_DG = relabel(UndirectedGraph(
        [
            0 0 1 0 1 1 1 0 0
            0 0 0 1 1 0 1 0 1
            1 0 0 1 0 1 1 0 0
            0 1 1 0 0 1 1 1 1
            1 1 0 0 0 1 1 1 1
            1 0 1 1 1 0 0 1 0
            1 1 1 1 1 0 0 0 0
            0 0 0 1 1 1 0 0 1
            0 1 0 1 1 0 0 1 0
        ]
    ), Dict(1:9 .=> [-1, (1:8)...]))
DT.clear_empty_features!(tri)
@test get_triangles(tri) == true_T
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
@test (get_graph ∘ get_graph)(tri) == true_DG

## Interior edge for a triangle that is on the boundary
i, j, r = 5, 3, 9
DT.split_edge!(tri, i, j, r)
DT.split_edge!(tri, j, i, r)
true_T = Set{NTuple{3,Int64}}([
    (1, 8, 4),
    (7, 8, 3),
    (6, 3, 1),
    (4, 6, 1),
    (4, 7, 5),
    (6, 2, 3),
    (8, 7, 4),
    (8, 1, 3),
    (9, 5, 7),
    (3, 9, 7),
    (9, 3, 2),
    (5, 9, 2)
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (1, 8) => 4, (8, 4) => 1, (4, 1) => 8,
        (7, 8) => 3, (8, 3) => 7, (3, 7) => 8,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (4, 7) => 5, (7, 5) => 4, (5, 4) => 7,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
        (8, 1) => 3, (1, 3) => 8, (3, 8) => 1,
        (7, 9) => 5, (9, 5) => 7, (5, 7) => 9,
        (3, 9) => 7, (9, 7) => 3, (7, 3) => 9,
        (9, 3) => 2, (3, 2) => 9, (2, 9) => 3,
        (5, 9) => 2, (9, 2) => 5, (2, 5) => 9,
        (4, 5) => DT.BoundaryIndex,
        (5, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex,
    ))
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(4, 6), (6, 3), (3, 8), (8, 4)]),
    2 => Set{NTuple{2,Int64}}([(5, 9), (9, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(2, 9), (9, 7), (7, 8), (8, 1), (1, 6), (6, 2)]),
    4 => Set{NTuple{2,Int64}}([(6, 1), (1, 8), (8, 7), (7, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 7), (7, 9), (9, 2)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
    7 => Set{NTuple{2,Int64}}([(9, 5), (5, 4), (4, 8), (8, 3), (3, 9)]),
    8 => Set{NTuple{2,Int64}}([(1, 3), (3, 7), (7, 4), (4, 1)]),
    9 => Set{NTuple{2,Int64}}([(2, 5), (5, 7), (7, 3), (3, 2)])
)
true_DG = relabel(UndirectedGraph(
        [
            0 0 1 0 1 1 1 0 0 0
            0 0 0 1 1 0 1 0 1 0
            1 0 0 1 0 1 1 0 0 1
            0 1 1 0 0 0 1 1 1 1
            1 1 0 0 0 1 1 1 1 0
            1 0 1 0 1 0 0 1 0 1
            1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 0 0 1 1
            0 1 0 1 1 0 0 1 0 0
            0 0 1 1 0 1 0 1 0 0
        ]
    ), Dict(1:10 .=> [-1, (1:9)...]))
DT.clear_empty_features!(tri)
@test get_triangles(tri) == true_T
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
@test (get_graph ∘ get_graph)(tri) == true_DG

## Boundary edge 
i, j, r = 5, 4, 10
DT.split_edge!(tri, i, j, r)
true_T = Set{NTuple{3,Int64}}([
    (1, 8, 4),
    (7, 8, 3),
    (6, 3, 1),
    (4, 6, 1),
    (6, 2, 3),
    (8, 7, 4),
    (8, 1, 3),
    (9, 5, 7),
    (3, 9, 7),
    (9, 3, 2),
    (5, 9, 2),
    (5, 10, 7),
    (10, 4, 7)
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (1, 8) => 4, (8, 4) => 1, (4, 1) => 8,
        (7, 8) => 3, (8, 3) => 7, (3, 7) => 8,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
        (8, 1) => 3, (1, 3) => 8, (3, 8) => 1,
        (7, 9) => 5, (9, 5) => 7, (5, 7) => 9,
        (3, 9) => 7, (9, 7) => 3, (7, 3) => 9,
        (9, 3) => 2, (3, 2) => 9, (2, 9) => 3,
        (5, 9) => 2, (9, 2) => 5, (2, 5) => 9,
        (5, 10) => 7, (10, 7) => 5, (7, 5) => 10,
        (10, 4) => 7, (4, 7) => 10, (7, 10) => 4,
        (4, 10) => DT.BoundaryIndex,
        (10, 5) => DT.BoundaryIndex,
        (5, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex,
    ))
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 10), (10, 5), (5, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(4, 6), (6, 3), (3, 8), (8, 4)]),
    2 => Set{NTuple{2,Int64}}([(5, 9), (9, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(2, 9), (9, 7), (7, 8), (8, 1), (1, 6), (6, 2)]),
    4 => Set{NTuple{2,Int64}}([(6, 1), (1, 8), (8, 7), (7, 10)]),
    5 => Set{NTuple{2,Int64}}([(10, 7), (7, 9), (9, 2)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
    7 => Set{NTuple{2,Int64}}([(5, 10), (10, 4), (4, 8), (8, 3), (3, 9), (9, 5)]),
    8 => Set{NTuple{2,Int64}}([(1, 3), (3, 7), (7, 4), (4, 1)]),
    9 => Set{NTuple{2,Int64}}([(2, 5), (5, 7), (7, 3), (3, 2)]),
    10 => Set{NTuple{2,Int64}}([(4, 7), (7, 5)])
)
true_DG = relabel(UndirectedGraph(
        [
            0 0 1 0 1 1 1 0 0 0 1
            0 0 0 1 1 0 1 0 1 0 0
            1 0 0 1 0 1 1 0 0 1 0
            0 1 1 0 0 0 1 1 1 1 0
            1 1 0 0 0 0 1 1 1 0 1
            1 0 1 0 0 0 0 1 0 1 1
            1 1 1 1 1 0 0 0 0 0 0
            0 0 0 1 1 1 0 0 1 1 1
            0 1 0 1 1 0 0 1 0 0 0
            0 0 1 1 0 1 0 1 0 0 0
            1 0 0 0 1 1 0 1 0 0 0
        ]
    ), Dict(1:11 .=> [-1, (1:10)...]))
DT.clear_empty_features!(tri)
@test get_triangles(tri) == true_T
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
@test (get_graph ∘ get_graph)(tri) == true_DG