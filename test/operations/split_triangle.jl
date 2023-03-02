using ..DelaunayTriangulation
const DT = DelaunayTriangulation
import SimpleGraphs: relabel, UndirectedGraph
using DataStructures

include("../helper_functions.jl")

tri = example_triangulation()
p9 = @SVector[-1.0, 1.0]
p10 = @SVector[4.0, 1.0]
push!(get_points(tri), p9, p10)
DT.add_triangle!(tri, 5, 2, 8)
DT.add_triangle!(tri, 6, 2, 3)

## Interior splitting 
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (4, 1, 5),
    (5, 2, 8),
    (6, 3, 1),
    (4, 6, 1),
    (1, 3, 7), (3, 5, 7), (5, 1, 7),
    (6, 2, 3),
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
        (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
        (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
        (5, 1) => 7, (1, 7) => 5, (7, 5) => 1,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (4, 5) => DT.BoundaryIndex,
        (5, 8) => DT.BoundaryIndex,
        (8, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex
    )
)
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(5, 4), (4, 6), (6, 3), (3, 7), (7, 5)]),
    2 => Set{NTuple{2,Int64}}([(8, 5), (5, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(5, 7), (7, 1), (1, 6), (6, 2), (2, 5)]),
    4 => Set{NTuple{2,Int64}}([(6, 1), (1, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 1), (1, 7), (7, 3), (3, 2), (2, 8)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
    7 => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 1)]),
    8 => Set{NTuple{2,Int64}}([(5, 2)])
)
true_DG = relabel(UndirectedGraph([
        0 0 1 0 1 1 1 0 1
        0 0 0 1 1 1 1 1 0
        1 0 0 1 0 1 1 0 1
        0 1 1 0 0 1 1 1 0
        1 1 0 0 0 1 1 0 0
        1 1 1 1 1 0 0 1 1
        1 1 1 1 1 0 0 0 0
        0 1 0 1 0 1 0 0 0
        1 0 1 0 0 1 0 0 0
    ]), Dict(1:9 .=> [-1, (1:8)...]))
DT.split_triangle!(tri, 1, 3, 5, 7)
DT.clear_empty_features!(tri)
@test get_triangles(tri) == true_T
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
@test (get_graph ∘ get_graph)(tri) == true_DG

## Splitting a triangle with one boundary edge 
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (4, 1, 5),
    (5, 2, 8),
    (6, 3, 1),
    (4, 6, 9), (6, 1, 9), (1, 4, 9),
    (1, 3, 7), (3, 5, 7), (5, 1, 7),
    (6, 2, 3),
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
        (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
        (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
        (5, 1) => 7, (1, 7) => 5, (7, 5) => 1,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (4, 6) => 9, (6, 9) => 4, (9, 4) => 6,
        (6, 1) => 9, (1, 9) => 6, (9, 6) => 1,
        (9, 1) => 4, (1, 4) => 9, (4, 9) => 1,
        (4, 5) => DT.BoundaryIndex,
        (5, 8) => DT.BoundaryIndex,
        (8, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex
    )
)
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(5, 4), (4, 9), (9, 6), (6, 3), (3, 7), (7, 5)]),
    2 => Set{NTuple{2,Int64}}([(8, 5), (5, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(5, 7), (7, 1), (1, 6), (6, 2), (2, 5)]),
    4 => Set{NTuple{2,Int64}}([(6, 9), (9, 1), (1, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 1), (1, 7), (7, 3), (3, 2), (2, 8)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 9), (9, 4)]),
    7 => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 1)]),
    8 => Set{NTuple{2,Int64}}([(5, 2)]),
    9 => Set{NTuple{2,Int64}}([(6, 1), (1, 4), (4, 6)])
)
true_DG = relabel(UndirectedGraph([
        0 0 1 0 1 1 1 0 1 0
        0 0 0 1 1 1 1 1 0 1
        1 0 0 1 0 1 1 0 1 0
        0 1 1 0 0 1 1 1 0 0
        1 1 0 0 0 1 1 0 0 1
        1 1 1 1 1 0 0 1 1 0
        1 1 1 1 1 0 0 0 0 1
        0 1 0 1 0 1 0 0 0 0
        1 0 1 0 0 1 0 0 0 0
        0 1 0 0 1 0 1 0 0 0
    ]), Dict(1:10 .=> [-1, (1:9)...]))
DT.split_triangle!(tri, 4, 6, 1, 9)
DT.clear_empty_features!(tri)
@test get_triangles(tri) == true_T
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
@test (get_graph ∘ get_graph)(tri) == true_DG

## Splitting two boundary edges 
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (4, 1, 5),
    (5, 2, 10), (2, 8, 10), (8, 5, 10),
    (6, 3, 1),
    (4, 6, 9), (6, 1, 9), (1, 4, 9),
    (1, 3, 7), (3, 5, 7), (5, 1, 7),
    (6, 2, 3),
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
        (5, 2) => 10, (2, 10) => 5, (10, 5) => 2,
        (2, 8) => 10, (8, 10) => 2, (10, 2) => 8,
        (8, 5) => 10, (5, 10) => 8, (10, 8) => 5,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
        (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
        (5, 1) => 7, (1, 7) => 5, (7, 5) => 1,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (4, 6) => 9, (6, 9) => 4, (9, 4) => 6,
        (6, 1) => 9, (1, 9) => 6, (9, 6) => 1,
        (9, 1) => 4, (1, 4) => 9, (4, 9) => 1,
        (4, 5) => DT.BoundaryIndex,
        (5, 8) => DT.BoundaryIndex,
        (8, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex
    )
)
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(5, 4), (4, 9), (9, 6), (6, 3), (3, 7), (7, 5)]),
    2 => Set{NTuple{2,Int64}}([(8, 10), (10, 5), (5, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(5, 7), (7, 1), (1, 6), (6, 2), (2, 5)]),
    4 => Set{NTuple{2,Int64}}([(6, 9), (9, 1), (1, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 1), (1, 7), (7, 3), (3, 2), (2, 10), (10, 8)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 9), (9, 4)]),
    7 => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 1)]),
    8 => Set{NTuple{2,Int64}}([(5, 10), (10, 2)]),
    9 => Set{NTuple{2,Int64}}([(6, 1), (1, 4), (4, 6)]),
    10 => Set{NTuple{2,Int64}}([(5, 2), (2, 8), (8, 5)])
)
true_DG = relabel(UndirectedGraph([
        0 0 1 0 1 1 1 0 1 0 0
        0 0 0 1 1 1 1 1 0 1 0
        1 0 0 1 0 1 1 0 1 0 1
        0 1 1 0 0 1 1 1 0 0 0
        1 1 0 0 0 1 1 0 0 1 0
        1 1 1 1 1 0 0 1 1 0 1
        1 1 1 1 1 0 0 0 0 1 0
        0 1 0 1 0 1 0 0 0 0 0
        1 0 1 0 0 1 0 0 0 0 1
        0 1 0 0 1 0 1 0 0 0 0
        0 0 1 0 0 1 0 0 1 0 0
    ]), Dict(1:11 .=> [-1, (1:10)...]))
DT.split_triangle!(tri, 5, 2, 8, 10)
DT.clear_empty_features!(tri)
@test get_triangles(tri) == true_T
@test (get_adjacent ∘ get_adjacent)(tri) == true_adj
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri) == true_adj2v
@test (get_graph ∘ get_graph)(tri) == true_DG