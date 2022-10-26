@testset "Deleting triangles" begin
    pts, T, DG, adj, adj2v = example_triangulation()
    p1, p2, p3, p4, p5, p6 = pts
    p7 = @SVector[2.0, 1.0]
    DT.add_point!(pts, p7)
    p8 = @SVector[5.0, 1.0]
    DT.add_point!(pts, p8)
    DT.add_triangle!(1, 3, 7, T, adj, adj2v, DG)
    DT.add_triangle!(5, 2, 8, T, adj, adj2v, DG)
    DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)

    ## Deleting an interior triangle 
    i, j, k = 1, 3, 7
    DT.delete_triangle!(i, j, k, T, adj, adj2v, DG)
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
        ), Dict(1:9 .=> 0:8))
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG

    ## Deleting a triangle with two boundary edges 
    for (i, j, k) in ((5, 2, 8), (2, 8, 5), (8, 5, 2))
        Tc, adjc, adj2vc, DGc = deepcopy(T), deepcopy(adj), deepcopy(adj2v), deepcopy(DG)
        DT.delete_triangle!(i, j, k, Tc, adjc, adj2vc, DGc)
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
                (1, 7) => DT.DefaultAdjacentValue,
                (7, 3) => DT.DefaultAdjacentValue,
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
            7 => Set{NTuple{2,Int64}}([]),
            8 => Set{NTuple{2,Int64}}([])
        )
        true_DG = relabel(UndirectedGraph(
                [
                    0 0 1 0 1 1 1 0 0
                    0 0 0 1 1 1 1 0 0
                    1 0 0 1 0 1 1 0 0
                    0 1 1 0 0 1 1 0 0
                    1 1 0 0 0 1 1 0 0
                    1 1 1 1 1 0 0 0 0
                    1 1 1 1 1 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                ]
            ), Dict(1:9 .=> 0:8))
        @test Tc == true_T
        @test adjacent(adjc) == true_adj
        @test adjacent2vertex(adj2vc) == true_adj2v
        @test graph(DGc) == true_DG
    end
    i, j, k = 5, 2, 8
    DT.delete_triangle!(i, j, k, T, adj, adj2v, DG)

    ## Deleting a triangle with a single boundary edge 
    for (i, j, k) in ((6, 2, 3), (2, 3, 6), (3, 6, 2))
        Tc, adjc, adj2vc, DGc = deepcopy(T), deepcopy(adj), deepcopy(adj2v), deepcopy(DG)
        DT.delete_triangle!(i, j, k, Tc, adjc, adj2vc, DGc)
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
                (1, 7) => DT.DefaultAdjacentValue,
                (7, 3) => DT.DefaultAdjacentValue,
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
            7 => Set{NTuple{2,Int64}}([]),
            8 => Set{NTuple{2,Int64}}([])
        )
        true_DG = relabel(UndirectedGraph(
                [
                    0 0 1 1 1 1 1 0 0
                    0 0 0 1 1 1 1 0 0
                    1 0 0 1 0 1 0 0 0
                    1 1 1 0 0 1 1 0 0
                    1 1 0 0 0 1 1 0 0
                    1 1 1 1 1 0 0 0 0
                    1 1 0 1 1 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                ]
            ), Dict(1:9 .=> 0:8))
        @test Tc == true_T
        @test adjacent(adjc) == true_adj
        @test adjacent2vertex(adj2vc) == true_adj2v
        @test graph(DGc) == true_DG
    end
    i, j, k = 6, 2, 3
    DT.delete_triangle!(i, j, k, T, adj, adj2v, DG)
end

@testset "Deleting the only triangle of a triangulation" begin
    pts, T, DG, adj, adj2v = example_empty_triangulation()
    p1, p2, p3 = pts
    DT.add_triangle!(1, 2, 3, T, adj, adj2v, DG)
    true_T = Set{NTuple{3,Int64}}([])
    true_adj = DefaultDict(DT.DefaultAdjacentValue, Dict())
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}(),
        1 => Set{NTuple{2,Int64}}(),
        2 => Set{NTuple{2,Int64}}(),
        3 => Set{NTuple{2,Int64}}()
    )
    true_DG = relabel(UndirectedGraph([0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]), Dict(1:4 .=> 0:3))
    DT.delete_triangle!(1, 2, 3, T, adj, adj2v, DG)
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
end