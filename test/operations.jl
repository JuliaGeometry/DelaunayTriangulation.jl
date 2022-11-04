@testset "Adding a triangle" begin
    pts, T, DG, adj, adj2v = example_triangulation()
    p1, p2, p3, p4, p5, p6 = pts
    p7 = @SVector[2.0, 1.0]
    DT.add_point!(pts, p7)
    p8 = @SVector[5.0, 1.0]
    DT.add_point!(pts, p8)

    ## Adding an interior triangle
    i, j, k = 1, 3, 7
    DT.add_triangle!(i, j, k, T, adj, adj2v, DG)
    true_T = Set{NTuple{3,Int64}}([
        (3, 2, 5),
        (1, 3, 7),
        (4, 1, 5),
        (6, 3, 1),
        (4, 6, 1),
        (5, 1, 3)
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (3, 5) => 1,
            (4, 5) => DT.BoundaryIndex, (5, 2) => DT.BoundaryIndex,
            (2, 3) => DT.BoundaryIndex, (3, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(3, 7), (3, 5), (6, 3), (5, 4), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3)]),
        3 => Set{NTuple{2,Int64}}([(2, 5), (5, 1), (7, 1), (1, 6)]),
        4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(3, 1), (1, 4)]),
        7 => Set{NTuple{2,Int64}}([(1, 3)])
    )
    true_DG = UndirectedGraph(
        [
            0 0 1 1 1 1 1 0
            0 0 0 1 1 1 1 1
            1 0 0 1 0 1 0 0
            1 1 1 0 0 1 1 1
            1 1 0 0 0 1 1 0
            1 1 1 1 1 0 0 0
            1 1 0 1 1 0 0 0
            0 1 0 1 0 0 0 0
        ]
    )
    true_DG = relabel(true_DG, Dict(1:8 .=> 0:7))
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG

    ## Adding a boundary triangle with one boundary edge
    for (i, j, k) in ((5, 2, 8), (2, 8, 5), (8, 5, 2))
        Tc, adjc, adj2vc, DGc = deepcopy(T), deepcopy(adj), deepcopy(adj2v), deepcopy(DG)
        DT.add_triangle!(i, j, k, Tc, adjc, adj2vc, DGc)
        true_T = Set{NTuple{3,Int64}}([
            (3, 2, 5),
            (1, 3, 7),
            (4, 1, 5),
            (6, 3, 1),
            (4, 6, 1),
            (5, 1, 3),
            (i, j, k)
        ])
        true_adj = DefaultDict(DT.DefaultAdjacentValue,
            Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
                (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
                (5, 1) => 3, (3, 5) => 1,
                (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
                (4, 5) => DT.BoundaryIndex, (5, 8) => DT.BoundaryIndex,
                (8, 2) => DT.BoundaryIndex, (2, 3) => DT.BoundaryIndex,
                (3, 6) => DT.BoundaryIndex, (6, 4) => DT.BoundaryIndex
            )
        )
        true_adj2v = Dict(
            DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 3), (3, 6), (6, 4)]),
            1 => Set{NTuple{2,Int64}}([(3, 7), (3, 5), (6, 3), (5, 4), (4, 6)]),
            2 => Set{NTuple{2,Int64}}([(5, 3), (8, 5)]),
            3 => Set{NTuple{2,Int64}}([(2, 5), (5, 1), (7, 1), (1, 6)]),
            4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
            5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2), (2, 8)]),
            6 => Set{NTuple{2,Int64}}([(3, 1), (1, 4)]),
            7 => Set{NTuple{2,Int64}}([(1, 3)]),
            8 => Set{NTuple{2,Int64}}([(5, 2)])
        )
        true_DG = relabel(UndirectedGraph(
                [
                    0 0 1 1 1 1 1 0 1
                    0 0 0 1 1 1 1 1 0
                    1 0 0 1 0 1 0 0 1
                    1 1 1 0 0 1 1 1 0
                    1 1 0 0 0 1 1 0 0
                    1 1 1 1 1 0 0 0 1
                    1 1 0 1 1 0 0 0 0
                    0 1 0 1 0 0 0 0 0
                    1 0 1 0 0 1 0 0 0
                ]
            ), Dict(1:9 .=> 0:8))
        @test Tc == true_T
        @test adjacent(adjc) == true_adj
        @test adjacent2vertex(adj2vc) == true_adj2v
        @test graph(DGc) == true_DG
    end
    DT.add_triangle!(5, 2, 8, T, adj, adj2v, DG) # Get an actual copy for later, test is above

    ## Adding a boundary triangle with two boundary edges 
    for (i, j, k) in ((3, 6, 2), (6, 2, 3), (2, 3, 6))
        Tc, adjc, adj2vc, DGc = deepcopy(T), deepcopy(adj), deepcopy(adj2v), deepcopy(DG)
        DT.add_triangle!(i, j, k, Tc, adjc, adj2vc, DGc)
        true_T = Set{NTuple{3,Int64}}([
            (3, 2, 5),
            (1, 3, 7),
            (4, 1, 5),
            (6, 3, 1),
            (4, 6, 1),
            (5, 1, 3),
            (5, 2, 8),
            (i, j, k)
        ])
        true_adj = DefaultDict(DT.DefaultAdjacentValue,
            Dict(
                (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
                (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
                (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
                (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
                (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
                (5, 1) => 3, (3, 5) => 1,
                (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
                (2, 3) => 6, (3, 6) => 2, (6, 2) => 3,
                (4, 5) => DT.BoundaryIndex, (5, 8) => DT.BoundaryIndex,
                (8, 2) => DT.BoundaryIndex, (2, 6) => DT.BoundaryIndex,
                (6, 4) => DT.BoundaryIndex
            )
        )
        true_adj2v = Dict(
            DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
            1 => Set{NTuple{2,Int64}}([(3, 7), (3, 5), (6, 3), (5, 4), (4, 6)]),
            2 => Set{NTuple{2,Int64}}([(5, 3), (8, 5), (3, 6)]),
            3 => Set{NTuple{2,Int64}}([(2, 5), (5, 1), (7, 1), (1, 6), (6, 2)]),
            4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
            5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2), (2, 8)]),
            6 => Set{NTuple{2,Int64}}([(3, 1), (1, 4), (2, 3)]),
            7 => Set{NTuple{2,Int64}}([(1, 3)]),
            8 => Set{NTuple{2,Int64}}([(5, 2)])
        )
        true_DG = relabel(UndirectedGraph(
                [
                    0 0 1 0 1 1 1 0 1
                    0 0 0 1 1 1 1 1 0
                    1 0 0 1 0 1 1 0 1
                    0 1 1 0 0 1 1 1 0
                    1 1 0 0 0 1 1 0 0
                    1 1 1 1 1 0 0 0 1
                    1 1 1 1 1 0 0 0 0
                    0 1 0 1 0 0 0 0 0
                    1 0 1 0 0 1 0 0 0
                ]
            ), Dict(1:9 .=> 0:8))
        @test Tc == true_T
        @test adjacent(adjc) == true_adj
        @test adjacent2vertex(adj2vc) == true_adj2v
        @test graph(DGc) == true_DG
    end
end

@testset "Adding into an empty triangulation" begin
    pts, T, DG, adj, adj2v = example_empty_triangulation()
    p1, p2, p3 = pts
    true_T = Set{NTuple{3,Int64}}([(1, 2, 3)])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict((1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (3, 2) => DT.BoundaryIndex, (2, 1) => DT.BoundaryIndex, (1, 3) => DT.BoundaryIndex))
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(3, 2), (2, 1), (1, 3)]),
        1 => Set{NTuple{2,Int64}}([(2, 3)]),
        2 => Set{NTuple{2,Int64}}([(3, 1)]),
        3 => Set{NTuple{2,Int64}}([(1, 2)])
    )
    true_DG = relabel(UndirectedGraph([
            0 1 1 1
            1 0 1 1
            1 1 0 1
            1 1 1 0
        ]), Dict(1:4 .=> 0:3))
    DT.add_triangle!(1, 2, 3, T, adj, adj2v, DG)
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
end

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

@testset "Flipping an edge" begin
    pts, T, DG, adj, adj2v = example_triangulation()
    hg = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(hg, (5, 1, 3), (6, 3, 1), (1, 4, 8), (4, 1, 7))
    p1, p2, p3, p4, p5, p6 = pts
    p7 = @SVector[0.0, 1.5]
    DT.add_point!(pts, p7)
    p8 = @SVector[-1.0, 1.0]
    DT.add_point!(pts, p8)
    DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)
    DT.split_triangle!(4, 6, 1, 8, T, adj, adj2v, DG)
    DT.split_triangle!(1, 5, 4, 7, T, adj, adj2v, DG)

    ## First flip 
    true_T = Set{NTuple{3,Int64}}([
        (5, 6, 3),
        (3, 2, 5),
        (4, 1, 7),
        (5, 4, 7),
        (5, 1, 6),
        (1, 5, 7),
        (6, 2, 3),
        (6, 1, 8),
        (4, 6, 8),
        (1, 4, 8)
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (5, 6) => 3, (6, 3) => 5, (3, 5) => 6,
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 7, (1, 7) => 4, (7, 4) => 1,
            (5, 4) => 7, (4, 7) => 5, (7, 5) => 4,
            (5, 1) => 6, (1, 6) => 5, (6, 5) => 1,
            (1, 5) => 7, (5, 7) => 1, (7, 1) => 5,
            (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
            (6, 1) => 8, (1, 8) => 6, (8, 6) => 1,
            (4, 6) => 8, (6, 8) => 4, (8, 4) => 6,
            (1, 4) => 8, (4, 8) => 1, (8, 1) => 4,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex,
            (2, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex,
            (3, 1) => DT.DefaultAdjacentValue
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(6, 5), (5, 7), (7, 4), (4, 8), (8, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3), (3, 6)]),
        3 => Set{NTuple{2,Int64}}([(5, 6), (6, 2), (2, 5)]),
        4 => Set{NTuple{2,Int64}}([(6, 8), (8, 1), (1, 7), (7, 5)]),
        5 => Set{NTuple{2,Int64}}([(4, 7), (7, 1), (1, 6), (6, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(2, 3), (3, 5), (5, 1), (1, 8), (8, 4)]),
        7 => Set{NTuple{2,Int64}}([(1, 5), (5, 4), (4, 1)]),
        8 => Set{NTuple{2,Int64}}([(1, 4), (4, 6), (6, 1)])
    )
    true_DG = relabel(UndirectedGraph(
            [
                0 0 1 0 1 1 1 0 0
                0 0 0 0 1 1 1 1 1
                1 0 0 1 0 1 1 0 0
                0 0 1 0 0 1 1 0 0
                1 1 0 0 0 1 1 1 1
                1 1 1 1 1 0 1 1 0
                1 1 1 1 1 1 0 0 1
                0 1 0 0 1 1 0 0 0
                0 1 0 0 1 0 1 0 0
            ]
        ), Dict(1:9 .=> 0:8))
    true_HG = DirectedGraph{NTuple{3,Int64}}()
    add!.(Ref(true_HG), ((5, 1, 3), (6, 3, 1), (5, 1, 6), (5, 6, 3), (1, 4, 8), (4, 1, 7)))
    add!.(Ref(true_HG), Ref((5, 1, 3)), ((5, 1, 6), (5, 6, 3)))
    add!.(Ref(true_HG), Ref((6, 3, 1)), ((5, 1, 6), (5, 6, 3)))
    DT.flip_edge!(1, 3, T, adj, adj2v, DG, hg)
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
    @test graph(hg) == true_HG

    ## Second flip 
    true_T = Set{NTuple{3,Int64}}([
        (5, 6, 3),
        (3, 2, 5),
        (5, 4, 7),
        (5, 1, 6),
        (1, 5, 7),
        (6, 2, 3),
        (6, 1, 8),
        (4, 6, 8),
        (8, 1, 7),
        (8, 7, 4)
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (5, 6) => 3, (6, 3) => 5, (3, 5) => 6,
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (5, 4) => 7, (4, 7) => 5, (7, 5) => 4,
            (5, 1) => 6, (1, 6) => 5, (6, 5) => 1,
            (1, 5) => 7, (5, 7) => 1, (7, 1) => 5,
            (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
            (6, 1) => 8, (1, 8) => 6, (8, 6) => 1,
            (4, 6) => 8, (6, 8) => 4, (8, 4) => 6,
            (8, 1) => 7, (1, 7) => 8, (7, 8) => 1,
            (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex,
            (2, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex,
            (3, 1) => DT.DefaultAdjacentValue,
            (4, 1) => DT.DefaultAdjacentValue
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(6, 5), (5, 7), (7, 8), (8, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3), (3, 6)]),
        3 => Set{NTuple{2,Int64}}([(5, 6), (6, 2), (2, 5)]),
        4 => Set{NTuple{2,Int64}}([(6, 8), (8, 7), (7, 5)]),
        5 => Set{NTuple{2,Int64}}([(4, 7), (7, 1), (1, 6), (6, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(2, 3), (3, 5), (5, 1), (1, 8), (8, 4)]),
        7 => Set{NTuple{2,Int64}}([(1, 5), (5, 4), (4, 8), (8, 1)]),
        8 => Set{NTuple{2,Int64}}([(1, 7), (7, 4), (4, 6), (6, 1)])
    )
    true_DG = relabel(UndirectedGraph(
            [
                0 0 1 0 1 1 1 0 0
                0 0 0 0 0 1 1 1 1
                1 0 0 1 0 1 1 0 0
                0 0 1 0 0 1 1 0 0
                1 0 0 0 0 1 1 1 1
                1 1 1 1 1 0 1 1 0
                1 1 1 1 1 1 0 0 1
                0 1 0 0 1 1 0 0 1
                0 1 0 0 1 0 1 1 0
            ]
        ), Dict(1:9 .=> 0:8))
    add!.(Ref(true_HG), ((1, 4, 8), (4, 1, 7), (8, 1, 7), (8, 7, 4)))
    add!.(Ref(true_HG), Ref((1, 4, 8)), ((8, 1, 7), (8, 7, 4)))
    add!.(Ref(true_HG), Ref((4, 1, 7)), ((8, 1, 7), (8, 7, 4)))
    DT.flip_edge!(1, 4, T, adj, adj2v, DG, hg)
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
    @test graph(hg) == true_HG

    ## Another example 
    pts, T, DG, adj, adj2v = example_triangulation()
    p1, p2, p3, p4, p5, p6 = pts
    p7 = @SVector[2.0, 1.0]
    DT.add_point!(pts, p7)
    DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)
    DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG)
    i, j = 5, 1
    r = 7
    e = DT.get_edge(adj, j, i)
    DT.flip_edge!(i, j, e, r, T, adj, adj2v, DG)
    true_T = Set{NTuple{3,Int64}}([
        (3, 2, 5),
        (1, 3, 7),
        (3, 5, 7),
        (6, 3, 1),
        (4, 6, 1),
        (6, 2, 3),
        (7, 5, 4),
        (7, 4, 1)
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
            (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
            (7, 5) => 4, (5, 4) => 7, (4, 7) => 5,
            (7, 4) => 1, (4, 1) => 7, (1, 7) => 4,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex,
            (2, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex,
            (1, 5) => DT.DefaultAdjacentValue
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(6, 3), (3, 7), (7, 4), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3), (3, 6)]),
        3 => Set{NTuple{2,Int64}}([(2, 5), (5, 7), (7, 1), (1, 6), (6, 2)]),
        4 => Set{NTuple{2,Int64}}([(6, 1), (1, 7), (7, 5)]),
        5 => Set{NTuple{2,Int64}}([(4, 7), (7, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
        7 => Set{NTuple{2,Int64}}([(3, 5), (5, 4), (4, 1), (1, 3)])
    )
    true_DG = relabel(UndirectedGraph([
            0 0 1 0 1 1 1 0
            0 0 0 1 1 0 1 1
            1 0 0 1 0 1 1 0
            0 1 1 0 0 1 1 1
            1 1 0 0 0 1 1 1
            1 0 1 1 1 0 0 1
            1 1 1 1 1 0 0 0
            0 1 0 1 1 1 0 0
        ]), Dict(1:8 .=> 0:7))
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
    @test all(DT.isoriented(T, pts) == 1 for T in T)
end

@testset "Splitting a triangle" begin
    pts, T, DG, adj, adj2v = example_triangulation()
    hg = DT.HistoryGraph{NTuple{3,Int64}}()
    p1, p2, p3, p4, p5, p6 = pts
    p7 = @SVector[2.0, 1.0]
    DT.add_point!(pts, p7)
    p8 = @SVector[5.0, 1.0]
    DT.add_point!(pts, p8)
    p9 = @SVector[-1.0, 1.0]
    DT.add_point!(pts, p9)
    p10 = @SVector[4.0, 1.0]
    DT.add_point!(pts, p10)
    DT.add_triangle!(5, 2, 8, T, adj, adj2v, DG)
    DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)

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
        ]), Dict(1:9 .=> 0:8))
    true_HG = DirectedGraph{NTuple{3,Int64}}()
    [add!(true_HG, (i, j, k)) for (i, j, k) in ((1, 3, 5), (1, 3, 7), (3, 5, 7), (5, 1, 7))]
    [add!(true_HG, (1, 3, 5), (i, j, k)) for (i, j, k) in ((1, 3, 7), (3, 5, 7), (5, 1, 7))]
    DT.add_triangle!(hg, (1, 3, 5))
    DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG, hg)
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
    @test graph(hg) == true_HG

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
        ]), Dict(1:10 .=> 0:9))
    [add!(true_HG, (i, j, k)) for (i, j, k) in ((4, 6, 1), (4, 6, 9), (6, 1, 9), (1, 4, 9))]
    [add!(true_HG, (4, 6, 1), (i, j, k)) for (i, j, k) in ((4, 6, 9), (6, 1, 9), (1, 4, 9))]
    DT.add_triangle!(hg, (4, 6, 1))
    DT.split_triangle!(4, 6, 1, 9, T, adj, adj2v, DG, hg)
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
    @test graph(hg) == true_HG

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
        ]), Dict(1:11 .=> 0:10))
    [add!(true_HG, (i, j, k)) for (i, j, k) in ((5, 2, 8), (5, 2, 10), (2, 8, 10), (8, 5, 10))]
    [add!(true_HG, (5, 2, 8), (i, j, k)) for (i, j, k) in ((5, 2, 10), (2, 8, 10), (8, 5, 10))]
    DT.add_triangle!(hg, (5, 2, 8))
    DT.split_triangle!(5, 2, 8, 10, T, adj, adj2v, DG, hg)
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
    @test graph(hg) == true_HG
end

@testset "Splitting an edge" begin
    pts, T, DG, adj, adj2v = example_triangulation()
    p1, p2, p3, p4, p5, p6 = pts
    p7 = @SVector[2.0, 1.0]
    DT.add_point!(pts, p7)
    DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)
    DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG)
    DT.flip_edge!(1, 5, T, adj, adj2v, DG)
    p8 = @SVector[1.0, 1.0]
    DT.add_point!(pts, p8)
    p9 = @SVector[2.5, 0.5]
    DT.add_point!(pts, p9)
    p10 = @SVector[1.5, 2.0]
    DT.add_point!(pts, p10)

    ## Interior edge 
    i, j, r = 1, 7, 8
    DT.split_edge!(i, j, r, T, adj, adj2v, DG)
    DT.split_edge!(j, i, r, T, adj, adj2v, DG)
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
            (5, 1) => DT.DefaultAdjacentValue,
            (1, 7) => DT.DefaultAdjacentValue
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
        ), Dict(1:9 .=> 0:8))
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG

    ## Interior edge for a triangle that is on the boundary
    i, j, r = 5, 3, 9
    DT.split_edge!(i, j, r, T, adj, adj2v, DG)
    DT.split_edge!(j, i, r, T, adj, adj2v, DG)
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
            (5, 1) => DT.DefaultAdjacentValue,
            (1, 7) => DT.DefaultAdjacentValue,
            (5, 3) => DT.DefaultAdjacentValue
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
        ), Dict(1:10 .=> 0:9))
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG


    ## Boundary edge 
    i, j, r = 5, 4, 10
    DT.split_edge!(i, j, r, T, adj, adj2v, DG)
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
            (5, 1) => DT.DefaultAdjacentValue,
            (1, 7) => DT.DefaultAdjacentValue,
            (5, 3) => DT.DefaultAdjacentValue
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
        ), Dict(1:11 .=> 0:10))
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
end

@testset "Legalising an edge" begin
    pts, T, DG, adj, adj2v = example_triangulation()
    hg = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(hg, T...)
    p1, p2, p3, p4, p5, p6 = pts
    p7 = @SVector[2.0, 1.0]
    DT.add_point!(pts, p7)
    DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)
    DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG)
    @test DT.islegal(1, 3, adj, pts)
    @test DT.islegal(3, 5, adj, pts)
    i, j, r = 5, 1, 7
    e = DT.get_edge(adj, j, i)
    DT.legalise_edge!(i, j, r, T, hg, adj, adj2v, DG, pts)
    true_T = Set{NTuple{3,Int64}}([
        (3, 2, 5),
        (1, 3, 7),
        (3, 5, 7),
        (6, 3, 1),
        (4, 6, 1),
        (6, 2, 3),
        (7, 5, 4),
        (7, 4, 1)
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
            (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
            (7, 5) => 4, (5, 4) => 7, (4, 7) => 5,
            (7, 4) => 1, (4, 1) => 7, (1, 7) => 4,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex,
            (2, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex,
            (1, 5) => DT.DefaultAdjacentValue
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(6, 3), (3, 7), (7, 4), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3), (3, 6)]),
        3 => Set{NTuple{2,Int64}}([(2, 5), (5, 7), (7, 1), (1, 6), (6, 2)]),
        4 => Set{NTuple{2,Int64}}([(6, 1), (1, 7), (7, 5)]),
        5 => Set{NTuple{2,Int64}}([(4, 7), (7, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
        7 => Set{NTuple{2,Int64}}([(3, 5), (5, 4), (4, 1), (1, 3)])
    )
    true_DG = relabel(UndirectedGraph([
            0 0 1 0 1 1 1 0
            0 0 0 1 1 0 1 1
            1 0 0 1 0 1 1 0
            0 1 1 0 0 1 1 1
            1 1 0 0 0 1 1 1
            1 0 1 1 1 0 0 1
            1 1 1 1 1 0 0 0
            0 1 0 1 1 1 0 0
        ]), Dict(1:8 .=> 0:7))
    @test T == true_T
    @test adjacent(adj) == true_adj
    @test adjacent2vertex(adj2v) == true_adj2v
    @test graph(DG) == true_DG
    @test all(DT.isoriented(T, pts) == 1 for T in T)
    for ((i, j), v) in adjacent(adj)
        if v ≠ DT.DefaultAdjacentValue
            @test DT.islegal(i, j, adj, pts)
        end
    end
end

@testset "Testing all the boundary addition cases" begin
    T = Set{NTuple{3,Int64}}()
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DG = DT.DelaunayGraph{Int64}()
    p1 = @SVector[0.0, 0.0]
    p2 = @SVector[1.0, 0.0]
    p3 = @SVector[0.0, 1.0]
    pts = [p1, p2, p3]
    DT.add_triangle!(1, 2, 3, T, adj, adj2v, DG; update_ghost_edges=true)
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
    true_DG = DT.DelaunayGraph{Int64}()
    DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 2, 3]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3]...)
    DT.add_neighbour!(true_DG, 3, [1, 2]...)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    p4 = @SVector[1.7, 1.7]
    push!(pts, p4)
    DT.add_triangle!(3, 2, 4, T, adj, adj2v, DG; update_ghost_edges=true)
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
    true_DG = DT.DelaunayGraph{Int64}()
    DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 3, 4, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4]...)
    DT.add_neighbour!(true_DG, 4, [2, 3]...)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    p5 = @SVector[1.0, 3.0]
    p6 = @SVector[3.0, 1.0]
    push!(pts, p5, p6)
    DT.add_triangle!(3, 4, 5, T, adj, adj2v, DG; update_ghost_edges=true)
    DT.add_triangle!(4, 2, 6, T, adj, adj2v, DG; update_ghost_edges=true)
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
    true_DG = DT.DelaunayGraph{Int64}()
    DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 3, 5, 4, 6, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4, 6]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4, 5]...)
    DT.add_neighbour!(true_DG, 4, [2, 3, 5, 6]...)
    DT.add_neighbour!(true_DG, 5, [3, 4]...)
    DT.add_neighbour!(true_DG, 6, [2, 4]...)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    DT.add_triangle!(5, 4, 6, T, adj, adj2v, DG; update_ghost_edges=true)
    true_T = Set{NTuple{3,Int64}}([
        (1, 2, 3),
        (3, 2, 4),
        (3, 4, 5),
        (4, 2, 6),
        (5, 4, 6),
        (2, 1, DT.BoundaryIndex),
        (1, 3, DT.BoundaryIndex),
        (3, 5, DT.BoundaryIndex),
        (5, 6, DT.BoundaryIndex),
        (6, 2, DT.BoundaryIndex)
    ])
    true_adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (3, 2) => 4, (2, 4) => 3, (4, 3) => 2,
            (3, 4) => 5, (4, 5) => 3, (5, 3) => 4,
            (4, 2) => 6, (2, 6) => 4, (6, 4) => 2,
            (5, 4) => 6, (4, 6) => 5, (6, 5) => 4,
            (2, 1) => DT.BoundaryIndex, (1, DT.BoundaryIndex) => 2, (DT.BoundaryIndex, 2) => 1,
            (1, 3) => DT.BoundaryIndex, (3, DT.BoundaryIndex) => 1, (DT.BoundaryIndex, 1) => 3,
            (3, 5) => DT.BoundaryIndex, (5, DT.BoundaryIndex) => 3, (DT.BoundaryIndex, 3) => 5,
            (5, 6) => DT.BoundaryIndex, (6, DT.BoundaryIndex) => 5, (DT.BoundaryIndex, 5) => 6,
            (6, 2) => DT.BoundaryIndex, (2, DT.BoundaryIndex) => 6, (DT.BoundaryIndex, 6) => 2
        )))
    true_adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.BoundaryIndex => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 6), (6, 2), (2, 1)]),
            1 => Set{NTuple{2,Int64}}([(2, 3), (3, DT.BoundaryIndex), (DT.BoundaryIndex, 2)]),
            2 => Set{NTuple{2,Int64}}([(3, 1), (4, 3), (6, 4), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 6)]),
            3 => Set{NTuple{2,Int64}}([(1, 2), (2, 4), (4, 5), (DT.BoundaryIndex, 1), (5, DT.BoundaryIndex)]),
            4 => Set{NTuple{2,Int64}}([(3, 2), (2, 6), (6, 5), (5, 3)]),
            5 => Set{NTuple{2,Int64}}([(3, 4), (4, 6), (6, DT.BoundaryIndex), (DT.BoundaryIndex, 3)]),
            6 => Set{NTuple{2,Int64}}([(4, 2), (5, 4), (DT.BoundaryIndex, 5), (2, DT.BoundaryIndex)])
        )
    )
    true_DG = DT.DelaunayGraph{Int64}()
    DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 3, 5, 6, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4, 6]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4, 5]...)
    DT.add_neighbour!(true_DG, 4, [2, 3, 5, 6]...)
    DT.add_neighbour!(true_DG, 5, [3, 4, 6]...)
    DT.add_neighbour!(true_DG, 6, [2, 4, 5]...)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    @test DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
end

@testset "Larger example" begin
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
    T = Set{NTuple{3,Int64}}()
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DG = DT.DelaunayGraph{Int64}()
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
        DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true)
    end
    true_T = Set{NTuple{3,Int64}}([
        (1, 6, 8),
        (1, 8, 9),
        (9, 8, 10),
        (10, 8, 11),
        (10, 11, 5),
        (11, 7, 5),
        (7, 4, 5),
        (7, 6, 4),
        (6, 3, 4),
        (6, 2, 3),
        (1, 2, 6),
        (8, 6, 7),
        (8, 7, 11),
        (2, 1, DT.BoundaryIndex),
        (1, 9, DT.BoundaryIndex),
        (9, 10, DT.BoundaryIndex),
        (10, 5, DT.BoundaryIndex),
        (5, 4, DT.BoundaryIndex),
        (4, 3, DT.BoundaryIndex),
        (3, 2, DT.BoundaryIndex)
    ])
    true_adj = DT.Adjacent(
        DefaultDict(
            DT.DefaultAdjacentValue,
            Dict(
                (1, 6) => 8, (6, 8) => 1, (8, 1) => 6,
                (1, 8) => 9, (8, 9) => 1, (9, 1) => 8,
                (9, 8) => 10, (8, 10) => 9, (10, 9) => 8,
                (10, 8) => 11, (8, 11) => 10, (11, 10) => 8,
                (10, 11) => 5, (11, 5) => 10, (5, 10) => 11,
                (11, 7) => 5, (7, 5) => 11, (5, 11) => 7,
                (7, 4) => 5, (4, 5) => 7, (5, 7) => 4,
                (7, 6) => 4, (6, 4) => 7, (4, 7) => 6,
                (6, 3) => 4, (3, 4) => 6, (4, 6) => 3,
                (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
                (1, 2) => 6, (2, 6) => 1, (6, 1) => 2,
                (8, 6) => 7, (6, 7) => 8, (7, 8) => 6,
                (8, 7) => 11, (7, 11) => 8, (11, 8) => 7,
                (2, 1) => DT.BoundaryIndex, (1, DT.BoundaryIndex) => 2, (DT.BoundaryIndex, 2) => 1,
                (1, 9) => DT.BoundaryIndex, (9, DT.BoundaryIndex) => 1, (DT.BoundaryIndex, 1) => 9,
                (9, 10) => DT.BoundaryIndex, (10, DT.BoundaryIndex) => 9, (DT.BoundaryIndex, 9) => 10,
                (10, 5) => DT.BoundaryIndex, (5, DT.BoundaryIndex) => 10, (DT.BoundaryIndex, 10) => 5,
                (5, 4) => DT.BoundaryIndex, (4, DT.BoundaryIndex) => 5, (DT.BoundaryIndex, 5) => 4,
                (4, 3) => DT.BoundaryIndex, (3, DT.BoundaryIndex) => 4, (DT.BoundaryIndex, 4) => 3,
                (3, 2) => DT.BoundaryIndex, (2, DT.BoundaryIndex) => 3, (DT.BoundaryIndex, 3) => 2
            )
        )
    )
    true_adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    for ij in edges(true_adj)
        DT.add_edge!(true_adj2v, DT.get_edge(true_adj, ij), ij)
    end
    true_DG = DT.DelaunayGraph{Int64}()
    DT.add_neighbour!(true_DG, DT.BoundaryIndex, 1, 9, 10, 5, 4, 3, 2)
    DT.add_neighbour!(true_DG, 1, 2, 6, 8, 9)
    DT.add_neighbour!(true_DG, 2, 1, 6, 3)
    DT.add_neighbour!(true_DG, 3, 2, 6, 4)
    DT.add_neighbour!(true_DG, 4, 6, 3, 7, 5)
    DT.add_neighbour!(true_DG, 5, 10, 11, 7, 4)
    DT.add_neighbour!(true_DG, 6, 1, 2, 3, 4, 7, 8)
    DT.add_neighbour!(true_DG, 7, 8, 6, 4, 5, 11)
    DT.add_neighbour!(true_DG, 8, 1, 6, 7, 11, 10, 9)
    DT.add_neighbour!(true_DG, 9, 1, 8, 10)
    DT.add_neighbour!(true_DG, 10, 9, 8, 11, 5)
    DT.add_neighbour!(true_DG, 11, 10, 8, 7, 5)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    @test DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
end

@testset "Testing all the boundary deletion cases" begin
    T = Set{NTuple{3,Int64}}()
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DG = DT.DelaunayGraph{Int64}()
    p1 = @SVector[0.0, 0.0]
    p2 = @SVector[1.0, 0.0]
    p3 = @SVector[0.0, 1.0]
    pts = [p1, p2, p3]
    DT.add_triangle!(1, 2, 3, T, adj, adj2v, DG; update_ghost_edges=true)
    p4 = @SVector[1.7, 1.7]
    push!(pts, p4)
    DT.add_triangle!(3, 2, 4, T, adj, adj2v, DG; update_ghost_edges=true)
    p5 = @SVector[1.0, 3.0]
    p6 = @SVector[3.0, 1.0]
    push!(pts, p5, p6)
    DT.add_triangle!(3, 4, 5, T, adj, adj2v, DG; update_ghost_edges=true)
    DT.add_triangle!(4, 2, 6, T, adj, adj2v, DG; update_ghost_edges=true)
    DT.add_triangle!(5, 4, 6, T, adj, adj2v, DG; update_ghost_edges=true)
    DT.delete_triangle!(5, 4, 6, T, adj, adj2v, DG; update_ghost_edges=true)
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
    true_DG = DT.DelaunayGraph{Int64}()
    DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 3, 5, 4, 6, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4, 6]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4, 5]...)
    DT.add_neighbour!(true_DG, 4, [2, 3, 5, 6]...)
    DT.add_neighbour!(true_DG, 5, [3, 4]...)
    DT.add_neighbour!(true_DG, 6, [2, 4]...)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    @test DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
    DT.delete_triangle!(4, 2, 6, T, adj, adj2v, DG; update_ghost_edges=true)
    DT.delete_triangle!(3, 4, 5, T, adj, adj2v, DG; update_ghost_edges=true)
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
    true_DG = DT.DelaunayGraph{Int64}()
    DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 3, 4, 2]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3, 4]...)
    DT.add_neighbour!(true_DG, 3, [1, 2, 4]...)
    DT.add_neighbour!(true_DG, 4, [2, 3]...)
    DT.clear_empty_points!(DG)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    @test DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
    DT.delete_triangle!(3, 2, 4, T, adj, adj2v, DG; update_ghost_edges=true)
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
    true_DG = DT.DelaunayGraph{Int64}()
    DT.add_neighbour!(true_DG, DT.BoundaryIndex, [1, 2, 3]...)
    DT.add_neighbour!(true_DG, 1, [2, 3]...)
    DT.add_neighbour!(true_DG, 2, [1, 3]...)
    DT.add_neighbour!(true_DG, 3, [1, 2]...)
    DT.clear_empty_points!(DG)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    @test DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
    DT.delete_triangle!(1, 2, 3, T, adj, adj2v, DG; update_ghost_edges=true)
    true_T = Set{NTuple{3,Int64}}()
    true_adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    true_adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    true_DG = DT.DelaunayGraph{Int64}()
    DT.clear_empty_points!(DG)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    @test DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
end

@testset "A larger example" begin
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
    T = Set{NTuple{3,Int64}}()
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DG = DT.DelaunayGraph{Int64}()
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
        DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true)
    end
    [DT.delete_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true) for (i, j, k) in (
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
    for ij in edges(true_adj)
        DT.add_edge!(true_adj2v, DT.get_edge(true_adj, ij), ij)
    end
    true_DG = DT.DelaunayGraph{Int64}()
    DT.add_neighbour!(true_DG, DT.BoundaryIndex, 1, 8, 10, 11, 5, 7, 4, 6, 1)
    DT.add_neighbour!(true_DG, 1, 6, 8)
    DT.add_neighbour!(true_DG, 4, 6, 7)
    DT.add_neighbour!(true_DG, 5, 7, 11)
    DT.add_neighbour!(true_DG, 6, 1, 8, 7, 4)
    DT.add_neighbour!(true_DG, 7, 4, 6, 8, 11, 5)
    DT.add_neighbour!(true_DG, 8, 1, 6, 7, 11, 10)
    DT.add_neighbour!(true_DG, 10, 8, 11)
    DT.add_neighbour!(true_DG, 11, 10, 8, 7, 5)
    DT.clear_empty_points!(DG)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, true_T, true_adj, true_adj2v, true_DG)
    @test DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
end

@testset "Can we add and delete ghost triangles from a ghost triangle?" begin
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
    DT.compute_centroid!(pts)
    T = Set{NTuple{3,Int64}}()
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DG = DT.DelaunayGraph{Int64}()
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
        DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true)
    end

    _T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
    DT.add_ghost_triangles!(_T, _adj, _adj2v, _DG)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)

    _T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
    DT.remove_ghost_triangles!(T, adj, adj2v, DG)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
end