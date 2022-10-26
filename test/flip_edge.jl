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
