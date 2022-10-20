
@testset "Testing all the boundary cases" begin
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