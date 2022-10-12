@testset "Can we identify boundary edges" begin
    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (1, 3) => 5, (3, 5) => 1,
            (4, 5) => DT.BoundaryIndex, (5, 2) => DT.BoundaryIndex,
            (2, 3) => DT.BoundaryIndex, (3, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex
        )
    ))
    adj2v = DT.Adjacent2Vertex(Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(5, 4), (3, 5), (6, 3), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3)]),
        3 => Set{NTuple{2,Int64}}([(1, 6), (5, 1), (2, 5)]),
        4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(1, 4), (3, 1)])
    ))
    for (k, S) in adjacent2vertex(adj2v)
        for (i, j) in S
            if k == DT.BoundaryIndex
                @test DT.is_boundary_edge((i, j), adj)
                @test DT.is_boundary_edge(i, j, adj)
            else
                @test !DT.is_boundary_edge((i, j), adj)
                @test !DT.is_boundary_edge(i, j, adj)
            end
        end
    end
end

@testset "Can we check edge validity?" begin 
    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (5, 1) => 3, (1, 3) => 5, (3, 5) => 1,
        (4, 5) => DT.BoundaryIndex, (5, 2) => DT.BoundaryIndex,
        (2, 3) => DT.BoundaryIndex, (3, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex
    )
))
    for (i, j) in edges(adj)
        @test DT.edge_exists(i, j, adj)
    end
    for _ in 1:1000
        i, j = abs.(rand(Int64, 3))
        @test !DT.edge_exists(i, j, adj)
    end
end

@testset "Is choose_uvw working correctly?" begin
    for _ in 1:1000
        i, j, k = abs.(rand(Int64, 3))
        @test DT.choose_uvw(true, false, false, i, j, k) == (i, j, k)
        @test DT.choose_uvw(false, true, false, i, j, k) == (j, k, i)
        @test DT.choose_uvw(false, false, true, i, j, k) == (k, i, j)
    end
end