function example_triangulation()
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    p4 = @SVector[-1.0, 2.0]
    p5 = @SVector[4.0, 2.0]
    p6 = @SVector[-2.0, -1.0]
    pts = [p1, p2, p3, p4, p5, p6]
    T = Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (3, 2, 5),
        (4, 1, 5),
        (4, 6, 1),
        (5, 1, 3)
    ])
    A = [
        0 0 1 1 1 1
        0 0 1 0 1 0
        1 1 0 0 1 1
        1 0 0 0 1 1
        1 1 1 1 0 0
        1 0 1 1 0 0
    ]
    DG = DT.DelaunayGraph(UndirectedGraph(A))
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
    all(all(DT.get_edge(adj, i, j) == k for (i, j) in DT.get_edge(adj2v, k)) for k in keys(adj2v.adjacent2vertex)) &&
        return pts, T, DG, adj, adj2v
end

function example_empty_triangulation()
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    pts = [p1, p2, p3]
    T = Set{NTuple{3,Int64}}([])
    A = zeros(Int64, 0, 0)
    DG = DT.DelaunayGraph(UndirectedGraph(A))
    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue, Dict{NTuple{2,Int64},Int64}()))
    adj2v = DT.Adjacent2Vertex(Dict(DT.BoundaryIndex => Set{NTuple{2,Int64}}()))
    return pts, T, DG, adj, adj2v
end
    