@testset "Triangulations with de Berg's method" begin
    for IntegerType in (Int64, Int32, Int16)
        # Standard 
        p1 = @SVector[-3.0, 3.0]
        p2 = @SVector[2.0, 4.0]
        p3 = @SVector[-1.0, 1.0]
        p4 = @SVector[3.0, -1.0]
        p5 = @SVector[6.0, 5.0]
        p6 = @SVector[8.0, -3.0]
        pts = [p1, p2, p3, p4, p5, p6]
        TriangleType = (NTuple{3,IntegerType})
        EdgeType = (NTuple{2,IntegerType})
        TrianglesType = Set{TriangleType}
        EdgesType = Set{EdgeType}
        Random.seed!(92994881)
        T, adj, adj2v, DG, HG = DT.triangulate_berg(pts;
            IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType)
        true_T = TrianglesType([
            (5, 1, 2),
            (3, 2, 1),
            (4, 2, 3),
            (4, 5, 2),
            (4, 6, 5)
        ])
        true_adj = DefaultDict(DT.DefaultAdjacentValue,
            Dict{EdgeType,IntegerType}(
                (5, 1) => 2, (1, 2) => 5, (2, 5) => 1,
                (1, 3) => 2, (3, 2) => 1, (2, 1) => 3,
                (4, 2) => 3, (2, 3) => 4, (3, 4) => 2,
                (4, 5) => 2, (5, 2) => 4, (2, 4) => 5,
                (4, 6) => 5, (6, 5) => 4, (5, 4) => 6,
                (((1, 5), (5, 6), (6, 4), (4, 3), (3, 1)) .=> DT.BoundaryIndex)...,
                (((1, -1), (3, -1), (-1, 2), (2, -2), (2, 6), (-3, 2), (3, 6)) .=> DT.DefaultAdjacentValue)...
            )
        )
        true_adj2v = Dict{IntegerType,EdgesType}(
            DT.BoundaryIndex => EdgesType([(1, 5), (5, 6), (6, 4), (4, 3), (3, 1)]),
            1 => EdgesType([(3, 2), (2, 5)]),
            2 => EdgesType([(4, 5), (5, 1), (1, 3), (3, 4)]),
            3 => EdgesType([(4, 2), (2, 1)]),
            4 => EdgesType([(6, 5), (5, 2), (2, 3)]),
            5 => EdgesType([(1, 2), (2, 4), (4, 6)]),
            6 => EdgesType([(5, 4)])
        )
        true_DG = UndirectedGraph(
            IntegerType[
                0 1 1 0 1 0
                1 0 1 1 1 0
                1 1 0 1 0 0
                0 1 1 0 1 1
                1 1 0 1 0 1
                0 0 0 1 1 0
            ]
        )
        @test T == true_T
        @test DT.adjacent(adj) == true_adj
        @test DT.adjacent2vertex(adj2v) == true_adj2v
        @test graph(DG) == true_DG
        @test DT.is_delaunay(adj, pts)
        @test DT.validate_triangulation(T, adj, adj2v, DG, pts)

        ## Triangulation with some collinear points
        p1 = @SVector[0.0, 1.0]
        p2 = @SVector[3.0, -1.0]
        p3 = @SVector[2.0, 0.0]
        p4 = @SVector[-1.0, 2.0]
        p5 = @SVector[4.0, 2.0]
        p6 = @SVector[-2.0, -1.0]
        p7 = @SVector[2.0, 1.0]
        p8 = @SVector[1.0, 1.0]
        p9 = @SVector[1.5, 2.0]
        p10 = @SVector[2.5, 0.5]
        pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
        Random.seed!(929947476781)
        T, adj, adj2v, DG, HG = DT.triangulate_berg(pts;
            IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType)
        @test DT.is_delaunay(adj, pts)
        @test DT.validate_triangulation(T, adj, adj2v, DG, pts)

        # Random triangulations
        for _ in 1:250
            n = rand(5:1000)
            x = rand(n)
            y = rand(n)
            pts = [(x, y) for (x, y) in zip(x, y)]
            T, adj, adj2v, DG, HG = DT.triangulate_berg(pts;
                IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType)
            @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
            @test all(DT.isoriented(T, pts) == 1 for T in T)
            num_pts = length(pts)
            num_eg = length(DT.graph(DG).E)
            num_tris = length(T)
            @test 1 == num_pts - num_eg + num_tris # Euler's formula
        end
    end
end