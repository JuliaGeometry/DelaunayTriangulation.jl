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
        (T, adj, adj2v, DG), HG = DT.triangulate_berg(pts;
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
                (((1, 5), (5, 6), (6, 4), (4, 3), (3, 1)) .=> DT.BoundaryIndex)...
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
        true_DG = relabel(UndirectedGraph(
                IntegerType[
                    0 1 0 1 1 1 1
                    1 0 1 1 0 1 0
                    0 1 0 1 1 1 0
                    1 1 1 0 1 0 0
                    1 0 1 1 0 1 1
                    1 1 1 0 1 0 1
                    1 0 0 0 1 1 0
                ]
            ), Dict(1:7 .=> 0:6))
        @test T == true_T
        @test DT.adjacent(adj) == true_adj
        @test DT.adjacent2vertex(adj2v) == true_adj2v
        @test graph(DG) == true_DG
        @test DT.is_delaunay(T, pts)
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
        tri, HG = DT.triangulate_berg(pts;
            IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType)
        T, adj, adj2v, DG = tri
        true_T = TrianglesType([
            (5, 10, 2),
            (5, 9, 7),
            (3, 6, 2),
            (1, 8, 9),
            (7, 9, 8),
            (7, 8, 3),
            (1, 3, 8),
            (1, 6, 3),
            (6, 1, 4),
            (10, 5, 7),
            (10, 3, 2),
            (3, 10, 7),
            (4, 1, 9)
        ]
        )
        true_adj = DefaultDict(DT.DefaultAdjacentValue,
            Dict{EdgeType,IntegerType}(
                (5, 10) => 2, (10, 2) => 5, (2, 5) => 10,
                (5, 9) => 7, (9, 7) => 5, (7, 5) => 9,
                (3, 6) => 2, (6, 2) => 3, (2, 3) => 6,
                (1, 8) => 9, (8, 9) => 1, (9, 1) => 8,
                (7, 9) => 8, (9, 8) => 7, (8, 7) => 9,
                (7, 8) => 3, (8, 3) => 7, (3, 7) => 8,
                (1, 3) => 8, (3, 8) => 1, (8, 1) => 3,
                (1, 6) => 3, (6, 3) => 1, (3, 1) => 6,
                (6, 1) => 4, (1, 4) => 6, (4, 6) => 1,
                (10, 5) => 7, (5, 7) => 10, (7, 10) => 5,
                (10, 3) => 2, (3, 2) => 10, (2, 10) => 3,
                (3, 10) => 7, (10, 7) => 3, (7, 3) => 10,
                (4, 1) => 9, (1, 9) => 4, (9, 4) => 1,
                (((4, 9), (9, 5), (5, 2), (2, 6), (6, 4)) .=> DT.BoundaryIndex)...
            )
        )
        true_adj2v = Dict{IntegerType,EdgesType}(
            DT.BoundaryIndex => EdgesType([(4, 9), (9, 5), (5, 2), (2, 6), (6, 4)]),
            1 => EdgesType([(4, 6), (6, 3), (3, 8), (8, 9), (9, 4)]),
            2 => EdgesType([(5, 10), (10, 3), (3, 6)]),
            3 => EdgesType([(10, 7), (7, 8), (8, 1), (1, 6), (6, 2), (2, 10)]),
            4 => EdgesType([(6, 1), (1, 9)]),
            5 => EdgesType([(9, 7), (7, 10), (10, 2)]),
            6 => EdgesType([(2, 3), (3, 1), (1, 4)]),
            7 => EdgesType([(9, 8), (8, 3), (3, 10), (10, 5), (5, 9)]),
            8 => EdgesType([(9, 1), (1, 3), (3, 7), (7, 9)]),
            9 => EdgesType([(4, 1), (1, 8), (8, 7), (7, 5)]),
            10 => EdgesType([(7, 3), (3, 2), (2, 5), (5, 7)])
        )
        true_DG = UndirectedGraph{Int64}()
        add!.(Ref(true_DG), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        add!.(Ref(true_DG), 0, [2, 4, 5, 6, 9])
        add!.(Ref(true_DG), 1, [4, 6, 3, 8, 9])
        add!.(Ref(true_DG), 2, [5, 10, 3, 6])
        add!.(Ref(true_DG), 3, [10, 7, 8, 1, 6, 2])
        add!.(Ref(true_DG), 4, [1, 6, 9])
        add!.(Ref(true_DG), 5, [9, 7, 10, 2])
        add!.(Ref(true_DG), 6, [2, 3, 1, 4])
        add!.(Ref(true_DG), 7, [9, 8, 3, 10, 5])
        add!.(Ref(true_DG), 8, [1, 9, 7, 3])
        add!.(Ref(true_DG), 9, [4, 1, 8, 7, 5])
        add!.(Ref(true_DG), 10, [7, 3, 2, 5])
        @test T == true_T
        @test adj.adjacent == true_adj
        @test adj2v.adjacent2vertex == true_adj2v
        @test DG.graph.N == true_DG.N
        @test DT.is_delaunay(T, pts)
        @test DT.validate_triangulation(T, adj, adj2v, DG, pts)

        ## Random triangulations
        for _ in 1:178
            n = rand(3:1000)
            x = rand(n)
            y = rand(n)
            pts = [(x, y) for (x, y) in zip(x, y)]
            (T, adj, adj2v, DG), HG = DT.triangulate_berg(pts;
                IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType)
            @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
            @test all(DT.isoriented(T, pts) == 1 for T in T)
            for (i, j) in adj2v.adjacent2vertex[DT.BoundaryIndex]
                @test i ≥ DT.FirstPointIndex && j ≥ DT.FirstPointIndex
            end
            num_pts = length(pts)
            num_eg = length(DT.graph(DG).E) - length(DT.get_neighbour(DG, DT.BoundaryIndex)) # don't include the boundary connections
            num_tris = length(T)
            @test 1 == num_pts - num_eg + num_tris # Euler's formula
        end

        ## Randomised insertion order 
        for _ in 1:32
            n = rand(3:1000)
            x = rand(n)
            y = rand(n)
            pts = [(x, y) for (x, y) in zip(x, y)]
            (T, adj, adj2v, DG), HG = DT.triangulate_berg(pts; randomise=false,
                IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType)
            @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
            @test all(DT.isoriented(T, pts) == 1 for T in T)
            for (i, j) in adj2v.adjacent2vertex[DT.BoundaryIndex]
                @test i ≥ DT.FirstPointIndex && j ≥ DT.FirstPointIndex
            end
            num_pts = length(pts)
            num_eg = length(DT.graph(DG).E) - length(DT.get_neighbour(DG, DT.BoundaryIndex))
            num_tris = length(T)
            @test 1 == num_pts - num_eg + num_tris # Euler's formula
            (_T, _adj, _adj2v, _DG), _HG = DT.triangulate_berg(pts; randomise=false,
                IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType)
            @test T == _T # Not only are the triangles equivalent, they are exactly equal 
            @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
        end

        ## Removing the bounding triangle 
        for _ in 1:32
            n = rand(3:1000)
            x = rand(n)
            y = rand(n)
            pts = [(x, y) for (x, y) in zip(x, y)]
            tri, HG = DT.triangulate_berg(pts;
                IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType)
            T, adj, adj2v, DG = tri
            @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
            @test all(DT.isoriented(T, pts) == 1 for T in T)
            for (i, j) in adj2v.adjacent2vertex[DT.BoundaryIndex]
                @test i ≥ DT.FirstPointIndex && j ≥ DT.FirstPointIndex
            end
            num_pts = length(pts)
            num_eg = length(DT.graph(DG).E) - length(DT.get_neighbour(DG, DT.BoundaryIndex))
            num_tris = length(T)
            @test 1 == num_pts - num_eg + num_tris # Euler's formula
            (_T, _adj, _adj2v, _DG), _HG = DT.triangulate_berg(pts; trim=false,
                IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType)
            @test DT.validate_triangulation(_T, _adj, _adj2v, _DG, pts)
            @test all(DT.isoriented(_T, pts) == 1 for _T in _T)
            @test !DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
            DT.remove_bounding_triangle!(_T, _adj, _adj2v, _DG)
            @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
            @test DT.validate_triangulation(_T, _adj, _adj2v, _DG, pts)
            @test all(DT.isoriented(_T, pts) == 1 for _T in _T)
            for (i, j) in adj2v.adjacent2vertex[DT.BoundaryIndex]
                @test i ≥ DT.FirstPointIndex && j ≥ DT.FirstPointIndex
            end
            num_pts = length(pts)
            num_eg = length(DT.graph(_DG).E) - length(DT.get_neighbour(DG, DT.BoundaryIndex))
            num_tris = length(_T)
            @test 1 == num_pts - num_eg + num_tris # Euler's formula
        end
    end
end

@testset "Matrix points" begin
    function DT._get_point(pts::AbstractMatrix, i)
        return @view pts[:, i]
    end
    function DT._eachindex(pts::AbstractMatrix)
        return axes(pts, 2)
    end
    for _ in 1:100
        pts = rand(2, 250)
        (T, adj, adj2v, DG), _ = DT.triangulate_berg(pts)
        pts2 = [pts[:, i] for i in axes(pts, 2)]
        (_T, _adj, _adj2v, _DG), _ = DT.triangulate_berg(pts2)
        @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    end
end

@testset "Testing point insertion with the Bowyer-Watson method" begin
    ## Existing triangulation with no bounding etc 
    p1 = (5.0, 6.0)
    p2 = (9.0, 6.0)
    p3 = (13.0, 5.0)
    p4 = (10.38, 0.0)
    p5 = (12.64, -1.69)
    p6 = (2.0, -2.0)
    p7 = (3.0, 4.0)
    p8 = (7.5, 3.53)
    p9 = (4.02, 1.85)
    p10 = (4.26, 0.0)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
    Random.seed!(928881)
    (T, adj, adj2v, DG), HG = DT.triangulate_berg(pts)
    p11 = (6.0, 2.5)
    push!(pts, p11)

    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 11)
    (_T, _adj, _adj2v, _DG), _HG = DT.triangulate_berg(pts)
    DT.clear_empty_keys!(adj)
    DT.clear_empty_keys!(_adj)
    @test DT.compare_triangle_sets(T, _T) &&
          adjacent(adj) == adjacent(_adj) &&
          adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
          graph(DG) == graph(_DG)

    p12 = (10.3, 2.85)
    push!(pts, p12)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 12)
    (_T, _adj, _adj2v, _DG), _HG = DT.triangulate_berg(pts)
    DT.clear_empty_keys!(adj)
    DT.clear_empty_keys!(_adj)
    @test DT.compare_triangle_sets(T, _T) &&
          adjacent(adj) == adjacent(_adj) &&
          adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
          graph(DG) == graph(_DG)

    p13 = (7.5, 3.5)
    push!(pts, p13)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 13)
    (_T, _adj, _adj2v, _DG), _HG = DT.triangulate_berg(pts)
    DT.clear_empty_keys!(adj)
    DT.clear_empty_keys!(_adj)
    @test DT.compare_triangle_sets(T, _T) &&
          adjacent(adj) == adjacent(_adj) &&
          adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
          graph(DG) == graph(_DG)

    ## Now do this for some circular points so that we can do some random tests - only testing for interior insertion here
    R = 10.0
    n = 1381
    pts = [(rand((-1, 1)) * R * rand(), rand((-1, 1)) * R * rand()) for _ in 1:n] # rand((-1, 1)) to get random signs
    pushfirst!(pts, (-11.0, -11.0), (11.0, -11.0), (11.0, 11.0), (-11.0, 11.0)) # bounding box to guarantee interior insertions only
    (T, adj, adj2v, DG), HG = @views DT.triangulate_berg(pts[1:7])
    (_T, _adj, _adj2v, _DG), _HG = @views DT.triangulate_berg(pts[1:7])
    n = length(pts)
    for i in 8:n
        DT.add_point_bowyer!(T, adj, adj2v, DG, pts, i)
        DT.add_point_berg!(_T, _adj, _adj2v, _DG, _HG, pts, i)
        @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    end
end

@testset "Adding a point outside of the triangulation" begin
    ## Small example
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
    p12 = @SVector[4.382, 3.2599]
    push!(pts, p12)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 12)
    (_T, _adj, _adj2v, _DG), _HG = DT.triangulate_berg(pts)
    DT.add_ghost_triangles!(_T, _adj, _adj2v, _DG)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    @test DT.compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)

    p13 = @SVector[-5.25303, 4.761]
    p14 = @SVector[-9.83801, 0.562]
    p15 = @SVector[-7.15986, -5.99]
    p16 = @SVector[4.79, 2.74]
    p17 = @SVector[3.77, 2.7689]
    push!(pts, p13, p14, p15, p16, p17)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 13)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 14)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 15)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 16)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 17)
    (_T, _adj, _adj2v, _DG), _HG = DT.triangulate_berg(pts)
    DT.add_ghost_triangles!(_T, _adj, _adj2v, _DG)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    @test DT.compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    @test DT.compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, pts)

    DT.remove_ghost_triangles!(T, adj, adj2v, DG)
    @test !DT.compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    @test !DT.compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, pts)
end

@testset "A small example" begin
    Random.seed!(929266666)
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
    p12 = @SVector[4.382, 3.2599]
    p13 = @SVector[-5.25303, 4.761]
    p14 = @SVector[-9.83801, 0.562]
    p15 = @SVector[-7.15986, -5.99]
    p16 = @SVector[4.79, 2.74]
    p17 = @SVector[3.77, 2.7689]
    p18 = @SVector[15.19, 15.40]
    p19 = @SVector[22.86, 9.93]
    p20 = @SVector[25.96, -4.457]
    p21 = @SVector[20.298, -13.444]
    p22 = @SVector[4.626, -14.2757]
    p23 = @SVector[-4.739, -14.313]
    p24 = @SVector[-7.3448, -14.50]
    p25 = @SVector[-10.06, -14.35]
    p26 = @SVector[-12.97, -13.63]
    p27 = @SVector[-15.46, -11.405]
    p28 = @SVector[-16.899, -4.419]
    p29 = @SVector[-15.766, 5.4368]
    p30 = @SVector[-11.31, 11.367]
    p31 = @SVector[-2.66, 14.95]
    p32 = @SVector[6.58, 15.746]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8,
        p9, p10, p11, p12, p13, p14, p15, p16, p17,
        p18, p19, p20, p21, p22, p23, p24, p25,
        p26, p27, p28, p29, p30, p31, p32]
    for _ in 1:500
        T, adj, adj2v, DG = DT.triangulate_bowyer(pts; trim=false)
        (_T, _adj, _adj2v, _DG), _HG = DT.triangulate_berg(pts)
        @test DT.compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    end
end

@testset "Larger random examples" begin
    for r in 1:137
        for IntegerType in (Int64, Int32, Int16)
            n = rand(3:1000)
            pts = rand(SVector{2,Float64}, n)
            T, adj, adj2v, DG = DT.triangulate_bowyer(pts; trim=false, IntegerType)
            @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
            _T, _adj, _adj2v, _DG = deepcopy(T), deepcopy(adj), deepcopy(adj2v), deepcopy(DG)
            DT.remove_ghost_triangles!(_T, _adj, _adj2v, _DG)
            @test DT.validate_triangulation(_T, _adj, _adj2v, _DG, pts)
            DT.add_ghost_triangles!(_T, _adj, _adj2v, _DG)
            DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
        end
    end
end

@testset "Matrix points" begin
    function DT._get_point(pts::AbstractMatrix, i)
        return @view pts[:, i]
    end
    function DT._eachindex(pts::AbstractMatrix)
        return axes(pts, 2)
    end
    for _ in 1:100
        pts = rand(2, 250)
        tri = DT.triangulate_bowyer(pts)
        T, adj, adj2v, DG = tri
        pts2 = [pts[:, i] for i in axes(pts, 2)]
        _T, _adj, _adj2v, _DG = DT.triangulate_bowyer(pts2)
        @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    end
end

@testset "Shifting the initial triangle" begin # Actually, none of this is needed - just use setdiff!. Oh well, here it is.
    pt_order = [1, 2, 3, 4, 5, 6, 7]
    pt_order_2 = deepcopy(pt_order)
    DT.select_valid_start_of_vector!(pt_order_2, Set{Int64}())
    @test pt_order_2 == pt_order

    pt_order_2 = deepcopy(pt_order)
    DT.select_valid_start_of_vector!(pt_order_2, Set{Int64}([1]))
    @test pt_order_2 == [2, 3, 4, 5, 6, 7, 1]

    pt_order_2 = deepcopy(pt_order)
    DT.select_valid_start_of_vector!(pt_order_2, Set{Int64}([2]))
    @test pt_order_2 == [1, 3, 4, 5, 6, 7, 2]

    pt_order_2 = deepcopy(pt_order)
    DT.select_valid_start_of_vector!(pt_order_2, Set{Int64}([3]))
    @test pt_order_2 == [1, 2, 4, 5, 6, 7, 3]

    pt_order_2 = deepcopy(pt_order)
    DT.select_valid_start_of_vector!(pt_order_2, Set{Int64}([1, 2]))
    @test pt_order_2 == [3, 4, 5, 6, 7, 1, 2]

    pt_order_2 = deepcopy(pt_order)
    DT.select_valid_start_of_vector!(pt_order_2, Set{Int64}([2, 3]))
    @test pt_order_2 == [1, 4, 5, 6, 7, 2, 3]

    pt_order_2 = deepcopy(pt_order)
    DT.select_valid_start_of_vector!(pt_order_2, Set{Int64}([3, 1]))
    @test pt_order_2 == [2, 4, 5, 6, 7, 1, 3]

    pt_order_2 = deepcopy(pt_order)
    DT.select_valid_start_of_vector!(pt_order_2, Set{Int64}([1, 2, 3]))
    @test pt_order_2 == [4, 5, 6, 7, 1, 2, 3]

    pt_order_2 = deepcopy(pt_order)
    DT.select_valid_start_of_vector!(pt_order_2, Set{Int64}([4, 5, 6, 7]))
    @test pt_order_2 == [1, 2, 3, 4, 5, 6, 7]

    for _ in 1:500
        n = rand(3:500)
        pt_order = collect(1:n)
        m = 0:(n-3)
        for _m in m
            skip_idx = Set{Int64}(rand(pt_order, _m))
            DT.select_valid_start_of_vector!(pt_order, skip_idx)
            @test pt_order[begin] ∉ skip_idx && pt_order[begin+1] ∉ skip_idx && pt_order[begin+2] ∉ skip_idx
        end
        shuffle!(pt_order)
        m = 0:(n-3)
        for _m in m
            skip_idx = Set{Int64}(rand(pt_order, _m))
            DT.select_valid_start_of_vector!(pt_order, skip_idx)
            @test pt_order[begin] ∉ skip_idx && pt_order[begin+1] ∉ skip_idx && pt_order[begin+2] ∉ skip_idx
        end
    end
end

@testset "Skipping points" begin
    for _ in 1:250
        n = rand(3:500)
        pts = rand(2, n)
        skip_pts = Set{Int64}(rand(1:n, n ÷ 4))
        T, adj, adj2v, DG = DT.triangulate_bowyer(pts; skip_pts)
        @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
        for i in skip_pts
            @test i ∉ DG.graph.V
        end
        (_T, _adj, _adj2v, _DG), _ = DT.triangulate_berg(pts; skip_pts, trim=false)
        @test DT.validate_triangulation(_T, _adj, _adj2v, _DG, pts)
        for i in skip_pts
            @test i ∉ DG.graph.V
        end

        T, adj, adj2v, DG = DT.triangulate_bowyer(pts; skip_pts, randomise=false)
        @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
        for i in skip_pts
            @test i ∉ DG.graph.V
        end
        (_T, _adj, _adj2v, _DG), _ = DT.triangulate_berg(pts; skip_pts, randomise=false, trim=false)
        @test DT.validate_triangulation(_T, _adj, _adj2v, _DG, pts)
        for i in skip_pts
            @test i ∉ DG.graph.V
        end
    end
end

@testset "Construction" begin
    a = 0.0
    b = 5.0
    c = 0.0
    d = 10.0
    Nˣ = 5
    Nʸ = 3
    (T, adj, adj2v, DG, pts), BN = triangulate_structured(a, b, c, d, Nˣ, Nʸ; return_boundary_types=true)
    @test length(edges(adj)) ÷ 2 == 30
    @test num_triangles(T) == 16
    @test pts == reduce(hcat, [[0, 0], [1.25, 0], [2.5, 0], [3.75, 0],
        [5, 0], [0, 5], [1.25, 5.0], [2.5, 5], [3.75, 5], [5, 5], [0, 10], [1.25, 10],
        [2.5, 10], [3.75, 10], [5, 10]])
    @test BN[1] == [1, 2, 3, 4, 5]
    @test BN[2] == [5, 10, 15]
    @test BN[3] == [15, 14, 13, 12, 11]
    @test BN[4] == [11, 6, 1]

    (T, adj, adj2v, DG, pts), BN = triangulate_structured(a, b, c, d, Nˣ, Nʸ; return_boundary_types=true, single_boundary=true)
    @test length(edges(adj)) ÷ 2 == 30
    @test num_triangles(T) == 16
    @test pts == reduce(hcat, [[0, 0], [1.25, 0], [2.5, 0], [3.75, 0],
        [5, 0], [0, 5], [1.25, 5.0], [2.5, 5], [3.75, 5], [5, 5], [0, 10], [1.25, 10],
        [2.5, 10], [3.75, 10], [5, 10]])
    @test BN == [1, 2, 3, 4, 5, 10, 15, 14, 13, 12, 11, 6]
end

@testset "Collinear" begin
    a = 0.0
    b = 10.0
    c = 0.0
    d = 10.0
    nx = 11
    ny = 11
    function DT._get_point(pts::AbstractMatrix, i)
        return @view pts[:, i]
    end
    function DT._eachindex(pts::AbstractMatrix)
        return axes(pts, 2)
    end
    Random.seed!(299756756791)
    (_T, _adj, _adj2v, _DG, pts), BN = triangulate_structured(a, b, c, d, nx, ny; return_boundary_types=true, single_boundary=true)

    T, adj, adj2v, DG = DT.triangulate_bowyer(pts; randomise=false)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(2820088771)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(20288888)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(1375560)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(17115)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(668591)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    for _ in 1:5000
        n = rand(1:2000000)
        Random.seed!(n)
        @show n
        T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
        e1 = DT.validate_triangulation(T, adj, adj2v, DG, pts)
        if !e1
            throw("") # Want to capture the seed. This won't always happen even if the triangulation fails, note, since we can have infinite loops when the point location fails.
        end
        @test e1
        @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]
    end

    a = 0.0
    b = 5.0
    c = 0.0
    d = 17.0
    nx = 67
    ny = 31
    (_T, _adj, _adj2v, _DG, pts), BN = triangulate_structured(a, b, c, d, nx, ny; return_boundary_types=true, single_boundary=true)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts; randomise=false)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(2820088771)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(20288888)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(1375560)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(17115)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]

    Random.seed!(668591)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
    @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]
    for _ in 1:60
        m = rand(1:2000000)
        Random.seed!(m)
        @show m
        T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
        e1 = DT.validate_triangulation(T, adj, adj2v, DG, pts)
        if !e1
            throw("") # Want to capture the seed. This won't always happen even if the triangulation fails, note, since we can have infinite loops when the point location fails.
        end
        @test e1
        @test DG.graph.N[DT.BoundaryIndex] == _DG.graph.N[DT.BoundaryIndex]
    end
end

@testset "Handling duplicate points" begin
    a = 0.0
    b = 1.0
    c = 0.0
    d = 0.0
    nx = 23
    ny = 28
    T, adj, adj2v, DG, pts = triangulate_structured(a, b, c, d, nx, ny)
    @test_throws "The points must all be unique." triangulate_bowyer(pts)
    @test_throws "The points must all be unique." triangulate_berg(pts)
end

@testset "Issue #21" begin
    Random.seed!(292888111)
    L = 2.0
    R = 1.0
    num_boundary_cells = 250
    num_interior_cells = 500

    ## The boundary 
    boundary_cells = [
        [[x, 0.0] for x in LinRange(0, L, num_boundary_cells ÷ 4)]...,
        [[L, y] for y in LinRange(0, L, num_boundary_cells ÷ 4)]...,
        [[x, L] for x in LinRange(L, 0, num_boundary_cells ÷ 4)]...,
        [[0.0, y] for y in LinRange(L, 0, num_boundary_cells ÷ 4)]...
    ]

    ## Generate the interior 
    x = L * rand(num_interior_cells)
    y = L * rand(num_interior_cells)
    interior_cells = [[x, y] for (x, y) in zip(x, y)]

    # Filter out the circle 
    bad_idx = Int64[]
    void_centre = [L / 2, L / 2]
    for i in eachindex(interior_cells)
        @views _x, _y = interior_cells[i]
        radsq = (_x - void_centre[1])^2 + (_y - void_centre[2])^2
        radsq < R^2 && push!(bad_idx, i)
    end
    deleteat!(interior_cells, bad_idx)

    ## Combine the boundary and interior cells
    cells = vcat(interior_cells, boundary_cells)
    interior_cell_idx = eachindex(interior_cells)
    boundary_cell_idx = lastindex(interior_cells) .+ eachindex(boundary_cells)
    cells = reduce(hcat, cells)
    cells = unique(cells; dims=2)

    for _ in 1:500
        (T, adj, adj2v, dg), HG = DT.triangulate_berg(cells)
        @test DT.validate_triangulation(T, adj, adj2v, dg, cells)
    end
end

@testset "Edges all exist" begin
    pts = rand(2, 500)
    T, adj, adj2v, DG = triangulate_bowyer(pts)
    for e in edges(adj)
        @test e ∈ edges(DG) || reverse(e) ∈ edges(DG)
    end
    for e in edges(DG)
        if DT.BoundaryIndex ∉ e
            @test e ∈ edges(adj)
        end
    end
end