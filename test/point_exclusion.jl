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
        _T, _adj, _adj2v, _DG = DT.triangulate_berg(pts; skip_pts)
        @test DT.validate_triangulation(_T, _adj, _adj2v, _DG, pts)
        for i in skip_pts
            @test i ∉ DG.graph.V
        end

        T, adj, adj2v, DG = DT.triangulate_bowyer(pts; skip_pts, randomise=false)
        @test DT.validate_triangulation(T, adj, adj2v, DG, pts)
        for i in skip_pts
            @test i ∉ DG.graph.V
        end
        _T, _adj, _adj2v, _DG = DT.triangulate_berg(pts; skip_pts, randomise=false)
        @test DT.validate_triangulation(_T, _adj, _adj2v, _DG, pts)
        for i in skip_pts
            @test i ∉ DG.graph.V
        end
    end
end