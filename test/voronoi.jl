@testset "Are we computing the circumcenters correctly?" begin
    for _ in 1:200
        pts = rand(SVector{2,Float64}, 200)
        T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
        for V in T
            i, j, k = indices(V)
            pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
            true_centxy, rad = circle_three_points(pᵢ, pⱼ, pₖ)
            true_centx, true_centy = true_centxy
            centx, centy = DT.circumcenter(V, pts)
            @test centx ≈ true_centx
            @test centy ≈ true_centy
        end
        cents, triangle_to_idx, idx_to_triangle = DT.circumcenters(T, pts)
        for V in T
            _i = triangle_to_idx[V]
            @test V == idx_to_triangle[_i]
            i, j, k = indices(V)
            pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
            true_centxy, rad = circle_three_points(pᵢ, pⱼ, pₖ)
            true_centx, true_centy = true_centxy
            centx, centy = cents[_i]
            @test centx ≈ true_centx
            @test centy ≈ true_centy
        end
    end
end

@testset "Voronoi polygons" begin
    Random.seed!(2992991)
    qa = (-1.4, 0.61)
    qb = (1.0, 1.0)
    qc = (0.0, -2.0)
    qd = (1.0, 3.0)
    qe = (0.0, 0.0)
    qf = (-0.56, 1.53)
    qg = (2.22, 1.43)
    qh = (2.46, -0.57)
    qi = (-0.68, -1.07)
    pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts; trim=true)
    cents, tri_to_idx, idx_to_tri = DT.circumcenters(T, pts)
    i = 5
    cell_idx = DT.get_voronoi_vertices(T, adj, adj2v, DG, i, tri_to_idx)
    slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx, adj, adj2v)
    shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
    circshift!(slow_cell_idx, shift_idx - 1)
    @test cell_idx == slow_cell_idx

    cell_idx = DT.get_voronoi_vertices(T, adj, adj2v, DG, i, tri_to_idx)
    slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx, adj, adj2v)
    shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
    circshift!(slow_cell_idx, shift_idx - 1)
    @test cell_idx == slow_cell_idx

    for i in eachindex(pts)
        cell_idx = DT.get_voronoi_vertices(T, adj, adj2v, DG, i, tri_to_idx)
        slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx, adj, adj2v)
        nonzero_cell_idx = findfirst(cell_idx .≠ DT.BoundaryIndex)
        circshift!(cell_idx, nonzero_cell_idx - 1)
        shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
        circshift!(slow_cell_idx, shift_idx - 1)
        @test cell_idx == slow_cell_idx
    end

    Random.seed!(2992991)
    qa = (-1.4, 0.61)
    qb = (1.0, 1.0)
    qc = (0.0, -2.0)
    qd = (1.0, 3.0)
    qe = (0.0, 0.0)
    qf = (-0.56, 1.53)
    qg = (2.22, 1.43)
    qh = (2.46, -0.57)
    qi = (-0.68, -1.07)
    pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts; trim=false)
    cents, tri_to_idx, idx_to_tri = DT.circumcenters(T, pts)
    i = 5
    cell_idx = DT.get_voronoi_vertices(T, adj, adj2v, DG, i, tri_to_idx)
    slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx, adj, adj2v)
    shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
    circshift!(slow_cell_idx, shift_idx - 1)
    @test cell_idx == slow_cell_idx

    cell_idx = DT.get_voronoi_vertices(T, adj, adj2v, DG, i, tri_to_idx)
    slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx, adj, adj2v)
    shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
    circshift!(slow_cell_idx, shift_idx - 1)
    @test cell_idx == slow_cell_idx

    for i in eachindex(pts)
        cell_idx = DT.get_voronoi_vertices(T, adj, adj2v, DG, i, tri_to_idx)
        slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx, adj, adj2v)
        nonzero_cell_idx = findfirst(cell_idx .≠ DT.BoundaryIndex)
        circshift!(cell_idx, nonzero_cell_idx - 1)
        shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
        circshift!(slow_cell_idx, shift_idx - 1)
        @test cell_idx == slow_cell_idx
    end
end

@testset "Voronoi tessellations" begin
    Random.seed!(2992991)
    qa = (-1.4, 0.61)
    qb = (1.0, 1.0)
    qc = (0.0, -2.0)
    qd = (1.0, 3.0)
    qe = (0.0, 0.0)
    qf = (-0.56, 1.53)
    qg = (2.22, 1.43)
    qh = (2.46, -0.57)
    qi = (-0.68, -1.07)
    pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    vorn = DT.voronoi(T, adj, adj2v, DG, pts)
    for (i, v) in vorn.idx_to_triangle
        @test vorn.triangle_to_idx[v] == i
    end
    for V in T
        i, j, k = indices(V)
        pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
        true_centxy, rad = circle_three_points(pᵢ, pⱼ, pₖ)
        true_centx, true_centy = true_centxy
        centx, centy = vorn.circumcenters[vorn.triangle_to_idx[V]]
        @test centx ≈ true_centx
        @test centy ≈ true_centy
    end
    for (i, poly) in vorn.polygons
        if DT.BoundaryIndex ∉ poly
            _pts = vorn.circumcenters[poly]
            A = DT.area(_pts)
            @test A ≥ 0.0
        end
    end

    for _ in 1:500
        pts = rand(SVector{2,Float64}, rand(3:1000))
        T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
        vorn = DT.voronoi(T, adj, adj2v, DG, pts)
        for (i, v) in vorn.idx_to_triangle
            @test vorn.triangle_to_idx[v] == i
        end
        for V in T
            i, j, k = indices(V)
            pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
            true_centxy, rad = circle_three_points(pᵢ, pⱼ, pₖ)
            true_centx, true_centy = true_centxy
            centx, centy = vorn.circumcenters[vorn.triangle_to_idx[V]]
            @test centx ≈ true_centx
            @test centy ≈ true_centy
        end
        for (i, poly) in vorn.polygons
            if DT.BoundaryIndex ∉ poly
                _pts = vorn.circumcenters[poly]
                A = DT.area(_pts)
                @test A ≥ 0.0
            end
        end
    end
end

@testset "Areas of Voronoi cells" begin
    for _ in 1:500
        pts = rand(SVector{2,Float64}, rand(3:500))
        T, adj, adj2v, DG = triangulate_bowyer(pts)
        vorn = voronoi(T, adj, adj2v, DG, pts)
        areas = DT.area(vorn)
        for (i, idx) in vorn.polygons
            if DT.BoundaryIndex ∈ idx
                @test areas[i] == Inf
            else
                _area = DT.area(vorn.circumcenters[idx])
                @test areas[i] ≈ _area
            end
            @test areas[i] ≈ DT.area(vorn, i)
        end
    end
end