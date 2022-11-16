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
    DT.add_ghost_triangles!(T, adj, adj2v, DG)
    cent_to_idx = Dict{Int64,Set{Int64}}()
    cell_idx = DT.get_voronoi_vertices(T, adj, DG, i, tri_to_idx, cent_to_idx)
    slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx)
    shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
    circshift!(slow_cell_idx, shift_idx - 1)
    @test cell_idx == slow_cell_idx

    cell_idx = DT.get_voronoi_vertices(T, adj, DG, i, tri_to_idx, cent_to_idx)
    slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx)
    shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
    circshift!(slow_cell_idx, shift_idx - 1)
    @test cell_idx == slow_cell_idx

    for i in eachindex(pts)
        cell_idx = DT.get_voronoi_vertices(T, adj, DG, i, tri_to_idx, cent_to_idx)
        slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx)
        nonzero_cell_idx = findfirst(cell_idx .≠ DT.BoundaryIndex)
        circshift!(cell_idx, nonzero_cell_idx - 1)
        shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
        circshift!(slow_cell_idx, shift_idx - 1)
        @test cell_idx == slow_cell_idx
    end
    DT.remove_ghost_triangles!(T, adj, adj2v, DG)

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
    cell_idx = DT.get_voronoi_vertices(T, adj, DG, i, tri_to_idx, cent_to_idx)
    slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx)
    shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
    circshift!(slow_cell_idx, shift_idx - 1)
    @test cell_idx == slow_cell_idx

    cell_idx = DT.get_voronoi_vertices(T, adj, DG, i, tri_to_idx, cent_to_idx)
    slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx)
    shift_idx = findfirst(slow_cell_idx .== cell_idx[1])
    circshift!(slow_cell_idx, shift_idx - 1)
    @test cell_idx == slow_cell_idx

    for i in eachindex(pts)
        cell_idx = DT.get_voronoi_vertices(T, adj, DG, i, tri_to_idx, cent_to_idx)
        slow_cell_idx = slow_get_voronoi_cells(T, DG, i, pts, tri_to_idx)
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
    for (i, v) in vorn.circumcenter_to_polygons
        for w in v
            poly = vorn.polygons[w]
            @test i ∈ poly
        end
    end

    for _ in 1:500
        pts = rand(SVector{2,Float64}, rand(3:5000))
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
        for (i, v) in vorn.circumcenter_to_polygons
            for w in v
                poly = vorn.polygons[w]
                @test i ∈ poly
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

@testset "Accessing triangles and circumcenters" begin
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
    for τ in T
        V = τ ∈ T ? τ : (DT.shift_triangle_1(τ) ∈ T ? DT.shift_triangle_1(τ) : DT.shift_triangle_2(τ))
        j = vorn.triangle_to_idx[V]
        @test DT.get_triangle_idx(vorn, τ) ==
              DT.get_triangle_idx(vorn, DT.shift_triangle_1(τ)) ==
              DT.get_triangle_idx(vorn, DT.shift_triangle_2(τ)) ==
              j
        @test DT.get_circumcenter(vorn, τ) ==
              DT.get_circumcenter(vorn, DT.shift_triangle_1(τ)) ==
              DT.get_circumcenter(vorn, DT.shift_triangle_2(τ)) ==
              DT.get_circumcenter(vorn, j) ==
              vorn.circumcenters[j]
    end

    for _ in 1:250
        pts = rand(SVector{2,Float64}, rand(3:500))
        T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
        vorn = DT.voronoi(T, adj, adj2v, DG, pts)
        for τ in T
            i, j, k = indices(τ)
            V = τ ∈ T ? τ : (DT.shift_triangle_1(τ) ∈ T ? DT.shift_triangle_1(τ) : DT.shift_triangle_2(τ))
            _j = vorn.triangle_to_idx[V]
            @test DT.get_triangle_idx(vorn, τ) ==
                  DT.get_triangle_idx(vorn, DT.shift_triangle_1(τ)) ==
                  DT.get_triangle_idx(vorn, DT.shift_triangle_2(τ)) ==
                  DT.get_triangle_idx(vorn, i, j, k) ==
                  DT.get_triangle_idx(vorn, j, k, i) ==
                  DT.get_triangle_idx(vorn, k, i, j) ==
                  _j
            @test DT.get_circumcenter(vorn, τ) ==
                  DT.get_circumcenter(vorn, DT.shift_triangle_1(τ)) ==
                  DT.get_circumcenter(vorn, DT.shift_triangle_2(τ)) ==
                  DT.get_circumcenter(vorn, i, j, k) ==
                  DT.get_circumcenter(vorn, j, k, i) ==
                  DT.get_circumcenter(vorn, k, i, j) ==
                  DT.get_circumcenter(vorn, _j) ==
                  vorn.circumcenters[_j]
        end
    end

    @test DT.triangle_type(vorn) == NTuple{3,Int64}
end

@testset "Unbounded polygons" begin
    Random.seed!(2992991)
    qa = (-1.4, 0.61)
    qb = (1.0, 1.0)
    qc = (-0.5, -1.5)
    qd = (1.0, 3.0)
    qe = (0.0, 0.0)
    qf = (-0.56, 1.23)
    qg = (2.22, 1.43)
    qh = (3.46, -0.2)
    qi = (-0.68, -1.07)
    pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]

    T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
    vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=false)

    j = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
    j = collect(keys(j))

    _j = j[1]
    polys = vorn.circumcenter_to_polygons[_j]
    polys = collect(polys)
    polys1 = vorn.polygons[9]
    polys2 = vorn.polygons[3]
    polys3 = vorn.polygons[1]
    @test DT.find_bounded_polygons(vorn, _j) == [9]

    _j = j[2]
    polys = vorn.circumcenter_to_polygons[_j]
    polys = collect(polys)
    polys1 = vorn.polygons[5]
    polys2 = vorn.polygons[8]
    polys3 = vorn.polygons[3]
    @test DT.find_bounded_polygons(vorn, _j) == [5]

    _j = j[3]
    polys = vorn.circumcenter_to_polygons[_j]
    polys = collect(polys)
    polys1 = vorn.polygons[4]
    polys2 = vorn.polygons[6]
    polys3 = vorn.polygons[1]
    @test DT.find_bounded_polygons(vorn, _j) == [6]

    @test DT.find_bounded_polygons(vorn, j[1]; find_first=true) == 9
    @test DT.find_bounded_polygons(vorn, j[2]; find_first=true) == 5
    @test DT.find_bounded_polygons(vorn, j[3]; find_first=true) == 6
end

@testset "Points outside the convex hull" begin
    Random.seed!(2992991)
    qa = (-1.4, 0.61)
    qb = (1.0, 1.0)
    qc = (-0.5, -1.5)
    qd = (1.0, 3.0)
    qe = (0.0, 0.0)
    qf = (-0.56, 1.23)
    qg = (2.22, 1.43)
    qh = (3.46, -0.2)
    qi = (-0.68, -1.07)
    pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]

    T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
    vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=false)
    j = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
    @test j == Dict{Int64,NTuple{2,Int64}}(2 => (3, 1), 10 => (8, 3), 11 => (1, 4))

    for _ in 1:670
        pts = rand(SVector{2,Float64}, rand(3:3000))
        T, adj, adj2v, _DG = DT.triangulate_bowyer(pts)
        vorn = DT.voronoi(T, adj, adj2v, _DG, pts)
        DT.add_ghost_triangles!(T, adj, adj2v, _DG)
        j = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
        DT.remove_ghost_triangles!(T, adj, adj2v, _DG)
        BN = convex_hull(_DG, pts)
        for (_j, (u, v)) in j
            @test DT.isinconvexhull(pts, BN, vorn.circumcenters[_j]) == -1
            bounded_polygon_idx = DT.find_bounded_polygons(vorn, _j; find_first=true)
            bounded_polygon = vorn.polygons[bounded_polygon_idx]
            _j_idx = findfirst(bounded_polygon .== _j) # Finds where _j is in the list of vertices for the bounded polygon
            # Find the circumcenters that connect with _j 
            next_j_idx = _j_idx == lastindex(bounded_polygon) ? firstindex(bounded_polygon) : _j_idx + 1
            prev_j_idx = _j_idx == firstindex(bounded_polygon) ? lastindex(bounded_polygon) : _j_idx - 1
            next_idx = _j_idx == lastindex(bounded_polygon) ? first(bounded_polygon) : bounded_polygon[next_j_idx]
            prev_idx = _j_idx == firstindex(bounded_polygon) ? last(bounded_polygon) : bounded_polygon[prev_j_idx]
            # Now get the triangle that is on the boundary, and rotate it so that its first two vertices give the boundary edge 

        end
    end

    for _ in 1:670
        pts = rand(SVector{2,Float64}, rand(3:3000))
        T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
        vorn = DT.voronoi(T, adj, adj2v, _DG, pts)
        j = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
        BN = convex_hull(_DG, pts)
        for _j in j
            @test DT.isinconvexhull(pts, BN, vorn.circumcenters[_j]) == -1
        end
    end
end

@testset "Bounded Voronoi tessellation" begin
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
    T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
    vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=false)
    @test !DT.polygon_is_bounded(vorn, 1)

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

    T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
    vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=true)
    BN = convex_hull(_DG, pts)
    areas = DT.area(vorn)
    for (i, idx) in vorn.polygons
        _area = DT.area(vorn.circumcenters[idx])
        @test areas[i] ≈ _area
        @test areas[i] ≈ DT.area(vorn, i)
    end
    for (i, v) in vorn.circumcenter_to_polygons
        for w in v
            poly = vorn.polygons[w]
            @test i ∈ poly
        end
    end

    @test vorn.polygons[1] == [13, 23, 3, 10, 1]
    @test vorn.polygons[2] == [5, 7, 16, 2, 12]
    @test vorn.polygons[3] == [22, 6, 9, 4, 15]
    @test vorn.polygons[4] == [12, 14, 17, 8, 5]
    @test vorn.polygons[5] == [9, 16, 7, 1, 10, 4]
    @test vorn.polygons[6] == [1, 7, 5, 8, 18, 13]
    @test vorn.polygons[7] == [2, 11, 19, 14, 12]
    @test vorn.polygons[8] == [6, 21, 11, 2, 16, 9]
    @test vorn.polygons[9] == [10, 3, 20, 15, 4]
    @test all(DT.polygon_is_bounded.(Ref(vorn), keys(vorn.polygons)))
    @test vorn.circumcenter_to_polygons[1] == Set{Int64}([1, 5, 6])
    @test vorn.circumcenter_to_polygons[2] == Set{Int64}([7, 2, 8])
    @test vorn.circumcenter_to_polygons[3] == Set{Int64}([1, 9])
    @test vorn.circumcenter_to_polygons[4] == Set{Int64}([3, 5, 9])
    @test vorn.circumcenter_to_polygons[5] == Set{Int64}([2, 6, 4])
    @test vorn.circumcenter_to_polygons[6] == Set{Int64}([3, 8])
    @test vorn.circumcenter_to_polygons[7] == Set{Int64}([2, 5, 6])
    @test vorn.circumcenter_to_polygons[8] == Set{Int64}([4, 6])
    @test vorn.circumcenter_to_polygons[9] == Set{Int64}([3, 5, 8])
    @test vorn.circumcenter_to_polygons[10] == Set{Int64}([1, 5, 9])
    @test vorn.circumcenter_to_polygons[11] == Set{Int64}([7, 8])
    @test vorn.circumcenter_to_polygons[12] == Set{Int64}([2, 4, 7])
    @test vorn.circumcenter_to_polygons[13] == Set{Int64}([1, 6])
    @test vorn.circumcenter_to_polygons[14] == Set{Int64}([4, 7])
    @test vorn.circumcenter_to_polygons[15] == Set{Int64}([3, 9])
    @test vorn.circumcenter_to_polygons[16] == Set{Int64}([5, 2, 8])
    @test vorn.circumcenter_to_polygons[17] == Set{Int64}([4])
    @test vorn.circumcenter_to_polygons[18] == Set{Int64}([6])
    @test vorn.circumcenter_to_polygons[19] == Set{Int64}([7])
    @test vorn.circumcenter_to_polygons[20] == Set{Int64}([9])
    @test vorn.circumcenter_to_polygons[21] == Set{Int64}([8])
    @test vorn.circumcenter_to_polygons[22] == Set{Int64}([3])

    @test DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG) == Set{Int64}([])
end

@testset "Adding a polygon" begin
    dict = Dict{Int64,Set{Int64}}()
    DT.add_polygon!(dict, 2, 1)
    @test dict[2] == Set{Int64}([1])
end