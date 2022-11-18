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

    @test DT.number_type(vorn) == Float64
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
    #=
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
            next_j_idx = _j_idx == lastindex(bounded_polygon) ? firstindex(bounded_polygon) : _j_idx + 1
            prev_j_idx = _j_idx == firstindex(bounded_polygon) ? lastindex(bounded_polygon) : _j_idx - 1
            next_idx = _j_idx == lastindex(bounded_polygon) ? first(bounded_polygon) : bounded_polygon[next_j_idx]
            prev_idx = _j_idx == firstindex(bounded_polygon) ? last(bounded_polygon) : bounded_polygon[prev_j_idx]
            pᵤ, pᵥ = get_point(vorn.generators, u, v)
            pᵣ = vorn.circumcenters[next_idx]
            pₛ = vorn.circumcenters[prev_idx]
            pⱼ = vorn.circumcenters[_j]
            @test ExactPredicates.meet(pᵤ, pᵥ, pⱼ, pᵣ) == ExactPredicates.meet(pᵤ, pᵥ,pⱼ,pₛ) == 1
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
    =#
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

    @test DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG) == Dict{Int64,NTuple{2,Int64}}([])
end

@testset "Adding a polygon" begin
    dict = Dict{Int64,Set{Int64}}()
    DT.add_polygon!(dict, 2, 1)
    @test dict[2] == Set{Int64}([1])
end

@testset "Finding neighbouring circumcenters" begin
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

    ngh = DT.get_neighbouring_circumcenters(vorn, 1, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pv, pw) == 1
    @test ngh == (12, 8, 11)

    ngh = DT.get_neighbouring_circumcenters(vorn, 2, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pw, pv) == 1 # orient order switched because BoundaryIndex was mapped to the centroid
    @test ngh == (4, 12, DT.BoundaryIndex)

    ngh = DT.get_neighbouring_circumcenters(vorn, 3, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pw, pv) == 1
    @test ngh == (14, 16, DT.BoundaryIndex)

    ngh = DT.get_neighbouring_circumcenters(vorn, 4, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pv, pw) == 1
    @test ngh == (2, 10, 12)

    ngh = DT.get_neighbouring_circumcenters(vorn, 5, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pv, pw) == 1
    @test ngh == (14, 11, 8)

    ngh = DT.get_neighbouring_circumcenters(vorn, 8, adj) # 6-7 skipped because those are ghost triangles
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pv, pw) == 1
    @test ngh == (16, 5, 1)

    ngh = DT.get_neighbouring_circumcenters(vorn, 10, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pw, pv) == 1
    @test ngh == (16, 4, DT.BoundaryIndex)

    ngh = DT.get_neighbouring_circumcenters(vorn, 11, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pw, pv) == 1
    @test ngh == (1, 5, DT.BoundaryIndex)

    ngh = DT.get_neighbouring_circumcenters(vorn, 12, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pv, pw) == 1
    @test ngh == (4, 1, 2)

    ngh = DT.get_neighbouring_circumcenters(vorn, 14, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pw, pv) == 1
    @test ngh == (5, 3, DT.BoundaryIndex)

    ngh = DT.get_neighbouring_circumcenters(vorn, 16, adj)
    pu, pv, pw = get_point(vorn.circumcenters, ngh...)
    @test DT.orient(pu, pv, pw) == 1
    @test ngh == (3, 8, 10)
end
#=
Random.seed!(29291919919)
fig = Figure()
ax = Axis(fig[1, 1])
mkp = voronoiplot!(ax, vorn, pts, _DG; markersize=0, strokewidth=0.5)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=0)
xlims!(ax, -4, 4)
ylims!(ax, -4, 4)
fig
=#
@testset "Intersection coordinates with only exterior circumcenters on unbounded polygons" begin
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
    exterior_circumcenters = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
    @test exterior_circumcenters == Dict{Int64,NTuple{2,Int64}}(2 => (3, 1), 10 => (8, 3), 11 => (1, 4))
    intersection_coordinates = DT.find_intersections_of_exterior_circumcenters_with_convex_hull(vorn, exterior_circumcenters, adj)
    p21, p22 = intersection_coordinates[2]
    p101, p102 = intersection_coordinates[10]
    p111, p112 = intersection_coordinates[11]
    p2, p10, p11 = get_point(vorn.circumcenters, 2, 10, 11)
    scatter!(ax, [p2, p21, p22], color=:red, markersize=14)
    scatter!(ax, [p10, p101, p102], color=:blue, markersize=14)
    scatter!(ax, [p11, p111, p112], color=:black, markersize=14)
    @test collect(p21) ≈ [-0.5914476760497487, -1.2856060039278114]
    @test collect(p22) ≈ [-1.041442472810534, -0.23061820263308141]
    @test collect(p101) ≈ [1.6906445661230798, -0.7808490060707061]
    @test collect(p102) ≈ [0.7595419847328246, -1.0865139949109415]
    @test collect(p111) ≈ [-1.026050660415118, 0.9823912173366114]
    @test collect(p112) ≈ [0.16233399796847348, 2.165824272976938]
    @test collect(p2) ≈ [-125.79166666666887, -53.695000000000924]
    @test collect(p10) ≈ [1.6557088846880907, -1.3852362948960304]
    @test collect(p11) ≈ [-5.73756928406467, 7.365739030023101]
end

@testset "Intersection coordinates with different intersection types" begin
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
    pts = rand(SVector{2,Float64}, 20)
    pts = Tuple.(pts)

    T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
    BN = convex_hull(_DG, pts)
    vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=false)

    exterior_circumcenters = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
    for j in keys(exterior_circumcenters)
        @test DT.isinconvexhull(pts, BN, vorn.circumcenters[j]) == -1
    end
    for j in setdiff(eachindex(vorn.circumcenters), keys(exterior_circumcenters))
        if !DT.is_ghost_triangle(vorn.idx_to_triangle[j]...)
            @test DT.isinconvexhull(pts, BN, vorn.circumcenters[j]) == 1
        end
    end
    j = collect(keys(exterior_circumcenters))
    @test j == [22, 33, 36, 9, 28, 8, 17]
    @test all(collect.(vorn.circumcenters[j]) .≈ [
        [0.5787302935270776, 1.2486600493531892],
        [0.5183479216513563, 2.2974283912607834],
        [0.9689117646555506, 0.8976718879670816],
        [-1.3594377951576813, -0.054169171986592435],
        [-0.12054429369233605, 0.7199135789107467],
        [1.9504373949703273, 0.5939190607059536],
        [2.989644794529847, 0.618812514220584]
    ])

    intersection_coordinates = DT.find_intersections_of_exterior_circumcenters_with_convex_hull(vorn, exterior_circumcenters, adj)
    @test length(intersection_coordinates) == 7
    @test collect.(collect(intersection_coordinates[22]))[[1, 2]] ≈ collect.(collect(((0.5599621186008643, 0.9867087257384604), (0.8304552790307996, 0.98804698791469), (NaN, NaN))))[[1, 2]]
    @test collect.(collect(intersection_coordinates[33]))[[1]] ≈ collect.(collect(((0.21857838170342844, 0.9850197331519431), (NaN, NaN), (NaN, NaN))))[[1]]
    @test collect.(collect(intersection_coordinates[36]))[[2, 3]] ≈ collect.(collect(((NaN, NaN), (0.9036487484104713, 0.9133790392147804), (0.9027980927121162, 0.8541809705789867))))[[2, 3]]
    @test collect.(collect(intersection_coordinates[9]))[[1, 2]] ≈ collect.(collect(((0.10517725105670811, 0.04678774762720798), (0.057615991175906774, 0.27704213616390855), (NaN, NaN))))[[1, 2]]
    @test collect.(collect(intersection_coordinates[28]))[[1, 2]] ≈ collect.(collect(((0.016653573001057704, 0.5840093935806543), (0.015305274150334884, 0.7606095401651458), (NaN, NaN))))[[1, 2]]
    @test collect.(collect(intersection_coordinates[8]))[[1, 2]] ≈ collect.(collect(((0.8991788214585335, 0.6023118906498992), (0.8954849455625512, 0.3452509913397629), (NaN, NaN))))[[1, 2]]
    @test all(isnan.((intersection_coordinates[22][3])))
    @test all(isnan.(collect(Iterators.flatten(intersection_coordinates[33][[2, 3]]))))
    @test all(isnan.(intersection_coordinates[36][1]))
    @test all(isnan.(intersection_coordinates[9][3]))
    @test all(isnan.(intersection_coordinates[28][3]))
    @test all(isnan.(intersection_coordinates[8][3]))
    @test all(isnan.(collect(Iterators.flatten(intersection_coordinates[17]))))
end

@testset "Testing that we can extract the correct list of modified polygons" begin
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
    exterior_circumcenters = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
    @test exterior_circumcenters == Dict{Int64,NTuple{2,Int64}}(2 => (3, 1), 10 => (8, 3), 11 => (1, 4))
    intersection_coordinates = DT.find_intersections_of_exterior_circumcenters_with_convex_hull(vorn, exterior_circumcenters, adj)
    exterior_circumcenters
    @test DT.get_modified_polygons(exterior_circumcenters) == Set{Int64}([1, 3, 4, 8])

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
    pts = rand(SVector{2,Float64}, 20)
    pts = Tuple.(pts)
    T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
    BN = convex_hull(_DG, pts)
    vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=false)
    exterior_circumcenters = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
    @test DT.get_modified_polygons(exterior_circumcenters) == Set{Int64}([14, 16, 3, 6, 11, 17])
end

@testset "Moving the exterior circumcenter to the intersection coordinates for a simple case" begin
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
    exterior_circumcenters = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
    @test exterior_circumcenters == Dict{Int64,NTuple{2,Int64}}(2 => (3, 1), 10 => (8, 3), 11 => (1, 4))
    intersection_coordinates = DT.find_intersections_of_exterior_circumcenters_with_convex_hull(vorn, exterior_circumcenters, adj)
    for j in keys(exterior_circumcenters)
        DT.move_exterior_circumcenter_to_convex_hull_intersection!(vorn, exterior_circumcenters, intersection_coordinates, j)
    end
    @test vorn.polygons[1] == [18, 12, 1, 21, DT.BoundaryIndex, DT.BoundaryIndex]
    @test vorn.polygons[3] == [DT.BoundaryIndex, 20, 4, 17, DT.BoundaryIndex]
    @test vorn.polygons[4] == [14, DT.BoundaryIndex, DT.BoundaryIndex, 22, 5]
    @test vorn.polygons[8] == [DT.BoundaryIndex, 3, 16, 19, DT.BoundaryIndex]
    @test collect.(vorn.circumcenters) ≈ collect.([
        (-0.6354382787597795, 0.4531744421906694)
        (-125.79166666666887, -53.695000000000924)
        (2.0368284558648586, 0.0039983345229595525)
        (0.272216494845361, -0.9240721649484538)
        (0.35048076923076915, 2.0)
        (-0.95, -0.445)
        (1.48, -0.85)
        (0.1769553072625698, 0.8230446927374302)
        (-0.19999999999999996, 1.805)
        (1.6557088846880907, -1.3852362948960304)
        (-5.73756928406467, 7.365739030023101)
        (-0.9085633626097867, -0.1736700125470515)
        (2.84, 0.615)
        (1.3333196721311475, 2.0)
        (1.61, 2.215)
        (1.695573770491803, -0.6955737704918031)
        (-0.5914476760497487, -1.2856060039278114)
        (-1.041442472810534, -0.23061820263308141)
        (1.6906445661230798, -0.7808490060707061)
        (0.7595419847328246, -1.0865139949109415)
        (-1.026050660415118, 0.9823912173366114)
        (0.16233399796847348, 2.165824272976938)
    ])

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
    exterior_circumcenters = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
    @test exterior_circumcenters == Dict{Int64,NTuple{2,Int64}}(2 => (3, 1), 10 => (8, 3), 11 => (1, 4))
    intersection_coordinates = DT.find_intersections_of_exterior_circumcenters_with_convex_hull(vorn, exterior_circumcenters, adj)
    num_ext = DT.move_exterior_circumcenter_to_convex_hull_intersection!(vorn, exterior_circumcenters, intersection_coordinates)
    @test vorn.polygons[1] == [18, 12, 1, 21, DT.BoundaryIndex, DT.BoundaryIndex]
    @test vorn.polygons[3] == [DT.BoundaryIndex, 20, 4, 17, DT.BoundaryIndex]
    @test vorn.polygons[4] == [14, DT.BoundaryIndex, DT.BoundaryIndex, 22, 5]
    @test vorn.polygons[8] == [DT.BoundaryIndex, 3, 16, 19, DT.BoundaryIndex]
    @test collect.(vorn.circumcenters) ≈ collect.([
        (-0.6354382787597795, 0.4531744421906694)
        (-125.79166666666887, -53.695000000000924)
        (2.0368284558648586, 0.0039983345229595525)
        (0.272216494845361, -0.9240721649484538)
        (0.35048076923076915, 2.0)
        (-0.95, -0.445)
        (1.48, -0.85)
        (0.1769553072625698, 0.8230446927374302)
        (-0.19999999999999996, 1.805)
        (1.6557088846880907, -1.3852362948960304)
        (-5.73756928406467, 7.365739030023101)
        (-0.9085633626097867, -0.1736700125470515)
        (2.84, 0.615)
        (1.3333196721311475, 2.0)
        (1.61, 2.215)
        (1.695573770491803, -0.6955737704918031)
        (-0.5914476760497487, -1.2856060039278114)
        (-1.041442472810534, -0.23061820263308141)
        (1.6906445661230798, -0.7808490060707061)
        (0.7595419847328246, -1.0865139949109415)
        (-1.026050660415118, 0.9823912173366114)
        (0.16233399796847348, 2.165824272976938)
    ])
    @test num_ext[4] == 1
    @test num_ext[8] == 1
    @test num_ext[3] == 2
    @test num_ext[1] == 2
end
#=
Random.seed!(29291919919)
fig = Figure()
ax = Axis(fig[1, 1])
mkp = voronoiplot!(ax, vorn, pts, _DG; markersize=0, strokewidth=0.5)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=0)
xlims!(ax, -2.0, 4.0)
ylims!(ax, -2.0, 4.0)
fig
=#

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
exterior_circumcenters = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)
@test exterior_circumcenters == Dict{Int64,NTuple{2,Int64}}(2 => (3, 1), 10 => (8, 3), 11 => (1, 4))
intersection_coordinates = DT.find_intersections_of_exterior_circumcenters_with_convex_hull(vorn, exterior_circumcenters, adj)
num_ext = DT.move_exterior_circumcenter_to_convex_hull_intersection!(vorn, exterior_circumcenters, intersection_coordinates)



added_generator_indices = Dict{Int64,Int64}()
modified_polygons = DT.get_modified_polygons(exterior_circumcenters)
DT.close_polygon_through_generator!(vorn, added_generator_indices, 8)


for (u, pu) in added_generator_indices
    @test pts[u] == vorn.circumcenters[pu]
end