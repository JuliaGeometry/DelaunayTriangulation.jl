@testset "A big first example" begin
    @testset "Cluttered" begin
        p0 = Float64[5, 5]
        p1 = Float64[4.5, 2.5]
        p2 = Float64[2.5, 1.5]
        p3 = Float64[3, 3.5]
        p4 = Float64[0, 2]
        p5 = Float64[1, 5]
        p6 = Float64[1, 3]
        p7 = Float64[4, -1]
        p8 = Float64[-1, 4]
        pts = [p0, p1, p2, p3, p4, p5, p6, p7, p8]

        # Initialisation 
        𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, root = DT.initialise_triangulation()

        # Add a point 
        r = 2
        pᵣ = pts[r]
        Tᵢⱼₖ, flag = DT.locate_triangle(𝒟, pts, pᵣ, root)
        _Tᵢⱼₖ, _flag = DT.locate_triangle(𝒟, pts, pᵣ)
        @test Tᵢⱼₖ == _Tᵢⱼₖ
        @test flag == _flag
        i, j, k = Tᵢⱼₖ
        𝒟𝒟 = deepcopy(𝒟)
        𝒜𝒜 = deepcopy(𝒜)
        𝒜⁻¹𝒜⁻¹ = deepcopy(𝒜⁻¹)
        𝒯𝒯 = deepcopy(𝒯)
        𝒱𝒱 = deepcopy(𝒟𝒢)
        push!(𝒯.triangles, DT.TriangleType((i, j, r)), DT.TriangleType((j, k, r)), DT.TriangleType((k, i, r)))
        delete!(𝒯, Tᵢⱼₖ)
        DT.add_triangle!(𝒟, DT.TriangleType((i, j, r)))
        DT.add_triangle!(𝒟, DT.TriangleType((j, k, r)))
        DT.add_triangle!(𝒟, DT.TriangleType((k, i, r)))
        DT.add_edge!(𝒟, DT.TriangleType((i, j, k)), DT.TriangleType((i, j, r)))
        DT.add_edge!(𝒟, DT.TriangleType((i, j, k)), DT.TriangleType((j, k, r)))
        DT.add_edge!(𝒟, DT.TriangleType((i, j, k)), DT.TriangleType((k, i, r)))
        DT.update_adjacent!(𝒜, DT.TriangleType((i, j, r)))
        DT.update_adjacent!(𝒜, DT.TriangleType((j, k, r)))
        DT.update_adjacent!(𝒜, DT.TriangleType((k, i, r)))
        DT.add_neighbour!(𝒟𝒢, 1, 2)
        DT.add_neighbour!(𝒟𝒢, -1, 2)
        DT.add_neighbour!(𝒟𝒢, 0, 2)
        DT.add_neighbour!(𝒟𝒢, 2, 1, -1, 0)
        DT.add_point!(𝒯𝒯, 𝒟𝒟, 𝒜𝒜, 𝒜⁻¹𝒜⁻¹, 𝒱𝒱, Tᵢⱼₖ, r)
        @test 𝒯.triangles == 𝒯𝒯.triangles
        @test 𝒟.graph == 𝒟𝒟.graph
        @test 𝒜.adjacent == 𝒜𝒜.adjacent
        @test DT.graph(𝒟𝒢) == DT.graph(𝒱𝒱)
        @test collect(𝒟.graph.N[(1, -1, 0)]) == DT.TriangleType[(1, -1, 2), (-1, 0, 2), (0, 1, 2)]
        @test collect(𝒟.graph.NN[(1, -1, 2)]) == DT.TriangleType[(1, -1, 0)]
        @test collect(𝒟.graph.NN[(-1, 0, 2)]) == DT.TriangleType[(1, -1, 0)]
        @test collect(𝒟.graph.NN[(0, 1, 2)]) == DT.TriangleType[(1, -1, 0)]
        @test !𝒟.graph.looped
        @test collect(𝒟.graph.V) == DT.TriangleType[(1, -1, 0), (1, -1, 2), (-1, 0, 2), (0, 1, 2)]
        @test 𝒜⁻¹𝒜⁻¹[0] == Set([(1, 2), (2, -1)])
        @test 𝒜⁻¹𝒜⁻¹[1] == Set([(2, 0), (-1, 2)])
        @test 𝒜⁻¹𝒜⁻¹[-1] == Set([(0, 2), (2, 1)])
        @test 𝒜⁻¹𝒜⁻¹[2] == Set([(1, -1), (0, 1), (-1, 0)])
        a1, a2, a3 = collect(𝒯.triangles)[1]
        b1, b2, b3 = collect(𝒯.triangles)[2]
        c1, c2, c3 = collect(𝒯.triangles)[3]
        @test DT.TriangleType((1, -1, 0)) ∉ 𝒯.triangles
        @test (a1, a2, a3) == (1, -1, 2)
        @test (b1, b2, b3) == (-1, 0, 2)
        @test (c1, c2, c3) == (0, 1, 2)
        @test 𝒜(-1, 1) == DT.BoundaryIdx
        @test 𝒜(0, -1) == DT.BoundaryIdx
        @test 𝒜(1, 0) == DT.BoundaryIdx
        @test 𝒜(1, -1) == 2
        @test 𝒜(-1, 0) == 2
        @test 𝒜(0, 1) == 2
        @test 𝒜(-1, 2) == 1
        @test 𝒜(2, -1) == 0
        @test 𝒜(2, 1) == -1
        @test 𝒜(1, 2) == 0
        @test 𝒜(2, 0) == 1
        @test 𝒜(0, 2) == -1
        @test length(𝒜.adjacent) == 12
        @test 𝒟𝒢(1) == Set([-1, 0, 2])
        @test 𝒟𝒢(-1) == Set([0, 1, 2])
        @test 𝒟𝒢(0) == Set([1, -1, 2])
        @test 𝒟𝒢(2) == Set([1, -1, 0])
        𝒜𝒜 = deepcopy(𝒜)
        𝒜⁻¹b = deepcopy(𝒜⁻¹𝒜⁻¹)
        𝒜⁻¹𝒜⁻¹ = deepcopy(𝒜⁻¹)
        𝒟𝒟 = deepcopy(𝒟)
        𝒯𝒯 = deepcopy(𝒯)
        𝒱𝒱 = deepcopy(𝒟𝒢)
        DT.legalise_edge!(𝒯𝒯, 𝒟𝒟, 𝒜𝒜, 𝒜⁻¹𝒜⁻¹, 𝒱𝒱, i, j, r, pts)
        DT.legalise_edge!(𝒯𝒯, 𝒟𝒟, 𝒜𝒜, 𝒜⁻¹𝒜⁻¹, 𝒱𝒱, j, k, r, pts)
        DT.legalise_edge!(𝒯𝒯, 𝒟𝒟, 𝒜𝒜, 𝒜⁻¹𝒜⁻¹, 𝒱𝒱, k, i, r, pts)
        @test 𝒯.triangles == 𝒯𝒯.triangles
        @test 𝒟.graph == 𝒟𝒟.graph
        @test 𝒜.adjacent == 𝒜𝒜.adjacent
        @test 𝒜⁻¹.adjacent2vertex == 𝒜⁻¹𝒜⁻¹.adjacent2vertex
        @test DT.graph(𝒟𝒢) == DT.graph(𝒱𝒱)
        𝒜⁻¹ = deepcopy(𝒜⁻¹b)

        # Add the next point
        r = 3
        pᵣ = pts[r]
        Tᵢⱼₖ, flag = DT.locate_triangle(𝒟, pts, pᵣ, root)
        _Tᵢⱼₖ, _flag = DT.locate_triangle(𝒟, pts, pᵣ)
        @test Tᵢⱼₖ == _Tᵢⱼₖ
        @test flag == _flag
        @test flag == 1
        @test Tᵢⱼₖ == DT.TriangleType((DT.LargeLeftIdx, DT.LargeRightIdx, 2))
        DT.add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, Tᵢⱼₖ, r)
        @test length(𝒜.adjacent) == 18
        @test 𝒜(-1, 1) == DT.BoundaryIdx
        @test 𝒜(1, -1) == 2
        @test 𝒜(-1, 0) == 3
        @test 𝒜(0, -1) == DT.BoundaryIdx
        @test 𝒜(0, 1) == 2
        @test 𝒜(1, 0) == DT.BoundaryIdx
        @test 𝒜(-1, 2) == 1
        @test 𝒜(2, -1) == 3
        @test 𝒜(2, 1) == -1
        @test 𝒜(1, 2) == 0
        @test 𝒜(2, 0) == 1
        @test 𝒜(0, 2) == 3
        @test 𝒜(2, 3) == 0
        @test 𝒜(3, 2) == -1
        @test 𝒜(3, -1) == 0
        @test 𝒜(-1, 3) == 2
        @test 𝒜(3, 0) == 2
        @test 𝒜(0, 3) == -1
        @test collect(𝒟.graph.N[(1, -1, 0)]) == DT.TriangleType[(1, -1, 2), (-1, 0, 2), (0, 1, 2)]
        @test collect(𝒟.graph.N[(-1, 0, 2)]) == DT.TriangleType[(-1, 0, 3), (0, 2, 3), (2, -1, 3)]
        @test collect(𝒟.graph.NN[(1, -1, 0)]) == []
        @test collect(𝒟.graph.NN[(1, -1, 2)]) == DT.TriangleType[(1, -1, 0)]
        @test collect(𝒟.graph.NN[(-1, 0, 2)]) == DT.TriangleType[(1, -1, 0)]
        @test collect(𝒟.graph.NN[(-1, 0, 3)]) == DT.TriangleType[(-1, 0, 2)]
        @test collect(𝒟.graph.NN[(0, 1, 2)]) == DT.TriangleType[(1, -1, 0)]
        @test collect(𝒟.graph.NN[(0, 2, 3)]) == DT.TriangleType[(-1, 0, 2)]
        @test collect(𝒟.graph.NN[(2, -1, 3)]) == DT.TriangleType[(-1, 0, 2)]
        @test 𝒜⁻¹[0] == Set([(1, 2), (2, 3), (3, -1)])
        @test length(𝒜⁻¹.adjacent2vertex) == 5
        @test 𝒜⁻¹[1] == Set([(-1, 2), (2, 0)])
        @test 𝒜⁻¹[2] == Set([(-1, 3), (0, 1), (3, 0), (1, -1)])
        @test 𝒜⁻¹[3] == Set([(-1, 0), (0, 2), (2, -1)])
        i, j, k = Tᵢⱼₖ
        @test DT.is_legal(i, j, 𝒜, pts)
        𝒜𝒜 = deepcopy(𝒜)
        𝒟𝒟 = deepcopy(𝒟)
        𝒯𝒯 = deepcopy(𝒯)
        DT.legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, i, j, r, pts)
        @test 𝒯.triangles == 𝒯𝒯.triangles
        @test 𝒟.graph == 𝒟𝒟.graph
        @test 𝒜.adjacent == 𝒜𝒜.adjacent
        @test !DT.is_legal(j, k, 𝒜, pts)
        @test 𝒜(1, 2) == DT.LargeRightIdx
        @test 𝒜(2, 1) == DT.LargeLeftIdx
        @test 𝒜(2, 3) == DT.LargeRightIdx
        @test 𝒜(3, 2) == DT.LargeLeftIdx
        @test 𝒜(DT.LargeLeftIdx, 3) == 2
        @test 𝒜(1, DT.LargeLeftIdx) == 2
        @test 𝒟𝒢(0) == Set([1, -1, 2, 3])
        @test 𝒟𝒢(1) == Set([-1, 0, 2])
        @test 𝒟𝒢(2) == Set([1, -1, 0, 3])
        @test 𝒟𝒢(3) == Set([-1, 0, 2])
        @test DT.is_legal(i, j, 𝒜, pts)
        𝒟𝒟 = deepcopy(𝒟)
        𝒜𝒜 = deepcopy(𝒜)
        𝒜⁻¹𝒜⁻¹ = deepcopy(𝒜⁻¹)
        𝒯𝒯 = deepcopy(𝒯)
        𝒱𝒱 = deepcopy(𝒟𝒢)
        DT.legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, i, j, r, pts)
        @test 𝒟𝒟.graph == 𝒟.graph
        @test 𝒜𝒜.adjacent == 𝒜.adjacent
        @test 𝒜⁻¹𝒜⁻¹.adjacent2vertex == 𝒜⁻¹.adjacent2vertex
        @test 𝒯𝒯.triangles == 𝒯.triangles
        @test DT.graph(𝒱𝒱) == DT.graph(𝒟𝒢)
        @test !DT.is_legal(j, k, 𝒜, pts)
        _i = 𝒜(k, j)
        DT.legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, j, k, r, pts)
        @test 𝒟𝒟.graph ≠ 𝒟.graph
        @test 𝒜𝒜.adjacent ≠ 𝒜.adjacent
        @test 𝒜⁻¹𝒜⁻¹.adjacent2vertex ≠ 𝒜⁻¹.adjacent2vertex
        @test 𝒯𝒯.triangles ≠ 𝒯.triangles
        @test DT.graph(𝒟𝒢) ≠ DT.graph(𝒱𝒱)
        @test DT.is_legal(_i, r, 𝒜, pts)
    end

    @testset "Cleaner run" begin
        p0 = Float64[5, 5]
        p1 = Float64[4.5, 2.5]
        p2 = Float64[2.5, 1.5]
        p3 = Float64[3, 3.5]
        p4 = Float64[0, 2]
        p5 = Float64[1, 5]
        p6 = Float64[1, 3]
        p7 = Float64[4, -1]
        p8 = Float64[-1, 4]
        pts = [p0, p1, p2, p3, p4, p5, p6, p7, p8]

        # Initialise 
        𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, root = DT.initialise_triangulation()

        # Add the second point 
        r = 2
        pᵣ = pts[r]
        𝒯ᵢⱼₖ, interior_flag = DT.locate_triangle(𝒟, pts, pᵣ, root)
        i, j, k = 𝒯ᵢⱼₖ
        DT.add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, 𝒯ᵢⱼₖ, r)
        @test DT.is_legal(i, j, 𝒜, pts)
        @test DT.is_legal(j, k, 𝒜, pts)
        @test DT.is_legal(k, i, 𝒜, pts)
        DT.legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, i, j, r, pts)
        DT.legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, j, k, r, pts)
        DT.legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, k, i, r, pts)
        @test 𝒯.triangles == Set(DT.TriangleType[(1, -1, 2), (-1, 0, 2), (0, 1, 2)])

        # Add the third point 
        r = 3
        pᵣ = pts[r]
        𝒯ᵢⱼₖ, interior_flag = DT.locate_triangle(𝒟, pts, pᵣ, root)
        i, j, k = 𝒯ᵢⱼₖ
        DT.add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, 𝒯ᵢⱼₖ, r)
        @test 𝒯.triangles == Set(DT.TriangleType[(1, -1, 2), (0, 1, 2), (-1, 0, 3), (0, 2, 3), (2, -1, 3)])
        @test 𝒜(1, -1) == 2
        @test 𝒜(-1, 2) == 1
        @test 𝒜(2, 1) == -1
        @test 𝒜(0, 1) == 2
        @test 𝒜(1, 2) == 0
        @test 𝒜(2, 0) == 1
        @test 𝒜(-1, 0) == 3
        @test 𝒜(0, 3) == -1
        @test 𝒜(3, -1) == 0
        @test 𝒜(0, 2) == 3
        @test 𝒜(2, 3) == 0
        @test 𝒜(3, 0) == 2
        @test 𝒜(2, -1) == 3
        @test 𝒜(-1, 3) == 2
        @test 𝒜(3, 2) == -1
        @test 𝒟𝒢(-1) == Set([0, 1, 2, 3])
        @test 𝒟𝒢(0) == Set([-1, 1, 2, 3])
        @test 𝒟𝒢(1) == Set([-1, 0, 2])
        @test 𝒟𝒢(2) == Set([-1, 0, 1, 3])
        @test 𝒟𝒢(3) == Set([-1, 0, 2])
        @test DT.is_legal(i, j, 𝒜, pts)
        DT.legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, i, j, r, pts)
        @test 𝒯.triangles == Set(DT.TriangleType[(1, -1, 2), (0, 1, 2), (-1, 0, 3), (0, 2, 3), (2, -1, 3)])
        @test !DT.is_legal(j, k, 𝒜, pts)
        @test 𝒜(j, k) == r
        @test 𝒜(k, j) == 1
        _i, _j, _k = j, k, 1
        DT.flip_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, _i, _j, _k, r)
        @test 𝒯.triangles == Set(DT.TriangleType[(1, -1, 2), (-1, 0, 3), (2, -1, 3), (3, 1, 2), (3, 0, 1)])
        @test DT.is_legal(_i, _k, 𝒜, pts)
        @test DT.is_legal(_k, _j, 𝒜, pts)
        @test DT.is_legal(j, 1, 𝒜, pts)
        @test DT.is_legal(1, k, 𝒜, pts)
        @test 𝒜(1, -1) == 2
        @test 𝒜(-1, 2) == 1
        @test 𝒜(2, 1) == -1
        @test 𝒜(-1, 0) == 3
        @test 𝒜(0, 3) == -1
        @test 𝒜(3, -1) == 0
        @test 𝒜(2, -1) == 3
        @test 𝒜(-1, 3) == 2
        @test 𝒜(3, 2) == -1
        @test 𝒜(3, 1) == 2
        @test 𝒜(1, 2) == 3
        @test 𝒜(2, 3) == 1
        @test 𝒜(3, 0) == 1
        @test 𝒜(0, 1) == 3
        @test 𝒜(1, 3) == 0
    end
end

## Another example 
p1 = [7.1693762282825, 19.9433794956457]
p2 = [13.407381993135, 17.4563314456065]
p3 = [4.7236716204664, 6.5423531089996]
pts = [p1, p2, p3]
𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, root = DT.initialise_triangulation()

r = 2
DT.add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, root, pts, r)
@test 𝒯.triangles == Set([(1, -1, 2), (-1, 0, 2), (0, 1, 2)])
@test 𝒜(1, -1) == 2
@test 𝒜(-1, 2) == 1
@test 𝒜(2, 1) == -1
@test 𝒜(-1, 0) == 2
@test 𝒜(0, 2) == -1
@test 𝒜(2, -1) == 0
@test 𝒜(0, 1) == 2
@test 𝒜(1, 2) == 0
@test 𝒜(2, 0) == 1
@test 𝒜⁻¹[-1] == Set([(2, 1), (0, 2)])
@test 𝒜⁻¹[0] == Set([(2, -1), (1, 2)])
@test 𝒜⁻¹[1] == Set([(-1, 2), (2, 0)])
@test 𝒜⁻¹[2] == Set([(1, -1), (-1, 0), (0, 1)])
@test 𝒟𝒢(-1) == Set([0, 1, 2])
@test 𝒟𝒢(0) == Set([-1, 1, 2])
@test 𝒟𝒢(1) == Set([-1, 0, 2])
@test 𝒟𝒢(2) == Set([-1, 0, 1])
@test DT.graph(𝒟).N[(1, -1, 0)] == Set([(1, -1, 2), (-1, 0, 2), (0, 1, 2)])
@test DT.graph(𝒟).N[(1, -1, 2)] == Set()
@test DT.graph(𝒟).N[(-1, 0, 2)] == Set()
@test DT.graph(𝒟).N[(0, 1, 2)] == Set()
@test DT.graph(𝒟).NN[(1, -1, 0)] == Set()
@test DT.graph(𝒟).NN[(1, -1, 2)] == Set([(1, -1, 0)])
@test DT.graph(𝒟).NN[(-1, 0, 2)] == Set([(1, -1, 0)])
@test DT.graph(𝒟).NN[(0, 1, 2)] == Set([(1, -1, 0)])


