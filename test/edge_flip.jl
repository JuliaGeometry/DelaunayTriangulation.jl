@testset "Can we correctly flip an edge?" begin
    I = [-2.0, -2.0]
    J = [-10.0, 2.0]
    K = [-4.0, 4.0]
    R = [-8.0, -3.0]
    IJKR = [I, J, K, R]
    i, j, k, r = 1, 2, 3, 4
    𝒟n = DT.HistoryDAG()
    Tᵢₖⱼ = DT.TriangleType((i, k, j))
    Tᵢⱼᵣ = DT.TriangleType((i, j, r))
    𝒯n = DT.Triangles(DT.TriangleType[Tᵢₖⱼ, Tᵢⱼᵣ])
    DT.add_triangle!(𝒟n, Tᵢₖⱼ)
    DT.add_triangle!(𝒟n, Tᵢⱼᵣ)
    𝒜n = DT.Adjacent()
    𝒜n[(i, j)] = r
    𝒜n[(j, r)] = i
    𝒜n[(r, i)] = j
    𝒜n[(i, k)] = j
    𝒜n[(k, j)] = i
    𝒜n[(j, i)] = k
    𝒜⁻¹n = DT.Adjacent2Vertex()
    𝒜⁻¹n[i] = Set([(k, j), (j, r)])
    𝒜⁻¹n[j] = Set([(r, i), (i, k)])
    𝒜⁻¹n[r] = Set([(i, j)])
    𝒜⁻¹n[k] = Set([(j, i)])
    𝒟𝒢n = DT.DelaunayGraph()
    DT.add_point!(𝒟𝒢n, i, j, k, r)
    DT.add_edge!(𝒟𝒢n, i, k, j, r)
    DT.add_edge!(𝒟𝒢n, j, r, i, k)
    DT.add_edge!(𝒟𝒢n, k, j, i)
    DT.add_edge!(𝒟𝒢n, r, i, j)
    𝒯𝒯 = deepcopy(𝒯n)
    𝒟𝒟 = deepcopy(𝒟n)
    𝒜𝒜 = deepcopy(𝒜n)
    𝒜⁻¹𝒜⁻¹ = deepcopy(𝒜⁻¹n)
    𝒟𝒢𝒟𝒢 = deepcopy(𝒟𝒢n)
    @test !DT.is_legal(i, j, 𝒜n, IJKR)
    𝒯𝒯n = deepcopy(𝒯n)
    𝒟𝒟n = deepcopy(𝒟n)
    𝒜𝒜n = deepcopy(𝒜n)
    𝒜⁻¹𝒜⁻¹n = deepcopy(𝒜⁻¹n)
    𝒟𝒢𝒟𝒢n = deepcopy(𝒟𝒢n)
    Tᵢₖⱼ = DT.TriangleType((i, k, j))
    Tᵢⱼᵣ = DT.TriangleType((i, j, r))
    DT.delete_triangle!(𝒯𝒯n, Tᵢₖⱼ, Tᵢⱼᵣ)
    @test isempty(𝒯𝒯n.triangles)
    DT.delete_edge!(𝒜𝒜n, i, j)
    @test 𝒜𝒜n(2, 4) == 1
    @test 𝒜𝒜n(3, 2) == 1
    @test 𝒜𝒜n(4, 1) == 2
    @test 𝒜𝒜n(1, 3) == 2
    @test (i, j) ∉ keys(𝒜𝒜n.adjacent)
    @test (j, i) ∉ keys(𝒜𝒜n.adjacent)
    DT.delete_neighbour!(𝒟𝒢𝒟𝒢n, i, j) #delete_neighbour!(𝒟𝒢, j, i)
    @test i ∉ 𝒟𝒢𝒟𝒢n[j]
    @test j ∉ 𝒟𝒢𝒟𝒢n[i]
    @test 𝒟𝒢𝒟𝒢n[j] == Set([4, 3])
    @test 𝒟𝒢𝒟𝒢n[i] == Set([4, 3])
    Tᵣₖⱼ = DT.TriangleType((r, k, j))
    Tᵣᵢₖ = DT.TriangleType((r, i, k))
    DT.add_triangle!(𝒯𝒯n, Tᵣₖⱼ, Tᵣᵢₖ)
    DT.add_triangle!(𝒟𝒟n, Tᵣₖⱼ, Tᵣᵢₖ)
    @test 𝒯𝒯n.triangles == Set([(4, 3, 2), (4, 1, 3)])
    DT.update_adjacent!(𝒜𝒜n, Tᵣₖⱼ)
    DT.update_adjacent!(𝒜𝒜n, Tᵣᵢₖ)
    @test 𝒜𝒜n(2, 4) == 3
    @test_throws KeyError 𝒜𝒜n(4, 2)
    @test 𝒜𝒜n(3, 2) == 4
    @test 𝒜𝒜n(4, 1) == 3
    @test 𝒜𝒜n(1, 3) == 4
    @test 𝒜𝒜n(4, 3) == 2
    @test 𝒜𝒜n(3, 4) == 1
    @test_throws KeyError 𝒜𝒜n(2, 3)
    @test_throws KeyError 𝒜𝒜n(1, 2)
    @test_throws KeyError 𝒜𝒜n(2, 1)
    DT.flip_edge!(𝒯n, 𝒟n, 𝒜n, 𝒜⁻¹n, 𝒟𝒢n, i, j, k, r)
    @test DT.is_legal(r, k, 𝒜n, IJKR)
    Tᵣₖⱼ = DT.TriangleType((r, k, j))
    Tᵣᵢₖ = DT.TriangleType((r, i, k))
    @test length(𝒜n.adjacent) == 6
    @test (i, j) ∉ keys(𝒜n.adjacent)
    @test (j, i) ∉ keys(𝒜n.adjacent)
    @test 𝒯n.triangles == Set(DT.TriangleType[Tᵣₖⱼ, Tᵣᵢₖ])
    @test 𝒜n(j, r) == k
    @test 𝒜n(r, k) == j
    @test 𝒜n(k, j) == r
    @test 𝒜n(r, i) == k
    @test 𝒜n(i, k) == r
    @test 𝒜n(k, r) == i
    @test 𝒜⁻¹n[k] == Set([(j, r), (r, i)])
    @test 𝒜⁻¹n[j] == Set([(r, k)])
    @test 𝒜⁻¹n[i] == Set([(k, r)])
    @test 𝒜⁻¹n[r] == Set([(k, j), (i, k)])
    @test 𝒟n.graph.N[(i, k, j)] == Set(DT.TriangleType[(r, k, j), (r, i, k)])
    @test 𝒟n.graph.N[(i, j, r)] == Set(DT.TriangleType[(r, k, j), (r, i, k)])
    @test all(==(1), orient.(𝒯n.triangles, Ref(IJKR)))
    @test 𝒟𝒢n(i) == Set([k, r])
    @test 𝒟𝒢n(j) == Set([r, k])
    @test 𝒟𝒢n(r) == Set([i, j, k])
    @test 𝒟𝒢n(k) == Set([j, i, r])
    DT.flip_edge!(𝒯n, 𝒟n, 𝒜n, 𝒜⁻¹n, 𝒟𝒢n, r, k, i, j) # This should go back to the original configuration
    @test 𝒯n.triangles == Set(DT.TriangleType[(2, 1, 3), (2, 4, 1)])
    @test 𝒜𝒜.adjacent == 𝒜n.adjacent
    @test all(𝒟𝒢𝒟𝒢(i) == 𝒟𝒢n(i) for i in 1:4)
    @test 𝒜⁻¹n.adjacent2vertex == 𝒜⁻¹𝒜⁻¹n.adjacent2vertex
    @test !DT.is_legal(i, j, 𝒜n, IJKR)
end