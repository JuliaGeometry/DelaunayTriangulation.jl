@testset "Can we correctly add and remove adjacencies?" begin
    𝒜n = DT.Adjacent()
    i, j, k, r = 1, 2, 3, 4
    𝒜n[(i, j)] = r
    𝒜n[(j, r)] = i
    𝒜n[(r, i)] = j
    𝒜n[(i, k)] = j
    𝒜n[(k, j)] = i
    𝒜n[(j, i)] = k
    𝒜⁻¹n = DT.Adjacent2Vertex()
    𝒜⁻¹n[1] = Set([(2, 4), (3, 2)])
    𝒜⁻¹n[2] = Set([(4, 1), (1, 3)])
    𝒜⁻¹n[3] = Set([(2, 1)])
    𝒜⁻¹n[4] = Set([(1, 2)])
    DT.delete_edge!(𝒜n, i, j)
    @test length(𝒜n.adjacent) == 4
    @test (i, j) ∉ keys(𝒜n.adjacent)
    @test (j, i) ∉ keys(𝒜n.adjacent)
    @test 𝒜n(j, r) == i
    @test 𝒜n(r, i) == j
    @test 𝒜n(i, k) == j
    @test 𝒜n(k, j) == i
    DT.delete_edge!(𝒜⁻¹n, 1, 2, 4)
    DT.delete_edge!(𝒜⁻¹n, 2, 1, 3)
    DT.delete_edge!(𝒜⁻¹n, 3, 2, 1)
    @test isempty(𝒜⁻¹n[k])
    @test 𝒜⁻¹n[2] == Set([(4, 1)])
    @test 𝒜⁻¹n[1] == Set([(3, 2)])
    DT.add_edge!(𝒜⁻¹n, r, 5, 7)
    @test 𝒜⁻¹n[4] == Set([(1, 2), (5, 7)])
end