@testset "Can we correctly initialise the triangulation structure?" begin
    𝒯 = DT.Triangles(Set(DT.TriangleType[(DT.LargeRightIdx + 1, DT.LargeLeftIdx, DT.LargeRightIdx)]))
    𝒟 = DT.HistoryDAG()
    𝒜 = DT.Adjacent()
    𝒜⁻¹ = DT.Adjacent2Vertex()
    𝒟𝒢 = DT.DelaunayGraph()
    DT.add_triangle!(𝒟, (DT.LargeRightIdx + 1, DT.LargeLeftIdx, DT.LargeRightIdx))
    root = DT.TriangleType((DT.LargeRightIdx + 1, DT.LargeLeftIdx, DT.LargeRightIdx))
    @test collect(keys(DT.graph(𝒟).V.dict))[1] == DT.TriangleType((DT.LargeRightIdx + 1, DT.LargeLeftIdx, DT.LargeRightIdx))
    𝒜[(DT.LargeRightIdx + 1, DT.LargeLeftIdx)] = DT.LargeRightIdx
    𝒜[(DT.LargeLeftIdx, DT.LargeRightIdx)] = DT.LargeRightIdx + 1
    𝒜[(DT.LargeRightIdx, DT.LargeRightIdx + 1)] = DT.LargeLeftIdx
    ℬ = DT.Adjacent()
    DT.update_adjacent!(ℬ, DT.TriangleType((DT.LargeRightIdx + 1, DT.LargeLeftIdx, DT.LargeRightIdx)))
    @test 𝒜.adjacent == ℬ.adjacent
    𝒜[(DT.LargeLeftIdx, DT.LargeRightIdx + 1)] = DT.BoundaryIdx
    𝒜[(DT.LargeRightIdx, DT.LargeLeftIdx)] = DT.BoundaryIdx
    𝒜[(DT.LargeRightIdx + 1, DT.LargeRightIdx)] = DT.BoundaryIdx
    𝒜⁻¹[0] = Set([(1, -1)])
    𝒜⁻¹[1] = Set([(-1, 0)])
    𝒜⁻¹[-1] = Set([(0, 1)])
    DT.add_neighbour!(𝒟𝒢, 1, -1, 0)
    DT.add_neighbour!(𝒟𝒢, -1, 0, 1)
    DT.add_neighbour!(𝒟𝒢, 0, 1, -1)
    T, D, A, A⁻¹, DG, _root = DT.initialise_triangulation()
    @test T.triangles == 𝒯.triangles
    @test DT.graph(𝒟) == DT.graph(D)
    @test 𝒜.adjacent == A.adjacent
    @test 𝒜⁻¹.adjacent2vertex == A⁻¹.adjacent2vertex
    @test DT.graph(𝒟𝒢) == DT.graph(DG)
    @test root == _root
end