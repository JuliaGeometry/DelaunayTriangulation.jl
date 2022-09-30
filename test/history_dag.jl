@testset "Can we correctly find the root of a graph?" begin
    𝒟 = DT.HistoryDAG()
    add!(DT.graph(𝒟), DT.TriangleType((1, 2, 3)))
    add!(DT.graph(𝒟), DT.TriangleType((4, 5, 6)))
    add!(DT.graph(𝒟), DT.TriangleType((7, 8, 9)))
    add!(DT.graph(𝒟), DT.TriangleType((10, 11, 12)))
    add!(DT.graph(𝒟), DT.TriangleType((13, 14, 15)))
    add!(DT.graph(𝒟), DT.TriangleType((1, 2, 3)), DT.TriangleType((4, 5, 6)))
    add!(DT.graph(𝒟), DT.TriangleType((1, 2, 3)), DT.TriangleType((7, 8, 9)))
    add!(DT.graph(𝒟), DT.TriangleType((7, 8, 9)), DT.TriangleType((10, 11, 12)))
    add!(DT.graph(𝒟), DT.TriangleType((7, 8, 9)), DT.TriangleType((4, 5, 6)))
    add!(DT.graph(𝒟), DT.TriangleType((4, 5, 6)), DT.TriangleType((13, 14, 15)))
    @test DT.find_root(𝒟; method=:brute) == DT.TriangleType((1, 2, 3))
    @test all(DT.find_root(𝒟; method=:rng) == DT.TriangleType((1, 2, 3)) for _ in 1:10)
end

@testset "Are we correctly handling identical triangles?" begin
    𝒟 = DT.HistoryDAG()
    DT.add_triangle!(𝒟, DT.TriangleType((1, 2, 3)))
    DT.add_triangle!(𝒟, DT.TriangleType((4, 5, 6)))
    DT.add_edge!(𝒟, DT.TriangleType((1, 2, 3)), DT.TriangleType((4, 5, 6)))
    𝒟true = deepcopy(𝒟)
    DT.add_edge!(𝒟, DT.TriangleType((2, 3, 1)), DT.TriangleType((4, 5, 6)))
    @test 𝒟true.graph == 𝒟.graph
    DT.add_edge!(𝒟, DT.TriangleType((2, 3, 1)), DT.TriangleType((5, 6, 4)))
    @test 𝒟true.graph == 𝒟.graph
    DT.add_edge!(𝒟, DT.TriangleType((1, 2, 3)), DT.TriangleType((6, 4, 5)))
    @test 𝒟true.graph == 𝒟.graph
    @test !DT.add_triangle!(𝒟, DT.TriangleType((3, 1, 2)))
    @test 𝒟true.graph == 𝒟.graph
end

@testset "Are we correctly adding multiple triangles and edges to the DAG?" begin
    𝒟 = DT.HistoryDAG()
    𝒟𝒟 = DT.HistoryDAG()
    T₁ = DT.TriangleType((1, 2, 3))
    T₂ = DT.TriangleType((4, 5, 6))
    T₃ = DT.TriangleType((7, 8, 9))
    DT.add_triangle!(𝒟, T₁, T₂, T₃)
    DT.add_triangle!(𝒟𝒟, T₁)
    DT.add_triangle!(𝒟𝒟, T₂)
    DT.add_triangle!(𝒟𝒟, T₃)
    @test DT.graph(𝒟) == DT.graph(𝒟𝒟)
    DT.add_triangle!(𝒟, DT.TriangleType((2, 3, 1))) # same triangle as T₁
    @test DT.graph(𝒟) == DT.graph(𝒟𝒟)
    DT.add_triangle!(𝒟, DT.TriangleType((2, 3, 1)), DT.TriangleType((9, 7, 8))) # same as T₁ and T₃
    @test DT.graph(𝒟) == DT.graph(𝒟𝒟)
    DT.add_edge!(𝒟, T₁, T₂, T₃)
    DT.add_edge!(𝒟𝒟, T₁, T₂)
    DT.add_edge!(𝒟𝒟, T₁, T₃)
    @test DT.graph(𝒟) == DT.graph(𝒟𝒟)
end