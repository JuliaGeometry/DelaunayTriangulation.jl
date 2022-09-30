@testset "Can we correctly construct and index the graph?" begin
    tvn = UndirectedGraph{Int64}()
    add!(tvn, 1)
    add!(tvn, 2)
    add!(tvn, 3)
    [add!(tvn, 1, i) for i in [4, 5, 6, 9, 3]]
    [add!(tvn, 2, i) for i in [11, 9, 8, 15, 16]]
    add!(tvn, 3, 1)
    𝒟𝒢 = DT.DelaunayGraph(tvn)
    @test DT.graph(𝒟𝒢) == tvn
    @test 𝒟𝒢(1) == Set([4, 5, 6, 9, 3])
    @test 𝒟𝒢(2) == Set([11, 9, 8, 15, 16])
    @test 𝒟𝒢(3) == Set([1])
    @test 𝒟𝒢[1] == Set([4, 5, 6, 9, 3])
    @test 𝒟𝒢[2] == Set([11, 9, 8, 15, 16])
    @test 𝒟𝒢[3] == Set([1])
    DT.add_point!(𝒟𝒢, 40)
    DT.add_edge!(𝒟𝒢, 40, 50)
    @test 𝒟𝒢(40) == Set([50])
    DT.add_edge!(𝒟𝒢, 60, 70)
    @test 𝒟𝒢[60] == Set([70])
    DT.add_point!(𝒟𝒢, 70, -1, -2, -3, -4, -5, -6, -7)
    @test all(has(DT.graph(𝒟𝒢), u) for u ∈ [70, -1, -2, -3, -4, -5, -6, -7])
    DT.add_edge!(𝒟𝒢, 70, 4, 5, 9, 3)
    @test 𝒟𝒢[70] == Set([4, 60, 5, 9, 3])
    @test DT.points(𝒟𝒢) == 𝒟𝒢.graph.V
    DT.add_neighbour!(𝒟𝒢, 1, 10)
    @test 𝒟𝒢(1) == Set([4, 5, 6, 9, 3, 10])
    DT.add_neighbour!(𝒟𝒢, 1, 10, 11, 12, 13)
    @test 𝒟𝒢(1) == Set([4, 5, 6, 9, 3, 10, 11, 12, 13])
    𝒟𝒢_empty = DT.DelaunayGraph()
    @test isempty(𝒟𝒢_empty.graph.V)
    @test isempty(DT.points(𝒟𝒢_empty))
    DT.add_neighbour!.(Ref(𝒟𝒢_empty), 2, [3, 4, 5, 1])
    @test 𝒟𝒢_empty[2] == Set([3, 4, 5, 1])
end

@testset "Can we remove points from a neighbourhood?" begin
    tvn = UndirectedGraph{Int64}()
    add!(tvn, 1)
    add!(tvn, 2)
    add!(tvn, 3)
    [add!(tvn, 1, i) for i in [4, 5, 6, 9, 3]]
    [add!(tvn, 2, i) for i in [11, 9, 8, 15, 16]]
    add!(tvn, 3, 1)
    𝒟𝒢 = DT.DelaunayGraph(tvn)
    DT.delete_neighbour!(𝒟𝒢, 1, 6)
    @test 𝒟𝒢[1] == Set([3, 4, 5, 9])
    @test 𝒟𝒢[2] == Set([11, 9, 8, 15, 16])
    @test 𝒟𝒢(3) == Set([1])
end

