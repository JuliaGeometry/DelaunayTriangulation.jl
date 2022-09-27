## Can we correctly find the root of a graph?
𝒟 = HistoryDAG()
add!(𝒟, TriangleType((1, 2, 3)))
add!(𝒟, TriangleType((4, 5, 6)))
add!(𝒟, TriangleType((7, 8, 9)))
add!(𝒟, TriangleType((10, 11, 12)))
add!(𝒟, TriangleType((13, 14, 15)))
add!(𝒟, TriangleType((1, 2, 3)), TriangleType((4, 5, 6)))
add!(𝒟, TriangleType((1, 2, 3)), TriangleType((7, 8, 9)))
add!(𝒟, TriangleType((7, 8, 9)), TriangleType((10, 11, 12)))
add!(𝒟, TriangleType((7, 8, 9)), TriangleType((4, 5, 6)))
add!(𝒟, TriangleType((4, 5, 6)), TriangleType((13, 14, 15)))
@test find_root(𝒟; method=:brute) == TriangleType((1, 2, 3))
@test all(find_root(𝒟; method=:rng) == TriangleType((1, 2, 3)) for _ in 1:10)

## Are we correctly handling identical triangles?
𝒟 = HistoryDAG()
add!(𝒟, TriangleType((1, 2, 3)))
add!(𝒟, TriangleType((4, 5, 6)))
add!(𝒟, TriangleType((1, 2, 3)), TriangleType((4, 5, 6)))
𝒟true = deepcopy(𝒟)
add!(𝒟, TriangleType((2, 3, 1)), TriangleType((4, 5, 6)))
@test 𝒟true.graph == 𝒟.graph
add!(𝒟, TriangleType((2, 3, 1)), TriangleType((5, 6, 4)))
@test 𝒟true.graph == 𝒟.graph
add!(𝒟, TriangleType((1, 2, 3)), TriangleType((6, 4, 5)))
@test 𝒟true.graph == 𝒟.graph
@test !add!(𝒟, TriangleType((3, 1, 2)))
@test 𝒟true.graph == 𝒟.graph

## Are we correctly adding multiple triangles and edges to the DAG?
𝒟 = HistoryDAG()
𝒟𝒟 = HistoryDAG()
T₁ = TriangleType((1, 2, 3))
T₂ = TriangleType((4, 5, 6))
T₃ = TriangleType((7, 8, 9))
add_triangle!(𝒟, T₁, T₂, T₃)
add!(𝒟𝒟, T₁)
add!(𝒟𝒟, T₂)
add!(𝒟𝒟, T₃)
@test graph(𝒟) == graph(𝒟𝒟)
add_triangle!(𝒟, TriangleType((2, 3, 1))) # same triangle as T₁
@test graph(𝒟) == graph(𝒟𝒟)
add_triangle!(𝒟, TriangleType((2, 3, 1)), TriangleType((9, 7, 8))) # same as T₁ and T₃
@test graph(𝒟) == graph(𝒟𝒟)
add_edge!(𝒟, T₁, T₂, T₃)
add!(𝒟𝒟, T₁, T₂)
add!(𝒟𝒟, T₁, T₃)
@test graph(𝒟) == graph(𝒟𝒟)