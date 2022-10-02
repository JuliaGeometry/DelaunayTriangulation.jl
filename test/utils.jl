############################################
##
## UTILITY FUNCTIONS 
##
############################################
@testset "Can we correctly find the root of a DAG?" begin
    𝒟 = DT.HistoryGraph()
    add!(DT.graph(𝒟), Triangle((1, 2, 3)))
    add!(DT.graph(𝒟), Triangle((4, 5, 6)))
    add!(DT.graph(𝒟), Triangle((7, 8, 9)))
    add!(DT.graph(𝒟), Triangle((10, 11, 12)))
    add!(DT.graph(𝒟), Triangle((13, 14, 15)))
    add!(DT.graph(𝒟), Triangle((1, 2, 3)), Triangle((4, 5, 6)))
    add!(DT.graph(𝒟), Triangle((1, 2, 3)), Triangle((7, 8, 9)))
    add!(DT.graph(𝒟), Triangle((7, 8, 9)), Triangle((10, 11, 12)))
    add!(DT.graph(𝒟), Triangle((7, 8, 9)), Triangle((4, 5, 6)))
    add!(DT.graph(𝒟), Triangle((4, 5, 6)), Triangle((13, 14, 15)))
    @test DT.find_root(𝒟; method=:brute) == Triangle((1, 2, 3))
    @test all(DT.find_root(𝒟; method=:rng) == Triangle((1, 2, 3)) for _ in 1:10)
end

@testset "Can we correctly compute whether a point is higher than another?" begin
    p = [2.0, 3.71] |> Point
    q = [2.0, 4.81] |> Point
    @test !DT.is_point_higher(p, q)
    q = [2.0, 3.60] |> Point
    @test DT.is_point_higher(p, q)
    q = [2.381, 3.71] |> Point
    @test DT.is_point_higher(p, q)
    p = [1.2999, 1.0] |> Point
    q = [1.981, 1.71] |> Point
    @test !DT.is_point_higher(p, q)
    @test DT.is_point_higher(q, p)
    p = [57.131, 4.0] |> Point
    q = [2.0, 3.1] |> Point
    @test DT.is_point_higher(p, q)
    @test DT.is_point_lower(q, p)
    p = [-2.31, 4.0] |> Point
    q = [5.0, 4.0] |> Point
    @test DT.is_point_higher(p, q)
    @test DT.is_point_lower(q, p)
end

@testset "Can we correctly sort points by height?" begin
    v = Vector{Point{Float64,Vector{Float64}}}(undef, 100)
    for _ in 1:500
        v .= [Point(rand(2)) for _ in 1:100]
        DT.partial_highest_point_sort!(v, 1)
        @test all(DT.is_point_higher(v[1], v[j]) for j in 2:lastindex(v))
    end
end

@testset "Can we count the number of negative values?" begin
    v = [-1, -2, 3, 5, 4.0]
    @test DT.num_less(0, v) == 2
    @test DT.num_less(7, v) == 5
    v = [2, 5, 10, 29.0]
    @test DT.num_less(0, v) == 0
    v = [2, 1, 5, 4]
    @test DT.num_less(10, v) == 4
    @test DT.num_less(5, v) == 3
end

@testset "Can we correctly identify a given triangulation as Delaunay?" begin
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    pts = Points(p1, p2, p3, p4)
    DTri = DT.initialise_triangulation(pts)
    T, HG, adj, adj2v, DG, root = triangles(DTri), history(DTri),
    adjacent(DTri), adjacent2vertex(DTri), graph(DTri), DT.root(DTri)
    r = 1
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    r = 2
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    r = 3
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test !DT.is_delaunay(DTri)
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    @test DT.is_delaunay(DTri)
end

@testset "Testing that we get the correct number of edges, points, and triangles" begin
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    pts = Points(p1, p2, p3, p4)
    DTri = triangulate(pts)
    @test num_triangles(DTri) == 2
    @test num_points(DTri) == 4
    @test num_edges(DTri) == 5
end