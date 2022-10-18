@testset "Testing point insertion with the Bowyer-Watson method" begin
    ## Existing triangulation with no bounding etc 
    p1 = (5.0, 6.0)
    p2 = (9.0, 6.0)
    p3 = (13.0, 5.0)
    p4 = (10.38, 0.0)
    p5 = (12.64, -1.69)
    p6 = (2.0, -2.0)
    p7 = (3.0, 4.0)
    p8 = (7.5, 3.53)
    p9 = (4.02, 1.85)
    p10 = (4.26, 0.0)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
    Random.seed!(928881)
    T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
    p11 = (6.0, 2.5)
    push!(pts, p11)

    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 11)
    _T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
    DT.clear_empty_keys!(adj, DG)
    DT.clear_empty_keys!(_adj, _DG)
    @test DT.compare_triangle_sets(T, _T) &&
          adjacent(adj) == adjacent(_adj) &&
          adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
          graph(DG) == graph(_DG)

    p12 = (10.3, 2.85)
    push!(pts, p12)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 12)
    _T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
    DT.clear_empty_keys!(adj, DG)
    DT.clear_empty_keys!(_adj, _DG)
    @test DT.compare_triangle_sets(T, _T) &&
          adjacent(adj) == adjacent(_adj) &&
          adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
          graph(DG) == graph(_DG)

    p13 = (7.5, 3.5)
    push!(pts, p13)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 13)
    _T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
    DT.clear_empty_keys!(adj, DG)
    DT.clear_empty_keys!(_adj, _DG)
    @test DT.compare_triangle_sets(T, _T) &&
          adjacent(adj) == adjacent(_adj) &&
          adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
          graph(DG) == graph(_DG)

    ## Now do this for some circular points so that we can do some random tests - only testing for interior insertion here
    R = 10.0
    n = 1381
    pts = [(rand((-1, 1)) * R * rand(), rand((-1, 1)) * R * rand()) for _ in 1:n] # rand((-1, 1)) to get random signs
    pushfirst!(pts, (-11.0, -11.0), (11.0, -11.0), (11.0, 11.0), (-11.0, 11.0)) # bounding box to guarantee interior insertions only
    T, adj, adj2v, DG, HG = @views DT.triangulate_berg(pts[1:7])
    _T, _adj, _adj2v, _DG, _HG = @views DT.triangulate_berg(pts[1:7])
    n = length(pts)
    for i in 8:n
        DT.add_point_bowyer!(T, adj, adj2v, DG, pts, i)
        DT.add_point_berg!(_T, _adj, _adj2v, _DG, _HG, pts, i)
        @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
    end
end