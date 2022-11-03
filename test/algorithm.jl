############################################
##
## UNCONSTRAINED INCREMENTAL TRIANGULATION
##
############################################
@testset "Initialisation" begin
    p1 = [2.7, 13.0]
    p2 = [-2.0, 3.0]
    p3 = [7.0, 2.0]
    p4 = [18.0, 3.5]
    pts = [p1, p2, p3, p4]
    Tri = UnconstrainedTriangulation(pts; method=:berg)
    @test adjacent(Tri) == Tri.adjacent
    @test adjacent2vertex(Tri) == Tri.adjacent2vertex
    @test graph(Tri) == Tri.graph
    @test pointlocation(Tri) == Tri.pointlocation
    @test triangles(Tri) == Tri.triangles
    @test points(Tri) == Tri.points
    @test triangles(Tri) == Set{NTuple{3,Int64}}([DT.BoundingTriangle])
    @test adjacent(Tri, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test adjacent(Tri, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test adjacent(Tri, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)) == DT.UpperBoundingIndex
    @test adjacent(Tri, (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri, (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.LowerRightBoundingIndex
    @test adjacent(Tri, (DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri, (DT.UpperBoundingIndex, DT.LowerRightBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri, (DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.BoundaryIndex
    @test adjacent2vertex(Tri, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)])
    @test adjacent2vertex(Tri, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)])
    @test adjacent2vertex(Tri, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)])
    @test adjacent2vertex(Tri, DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (DT.UpperBoundingIndex, DT.LowerRightBoundingIndex),
        (DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test neighbours(Tri, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex])
    @test pointlocation(Tri).graph.N[DT.BoundingTriangle] == pointlocation(Tri).graph.N[NTuple{3,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))] == Set()
    @test points(Tri) == pts
    @test adjacent(Tri) isa DT.Adjacent{Int64,NTuple{2,Int64}}
    @test adjacent2vertex(Tri) isa DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}
    @test graph(Tri) isa DT.DelaunayGraph{Int64}
    @test pointlocation(Tri) isa DT.HistoryGraph{NTuple{3,Int64}}
    @test triangles(Tri) isa Set{NTuple{3,Int64}}
    @test points(Tri) isa Vector{Vector{Float64}}
    @test Tri isa UnconstrainedTriangulation{NTuple{3,Int64},NTuple{2,Int64},Set{NTuple{3,Int64}},Set{NTuple{2,Int64}},Vector{Vector{Float64}},Int64,DT.HistoryGraph{NTuple{3,Int64}}}
    @test typeof(Tri) <: DT.AbstractTriangulation
    @test typeof(Tri) <: DT.AbstractUnconstrainedTriangulation

    Tri3 = UnconstrainedTriangulation(pts; IntegerType=Int16, method=:berg)
    @test adjacent(Tri3) == Tri3.adjacent
    @test adjacent2vertex(Tri3) == Tri3.adjacent2vertex
    @test graph(Tri3) == Tri3.graph
    @test pointlocation(Tri3) == Tri3.pointlocation
    @test triangles(Tri3) == Tri3.triangles
    @test points(Tri3) == Tri3.points
    @test triangles(Tri3) == Set{NTuple{3,Int16}}([NTuple{3,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))])
    @test adjacent(Tri3, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test adjacent(Tri3, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri3, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test adjacent(Tri3, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex))) == DT.UpperBoundingIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))) == DT.LowerRightBoundingIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex))) == DT.BoundaryIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.UpperBoundingIndex, DT.LowerRightBoundingIndex))) == DT.BoundaryIndex
    @test adjacent(Tri3, NTuple{2,Int16}((DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex))) == DT.BoundaryIndex
    @test adjacent2vertex(Tri3, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int16}}([NTuple{2,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))])
    @test adjacent2vertex(Tri3, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int16}}([NTuple{2,Int16}((DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))])
    @test adjacent2vertex(Tri3, DT.UpperBoundingIndex) == Set{NTuple{2,Int16}}([NTuple{2,Int16}((DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex))])
    @test neighbours(Tri3, DT.LowerLeftBoundingIndex) == Set{Int16}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri3, DT.LowerRightBoundingIndex) == Set{Int16}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri3, DT.UpperBoundingIndex) == Set{Int16}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex])
    @test pointlocation(Tri3).graph.N[NTuple{3,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))] == pointlocation(Tri3).graph.N[NTuple{3,Int16}((DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex))] == Set()
    @test points(Tri3) == pts
    @test adjacent(Tri3) isa DT.Adjacent{Int16,NTuple{2,Int16}}
    @test adjacent2vertex(Tri3) isa DT.Adjacent2Vertex{Int16,Set{NTuple{2,Int16}},NTuple{2,Int16}}
    @test graph(Tri3) isa DT.DelaunayGraph{Int16}
    @test pointlocation(Tri3) isa DT.HistoryGraph{NTuple{3,Int16}}
    @test triangles(Tri3) isa Set{NTuple{3,Int16}}
    @test points(Tri3) isa Vector{Vector{Float64}}
    @test Tri3 isa UnconstrainedTriangulation{NTuple{3,Int16},NTuple{2,Int16},Set{NTuple{3,Int16}},Set{NTuple{2,Int16}},Vector{Vector{Float64}},Int16,DT.HistoryGraph{NTuple{3,Int16}}}
    @test typeof(Tri) <: DT.AbstractTriangulation
    @test typeof(Tri) <: DT.AbstractUnconstrainedTriangulation
end

@testset "Testing that the point location data structure can correctly locate triangles" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)

    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))

    p = (33.373, 15.2287)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (-31.0689, 52.90257)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = (3.63, 1.679)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    @test flag == 1

    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))

    p = (3.63, 1.679)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = (27.706, 0.968)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = (-13.6689, 1.3567)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1)
    @test flag == 1
    p = (-3.56804, 1.745279)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1)
    @test flag == 1
    p = (5.95, -2.91669)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = (32.9507, -4.2764)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = (41.4976, 46.81)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (10.0, 10.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (-33.48, 23.11)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = (-10.0, 10.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1

    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))

    p = (-10.0, 10.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = (-36.59, 13.594)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = (35.86, 34.379)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (15.66, 7.766)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = (5.173, -3.305)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = (-6.09, -3.305)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = (14.57, 0.48686)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = (9.89, 2.198)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = (-5.735, 3.11)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = (-3.7957, 3.11)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = (-11.21, 2.54)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = (-3.68, 1.057)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3)
    @test flag == 1
    p = (0.0, 0.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test flag == 0
    p = (0.916, 1.408)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)
    @test flag == 1
    p = (0.0, 2.0)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)
    @test flag == 1
    p = (2.5057, 2.8986)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts), DT.BoundingTriangle)
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)
    @test flag == 1
end

@testset "Testing that we can correctly add a point into a triangulation" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    # Building an example triangulation
    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    T = Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    ])
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 2, 3)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 3, 1)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 2, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 1, 2)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, 1, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, 2, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.LowerLeftBoundingIndex, 3)
    DT.add_edge!(adj2v, 1, 3, 2)
    DT.add_edge!(adj2v, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 2, 1, 3)
    DT.add_edge!(adj2v, 2, 3, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerRightBoundingIndex, 1)
    DT.add_edge!(adj2v, 3, DT.LowerLeftBoundingIndex, 2)
    DT.add_edge!(adj2v, 3, 1, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 3, 2, 1)
    DG = DT.DelaunayGraph{Int64}()
    DT.add_point!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1, 2, 3)
    DT.add_neighbour!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2, 3, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    DT.add_neighbour!(DG, 1, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 2, 3)
    DT.add_neighbour!(DG, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3)
    DT.add_neighbour!(DG, 3, DT.LowerLeftBoundingIndex, 1, 2)

    # Now do the actual test 
    DT.split_triangle!(T, HG, adj, adj2v, DG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1), 4)
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4) ∈ T
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test DT.get_edge(adj, 3, 1) == 4
    @test DT.get_edge(adj, 1, 4) == 3
    @test DT.get_edge(adj, 4, 3) == 1
    @test DT.get_edge(adj, 4, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == 4
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 4) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 2, 1) == 3
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 3, 2) == 1
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 3),
        (3, 4),
        (4, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (1, 2),
        (2, DT.LowerLeftBoundingIndex),
        (DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (2, DT.LowerRightBoundingIndex),
        (3, 2),
        (4, 3),
        (DT.LowerLeftBoundingIndex, 4),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (3, DT.LowerLeftBoundingIndex),
        (1, 3),
        (DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 2),
        (4, DT.LowerLeftBoundingIndex),
        (1, 4),
        (2, 1)
    ])
    @test DT.get_edge(adj2v, 4) == Set{NTuple{2,Int64}}([
        (1, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, 3),
        (3, 1)
    ])
    @test_throws KeyError DT.get_edge(adj2v, 5)
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, 4, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 4, 3, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 3, 1])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 1, 2, DT.LowerLeftBoundingIndex])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([1, DT.LowerLeftBoundingIndex, 1, 3])
    @test !DT.islegal(3, 1, adj, pts)
end

@testset "Testing that we can correctly flip an edge" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    # Building an example triangulation
    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    T = Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    ])
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 2, 3)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 3, 1)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 2, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 1, 2)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, 1, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, 2, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.LowerLeftBoundingIndex, 3)
    DT.add_edge!(adj2v, 1, 3, 2)
    DT.add_edge!(adj2v, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 2, 1, 3)
    DT.add_edge!(adj2v, 2, 3, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerRightBoundingIndex, 1)
    DT.add_edge!(adj2v, 3, DT.LowerLeftBoundingIndex, 2)
    DT.add_edge!(adj2v, 3, 1, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 3, 2, 1)
    DG = DT.DelaunayGraph{Int64}()
    DT.add_point!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1, 2, 3)
    DT.add_neighbour!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2, 3, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    DT.add_neighbour!(DG, 1, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 2, 3)
    DT.add_neighbour!(DG, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3)
    DT.add_neighbour!(DG, 3, DT.LowerLeftBoundingIndex, 1, 2)
    DT.split_triangle!(T, HG, adj, adj2v, DG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    DT.flip_edge!(T, HG, adj, adj2v, DG, i, j, k, r)
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4) ∈ T
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    ])
    p = (1.88, 2.5834)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)
    @test flag == 1
    p = (-0.9802, 2.5834)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    @test flag == 1
    p = (-3.642, 2.8615)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)
    @test flag == 1
    p = (4.3036, -3.0977)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)] || in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4), DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)] || in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4), DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test DT.get_edge(adj, 3, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 4) == 2
    @test DT.get_edge(adj, 4, 3) == 2
    @test DT.get_edge(adj, 4, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == 4
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 4) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 2, 1) == 4
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 3, 2) == 4
    @test DT.get_edge(adj, 4, 2) == 1
    @test DT.get_edge(adj, 2, 4) == 3
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 3),
        (3, 4),
        (4, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (1, 2),
        (2, DT.LowerLeftBoundingIndex),
        (DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (2, DT.LowerRightBoundingIndex),
        (4, 2),
        (DT.LowerLeftBoundingIndex, 4),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (3, DT.LowerLeftBoundingIndex),
        (1, 4),
        (DT.LowerRightBoundingIndex, 1),
        (4, 3)
    ])
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 2),
        (4, DT.LowerLeftBoundingIndex),
        (2, 4)
    ])
    @test DT.get_edge(adj2v, 4) == Set{NTuple{2,Int64}}([
        (1, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, 3),
        (3, 2),
        (2, 1)
    ])
    @test_throws KeyError DT.get_edge(adj2v, 5)
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, 4, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 4, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 3, 1, 4])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 2, DT.LowerLeftBoundingIndex])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([1, DT.LowerLeftBoundingIndex, 1, 3, 2])
    @test DT.islegal(r, k, adj, pts)
end

@testset "Testing legalise edge" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    # Building an example triangulation
    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    T = Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    ])
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 2, 3)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 3, 1)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 2, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 1, 2)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, 1, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, 2, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.LowerLeftBoundingIndex, 3)
    DT.add_edge!(adj2v, 1, 3, 2)
    DT.add_edge!(adj2v, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 2, 1, 3)
    DT.add_edge!(adj2v, 2, 3, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerRightBoundingIndex, 1)
    DT.add_edge!(adj2v, 3, DT.LowerLeftBoundingIndex, 2)
    DT.add_edge!(adj2v, 3, 1, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 3, 2, 1)
    DG = DT.DelaunayGraph()
    DT.add_point!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1, 2, 3)
    DT.add_neighbour!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2, 3, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    DT.add_neighbour!(DG, 1, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 2, 3)
    DT.add_neighbour!(DG, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3)
    DT.add_neighbour!(DG, 3, DT.LowerLeftBoundingIndex, 1, 2)
    DT.split_triangle!(T, HG, adj, adj2v, DG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    @test !DT.islegal(i, j, adj, pts)
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4) ∉ T
    @test DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2) ∈ T
    @test DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4) ∈ T
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    ])
    @test out_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    ])
    p = (1.88, 2.5834)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)
    @test flag == 1
    p = (-0.9802, 2.5834)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)
    @test flag == 1
    p = (-3.642, 2.8615)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)
    @test flag == 1
    p = (4.3036, -3.0977)
    push!(pts, p)
    intri, flag = DT.locate_triangle(HG, pts, length(pts))
    @test intri == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 4)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)] || in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4), DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)]
    @test in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1), DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4)] || in_neighbors(HG, DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)) == [DT.construct_triangle(NTuple{3,Int64}, 3, 1, 4), DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test DT.get_edge(adj, 3, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 4) == 2
    @test DT.get_edge(adj, 4, 3) == 2
    @test DT.get_edge(adj, 4, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == 4
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 4) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 2, 1) == 4
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 3, 2) == 4
    @test DT.get_edge(adj, 4, 2) == 1
    @test DT.get_edge(adj, 2, 4) == 3
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 3),
        (3, 4),
        (4, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (1, 2),
        (2, DT.LowerLeftBoundingIndex),
        (DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (2, DT.LowerRightBoundingIndex),
        (4, 2),
        (DT.LowerLeftBoundingIndex, 4),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (3, DT.LowerLeftBoundingIndex),
        (1, 4),
        (DT.LowerRightBoundingIndex, 1),
        (4, 3)
    ])
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 2),
        (4, DT.LowerLeftBoundingIndex),
        (2, 4)
    ])
    @test DT.get_edge(adj2v, 4) == Set{NTuple{2,Int64}}([
        (1, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, 3),
        (3, 2),
        (2, 1)
    ])
    @test_throws KeyError DT.get_edge(adj2v, 5)
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, 4, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 4, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 3, 1, 4])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 2, DT.LowerLeftBoundingIndex])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([1, DT.LowerLeftBoundingIndex, 1, 3, 2])
    @test DT.islegal(r, k, adj, pts)
end

@testset "Testing that we can remove the bounding triangle" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    p5 = (2.0, 3.0)
    pts = [p1, p2, p3, p4, p5]

    # Building an example triangulation
    HG = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    DT.add_edge!(HG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1))
    T = Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 2, DT.LowerRightBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 2, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1)
    ])
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 2, 3)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 3, 1)
    DT.add_edge!(adj2v, DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 2, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, 1, 2)
    DT.add_edge!(adj2v, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)
    DT.add_edge!(adj2v, DT.UpperBoundingIndex, 1, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, 2, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.LowerLeftBoundingIndex, 3)
    DT.add_edge!(adj2v, 1, 3, 2)
    DT.add_edge!(adj2v, 1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    DT.add_edge!(adj2v, 1, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    DT.add_edge!(adj2v, 2, 1, 3)
    DT.add_edge!(adj2v, 2, 3, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 2, DT.LowerRightBoundingIndex, 1)
    DT.add_edge!(adj2v, 3, DT.LowerLeftBoundingIndex, 2)
    DT.add_edge!(adj2v, 3, 1, DT.LowerLeftBoundingIndex)
    DT.add_edge!(adj2v, 3, 2, 1)
    DG = DT.DelaunayGraph{Int64}()
    DT.add_point!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1, 2, 3)
    DT.add_neighbour!(DG, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2, 3, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex)
    DT.add_neighbour!(DG, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    DT.add_neighbour!(DG, 1, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 2, 3)
    DT.add_neighbour!(DG, 2, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3)
    DT.add_neighbour!(DG, 3, DT.LowerLeftBoundingIndex, 1, 2)
    DT.split_triangle!(T, HG, adj, adj2v, DG, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    DT.remove_bounding_triangle!(T, adj, adj2v, DG)
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, 4, 3, 2),
        DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)
    ])
    @test DT.get_edge(adj, 3, 2) == 4
    @test DT.get_edge(adj, 2, 3) == DT.BoundaryIndex
    @test DT.get_edge(adj, 2, 1) == 4
    @test DT.get_edge(adj, 1, 2) == DT.BoundaryIndex
    @test DT.get_edge(adj, 1, 4) == 2
    @test DT.get_edge(adj, 4, 1) == DT.BoundaryIndex
    @test DT.get_edge(adj, 4, 3) == 2
    @test DT.get_edge(adj, 3, 4) == DT.BoundaryIndex
    @test DT.get_edge(adj, 4, 2) == 1
    @test DT.get_edge(adj, 2, 4) == 3
    @test length(adj.adjacent) == 10
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (2, 4)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (4, 2)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (1, 4),
        (4, 3)
    ])
    @test DT.get_edge(adj2v, 4) == Set{NTuple{2,Int64}}([
        (3, 2),
        (2, 1)
    ])
    @test DT.get_edge(adj2v, DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (1, 2),
        (4, 1),
        (3, 4),
        (2, 3)
    ])
    @test length(adj2v.adjacent2vertex) == 5
    @test DT.get_neighbour(DG, 1) == Set{Int64}([4, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([1, 4, 3])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 2])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([2, 3, 1])
end

@testset "Testing the steps for a small problem" begin
    p1 = (5.0, 5.0)
    p2 = (1.0, -1.0)
    p3 = (-2.0, 2.0)
    p4 = (-1.0, 4.0)
    pts = [p1, p2, p3, p4]

    DTri = DT.UnconstrainedTriangulation(pts; method=:berg)
    T, hg, adj, adj2v, DG, root = triangles(DTri), pointlocation(DTri),
    adjacent(DTri), adjacent2vertex(DTri), graph(DTri), DT.BoundingTriangle

    r = 1
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(hg, pts, r, root)
    @test Tᵢⱼₖ == DT.BoundingTriangle
    @test interior_flag == 1
    i, j, k = Tᵢⱼₖ
    @test i == DT.LowerRightBoundingIndex
    @test j == DT.UpperBoundingIndex
    @test k == DT.LowerLeftBoundingIndex
    DT.split_triangle!(T, hg, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 1) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, 1) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.UpperBoundingIndex, 1),
        (1, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, 1])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test out_neighbors(hg, DT.BoundingTriangle) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == [DT.BoundingTriangle]

    r = 2
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(hg, pts, r, root)
    i, j, k = Tᵢⱼₖ
    @test Tᵢⱼₖ == DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    @test interior_flag == 1
    DT.split_triangle!(T, hg, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, 1, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 2)
    ])
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 2) == 1
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 1) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, 2) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, 1) == 2
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, 1) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 1, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 2, 1) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 2, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 2, DT.LowerRightBoundingIndex) == 1
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 1),
        (1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.UpperBoundingIndex, 1),
        (1, 2),
        (2, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 2),
        (2, DT.LowerRightBoundingIndex),
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 1),
        (1, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test out_neighbors(hg, DT.BoundingTriangle) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test out_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, 1, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 2)
    ])
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, 1, 2)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 2)) == [DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test DT.islegal(i, j, adj, pts)
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)

    r = 3
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(hg, pts, r, root)
    i, j, k = Tᵢⱼₖ
    DT.split_triangle!(T, hg, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test !DT.islegal(i, j, adj, pts)
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)
    DT.legalise_edge!(T, hg, adj, adj2v, DG, i, j, r, pts)
    @test Tᵢⱼₖ == DT.construct_triangle(NTuple{3,Int64}, 1, DT.LowerLeftBoundingIndex, 2)
    @test interior_flag == 1
    @test T == Set{NTuple{3,Int64}}([
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        DT.construct_triangle(NTuple{3,Int64}, 3, 1, DT.UpperBoundingIndex),
        DT.construct_triangle(NTuple{3,Int64}, 2, 1, 3),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerRightBoundingIndex, 1, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        DT.construct_triangle(NTuple{3,Int64}, DT.LowerLeftBoundingIndex, 2, 3),
        DT.construct_triangle(NTuple{3,Int64}, 3, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == 2
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 2) == 3
    @test DT.get_edge(adj, 2, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == 3
    @test DT.get_edge(adj, 3, 1) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, 1, 3) == 2
    @test DT.get_edge(adj, 2, 1) == 3
    @test DT.get_edge(adj, 1, 2) == DT.LowerRightBoundingIndex
    @test DT.get_edge(adj, 1, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, 1) == 2
    @test DT.get_edge(adj, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == 1
    @test DT.get_edge(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, 2),
        (2, 3),
        (3, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.UpperBoundingIndex, 1)
        (1, 2)
        (2, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{NTuple{2,Int64}}([
        (DT.LowerLeftBoundingIndex, 3),
        (3, 1),
        (1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{NTuple{2,Int64}}([
        (DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        (3, 2),
        (2, DT.LowerRightBoundingIndex),
        (DT.UpperBoundingIndex, 3)
    ])
    @test DT.get_edge(adj2v, 2) == Set{NTuple{2,Int64}}([
        (1, 3),
        (3, DT.LowerLeftBoundingIndex),
        (DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        (DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([
        (1, DT.UpperBoundingIndex),
        (DT.LowerLeftBoundingIndex, 2),
        (DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        (2, 1)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.UpperBoundingIndex, 3, DT.LowerRightBoundingIndex, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([1, 2, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)

    DTri = triangulate(pts; method=:berg)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false, method=:berg)
    @test DT.is_delaunay(DTri)
end

@testset "Random triangulations" begin
    p1 = (-6.88, 3.61)
    p2 = (-6.08, -0.43)
    p3 = (-0.3, 2.01)
    p4 = (5.1, -1.27)
    p5 = (6.18, 1.87)
    p6 = (3.08, 4.43)
    p7 = (-1.34, 4.83)
    p8 = (-1.68, -0.77)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    DTri = triangulate(pts; method=:berg)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false, method=:berg)
    @test DT.is_delaunay(DTri)
    pts = [[-6.88, 3.61], [-6.08, -0.43], [-0.3, 2.01], [5.1, -1.27], [6.18, 1.87], [3.08, 4.43], [-1.34, 4.83], [-1.68, -0.77]]
    DTri = triangulate(pts)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false, method=:berg)
    @test DT.is_delaunay(DTri)

    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    DTri = triangulate(pts; randomise=false, trim=false, method=:berg)
    DTri2 = triangulate([p1, p2, p3, p4, p5, p6, p7]; randomise=false, trim=false, method=:berg)
    add_point!(DTri2, (-1.68, -0.77))
    @test triangles(DTri2) == triangles(DTri)

    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    DTri = triangulate(pts; randomise=false, trim=true, method=:berg)
    DTri2 = triangulate([p1, p2, p3, p4, p5, p6, p7]; randomise=false, trim=false, method=:berg)
    add_point!(DTri2, (-1.68, -0.77))
    DT.remove_bounding_triangle!(DTri2)
    @test triangles(DTri2) == triangles(DTri)

    for _ in 1:250
        x = rand(100)
        y = rand(100)
        pts = [(x, y) for (x, y) in zip(x, y)]
        DTri = triangulate(pts; method=:berg)
        @test DT.is_delaunay(DTri)
        @test all(DT.isoriented(T, pts) == 1 for T in triangles(DTri))
        DTri = triangulate(pts; method=:berg, trim=false)
        @test DT.is_delaunay(DTri)
        @test all(DT.isoriented(T, pts) == 1 for T in triangles(DTri))
        @test 1 == num_points(DTri) - num_edges(DTri) + num_triangles(DTri) # Euler's formula
    end
end

import CairoMakie: poly!, Figure, Axis, lines!, save
fig = Figure()
ax = Axis(fig[1, 1])
x = rand(100)
y = rand(100)
pts = [(x, y) for (x, y) in zip(x, y)]
DTri = triangulate(pts; method=:berg)
Tmat = zeros(Int64, num_triangles(DTri), 3)
for (i, T) in enumerate(triangles(DTri))
    Tmat[i, :] = [geti(T), getj(T), getk(T)]
end
pmat = zeros(num_points(DTri), 2)
for (i, p) in enumerate(points(DTri))
    pmat[i, :] = [getx(p), gety(p)]
end
poly!(ax, pmat, Tmat, strokewidth=2)
for (i, j) in adjacent2vertex(DTri, DelaunayTriangulation.BoundaryIndex)
    p = DT.get_point(DTri, i)
    q = DT.get_point(DTri, j)
    lines!(ax, [getx(p), getx(q)], [gety(p), gety(q)], color=:red, linewidth=5)
end
save("figures/test_triangulation.png", fig)

@testset "Custom integer types" begin
    p1 = (-6.88, 3.61)
    p2 = (-6.08, -0.43)
    p3 = (-0.3, 2.01)
    p4 = (5.1, -1.27)
    p5 = (6.18, 1.87)
    p6 = (3.08, 4.43)
    p7 = (-1.34, 4.83)
    p8 = (-1.68, -0.77)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    DTri = triangulate(pts; IntegerType=Int16, randomise=false, method=:berg)
    @test points(DTri) == pts
    @test adjacent(DTri) isa DT.Adjacent{Int16,NTuple{2,Int16}}
    @test adjacent2vertex(DTri) isa DT.Adjacent2Vertex{Int16,Set{NTuple{2,Int16}},NTuple{2,Int16}}
    @test graph(DTri) isa DT.DelaunayGraph{Int16}
    @test pointlocation(DTri) isa DT.HistoryGraph{NTuple{3,Int16}}
    @test triangles(DTri) isa Set{NTuple{3,Int16}}
    @test points(DTri) isa Vector{NTuple{2,Float64}}
    @test DTri isa UnconstrainedTriangulation{NTuple{3,Int16},NTuple{2,Int16},Set{NTuple{3,Int16}},
        Set{NTuple{2,Int16}},Vector{NTuple{2,Float64}},Int16}
    @test DT.is_delaunay(DTri)
    @test DT.get_point(DTri, 1) == p1
    @test DT.get_point(DTri, 2) == p2
    @test DT.get_point(DTri, 3) == p3
    @test DT.get_point(DTri, 4) == p4
    @test DT.get_point(DTri, 5) == p5
    @test DT.get_point(DTri, 6) == p6
    @test DT.get_point(DTri, 7) == p7
    @test DT.get_point(DTri, 8) == p8
end

############################################
##
## UTILITIES
##
############################################ 
@testset "Can we correctly find the root of a DAG?" begin
    D = DT.HistoryGraph{NTuple{3,Int64}}()
    add!(DT.graph(D), (1, 2, 3))
    add!(DT.graph(D), (4, 5, 6))
    add!(DT.graph(D), (7, 8, 9))
    add!(DT.graph(D), (10, 11, 12))
    add!(DT.graph(D), (13, 14, 15))
    add!(DT.graph(D), (1, 2, 3), (4, 5, 6))
    add!(DT.graph(D), (1, 2, 3), (7, 8, 9))
    add!(DT.graph(D), (7, 8, 9), (10, 11, 12))
    add!(DT.graph(D), (7, 8, 9), (4, 5, 6))
    add!(DT.graph(D), (4, 5, 6), (13, 14, 15))
    @test DT.find_root(D; method=:brute) == (1, 2, 3)
    @test all(DT.find_root(D; method=:rng) == (1, 2, 3) for _ in 1:10)
end

############################################
##
## DELETING A TRIANGLE 
##
############################################
@testset "Deleting a triangle" begin
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    p4 = @SVector[-1.0, 2.0]
    p5 = @SVector[4.0, 2.0]
    p6 = @SVector[-2.0, -1.0]
    pts = [p1, p2, p3, p4, p5, p6]
    tri = triangulate(pts; randomise=false)
    p7 = @SVector[2.0, 1.0]
    newtri = deepcopy(tri)
    push!(points(newtri), p7)
    DT.split_triangle!(triangles(newtri),
        pointlocation(newtri), adjacent(newtri),
        adjacent2vertex(newtri), graph(newtri),
        (1, 3, 5), 7)
    newtri2 = deepcopy(newtri)
    adj = adjacent(newtri2)
    DT.delete_triangle!(triangles(newtri2),
        adjacent(newtri2),
        adjacent2vertex(newtri2), graph(newtri2),
        1, 3, 7)
    @test triangles(newtri2) == Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (6, 2, 3),
        (5, 3, 2),
        (5, 1, 7),
        (3, 5, 7),
        (6, 1, 4),
        (5, 4, 1)
    ])
    @test length(triangles(newtri2)) == 7
    @test (1, 3, 7) ∉ triangles(newtri2) &&
          (3, 7, 1) ∉ triangles(newtri2) &&
          (7, 1, 3) ∉ triangles(newtri2)
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == 2
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 3, 2) == 5
    @test DT.get_edge(adjacent(newtri2), 2, 3) == 6
    @test DT.get_edge(adjacent(newtri2), 6, 2) == 3
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 5) == 3
    @test DT.get_edge(adjacent(newtri2), 5, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 3) == 2
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.is_boundary_edge(2, 6, adjacent(newtri2))
    @test DT.is_boundary_edge(6, 4, adjacent(newtri2))
    @test DT.is_boundary_edge(4, 5, adjacent(newtri2))
    @test DT.is_boundary_edge(5, 2, adjacent(newtri2))
    @test DT.get_edge(adjacent2vertex(newtri2), 1) == Set{NTuple{2,Int64}}([
        (7, 5),
        (4, 6),
        (6, 3),
        (5, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 2) == Set{NTuple{2,Int64}}([
        (5, 3),
        (3, 6)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 3) == Set{NTuple{2,Int64}}([
        (5, 7),
        (1, 6),
        (2, 5),
        (6, 2)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 4) == Set{NTuple{2,Int64}}([
        (6, 1),
        (1, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 5) == Set{NTuple{2,Int64}}([
        (4, 1),
        (1, 7),
        (7, 3),
        (3, 2)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 6) == Set{NTuple{2,Int64}}([
        (3, 1),
        (2, 3),
        (1, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 7) == Set{NTuple{2,Int64}}([
        (5, 1),
        (3, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (4, 5),
        (5, 2),
        (2, 6),
        (6, 4)
    ])
    @test (1, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (3, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (7, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test DT.get_neighbour(graph(newtri2), 1) == Set{Int64}([
        4, 5, 7, 3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 2) == Set{Int64}([
        3, 6, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 3) == Set{Int64}([
        7, 1, 5, 6, 2
    ])
    @test DT.get_neighbour(graph(newtri2), 4) == Set{Int64}([
        6, 1, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 5) == Set{Int64}([
        7, 1, 4, 2, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 6) == Set{Int64}([
        4, 1, 3, 2
    ])
    @test DT.get_neighbour(graph(newtri2), 7) == Set{Int64}([
        5, 3, 1
    ])
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 5, 1) == 7
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == 2
    @test DT.get_edge(adjacent(newtri2), 3, 2) == 5
    @test DT.get_edge(adjacent(newtri2), 2, 3) == 6
    @test DT.get_edge(adjacent(newtri2), 7, 5) == 1
    @test DT.get_edge(adjacent(newtri2), 5, 7) == 3
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 5, 3) == 2
    @test DT.get_edge(adjacent(newtri2), 2, 5) == 3
    @test DT.get_edge(adjacent(newtri2), 5, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 2) == 3
    @test DT.get_edge(adjacent(newtri2), 2, 5) == 3
    @test DT.get_edge(adjacent(newtri2), 5, 2) == DT.BoundaryIndex

    DT.delete_triangle!(triangles(newtri2),
        adjacent(newtri2),
        adjacent2vertex(newtri2), graph(newtri2),
        3, 2, 5)

    @test triangles(newtri2) == Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (6, 2, 3),
        (5, 1, 7),
        (3, 5, 7),
        (6, 1, 4),
        (5, 4, 1)
    ])
    @test length(triangles(newtri2)) == 6
    @test (1, 3, 7) ∉ triangles(newtri2) &&
          (3, 7, 1) ∉ triangles(newtri2) &&
          (7, 1, 3) ∉ triangles(newtri2) &&
          (3, 2, 5) ∉ triangles(newtri2) &&
          (2, 5, 3) ∉ triangles(newtri2) &&
          (5, 3, 2) ∉ triangles(newtri2)
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == 2
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 3, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 3) == 6
    @test DT.get_edge(adjacent(newtri2), 6, 2) == 3
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.is_boundary_edge(2, 6, adjacent(newtri2))
    @test DT.is_boundary_edge(6, 4, adjacent(newtri2))
    @test DT.is_boundary_edge(4, 5, adjacent(newtri2))
    @test DT.is_boundary_edge(5, 3, adjacent(newtri2))
    @test DT.is_boundary_edge(3, 2, adjacent(newtri2))
    @test DT.get_edge(adjacent2vertex(newtri2), 1) == Set{NTuple{2,Int64}}([
        (7, 5),
        (4, 6),
        (6, 3),
        (5, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 2) == Set{NTuple{2,Int64}}([
        (3, 6)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 3) == Set{NTuple{2,Int64}}([
        (5, 7),
        (1, 6),
        (6, 2)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 4) == Set{NTuple{2,Int64}}([
        (6, 1),
        (1, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 5) == Set{NTuple{2,Int64}}([
        (4, 1),
        (1, 7),
        (7, 3),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 6) == Set{NTuple{2,Int64}}([
        (3, 1),
        (2, 3),
        (1, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 7) == Set{NTuple{2,Int64}}([
        (5, 1),
        (3, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (4, 5),
        (2, 6),
        (6, 4),
        (5, 3),
        (3, 2)
    ])
    @test (1, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (3, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (7, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test (3, 2) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (2, 5) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (5, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_neighbour(graph(newtri2), 1) == Set{Int64}([
        4, 5, 7, 3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 2) == Set{Int64}([
        3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 3) == Set{Int64}([
        7, 1, 5, 6, 2
    ])
    @test DT.get_neighbour(graph(newtri2), 4) == Set{Int64}([
        6, 1, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 5) == Set{Int64}([
        7, 1, 4, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 6) == Set{Int64}([
        4, 1, 3, 2
    ])
    @test DT.get_neighbour(graph(newtri2), 7) == Set{Int64}([
        5, 3, 1
    ])
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 5, 1) == 7
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 2, 3) == 6
    @test DT.get_edge(adjacent(newtri2), 3, 2) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 2, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 2) == 3
    @test DT.get_edge(adjacent(newtri2), 3, 6) == 2
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 7) == 3
    @test DT.get_edge(adjacent(newtri2), 7, 5) == 1

    DT.delete_triangle!(triangles(newtri2),
        adjacent(newtri2),
        adjacent2vertex(newtri2), graph(newtri2),
        3, 6, 2)

    @test triangles(newtri2) == Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (5, 1, 7),
        (3, 5, 7),
        (6, 1, 4),
        (5, 4, 1)
    ])
    @test length(triangles(newtri2)) == 5
    @test (1, 3, 7) ∉ triangles(newtri2) &&
          (3, 7, 1) ∉ triangles(newtri2) &&
          (7, 1, 3) ∉ triangles(newtri2) &&
          (3, 2, 5) ∉ triangles(newtri2) &&
          (2, 5, 3) ∉ triangles(newtri2) &&
          (5, 3, 2) ∉ triangles(newtri2) &&
          (3, 6, 2) ∉ triangles(newtri2) &&
          (6, 2, 3) ∉ triangles(newtri2) &&
          (2, 3, 6) ∉ triangles(newtri2)
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test !DT.is_valid_edge(7, 1, adjacent(newtri2))
    @test !DT.is_valid_edge(1, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(3, 7, adjacent(newtri2))
    @test DT.is_valid_edge(5, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(2, 5, adjacent(newtri2))
    @test !DT.is_valid_edge(3, 2, adjacent(newtri2))
    @test DT.is_valid_edge(3, 6, adjacent(newtri2))
    @test !DT.is_valid_edge(6, 2, adjacent(newtri2))
    @test !DT.is_valid_edge(2, 3, adjacent(newtri2))
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue
    @test DT.is_boundary_edge(6, 4, adjacent(newtri2))
    @test DT.is_boundary_edge(4, 5, adjacent(newtri2))
    @test DT.is_boundary_edge(5, 3, adjacent(newtri2))
    @test DT.is_boundary_edge(3, 6, adjacent(newtri2))
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent2vertex(newtri2), 1) == Set{NTuple{2,Int64}}([
        (7, 5),
        (4, 6),
        (6, 3),
        (5, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 2) == Set{NTuple{2,Int64}}([
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 3) == Set{NTuple{2,Int64}}([
        (5, 7),
        (1, 6),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 4) == Set{NTuple{2,Int64}}([
        (6, 1),
        (1, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 5) == Set{NTuple{2,Int64}}([
        (4, 1),
        (1, 7),
        (7, 3),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 6) == Set{NTuple{2,Int64}}([
        (3, 1),
        (1, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 7) == Set{NTuple{2,Int64}}([
        (5, 1),
        (3, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (4, 5),
        (6, 4),
        (5, 3),
        (3, 6)
    ])
    @test (1, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (3, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (7, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test (3, 2) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (2, 5) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (5, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test DT.get_neighbour(graph(newtri2), 1) == Set{Int64}([
        4, 5, 7, 3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 2) == Set{Int64}([
    ])
    @test DT.get_neighbour(graph(newtri2), 3) == Set{Int64}([
        7, 1, 5, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 4) == Set{Int64}([
        6, 1, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 5) == Set{Int64}([
        7, 1, 4, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 6) == Set{Int64}([
        4, 1, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 7) == Set{Int64}([
        5, 3, 1
    ])
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 7) == 5
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 5, 1) == 7
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 5, 7) == 3
    @test DT.get_edge(adjacent(newtri2), 7, 5) == 1
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue

    DT.delete_triangle!(triangles(newtri2),
        adjacent(newtri2),
        adjacent2vertex(newtri2), graph(newtri2),
        5, 1, 7)

    @test triangles(newtri2) == Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (3, 5, 7),
        (6, 1, 4),
        (5, 4, 1)
    ])
    @test length(triangles(newtri2)) == 4
    @test (1, 3, 7) ∉ triangles(newtri2) &&
          (3, 7, 1) ∉ triangles(newtri2) &&
          (7, 1, 3) ∉ triangles(newtri2) &&
          (3, 2, 5) ∉ triangles(newtri2) &&
          (2, 5, 3) ∉ triangles(newtri2) &&
          (5, 3, 2) ∉ triangles(newtri2) &&
          (3, 6, 2) ∉ triangles(newtri2) &&
          (6, 2, 3) ∉ triangles(newtri2) &&
          (2, 3, 6) ∉ triangles(newtri2) &&
          (5, 1, 7) ∉ triangles(newtri2) &&
          (1, 7, 5) ∉ triangles(newtri2) &&
          (7, 5, 1) ∉ triangles(newtri2)
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 1) == 6
    @test !DT.is_valid_edge(7, 1, adjacent(newtri2))
    @test !DT.is_valid_edge(1, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(3, 7, adjacent(newtri2))
    @test DT.is_valid_edge(5, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(2, 5, adjacent(newtri2))
    @test !DT.is_valid_edge(3, 2, adjacent(newtri2))
    @test DT.is_valid_edge(3, 6, adjacent(newtri2))
    @test !DT.is_valid_edge(6, 2, adjacent(newtri2))
    @test !DT.is_valid_edge(2, 3, adjacent(newtri2))
    @test !DT.is_valid_edge(5, 1, adjacent(newtri2))
    @test !DT.is_valid_edge(1, 7, adjacent(newtri2))
    @test !DT.is_valid_edge(7, 5, adjacent(newtri2))
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 6) == DT.DefaultAdjacentValue
    @test DT.is_boundary_edge(6, 4, adjacent(newtri2))
    @test DT.is_boundary_edge(4, 5, adjacent(newtri2))
    @test DT.is_boundary_edge(5, 3, adjacent(newtri2))
    @test DT.is_boundary_edge(3, 6, adjacent(newtri2))
    @test DT.get_edge(adjacent2vertex(newtri2), 1) == Set{NTuple{2,Int64}}([
        (4, 6),
        (6, 3),
        (5, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 2) == Set{NTuple{2,Int64}}([
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 3) == Set{NTuple{2,Int64}}([
        (5, 7),
        (1, 6),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 4) == Set{NTuple{2,Int64}}([
        (6, 1),
        (1, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 5) == Set{NTuple{2,Int64}}([
        (4, 1),
        (7, 3),
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 6) == Set{NTuple{2,Int64}}([
        (3, 1),
        (1, 4)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), 7) == Set{NTuple{2,Int64}}([
        (3, 5)
    ])
    @test DT.get_edge(adjacent2vertex(newtri2), DT.BoundaryIndex) == Set{NTuple{2,Int64}}([
        (4, 5),
        (6, 4),
        (5, 3),
        (3, 6)
    ])
    @test (1, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (3, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (7, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test (3, 2) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (2, 5) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test (5, 3) ∉ DT.get_edge(adjacent2vertex(newtri2), 3)
    @test (5, 1) ∉ DT.get_edge(adjacent2vertex(newtri2), 7)
    @test (1, 7) ∉ DT.get_edge(adjacent2vertex(newtri2), 5)
    @test (7, 5) ∉ DT.get_edge(adjacent2vertex(newtri2), 1)
    @test DT.get_neighbour(graph(newtri2), 1) == Set{Int64}([
        4, 5, 3, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 2) == Set{Int64}([
    ])
    @test DT.get_neighbour(graph(newtri2), 3) == Set{Int64}([
        7, 1, 5, 6
    ])
    @test DT.get_neighbour(graph(newtri2), 4) == Set{Int64}([
        6, 1, 5
    ])
    @test DT.get_neighbour(graph(newtri2), 5) == Set{Int64}([
        7, 1, 4, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 6) == Set{Int64}([
        4, 1, 3
    ])
    @test DT.get_neighbour(graph(newtri2), 7) == Set{Int64}([
        5, 3
    ])
    @test DT.get_edge(adjacent(newtri2), 1, 4) == 6
    @test DT.get_edge(adjacent(newtri2), 1, 6) == 3
    @test DT.get_edge(adjacent(newtri2), 1, 5) == 4
    @test DT.get_edge(adjacent(newtri2), 4, 1) == 5
    @test DT.get_edge(adjacent(newtri2), 6, 1) == 4
    @test DT.get_edge(adjacent(newtri2), 3, 1) == 6
    @test DT.get_edge(adjacent(newtri2), 3, 5) == 7
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adjacent(newtri2), 3, 6) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 3) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 7, 3) == 5
    @test DT.get_edge(adjacent(newtri2), 6, 3) == 1
    @test DT.get_edge(adjacent(newtri2), 4, 6) == 1
    @test DT.get_edge(adjacent(newtri2), 6, 4) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 4, 5) == DT.BoundaryIndex
    @test DT.get_edge(adjacent(newtri2), 5, 4) == 1
    @test DT.get_edge(adjacent(newtri2), 5, 7) == 3
    @test DT.get_edge(adj, 1, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 3, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 5) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 2) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 2, 3) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 5, 1) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 1, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 7, 5) == DT.DefaultAdjacentValue
end

############################################
##
## ADDING A POINT INTO A TRIANGULATION WITH BOWYER-WATSON 
##
############################################
#=
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
tri = triangulate(pts; randomise=false)
T, adj, adj2v, DG = triangles(tri), adjacent(tri),
adjacent2vertex(tri), graph(tri)
p11 = (6.0, 2.5)
push!(pts, p11)
@test DT.locate_triangle(T, pts, 11) == (10, 8, 9)
DT.add_point_bowyer!(T, adj,adj2v,DG,pts,11)
=#