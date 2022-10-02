############################################
##
## MAIN TRIANGULATION FUNCTIONS
##
############################################
@testset "Initialisation" begin
    p1 = Point(2.7, 13.0)
    p2 = Point(-2.0, 3.0)
    p3 = Point(17.0, 2.0)
    p4 = Point(18.0, 3.5)
    pts = Points(p1, p2, p3, p4)
    Tri = DT.initialise_triangulation(pts)
    @test adjacent(Tri) == Tri.adjacent
    @test adjacent2vertex(Tri) == Tri.adjacent2vertex
    @test graph(Tri) == Tri.graph
    @test history(Tri) == Tri.history
    @test triangles(Tri) == Tri.triangles
    @test points(Tri) == Tri.points
    @test DT.root(Tri) == Tri.root
    @test DT.root(Tri) == DT.BoundingTriangle
    @test triangles(Tri).triangles == Set{Triangle{Int64}}([DT.BoundingTriangle])
    @test adjacent(Tri, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test adjacent(Tri, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test adjacent(Tri, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri, Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)) == DT.UpperBoundingIndex
    @test adjacent(Tri, Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri, Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.LowerRightBoundingIndex
    @test adjacent(Tri, Edge{Int64}(DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri, Edge{Int64}(DT.UpperBoundingIndex, DT.LowerRightBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri, Edge{Int64}(DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.BoundaryIndex
    @test adjacent2vertex(Tri, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)])
    @test adjacent2vertex(Tri, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)])
    @test adjacent2vertex(Tri, DT.UpperBoundingIndex) == Set{Edge{Int64}}([Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)])
    @test neighbours(Tri, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex])
    @test history(Tri).graph.N[DT.BoundingTriangle] == history(Tri).graph.N[DT.root(Tri)] == Set()
    @test points(Tri) == pts
    @test adjacent(Tri) isa DT.Adjacent{Int64,Edge{Int64}}
    @test adjacent2vertex(Tri) isa DT.Adjacent2Vertex{Int64,Edge{Int64}}
    @test graph(Tri) isa DT.DelaunayGraph{Int64}
    @test history(Tri) isa DT.HistoryDAG{Int64,Triangle{Int64}}
    @test triangles(Tri) isa Triangles{Int64,Triangle{Int64}}
    @test points(Tri) isa Points{Float64,Vector{Float64},Point{Float64,Vector{Float64}}}
    @test DT.root(Tri) isa Triangle{Int64}
    @test Tri isa Triangulation{typeof(adjacent(Tri)),typeof(adjacent2vertex(Tri)),
        typeof(graph(Tri)),typeof(history(Tri)),typeof(triangles(Tri)),
        typeof(points(Tri)),typeof(DT.root(Tri))}

    pts = [p1, p2, p3, p4]
    Tri2 = DT.initialise_triangulation(pts)
    @test Tri2.adjacent.adjacent == Tri.adjacent.adjacent
    @test Tri2.points.points == Tri.points.points
    @test Tri2.adjacent2vertex.adjacent2vertex == Tri.adjacent2vertex.adjacent2vertex
    @test Tri2.graph.graph == Tri.graph.graph
    @test Tri2.history.graph == Tri.history.graph
    @test Tri2.triangles.triangles == Tri.triangles.triangles
    @test Tri2.root == Tri.root

    Tri3 = DT.initialise_triangulation(pts; IntegerType=Int16)
    @test adjacent(Tri3) == Tri3.adjacent
    @test adjacent2vertex(Tri3) == Tri3.adjacent2vertex
    @test graph(Tri3) == Tri3.graph
    @test history(Tri3) == Tri3.history
    @test triangles(Tri3) == Tri3.triangles
    @test points(Tri3) == Tri3.points
    @test DT.root(Tri3) == Tri3.root
    @test DT.root(Tri3) == Triangle{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    @test triangles(Tri3).triangles == Set{Triangle{Int16}}([Triangle{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)])
    @test adjacent(Tri3, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex) == DT.UpperBoundingIndex
    @test adjacent(Tri3, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri3, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex) == DT.LowerRightBoundingIndex
    @test adjacent(Tri3, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex) == DT.BoundaryIndex
    @test adjacent(Tri3, Edge{Int16}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)) == DT.UpperBoundingIndex
    @test adjacent(Tri3, Edge{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)) == DT.LowerLeftBoundingIndex
    @test adjacent(Tri3, Edge{Int16}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.LowerRightBoundingIndex
    @test adjacent(Tri3, Edge{Int16}(DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri3, Edge{Int16}(DT.UpperBoundingIndex, DT.LowerRightBoundingIndex)) == DT.BoundaryIndex
    @test adjacent(Tri3, Edge{Int16}(DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex)) == DT.BoundaryIndex
    @test adjacent2vertex(Tri3, DT.LowerLeftBoundingIndex) == Set{Edge{Int16}}([Edge{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)])
    @test adjacent2vertex(Tri3, DT.LowerRightBoundingIndex) == Set{Edge{Int16}}([Edge{Int16}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)])
    @test adjacent2vertex(Tri3, DT.UpperBoundingIndex) == Set{Edge{Int16}}([Edge{Int16}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)])
    @test neighbours(Tri3, DT.LowerLeftBoundingIndex) == Set{Int16}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri3, DT.LowerRightBoundingIndex) == Set{Int16}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test neighbours(Tri3, DT.UpperBoundingIndex) == Set{Int16}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex])
    @test history(Tri3).graph.N[Triangle{Int16}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)] == history(Tri3).graph.N[DT.root(Tri3)] == Set()
    @test points(Tri3).points == Points(pts).points
    @test adjacent(Tri3) isa DT.Adjacent{Int16,Edge{Int16}}
    @test adjacent2vertex(Tri3) isa DT.Adjacent2Vertex{Int16,Edge{Int16}}
    @test graph(Tri3) isa DT.DelaunayGraph{Int16}
    @test history(Tri3) isa DT.HistoryDAG{Int16,Triangle{Int16}}
    @test triangles(Tri3) isa Triangles{Int16,Triangle{Int16}}
    @test points(Tri3) isa Points{Float64,Vector{Float64},Point{Float64,Vector{Float64}}}
    @test DT.root(Tri3) isa Triangle{Int16}
    @test Tri3 isa Triangulation{typeof(adjacent(Tri3)),typeof(adjacent2vertex(Tri3)),
        typeof(graph(Tri3)),typeof(history(Tri3)),typeof(triangles(Tri3)),
        typeof(points(Tri3)),typeof(DT.root(Tri3))}
end

@testset "Testing that the point location data structure can correctly locate triangles" begin
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)

    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))

    p = Point(33.373, 15.2287)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(-31.0689, 52.90257)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(3.63, 1.679)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    @test flag == 1

    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))

    p = Point(3.63, 1.679)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = Point(27.706, 0.968)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = Point(-13.6689, 1.3567)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 2, 1)
    @test flag == 1
    p = Point(-3.56804, 1.745279)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 2, 1)
    @test flag == 1
    p = Point(5.95, -2.91669)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = Point(32.9507, -4.2764)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = Point(41.4976, 46.81)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(10.0, 10.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(-33.48, 23.11)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(-10.0, 10.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1

    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))

    p = Point(-10.0, 10.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(-36.59, 13.594)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(35.86, 34.379)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(15.66, 7.766)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    @test flag == 1
    p = Point(5.173, -3.305)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = Point(-6.09, -3.305)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    p = Point(14.57, 0.48686)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = Point(9.89, 2.198)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(2, DT.LowerRightBoundingIndex, 1)
    @test flag == 1
    p = Point(-5.735, 3.11)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = Point(-3.7957, 3.11)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = Point(-11.21, 2.54)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    @test flag == 1
    p = Point(-3.68, 1.057)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 2, 3)
    @test flag == 1
    p = Point(0.0, 0.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test flag == 0
    p = Point(0.916, 1.408)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(3, 2, 1)
    @test flag == 1
    p = Point(0.0, 2.0)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(3, 2, 1)
    @test flag == 1
    p = Point(2.5057, 2.8986)
    intri, flag = DT.locate_triangle(HG, pts, p, DT.BoundingTriangle)
    @test intri == Triangle(3, 2, 1)
    @test flag == 1
end

@testset "Testing that we can correctly add a point into a triangulation" begin
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    # Building an example triangulation
    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    T = Triangles(
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    )
    adj = DT.Adjacent()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex()
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

    # Now do the actual test 
    DT.add_point!(T, HG, adj, adj2v, DG, Triangle(DT.LowerLeftBoundingIndex, 3, 1), 4)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 1) ∉ triangles(T)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 4) ∈ triangles(T)
    @test Triangle(3, 1, 4) ∈ triangles(T)
    @test Triangle(1, DT.LowerLeftBoundingIndex, 4) ∈ triangles(T)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(3, 1, 4),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 1)) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(3, 1, 4),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(3, 1, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(1, DT.LowerLeftBoundingIndex, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
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
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 3),
        Edge{Int64}(3, 4),
        Edge{Int64}(4, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(1, 2),
        Edge{Int64}(2, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(3, 2),
        Edge{Int64}(4, 3),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 4),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(3, DT.LowerLeftBoundingIndex),
        Edge{Int64}(1, 3),
        Edge{Int64}(DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(4, DT.LowerLeftBoundingIndex),
        Edge{Int64}(1, 4),
        Edge{Int64}(2, 1)
    ])
    @test DT.get_edge(adj2v, 4) == Set{Edge{Int64}}([
        Edge{Int64}(1, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 3),
        Edge{Int64}(3, 1)
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
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    # Building an example triangulation
    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    T = Triangles(
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    )
    adj = DT.Adjacent()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex()
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
    DT.add_point!(T, HG, adj, adj2v, DG, Triangle(DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    DT.flip_edge!(T, HG, adj, adj2v, DG, i, j, k, r)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 1) ∉ triangles(T)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 4) ∈ triangles(T)
    @test Triangle(3, 1, 4) ∉ triangles(T)
    @test Triangle(4, 3, 2) ∈ triangles(T)
    @test Triangle(1, DT.LowerLeftBoundingIndex, 4) ∈ triangles(T)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(4, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(4, 3, 2),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 1)) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(3, 1, 4),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(3, 1, 4)) == Set{Triangle{Int64}}([
        Triangle(4, 2, 1), Triangle(4, 3, 2)
    ])
    @test out_neighbors(HG, Triangle(3, 2, 1)) == Set{Triangle{Int64}}([
        Triangle(4, 2, 1), Triangle(4, 3, 2)
    ])
    p = Point(1.88, 2.5834)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(4, 2, 1)
    @test flag == 1
    p = Point(-0.9802, 2.5834)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(4, 3, 2)
    @test flag == 1
    p = Point(-3.642, 2.8615)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 4)
    @test flag == 1
    p = Point(4.3036, -3.0977)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(3, 1, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(1, DT.LowerLeftBoundingIndex, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(4, 3, 2)) == [Triangle(3, 2, 1), Triangle(3, 1, 4)] || in_neighbors(HG, Triangle(4, 3, 2)) == [Triangle(3, 1, 4), Triangle(3, 2, 1)]
    @test in_neighbors(HG, Triangle(4, 2, 1)) == [Triangle(3, 2, 1), Triangle(3, 1, 4)] || in_neighbors(HG, Triangle(4, 2, 1)) == [Triangle(3, 1, 4), Triangle(3, 2, 1)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test_throws KeyError DT.get_edge(adj, 3, 1) == 4
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
    @test_throws KeyError DT.get_edge(adj, 1, 3)
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 3),
        Edge{Int64}(3, 4),
        Edge{Int64}(4, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(1, 2),
        Edge{Int64}(2, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(4, 2),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 4),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(3, DT.LowerLeftBoundingIndex),
        Edge{Int64}(1, 4),
        Edge{Int64}(DT.LowerRightBoundingIndex, 1),
        Edge{Int64}(4, 3)
    ])
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(4, DT.LowerLeftBoundingIndex),
        Edge{Int64}(2, 4)
    ])
    @test DT.get_edge(adj2v, 4) == Set{Edge{Int64}}([
        Edge{Int64}(1, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 3),
        Edge{Int64}(3, 2),
        Edge{Int64}(2, 1)
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
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    # Building an example triangulation
    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    T = Triangles(
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    )
    adj = DT.Adjacent()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex()
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
    DT.add_point!(T, HG, adj, adj2v, DG, Triangle(DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    @test !DT.islegal(i, j, adj, pts)
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 1) ∉ triangles(T)
    @test Triangle(DT.LowerLeftBoundingIndex, 3, 4) ∈ triangles(T)
    @test Triangle(3, 1, 4) ∉ triangles(T)
    @test Triangle(4, 3, 2) ∈ triangles(T)
    @test Triangle(1, DT.LowerLeftBoundingIndex, 4) ∈ triangles(T)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(4, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(4, 3, 2),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 1)) == Set{Triangle{Int64}}([
        Triangle(DT.LowerLeftBoundingIndex, 3, 4),
        Triangle(3, 1, 4),
        Triangle(1, DT.LowerLeftBoundingIndex, 4)
    ])
    @test out_neighbors(HG, Triangle(3, 1, 4)) == Set{Triangle{Int64}}([
        Triangle(4, 2, 1), Triangle(4, 3, 2)
    ])
    @test out_neighbors(HG, Triangle(3, 2, 1)) == Set{Triangle{Int64}}([
        Triangle(4, 2, 1), Triangle(4, 3, 2)
    ])
    p = Point(1.88, 2.5834)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(4, 2, 1)
    @test flag == 1
    p = Point(-0.9802, 2.5834)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(4, 3, 2)
    @test flag == 1
    p = Point(-3.642, 2.8615)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, 3, 4)
    @test flag == 1
    p = Point(4.3036, -3.0977)
    intri, flag = DT.locate_triangle(HG, pts, p)
    @test intri == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)
    @test flag == 1
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, 3, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(3, 1, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(1, DT.LowerLeftBoundingIndex, 4)) == [Triangle(DT.LowerLeftBoundingIndex, 3, 1)]
    @test in_neighbors(HG, Triangle(4, 3, 2)) == [Triangle(3, 2, 1), Triangle(3, 1, 4)] || in_neighbors(HG, Triangle(4, 3, 2)) == [Triangle(3, 1, 4), Triangle(3, 2, 1)]
    @test in_neighbors(HG, Triangle(4, 2, 1)) == [Triangle(3, 2, 1), Triangle(3, 1, 4)] || in_neighbors(HG, Triangle(4, 2, 1)) == [Triangle(3, 1, 4), Triangle(3, 2, 1)]
    @test DT.get_edge(adj, DT.LowerLeftBoundingIndex, 3) == 4
    @test DT.get_edge(adj, 3, DT.LowerLeftBoundingIndex) == 2
    @test DT.get_edge(adj, 3, 4) == DT.LowerLeftBoundingIndex
    @test DT.get_edge(adj, 4, DT.LowerLeftBoundingIndex) == 3
    @test_throws KeyError DT.get_edge(adj, 3, 1) == 4
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
    @test_throws KeyError DT.get_edge(adj, 1, 3)
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 3),
        Edge{Int64}(3, 4),
        Edge{Int64}(4, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(1, 2),
        Edge{Int64}(2, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(4, 2),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 4),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(3, DT.LowerLeftBoundingIndex),
        Edge{Int64}(1, 4),
        Edge{Int64}(DT.LowerRightBoundingIndex, 1),
        Edge{Int64}(4, 3)
    ])
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(4, DT.LowerLeftBoundingIndex),
        Edge{Int64}(2, 4)
    ])
    @test DT.get_edge(adj2v, 4) == Set{Edge{Int64}}([
        Edge{Int64}(1, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 3),
        Edge{Int64}(3, 2),
        Edge{Int64}(2, 1)
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
    p1 = Point(5.0, 5.0)
    p2 = Point(1.0, -1.0)
    p3 = Point(-2.0, 2.0)
    p4 = Point(-1.0, 4.0)
    p5 = Point(2.0, 3.0)
    pts = Points(p1, p2, p3, p4, p5)

    # Building an example triangulation
    HG = DT.HistoryDAG()
    DT.add_triangle!(HG, DT.BoundingTriangle)
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_edge!(HG, DT.BoundingTriangle, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(2, DT.LowerRightBoundingIndex, 1))
    DT.add_triangle!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    DT.add_edge!(HG, Triangle(DT.LowerLeftBoundingIndex, 2, 1),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1),
        Triangle(3, 2, 1))
    T = Triangles(
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(2, DT.LowerRightBoundingIndex, 1),
        Triangle(3, 2, 1),
        Triangle(1, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 1, DT.UpperBoundingIndex),
        Triangle(DT.LowerLeftBoundingIndex, 3, 1)
    )
    adj = DT.Adjacent()
    [DT.add_triangle!(adj, T) for T in T]
    DT.add_edge!(adj, DT.LowerRightBoundingIndex, DT.LowerLeftBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.UpperBoundingIndex, DT.LowerRightBoundingIndex, DT.BoundaryIndex)
    DT.add_edge!(adj, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, DT.BoundaryIndex)
    adj2v = DT.Adjacent2Vertex()
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
    DT.add_point!(T, HG, adj, adj2v, DG, Triangle(DT.LowerLeftBoundingIndex, 3, 1), 4)
    i, j = 3, 1
    r, k = 4, 2
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    DT.remove_bounding_triangle!(T, adj, adj2v, DG)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(4, 3, 2),
        Triangle(4, 2, 1)
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
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(2, 4)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(4, 2)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(1, 4),
        Edge{Int64}(4, 3)
    ])
    @test DT.get_edge(adj2v, 4) == Set{Edge{Int64}}([
        Edge{Int64}(3, 2),
        Edge{Int64}(2, 1)
    ])
    @test DT.get_edge(adj2v, DT.BoundaryIndex) == Set{Edge{Int64}}([
        Edge{Int64}(1, 2),
        Edge{Int64}(4, 1),
        Edge{Int64}(3, 4),
        Edge{Int64}(2, 3)
    ])
    @test length(adj2v.adjacent2vertex) == 5
    @test DT.get_neighbour(DG, 1) == Set{Int64}([4, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([1, 4, 3])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([4, 2])
    @test DT.get_neighbour(DG, 4) == Set{Int64}([2, 3, 1])
end

@testset "Testing the steps for a small problem" begin
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
    @test Tᵢⱼₖ == DT.BoundingTriangle
    @test interior_flag == 1
    i, j, k = Tᵢⱼₖ
    @test i == DT.LowerRightBoundingIndex
    @test j == DT.UpperBoundingIndex
    @test k == DT.LowerLeftBoundingIndex
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
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
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.UpperBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex, 1])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, DT.UpperBoundingIndex])
    @test out_neighbors(HG, DT.BoundingTriangle) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test in_neighbors(HG, Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == [DT.BoundingTriangle]

    r = 2
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    @test Tᵢⱼₖ == Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    @test interior_flag == 1
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(DT.LowerRightBoundingIndex, 1, 2),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        Triangle(1, DT.LowerLeftBoundingIndex, 2)
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
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 1),
        Edge{Int64}(1, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.UpperBoundingIndex, 1),
        Edge{Int64}(1, 2),
        Edge{Int64}(2, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 1),
        Edge{Int64}(1, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1])
    @test out_neighbors(HG, DT.BoundingTriangle) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)
    ])
    @test in_neighbors(HG, Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, Triangle(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == [DT.BoundingTriangle]
    @test out_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, 1, 2),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(1, DT.LowerLeftBoundingIndex, 2)
    ])
    @test in_neighbors(HG, Triangle(DT.LowerRightBoundingIndex, 1, 2)) == [Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test in_neighbors(HG, Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2)) == [Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test in_neighbors(HG, Triangle(1, DT.LowerLeftBoundingIndex, 2)) == [Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1)]
    @test DT.islegal(i, j, adj, pts)
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)

    r = 3
    pᵣ = DT.get_point(pts, r)
    Tᵢⱼₖ, interior_flag = DT.locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
    @test !DT.islegal(i, j, adj, pts)
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)
    DT.legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    @test Tᵢⱼₖ == Triangle(1, DT.LowerLeftBoundingIndex, 2)
    @test interior_flag == 1
    @test triangles(T) == Set{Triangle{Int64}}([
        Triangle(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex, 1),
        Triangle(3, 1, DT.UpperBoundingIndex),
        Triangle(2, 1, 3),
        Triangle(DT.LowerRightBoundingIndex, 1, 2),
        Triangle(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 2),
        Triangle(DT.LowerLeftBoundingIndex, 2, 3),
        Triangle(3, DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex)
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
    @test DT.get_edge(adj2v, DT.LowerLeftBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, 2),
        Edge{Int64}(2, 3),
        Edge{Int64}(3, DT.UpperBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.LowerRightBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.UpperBoundingIndex, 1)
        Edge{Int64}(1, 2)
        Edge{Int64}(2, DT.LowerLeftBoundingIndex)
    ])
    @test DT.get_edge(adj2v, DT.UpperBoundingIndex) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerLeftBoundingIndex, 3),
        Edge{Int64}(3, 1),
        Edge{Int64}(1, DT.LowerRightBoundingIndex)
    ])
    @test DT.get_edge(adj2v, 1) == Set{Edge{Int64}}([
        Edge{Int64}(DT.LowerRightBoundingIndex, DT.UpperBoundingIndex),
        Edge{Int64}(3, 2),
        Edge{Int64}(2, DT.LowerRightBoundingIndex),
        Edge{Int64}(DT.UpperBoundingIndex, 3)
    ])
    @test DT.get_edge(adj2v, 2) == Set{Edge{Int64}}([
        Edge{Int64}(1, 3),
        Edge{Int64}(3, DT.LowerLeftBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex),
        Edge{Int64}(DT.LowerRightBoundingIndex, 1)
    ])
    @test DT.get_edge(adj2v, 3) == Set{Edge{Int64}}([
        Edge{Int64}(1, DT.UpperBoundingIndex),
        Edge{Int64}(DT.LowerLeftBoundingIndex, 2),
        Edge{Int64}(DT.UpperBoundingIndex, DT.LowerLeftBoundingIndex),
        Edge{Int64}(2, 1)
    ])
    @test DT.get_neighbour(DG, DT.LowerLeftBoundingIndex) == Set{Int64}([DT.LowerRightBoundingIndex, 2, 3, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.LowerRightBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, 2, 1, DT.UpperBoundingIndex])
    @test DT.get_neighbour(DG, DT.UpperBoundingIndex) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3])
    @test DT.get_neighbour(DG, 1) == Set{Int64}([DT.UpperBoundingIndex, 3, DT.LowerRightBoundingIndex, 2])
    @test DT.get_neighbour(DG, 2) == Set{Int64}([DT.LowerLeftBoundingIndex, DT.LowerRightBoundingIndex, 1, 3])
    @test DT.get_neighbour(DG, 3) == Set{Int64}([1, 2, DT.LowerLeftBoundingIndex, DT.UpperBoundingIndex])
    @test DT.islegal(j, k, adj, pts)
    @test DT.islegal(k, i, adj, pts)

    DTri = triangulate(pts)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false)
    @test DT.is_delaunay(DTri)
end

@testset "Random triangulations" begin
    p1 = Point(-6.88, 3.61)
    p2 = Point(-6.08, -0.43)
    p3 = Point(-0.3, 2.01)
    p4 = Point(5.1, -1.27)
    p5 = Point(6.18, 1.87)
    p6 = Point(3.08, 4.43)
    p7 = Point(-1.34, 4.83)
    p8 = Point(-1.68, -0.77)
    pts = Points(p1, p2, p3, p4, p5, p6, p7, p8)
    DTri = triangulate(pts)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false)
    @test DT.is_delaunay(DTri)
    pts = [[-6.88, 3.61], [-6.08, -0.43], [-0.3, 2.01], [5.1, -1.27], [6.18, 1.87], [3.08, 4.43], [-1.34, 4.83], [-1.68], [-0.77]]
    DTri = triangulate(pts)
    @test DT.is_delaunay(DTri)
    DTri = triangulate(pts; trim=false)
    @test DT.is_delaunay(DTri)

    pts = Points(p1, p2, p3, p4, p5, p6, p7, p8)
    DTri = triangulate(pts; shuffle_pts=false, trim=false)
    DTri2 = triangulate(Points(p1, p2, p3, p4, p5, p6, p7); shuffle_pts=false, trim=false)
    add_point!(DTri2, [-1.68, -0.77])
    @test triangles(DTri2).triangles == triangles(DTri).triangles

    pts = Points(p1, p2, p3, p4, p5, p6, p7, p8)
    DTri = triangulate(pts; shuffle_pts=false, trim=true)
    DTri2 = triangulate(Points(p1, p2, p3, p4, p5, p6, p7); shuffle_pts=false, trim=false)
    add_point!(DTri2, [-1.68, -0.77])
    DT.remove_bounding_triangle!(DTri2)
    @test triangles(DTri2).triangles == triangles(DTri).triangles

    for _ in 1:10000
        x = rand(100)
        y = rand(100)
        pts = Points([Point(x, y) for (x, y) in zip(x, y)])
        DTri = triangulate(pts)
        @test DT.is_delaunay(DTri)
        DTri = triangulate(pts; trim=false)
        @test DT.is_delaunay(DTri)
        @test 1 == num_points(DTri) - num_edges(DTri) + num_triangles(DTri) # Euler's formula
    end
end

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1])
x = rand(100)
y = rand(100)
pts = Points([Point(x, y) for (x, y) in zip(x, y)])
DTri = triangulate(pts)
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
    p1 = Point(-6.88, 3.61)
    p2 = Point(-6.08, -0.43)
    p3 = Point(-0.3, 2.01)
    p4 = Point(5.1, -1.27)
    p5 = Point(6.18, 1.87)
    p6 = Point(3.08, 4.43)
    p7 = Point(-1.34, 4.83)
    p8 = Point(-1.68, -0.77)
    pts = Points(p1, p2, p3, p4, p5, p6, p7, p8)
    DTri = triangulate(pts; IntegerType=Int16, shuffle_pts=false)
    @test points(DTri).points == Points(pts).points
    @test adjacent(DTri) isa DT.Adjacent{Int16,Edge{Int16}}
    @test adjacent2vertex(DTri) isa DT.Adjacent2Vertex{Int16,Edge{Int16}}
    @test graph(DTri) isa DT.DelaunayGraph{Int16}
    @test history(DTri) isa DT.HistoryDAG{Int16,Triangle{Int16}}
    @test triangles(DTri) isa Triangles{Int16,Triangle{Int16}}
    @test points(DTri) isa Points{Float64,Vector{Float64},Point{Float64,Vector{Float64}}}
    @test DT.root(DTri) isa Triangle{Int16}
    @test DTri isa Triangulation{typeof(adjacent(DTri)),typeof(adjacent2vertex(DTri)),
        typeof(graph(DTri)),typeof(history(DTri)),typeof(triangles(DTri)),
        typeof(points(DTri)),typeof(DT.root(DTri))}
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