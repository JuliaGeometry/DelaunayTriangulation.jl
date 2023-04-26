using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StructEquality
using Setfield

@struct_equal Triangulation
@struct_equal DT.Adjacent
@struct_equal DT.Adjacent2Vertex
@struct_equal DT.Graph
include("../helper_functions.jl")

global pts = rand(2, 500)
global tri = Triangulation(pts; IntegerType=Int32)

throw_f = expr -> @static if VERSION < v"1.8"
      ErrorException(expr)
  else
      expr
  end

@testset "Initialising a triangulation" begin
      @test tri == Triangulation(pts,
            Set{NTuple{3,Int32}}(),
            DT.Adjacent{Int32,NTuple{2,Int32}}(),
            DT.Adjacent2Vertex{Int32,Set{NTuple{2,Int32}},NTuple{2,Int32}}(),
            DT.Graph{Int32}(),
            Int32[],
            Dict{Tuple{Int32,Int32},Tuple{Vector{Int64},Int64}}(),
            OrderedDict{Int32,Vector{Int32}}(Int32(DT.BoundaryIndex) => Int32[]),
            OrderedDict{Int32,UnitRange{Int32}}(-1 => -1:-1),
            Set{NTuple{2,Int32}}(),
            Set{NTuple{2,Int32}}(),
            ConvexHull(pts, Int32[]),
            DT.get_empty_representative_points())
end

_x, _y = complicated_geometry()
global x = _x
global y = _y
boundary_nodes, points = convert_boundary_points_to_indices(x, y)
global tri = triangulate(points; boundary_nodes, delete_ghosts=false)
A = get_total_area(tri)
refine!(tri; max_area=1e-1A)
boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
global tri_2 = triangulate(points; boundary_nodes, delete_ghosts=false)
A = get_total_area(tri_2)
refine!(tri_2; max_area=1e-1A)

boundary_nodes, points = convert_boundary_points_to_indices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
global tri_3 = triangulate(points; boundary_nodes, delete_ghosts=false)
A = get_total_area(tri_3)
refine!(tri_3; max_area=1e-1A)
boundary_nodes, points = convert_boundary_points_to_indices(reverse(reverse.(x[2])), reverse(reverse.(y[2])))
global tri_4 = triangulate(points; boundary_nodes, delete_ghosts=false)
A = get_total_area(tri_4)
refine!(tri_4; max_area=1e-1A)

@testset "Triangulation getters" begin
      @test DT.get_points(tri) == tri.points
      @test DT.get_triangles(tri) == tri.triangles
      @test DT.get_adjacent(tri) == tri.adjacent
      @test DT.get_adjacent2vertex(tri) == tri.adjacent2vertex
      @test DT.get_graph(tri) == tri.graph
      @test DT.get_boundary_nodes(tri) == tri.boundary_nodes
      @test DT.get_boundary_map(tri) == tri.boundary_map
      @test DT.get_constrained_edges(tri) == tri.constrained_edges
      @test DT.get_convex_hull(tri) == tri.convex_hull
      @test compare_edge_vectors(collect(DT.get_all_constrained_edges(tri)), collect(tri.all_constrained_edges))
      @test compare_edge_vectors(collect(DT.get_all_constrained_edges(tri)), collect(DT.merge_constrained_edges(get_boundary_map(tri), get_boundary_nodes(tri), get_constrained_edges(tri))))
      @test DT.get_boundary_edge_map(tri) == tri.boundary_edge_map == DT.construct_boundary_edge_map(get_boundary_nodes(tri))
      @inferred DT.get_points(tri)
      @inferred DT.get_triangles(tri)
      @inferred DT.get_adjacent(tri)
      @inferred DT.get_adjacent2vertex(tri)
      @inferred DT.get_graph(tri)
      @inferred DT.get_boundary_nodes(tri)
      @inferred DT.get_boundary_map(tri)
      @inferred DT.get_constrained_edges(tri)
      @inferred DT.get_convex_hull(tri)
      @inferred DT.get_all_constrained_edges(tri)
      @inferred DT.get_boundary_edge_map(tri)
end

@testset "Forwarded methods" begin
      @testset "Adjacent" begin
            @test DT.get_adjacent(tri, 92, -6) ==
                  DT.get_adjacent(tri, 92, -6; check_existence=Val(true)) ==
                  tri.adjacent.adjacent[(92, -6)]
            @test DT.get_adjacent(tri, (723, 1356)) ==
                  DT.get_adjacent(tri, (723, 1356); check_existence=Val(true)) ==
                  tri.adjacent.adjacent[(723, 1356)]
            @inferred DT.get_adjacent(tri, (723, 1356))
            @inferred DT.get_adjacent(tri, (723, 1356); check_existence=Val(true))
            DT.add_adjacent!(tri, 101117, 20311, 5)
            DT.add_adjacent!(tri, (27311, 50511), 10)
            @test DT.get_adjacent(tri, 101117, 20311) == 5
            @inferred DT.get_adjacent(tri, 101117, 20311)
            @inferred DT.get_adjacent(tri, 101117, 20311; check_existence=Val(true))
            @test DT.get_adjacent(tri, 27311, 50511) == 10
            DT.delete_adjacent!(tri, 101117, 20311)
            DT.delete_adjacent!(tri, (27311, 50511))
            @test DT.get_adjacent(tri, 101117, 20311) == DT.DefaultAdjacentValue
            @test DT.get_adjacent(tri, 27311, 50511) == DT.DefaultAdjacentValue
            @test DT.get_adjacent(tri, 711, DT.BoundaryIndex) == DT.DefaultAdjacentValue
      end

      @testset "Adjacent2Vertex" begin
            @test DT.get_adjacent2vertex(tri, -6) == tri.adjacent2vertex.adjacent2vertex[-6]
            @inferred DT.get_adjacent2vertex(tri, -6)
            DT.add_adjacent2vertex!(tri, 10005, 23, 27)
            DT.add_adjacent2vertex!(tri, 19988, (50, 171))
            @test DT.contains_edge(23, 27, tri.adjacent2vertex.adjacent2vertex[10005])
            @test DT.contains_edge(50, 171, tri.adjacent2vertex.adjacent2vertex[19988])
            @inferred DT.contains_edge(50, 171, tri.adjacent2vertex.adjacent2vertex[19988])
            DT.delete_adjacent2vertex!(tri, 10005, 23, 27)
            DT.delete_adjacent2vertex!(tri, 19988, (50, 171))
            @test !DT.contains_edge(23, 27, tri.adjacent2vertex.adjacent2vertex[10005])
            @test !DT.contains_edge(50, 171, tri.adjacent2vertex.adjacent2vertex[19988])
            @test 6 ∈ keys(tri.adjacent2vertex.adjacent2vertex)
            DT.delete_adjacent2vertex!(tri, -6)
            @test -6 ∉ keys(tri.adjacent2vertex.adjacent2vertex)
      end

      @testset "Graph" begin
            @test DT.get_edges(tri) == tri.graph.graph.E == DT.get_edges(tri.graph)
            @inferred DT.get_edges(tri)
            @test DT.get_neighbours(tri, 3) == tri.graph.graph.N[3]
            @inferred DT.get_neighbours(tri, 3)
            DT.add_vertex!(tri, 19998, 23721)
            DT.add_vertex!(tri, 28371)
            @test all(∈(tri.graph.graph.V), (19998, 23721, 28371))
            DT.add_neighbour!(tri, 28371, 50912)
            DT.add_neighbour!(tri, 28371, 271, 501)
            @test all(∈(tri.graph.graph.N[28371]), (50912, 271, 501))
            DT.delete_neighbour!(tri, 28371, 50912)
            DT.delete_neighbour!(tri, 28371, 271, 501)
            @test all(∉(tri.graph.graph.N[28371]), (50912, 271, 501))
            DT.delete_vertex!(tri, 19998)
            DT.delete_vertex!(tri, 28371, 3)
            @test all(∉(tri.graph.graph.V), (19998, 28371, 3))
            DT.delete_boundary_vertices_from_graph!(tri)
            @test all(∉(tri.graph.graph.V), -11:-1)
            @test DT.get_neighbours(tri) == tri.graph.graph.N
            @test DT.get_vertices(tri) == tri.graph.graph.V
            @test get_vertices(tri) == each_vertex(tri)
            @test num_vertices(tri) == length(get_vertices(tri))
            @test all((DT.num_neighbours(tri, u) == DT.num_neighbours(get_graph(tri), u) for u in DT.get_vertices(tri)))
      end

      @testset "ConvexHull" begin
            @test DT.get_convex_hull(tri) == tri.convex_hull
            @test DT.get_convex_hull_indices(tri) == tri.convex_hull.indices
            _tri = triangulate_rectangle(0.0, 2.0, 5.0, 7.3, 5, 15; single_boundary=true)
            points = get_points(_tri)
            ch = convex_hull(points)
            @test ch.indices == _tri.boundary_nodes
            _indices = deepcopy(ch.indices)
            __indices = DT.get_convex_hull_indices(_tri)
            empty!(__indices)
            convex_hull!(_tri; reconstruct=false)
            @test length(__indices) == length(_indices)
            unique!(__indices)
            unique!(ch.indices)
            shift = findfirst(ch.indices .== first(__indices))
            ch.indices .= circshift(ch.indices, 1 - shift)
            @test __indices == ch.indices
            convex_hull!(_tri)
            @test length(__indices) == length(_indices)
            unique!(__indices)
            shift = findfirst(ch.indices .== first(__indices))
            ch.indices .= circshift(ch.indices, 1 - shift)
            @test __indices == ch.indices
            delete_ghost_triangles!(_tri)
            convex_hull!(_tri; reconstruct=false)
            @test length(__indices) == length(_indices)
            unique!(__indices)
            unique!(ch.indices)
            shift = findfirst(ch.indices .== first(__indices))
            ch.indices .= circshift(ch.indices, 1 - shift)
            @test __indices == ch.indices
            convex_hull!(_tri)
            @test length(__indices) == length(_indices)
            unique!(__indices)
            shift = findfirst(ch.indices .== first(__indices))
            ch.indices .= circshift(ch.indices, 1 - shift)
            @test __indices == ch.indices
            @test !DT.has_ghost_triangles(_tri)
      end

      @testset "Boundary nodes" begin
            @test DT.has_multiple_curves(tri)
            @inferred DT.has_multiple_curves(tri)
            @test !DT.has_multiple_curves(tri_2)
            @test DT.has_multiple_segments(tri)
            @test DT.has_multiple_segments(tri_2)
            @test !DT.has_multiple_segments(tri_3)
            @inferred DT.has_multiple_segments(tri_2)
            @test DT.num_curves(tri) == 5
            @inferred DT.num_curves(tri)
            a, b = 0.0, 5.0
            c, d = 3.0, 7.0
            nx = 12
            ny = 15
            @test DT.num_curves(triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 10, 10; add_ghost_triangles=true, single_boundary=false)) == 1
            @test_throws throw_f("The num_segments function has not been defined for the type Vector{Vector{Vector{Int64}}}.") DT.num_segments(tri)
            @test DT.num_segments(tri_2) == 4
            @inferred DT.num_segments(tri_2)
            @test DT.get_boundary_nodes(tri, 1) == tri.boundary_nodes[1]
            @test DT.get_boundary_nodes(tri, 1, 3) == tri.boundary_nodes[1][3]
            @test DT.get_boundary_nodes(tri, (5, 1)) == tri.boundary_nodes[5][1]
            @test DT.get_boundary_nodes(tri_2, 1) == tri_2.boundary_nodes[1]
            @test DT.get_boundary_nodes(tri_2, 3) == tri_2.boundary_nodes[3]
            @test DT.get_boundary_nodes(tri_3) == tri_3.boundary_nodes
            @test DT.get_boundary_nodes(tri_2, 3) == tri_2.boundary_nodes[3]
            @test DT.map_boundary_index(tri, -2) == tri.boundary_map[-2] == (1, 2)
            @test DT.map_boundary_index(tri_2, -3) == tri_2.boundary_map[-3] == 3
            @test DT.map_boundary_index(tri_3, -1) == tri_3.boundary_nodes ==
                  DT.map_boundary_index(tri_3, -1)
            @inferred DT.map_boundary_index(tri_3, -1)
            @test DT.get_curve_index(tri, -1) == 1
            @test DT.get_curve_index(tri, -3) == 1
            @test DT.get_curve_index(tri, -5) == 2
            @test DT.get_curve_index(tri, -11) == 5
            @test DT.get_curve_index(tri_2, -3) == 1
            @test DT.get_curve_index(tri_2, -2) == 1
            @test DT.get_curve_index(tri_3, -1) == 1
            @inferred DT.get_curve_index(tri_3, -1)
            @test DT.num_outer_boundary_segments(tri) == 4
            @test DT.num_outer_boundary_segments(tri_2) == 4
            @test DT.num_outer_boundary_segments(tri_3) == 1
            @test DT.num_outer_boundary_segments(tri_4) == 1
            @inferred DT.num_outer_boundary_segments(tri_4)
      end

      @testset "Triangles" begin
            rng = StableRNG(9882881)
            boundary_nodes, points = convert_boundary_points_to_indices(x, y)
            tri = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
            A = get_total_area(tri)
            refine!(tri; max_area=1e-1A, rng)
            @test DT.triangle_type(tri) == NTuple{3,Int64}
            @inferred DT.triangle_type(tri)
            @test DT.num_triangles(tri) == length(tri.triangles)
            @test DT.each_triangle(tri) == tri.triangles
            @inferred DT.num_triangles(tri)
            @test DT.contains_triangle(tri, (68, 67, -6))[2]
            @inferred DT.contains_triangle(tri, (68, 67, -6))
            @test DT.contains_triangle(tri, (317, 24, 2))[2]
            @test !DT.contains_triangle(tri, (1, 1, 5))[2]
            @test !DT.contains_triangle(tri, 1, 5, 5)[2]
            @test DT.contains_triangle(tri, 199, 185, 170)[2]
            @inferred DT.contains_triangle(tri, 3, 140, 1126)
            @test DT.construct_positively_oriented_triangle(tri, 426, 394, 188) == (426, 394, 188)
            @test DT.construct_positively_oriented_triangle(tri, 188, 394, 426) == (394, 188, 426)
            @inferred DT.construct_positively_oriented_triangle(tri, 188, 394, 426)
            _solid_itr = each_solid_triangle(tri)
            @test DelaunayTriangulation.each_triangle(_solid_itr) == _solid_itr
            @test DT.initialise_triangles(typeof(_solid_itr)) == Set{NTuple{3,Int64}}() && eltype(DT.initialise_triangles(typeof(_solid_itr))) == NTuple{3,Int64}
            @test Base.IteratorSize(_solid_itr) == Base.HasLength()
            @test Base.IteratorEltype(_solid_itr) == Base.HasEltype()
            @test Base.eltype(_solid_itr) == NTuple{3,Int64}
            @test each_solid_triangle(tri) isa DT.EachSolidTriangle
            _solid_tri = collect(_solid_itr)
            @inferred collect(_solid_itr)
            @test all(!DT.is_ghost_triangle, _solid_tri)
            _ghost_itr = each_ghost_triangle(tri)
            @test DT.initialise_triangles(typeof(_ghost_itr)) == Set{NTuple{3,Int64}}() && eltype(DT.initialise_triangles(typeof(_ghost_itr))) == NTuple{3,Int64}
            @test Base.IteratorSize(_ghost_itr) == Base.HasLength()
            @test Base.IteratorEltype(_ghost_itr) == Base.HasEltype()
            @test Base.eltype(_ghost_itr) == NTuple{3,Int64}
            @test each_ghost_triangle(tri) isa DT.EachGhostTriangle
            _ghost_tri = collect(_ghost_itr)
            @inferred collect(_ghost_itr)
            @test DelaunayTriangulation.each_triangle(_ghost_itr) == _ghost_itr
            @test all(DT.is_ghost_triangle, _ghost_tri)
            @test length(_ghost_tri) + length(_solid_tri) == num_triangles(tri)
            @test length(_solid_tri) == DT.num_solid_triangles(tri) == length(_solid_itr)
            @test length(_ghost_tri) == DT.num_ghost_triangles(tri) == length(_ghost_itr)
            boundary_nodes, points = convert_boundary_points_to_indices(x, y)
            rng = StableRNG(9882881)
            ___tri = triangulate(points; boundary_nodes, delete_ghosts=false, rng)
            A = get_total_area(___tri)
            refine!(___tri; max_area=1e-1A, rng)
            DT.delete_ghost_triangles!(___tri)
            @test collect(each_triangle(___tri)) == collect(each_solid_triangle(___tri))
            @test length(collect(each_ghost_triangle(___tri))) == 0
            @test sort(collect(filter(!DT.is_ghost_triangle, each_triangle(___tri)))) == sort(collect(each_solid_triangle(___tri)))
            @test sort(collect(filter(DT.is_ghost_triangle, each_triangle(___tri)))) == sort(collect(each_ghost_triangle(___tri)))
            @test DelaunayTriangulation.sort_triangles(each_solid_triangle(tri)) == DelaunayTriangulation.sort_triangles(get_triangles(___tri))
            @test DelaunayTriangulation.triangle_type(_ghost_itr) == DelaunayTriangulation.triangle_type(_ghost_itr)
            @test DelaunayTriangulation.triangle_type(_solid_itr) == DelaunayTriangulation.triangle_type(_solid_itr)
      end

      @testset "Edges" begin
            rng = StableRNG(998871)
            boundary_nodes, points = convert_boundary_points_to_indices(x, y)
            tri = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
            A = get_total_area(tri)
            refine!(tri; max_area=1e-1A, rng)
            @test DT.edge_type(tri) == NTuple{2,Int64}
            @inferred DT.edge_type(tri)
            @test DT.num_edges(tri) == length(tri.graph.graph.E)
            @inferred DT.num_edges(tri)
            @test DT.each_edge(tri) == tri.graph.graph.E
            @inferred DT.each_edge(tri)
            _solid_itr = each_solid_edge(tri)
            @test DelaunayTriangulation.each_edge(_solid_itr) == _solid_itr
            @test DT.initialise_edges(typeof(_solid_itr)) == Set{NTuple{2,Int64}}() && eltype(DT.initialise_edges(typeof(_solid_itr))) == NTuple{2,Int64}
            @test Base.IteratorSize(_solid_itr) == Base.HasLength()
            @test Base.IteratorEltype(_solid_itr) == Base.HasEltype()
            @test Base.eltype(_solid_itr) == NTuple{2,Int64}
            @test each_solid_edge(tri) isa DT.EachSolidEdge
            _solid_tri = collect(_solid_itr)
            @inferred collect(_solid_itr)
            @test all(!DT.is_ghost_edge, _solid_tri)
            _ghost_itr = each_ghost_edge(tri)
            @test DT.initialise_edges(typeof(_ghost_itr)) == Set{NTuple{2,Int64}}() && eltype(DT.initialise_edges(typeof(_ghost_itr))) == NTuple{2,Int64}
            @test Base.IteratorSize(_ghost_itr) == Base.HasLength()
            @test Base.IteratorEltype(_ghost_itr) == Base.HasEltype()
            @test Base.eltype(_ghost_itr) == NTuple{2,Int64}
            @test each_ghost_edge(tri) isa DT.EachGhostEdge
            _ghost_tri = collect(_ghost_itr)
            @inferred collect(_ghost_itr)
            @test DelaunayTriangulation.each_edge(_ghost_itr) == _ghost_itr
            @test all(DT.is_ghost_edge, _ghost_tri)
            @test length(_ghost_tri) + length(_solid_tri) == num_edges(tri)
            rng = StableRNG(998871)
            boundary_nodes, points = convert_boundary_points_to_indices(x, y)
            ___tri = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
            A = get_total_area(___tri)
            refine!(___tri; max_area=1e-1A, rng)
            DT.delete_ghost_triangles!(___tri)
            @test sort(collect(filter(!DT.is_ghost_edge, each_edge(___tri)))) == sort(collect(each_solid_edge(___tri)))
            @test sort(collect(filter(DT.is_ghost_edge, each_edge(___tri)))) == sort(collect(each_ghost_edge(___tri)))
            @test length(collect(each_ghost_edge(___tri))) == num_edges(___tri) .- length(sort(collect(filter(!DT.is_ghost_edge, each_edge(___tri)))))
            @test length(collect(each_ghost_edge(___tri))) == DT.num_ghost_edges(___tri) == length(_ghost_itr)
            @test length(collect(each_solid_edge(___tri))) == DT.num_solid_edges(___tri) == length(_solid_itr)
            @test DelaunayTriangulation.edge_type(_ghost_itr) == DelaunayTriangulation.edge_type(_ghost_itr)
            @test DelaunayTriangulation.edge_type(_solid_itr) == DelaunayTriangulation.edge_type(_solid_itr)
      end

      @testset "Points" begin
            @test DT.get_point(tri, 2) == Tuple(tri.points[2])
            @inferred DT.get_point(tri, 2)
            @test DT.get_point(tri, 17) == Tuple(tri.points[17])
            @test DT.each_point_index(tri) == 1:length(tri.points)
            @inferred DT.each_point_index(tri)
            @test DT.each_point(tri) == tri.points
            @inferred DT.each_point(tri)
            @test DT.num_points(tri) == length(tri.points)
            @inferred DT.num_points(tri)
            @test DT.get_point(tri, 2, 5, 6, 10) ==
                  ntuple(i -> Tuple(tri.points[(2, 5, 6, 10)[i]]), 4)
            @inferred DT.get_point(tri, 2, 5, 6, 10)
            rep = DT.get_representative_point_list(tri)
            rep[1] = DT.RepresentativeCoordinates(0.5, 0.3, 2)
            rep[2] = DT.RepresentativeCoordinates(2.5, 7.3, 7)
            rep[3] = DT.RepresentativeCoordinates(6.5, -0.6, 13)
            rep[4] = DT.RepresentativeCoordinates(6.5, -0.66234, 13)
            rep[5] = DT.RepresentativeCoordinates(6.534234, -0.6, 13)
            @test DT.get_point(tri, -1) == DT.getxy(rep[1])
            @test DT.get_point(tri, -2) == DT.getxy(rep[1])
            @test DT.get_point(tri, -3) == DT.getxy(rep[1])
            @test DT.get_point(tri, -4) == DT.getxy(rep[1])
            @test DT.get_point(tri, -5) == DT.getxy(rep[2])
            @test DT.get_point(tri, -6) == DT.getxy(rep[3])
            @test DT.get_point(tri, -7) == DT.getxy(rep[4])
            @test DT.get_point(tri, -8) == DT.getxy(rep[4])
            @test DT.get_point(tri, -9) == DT.getxy(rep[4])
            @test DT.get_point(tri, -10) == DT.getxy(rep[4])
            @test DT.get_point(tri, -11) == DT.getxy(rep[5])
            rep = DT.get_representative_point_list(tri_2)
            @test DT.get_point(tri_2, -1) == DT.getxy(rep[1])
            @inferred DT.get_point(tri_2, -1)
            @test DT.get_point(tri_2, -2) == DT.getxy(rep[1])
            @test DT.get_point(tri_2, -3) == DT.getxy(rep[1])
            @test DT.get_point(tri_2, -4) == DT.getxy(rep[1])
            @test DT.get_point(tri_3, -1) == DT.getxy(DT.get_representative_point_coordinates(tri_3, 1))
            @test_throws KeyError DT.get_point(tri, -12)
            @test_throws KeyError DT.get_point(tri_2, -5)
            @test_throws KeyError DT.get_point(tri_3, -2)
            @test get_point(tri, tri.points[2]) == tri.points[2]
            @inferred get_point(tri, tri.points[2])
            @test collect(DT.all_boundary_indices(tri)) == [DT.BoundaryIndex,
                  DT.BoundaryIndex - 1,
                  DT.BoundaryIndex - 2,
                  DT.BoundaryIndex - 3,
                  DT.BoundaryIndex - 4,
                  DT.BoundaryIndex - 5,
                  DT.BoundaryIndex - 6,
                  DT.BoundaryIndex - 7,
                  DT.BoundaryIndex - 8,
                  DT.BoundaryIndex - 9,
                  DT.BoundaryIndex - 10]
            rng = StableRNG(887271)
            boundary_nodes, points = convert_boundary_points_to_indices(x, y)
            _triq = triangulate(points; boundary_nodes, rng)
            A = get_total_area(_triq)
            refine!(_triq; max_area=1e-1A, rng)
            _solid_itr = each_solid_vertex(_triq)
            @test DelaunayTriangulation.each_vertex(_solid_itr) == _solid_itr
            @test Base.IteratorSize(_solid_itr) == Base.HasLength()
            @test Base.IteratorEltype(_solid_itr) == Base.HasEltype()
            @test Base.eltype(_solid_itr) == Int64
            @test each_solid_vertex(_triq) isa DT.EachSolidVertex
            _solid_tri = collect(_solid_itr)
            @inferred collect(_solid_itr)
            @test sort(_solid_tri) == sort(collect(1:num_points(_triq)))
            @test all(!DT.is_boundary_index, _solid_tri)
            _ghost_itr = each_ghost_vertex(_triq)
            @test Base.IteratorSize(_ghost_itr) == Base.HasLength()
            @test Base.IteratorEltype(_ghost_itr) == Base.HasEltype()
            @test Base.eltype(_ghost_itr) == Int64
            @test each_ghost_vertex(_triq) isa DT.EachGhostVertex
            _ghost_tri = collect(_ghost_itr)
            @inferred collect(_ghost_itr)
            @test DelaunayTriangulation.each_vertex(_ghost_itr) == _ghost_itr
            @test all(DT.is_boundary_index, _ghost_tri)
            @test _ghost_tri == [DT.BoundaryIndex,
                  DT.BoundaryIndex - 1,
                  DT.BoundaryIndex - 2,
                  DT.BoundaryIndex - 3,
                  DT.BoundaryIndex - 4,
                  DT.BoundaryIndex - 5,
                  DT.BoundaryIndex - 6,
                  DT.BoundaryIndex - 7,
                  DT.BoundaryIndex - 8,
                  DT.BoundaryIndex - 9,
                  DT.BoundaryIndex - 10]
            @test length(_ghost_tri) == length(_ghost_itr) == sum(<(0), keys(_triq.adjacent2vertex.adjacent2vertex))
            @test length(_solid_tri) == length(_solid_itr) == num_points(_triq)
            rng = StableRNG(887271)
            ___tri = triangulate(points; boundary_nodes, rng)
            A = get_total_area(___tri)
            refine!(___tri; max_area=1e-1A, rng)
            DT.delete_ghost_triangles!(___tri)
            @test sort(collect(filter(!DT.is_boundary_index, each_vertex(___tri)))) == sort(collect(each_solid_vertex(___tri)))
            @test sort(collect(filter(DT.is_boundary_index, each_vertex(___tri)))) == sort(collect(each_ghost_vertex(___tri)))
            @test length(collect(each_ghost_vertex(___tri))) == num_vertices(___tri) .- length(sort(collect(filter(!DT.is_boundary_index, each_vertex(___tri)))))
            tri1 = Triangulation([[1.0, 2.0], [3.0, 4.0]])
            DT.push_point!(tri1, 13.7, 5.0)
            @test get_points(tri1) == [[1.0, 2.0], [3.0, 4.0], [13.7, 5.0]]
            @test get_point(tri1, 3) == (13.7, 5.0)
            tri1 = Triangulation([(1.0, 2.0), (3.0, 4.0)])
            DT.push_point!(tri1, 13.7, 5.0)
            @test get_points(tri1) == [(1.0, 2.0), (3.0, 4.0), (13.7, 5.0)]
            @test get_point(tri1, 3) == (13.7, 5.0)
            DT.push_point!(tri1, (19.0, 17.05))
            @test get_points(tri1) == [(1.0, 2.0), (3.0, 4.0), (13.7, 5.0), (19.0, 17.05)]
            DT.pop_point!(tri1)
            @test get_points(tri1) == [(1.0, 2.0), (3.0, 4.0), (13.7, 5.0)]
            DT.pop_point!(tri1)
            @test get_points(tri1) == [(1.0, 2.0), (3.0, 4.0)]
            DT.set_point!(tri1, 2, 13.7, 5.0)
            @test get_points(tri1) == [(1.0, 2.0), (13.7, 5.0)]
            DT.set_point!(tri1, 2, (19.0, 17.05))
            @test get_points(tri1) == [(1.0, 2.0), (19.0, 17.05)]
      end

      @testset "Miscellaneous" begin
            @test DT.integer_type(tri) == Int64
            @test DT.number_type(tri) == Float64
            @inferred DT.integer_type(tri)
            @inferred DT.number_type(tri)
            DT.clear_empty_features!(tri)
            clean_tri = deepcopy(tri)
            get_adjacent(tri, 17, 32)
            get_adjacent(tri, 58, 37)
            DT.add_adjacent2vertex!(tri, 3, 58, 60)
            DT.add_neighbour!(tri, 5, 171)
            @test clean_tri ≠ tri
            DT.delete_adjacent2vertex!(tri, 3, 58, 60)
            DT.delete_neighbour!(tri, 5, 171)
            @test clean_tri ≠ tri
            DT.clear_empty_features!(tri)
            @test clean_tri == tri
            _tri = triangulate_rectangle(0.0, 10.0, 0.0, 20.0, 11, 21)
            T = (51, 41, 52)
            ℓ = (7.0, 3.5)
            @test DT.find_edge(_tri, T, ℓ) == (41, 52)
            ℓ = (6.5, 4.0)
            @test DT.find_edge(_tri, T, ℓ) == (52, 51)
            ℓ = (6.5, 3.5)
            @test DT.find_edge(_tri, T, ℓ) == (51, 41)
            @inferred DT.find_edge(_tri, T, ℓ)
            push!(get_all_constrained_edges(_tri), (1, 2))
            @test DT.is_constrained(_tri)
            empty!(get_all_constrained_edges(_tri))
            @test !DT.is_constrained(_tri)
            @test collect(DT.all_boundary_indices(_tri)) ==
                  [DT.BoundaryIndex, DT.BoundaryIndex - 1, DT.BoundaryIndex - 2, DT.BoundaryIndex - 3]
      end
end

@testset "merge_constrained_edges" begin
      all_bn = get_boundary_nodes(tri)
      i = rand(1:100000, 50)
      j = rand(1:100000, 50)
      all_ce = Set(((i, j) for (i, j) in zip(i, j)))
      bn_map = get_boundary_map(tri)
      bn1 = all_bn[1]
      bn11 = bn1[1]
      bn12 = bn1[2]
      bn13 = bn1[3]
      bn14 = bn1[4]
      bn2 = all_bn[2][1]
      bn3 = all_bn[3][1]
      bn4 = all_bn[4]
      bn41 = bn4[1]
      bn42 = bn4[2]
      bn43 = bn4[3]
      bn44 = bn4[4]
      bn5 = all_bn[5][1]
      e11 = Set(((bn11[i], bn11[i+1]) for i in 1:(length(bn11)-1)))
      e12 = Set(((bn12[i], bn12[i+1]) for i in 1:(length(bn12)-1)))
      e13 = Set(((bn13[i], bn13[i+1]) for i in 1:(length(bn13)-1)))
      e14 = Set(((bn14[i], bn14[i+1]) for i in 1:(length(bn14)-1)))
      e2 = Set(((bn2[i], bn2[i+1]) for i in 1:(length(bn2)-1)))
      e3 = Set(((bn3[i], bn3[i+1]) for i in 1:(length(bn3)-1)))
      e41 = Set(((bn41[i], bn41[i+1]) for i in 1:(length(bn41)-1)))
      e42 = Set(((bn42[i], bn42[i+1]) for i in 1:(length(bn42)-1)))
      e43 = Set(((bn43[i], bn43[i+1]) for i in 1:(length(bn43)-1)))
      e44 = Set(((bn44[i], bn44[i+1]) for i in 1:(length(bn44)-1)))
      e5 = Set(((bn5[i], bn5[i+1]) for i in 1:(length(bn5)-1)))
      ace = Set{NTuple{2,Int64}}()
      for es in (e11, e12, e13, e14, e2, e3, e41, e42, e43, e44, e5)
            for e in es
                  push!(ace, e)
            end
      end
      for e in all_ce
            push!(ace, e)
      end
      @test ace == DT.merge_constrained_edges(bn_map, all_bn, all_ce)

      all_bn = get_boundary_nodes(tri_2)
      i = rand(1:100000, 50)
      j = rand(1:100000, 50)
      all_ce = Set(((i, j) for (i, j) in zip(i, j)))
      bn_map = get_boundary_map(tri_2)
      bn1 = all_bn[1]
      bn2 = all_bn[2]
      bn3 = all_bn[3]
      bn4 = all_bn[4]
      e1 = Set(((bn1[i], bn1[i+1]) for i in 1:(length(bn1)-1)))
      e2 = Set(((bn2[i], bn2[i+1]) for i in 1:(length(bn2)-1)))
      e3 = Set(((bn3[i], bn3[i+1]) for i in 1:(length(bn3)-1)))
      e4 = Set(((bn4[i], bn4[i+1]) for i in 1:(length(bn4)-1)))
      ace = Set{NTuple{2,Int64}}()
      for es in (e1, e2, e3, e4)
            for e in es
                  push!(ace, e)
            end
      end
      for e in all_ce
            push!(ace, e)
      end
      @test ace == DT.merge_constrained_edges(bn_map, all_bn, all_ce)

      all_bn = get_boundary_nodes(tri_3)
      i = rand(1:100000, 50)
      j = rand(1:100000, 50)
      all_ce = Set(((i, j) for (i, j) in zip(i, j)))
      bn_map = get_boundary_map(tri_3)
      e = Set(((all_bn[i], all_bn[i+1]) for i in 1:(length(all_bn)-1)))
      ace = Set{NTuple{2,Int64}}()
      for e in e
            push!(ace, e)
      end
      for e in all_ce
            push!(ace, e)
      end
      @test ace == DT.merge_constrained_edges(bn_map, all_bn, all_ce)

      all_bn = get_boundary_nodes(tri_4)[1]
      i = rand(1:100000, 50)
      j = rand(1:100000, 50)
      all_ce = Set(((i, j) for (i, j) in zip(i, j)))
      bn_map = get_boundary_map(tri_4)
      e = Set(((all_bn[i], all_bn[i+1]) for i in 1:(length(all_bn)-1)))
      ace = Set{NTuple{2,Int64}}()
      for e in e
            push!(ace, e)
      end
      for e in all_ce
            push!(ace, e)
      end
      @test ace == DT.merge_constrained_edges(bn_map, [all_bn], all_ce)
end

@testset "sort_edge_by_degree" begin
      tri = triangulate(rand(2, 500); delete_ghosts=false)
      for e in each_edge(tri)
            new_e = DT.sort_edge_by_degree(tri, e)
            d1 = DT.num_neighbours(tri, e[1])
            d2 = DT.num_neighbours(tri, e[2])
            if d1 ≤ d2
                  @test new_e == e
            else
                  @test new_e == (e[2], e[1])
            end
      end
end

@testset "triangle_line_segment_intersection" begin
      n = 60
      for _ in 1:10
            n += rand(1:125)
            tri = triangulate(12randn(2, n))
            for qi in each_solid_vertex(tri)
                  for k in each_solid_vertex(tri)
                        q = get_point(tri, qi)
                        history = DT.PointLocationHistory{NTuple{3,Int64},NTuple{2,Int64},Int64}()
                        jump_and_march(tri, q;
                              k,
                              store_history=true,
                              history)
                        visited_triangles = history.triangles
                        collinear_segments = history.collinear_segments
                        @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T)), visited_triangles)
                        @test all(!DT.is_none, [DT.triangle_line_segment_intersection(tri, T, (qi, k)) for T in visited_triangles])
                        @test allunique(visited_triangles)
                        if !isempty(collinear_segments)
                              @test all(E -> DT.is_collinear(DT.point_position_relative_to_line(tri, qi, k, E[1])), collinear_segments)
                              @test all(E -> DT.is_collinear(DT.point_position_relative_to_line(tri, qi, k, E[2])), collinear_segments)
                        end
                  end
            end
      end
end

@testset "point_closest_to_line" begin
      tri = fixed_shewchuk_example_constrained()
      i, j = 2, 7
      u, v = 9, 8
      @test DT.is_closer(DT.point_closest_to_line(tri, i, j, u, v))
      @test DT.is_further(DT.point_closest_to_line(tri, i, j, v, u))
      u, v = 3, 11
      @test DT.is_closer(DT.point_closest_to_line(tri, i, j, u, v))
      @test DT.is_further(DT.point_closest_to_line(tri, i, j, v, u))
      i, j = 7, 2
      u, v = 6, 1
      @test DT.is_closer(DT.point_closest_to_line(tri, i, j, u, v))
      @test DT.is_further(DT.point_closest_to_line(tri, i, j, v, u))
      u, v = 2, 6
      @test DT.is_closer(DT.point_closest_to_line(tri, i, j, u, v))
      @test DT.is_further(DT.point_closest_to_line(tri, i, j, v, u))
      u, v = 7, 4
      @test DT.is_closer(DT.point_closest_to_line(tri, i, j, u, v))
      @test DT.is_further(DT.point_closest_to_line(tri, i, j, v, u))
end

rng = StableRNG(91928281)
pts = [(rand(rng), rand(rng)) for _ in 1:50]
bnd_pts = [(0.3cos(θ), 0.3sin(θ)) .+ 0.5 for θ in LinRange(0, 2π - 1 / 250, 25)]
bnd_id = [(51:75)..., 51]
append!(pts, bnd_pts)
global tric = triangulate(pts; boundary_nodes=bnd_id, rng)

@testset "each_constrained_edge" begin
      @test each_constrained_edge(tric) == each_edge(get_all_constrained_edges(tric))
end

@testset "contains_constrained_edge" begin
      @test !DT.contains_constrained_edge(tric, 12, 17)
      @test DT.contains_constrained_edge(tric, 69, 70)
      @test DT.contains_constrained_edge(tric, 70, 69)
      @test !DT.contains_constrained_edge(tric, 32, 41)
      @test !DT.contains_constrained_edge(tric, 45, 38)
      @test DT.contains_constrained_edge(tric, 63, 64)
      @test !DT.contains_constrained_edge(tric, (45, 38))
      @test !DT.contains_constrained_edge(tric, 26, 22)
      @test DT.contains_constrained_edge(tric, 64, 65)
      @test DT.contains_constrained_edge(tric, 55, 54)
      @test DT.contains_constrained_edge(tric, 58, 57)
      @test DT.contains_constrained_edge(tric, 59, 60)
      @test !DT.contains_constrained_edge(tric, 30, 70)
      @test !DT.contains_constrained_edge(tric, 56, 37)
      @test DT.contains_constrained_edge(tric, 73, 74)
end

@testset "get_all_boundary_nodes" begin
      x, y = complicated_geometry()
      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x, y)
      tri = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri)
      refine!(tri, rng=rng, max_area=1e-1A)
      all_bn = DT.get_all_boundary_nodes(tri)
      @test all_bn == Set(reduce(vcat, reduce(vcat, get_boundary_nodes(tri))))
      tri2, label_map, index_map = simple_geometry()
      all_bn = DT.get_all_boundary_nodes(tri2)
      @test all_bn == Set(reduce(vcat, reduce(vcat, get_boundary_nodes(tri2))))
      tri3 = triangulate_rectangle(0, 1, 0, 1, 50, 50; add_ghost_triangles=true, single_boundary=false)
      all_bn = DT.get_all_boundary_nodes(tri3)
      @test all_bn == Set(reduce(vcat, reduce(vcat, get_boundary_nodes(tri3))))
      tri4 = triangulate_rectangle(0, 1, 0, 1, 50, 50; add_ghost_triangles=true, single_boundary=true)
      all_bn = DT.get_all_boundary_nodes(tri4)
      @test all_bn == Set(reduce(vcat, reduce(vcat, get_boundary_nodes(tri4))))
      tri = triangulate(rand(2, 50))
      @test isempty(DT.get_all_boundary_nodes(tri))
end

@testset "get_boundary_edge_map" begin
      x, y = complicated_geometry()
      x, y = complicated_geometry()
      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x, y)
      tri = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri)
      refine!(tri, rng=rng, max_area=1e-1A)
      @test DT.get_boundary_edge_map(tri, 111, 221) == ((4, 1), 1)
      @test DT.get_boundary_edge_map(tri, 49, 50) == ((2, 1), 37)
      @test DT.get_boundary_edge_map(tri, (1, 200)) == ((1, 1), 1)
      @test DT.get_boundary_edge_map(tri, (7, 8)) == ((1, 3), 1)
end

@testset "split_boundary_edge!" begin
      x, y = complicated_geometry()
      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x, y)
      tri_1 = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri_1)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
      tri_2 = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri_2)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
      tri_3 = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri_3)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(reverse(reverse.(x[2])), reverse(reverse.(y[2])))
      tri_4 = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri_4)

      DT.split_boundary_edge!(tri_1, (21, 22), 500)
      @test tri_1.boundary_nodes[2][1][1:11] == [13, 14, 15, 16, 17, 18, 19, 20, 21, 500, 22]
      @test DT.get_boundary_edge_map(tri_1, 21, 500) == ((2, 1), 9)
      @test DT.get_boundary_edge_map(tri_1, 500, 22) == ((2, 1), 10)
      @test_throws KeyError DT.get_boundary_edge_map(tri_1, 21, 22)
      DT.split_boundary_edge!(tri_1, (7, 8), 5000)
      @test tri_1.boundary_nodes[1][3] == [7, 5000, 8, 9, 10]
      @test DT.get_boundary_edge_map(tri_1, 7, 5000) == ((1, 3), 1)
      @test DT.get_boundary_edge_map(tri_1, 5000, 8) == ((1, 3), 2)
      @test_throws KeyError DT.get_boundary_edge_map(tri_1, 7, 8)

      DT.split_boundary_edge!(tri_2, 8, 9, 300)
      @test tri_2.boundary_nodes[3] == [7, 8, 300, 9, 10]
      @test DT.get_boundary_edge_map(tri_2, 8, 300) == (3, 2)
      @test DT.get_boundary_edge_map(tri_2, 300, 9) == (3, 3)
      @test_throws KeyError DT.get_boundary_edge_map(tri_2, 8, 9)

      DT.split_boundary_edge!(tri_3, 3, 4, 5000)
      @test tri_3.boundary_nodes == [1, 2, 3, 5000, 4, 1]
      @test DT.get_boundary_edge_map(tri_3, 3, 5000) == (tri_3.boundary_nodes, 3)
      @test DT.get_boundary_edge_map(tri_3, 5000, 4) == (tri_3.boundary_nodes, 4)
      @test_throws KeyError DT.get_boundary_edge_map(tri_3, 3, 4)

      DT.split_boundary_edge!(tri_4, 6, 7, 1200)
      @test DT.get_boundary_edge_map(tri_4, 6, 1200) == (1, 6)
      @test DT.get_boundary_edge_map(tri_4, 1200, 7) == (1, 7)
      @test_throws KeyError DT.get_boundary_edge_map(tri_4, 6, 7)
end

@testset "merge_boundary_edge!" begin
      x, y = complicated_geometry()
      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x, y)
      tri_1 = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri_1)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
      tri_2 = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri_2)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
      tri_3 = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri_3)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(reverse(reverse.(x[2])), reverse(reverse.(y[2])))
      tri_4 = triangulate(points; rng, boundary_nodes)
      A = get_total_area(tri_4)

      orig_bn = deepcopy(get_boundary_nodes(tri_1))
      orig_bnn = deepcopy(get_boundary_edge_map(tri_1))
      DT.split_boundary_edge!(tri_1, (21, 22), 170)
      @test get_boundary_nodes(tri_1) ≠ orig_bn
      @test get_boundary_edge_map(tri_1) ≠ orig_bnn
      DT.merge_boundary_edge!(tri_1, (21, 22), 170)
      @test get_boundary_nodes(tri_1) == orig_bn
      @test get_boundary_edge_map(tri_1) == orig_bnn
      @test_throws KeyError DT.get_boundary_edge_map(tri_1, 21, 170)
      @test_throws KeyError DT.get_boundary_edge_map(tri_1, 170, 22)

      orig_bn = deepcopy(get_boundary_nodes(tri_1))
      orig_bnn = deepcopy(get_boundary_edge_map(tri_1))
      DT.split_boundary_edge!(tri_1, (7, 8), 5000)
      @test get_boundary_nodes(tri_1) ≠ orig_bn
      @test get_boundary_edge_map(tri_1) ≠ orig_bnn
      DT.merge_boundary_edge!(tri_1, (7, 8), 5000)
      @test get_boundary_nodes(tri_1) == orig_bn
      @test get_boundary_edge_map(tri_1) == orig_bnn
      @test_throws KeyError DT.get_boundary_edge_map(tri_1, 7, 5000)
      @test_throws KeyError DT.get_boundary_edge_map(tri_1, 5000, 8)

      orig_bn = deepcopy(get_boundary_nodes(tri_2))
      orig_bnn = deepcopy(get_boundary_edge_map(tri_2))
      DT.split_boundary_edge!(tri_2, 8, 9, 8182)
      @test get_boundary_nodes(tri_2) ≠ orig_bn
      @test get_boundary_edge_map(tri_2) ≠ orig_bnn
      DT.merge_boundary_edge!(tri_2, 8, 9, 8182)
      @test get_boundary_nodes(tri_2) == orig_bn
      @test get_boundary_edge_map(tri_2) == orig_bnn
      @test_throws KeyError DT.get_boundary_edge_map(tri_2, 8, 8182)
      @test_throws KeyError DT.get_boundary_edge_map(tri_2, 8182, 9)

      orig_bn = deepcopy(get_boundary_nodes(tri_3))
      orig_bnn = deepcopy(get_boundary_edge_map(tri_3))
      DT.split_boundary_edge!(tri_3, 3, 4, 18289)
      @test get_boundary_nodes(tri_3) ≠ orig_bn
      @test get_boundary_edge_map(tri_3) ≠ orig_bnn
      DT.merge_boundary_edge!(tri_3, 3, 4, 18289)
      @test get_boundary_nodes(tri_3) == orig_bn
      @test get_boundary_edge_map(tri_3) == orig_bnn
      @test_throws KeyError DT.get_boundary_edge_map(tri_3, 3, 18289)
      @test_throws KeyError DT.get_boundary_edge_map(tri_3, 18289, 4)

      orig_bn = deepcopy(get_boundary_nodes(tri_4))
      orig_bnn = deepcopy(get_boundary_edge_map(tri_4))
      DT.split_boundary_edge!(tri_4, 6, 7, 1200)
      @test get_boundary_nodes(tri_4) ≠ orig_bn
      @test get_boundary_edge_map(tri_4) ≠ orig_bnn
      DT.merge_boundary_edge!(tri_4, 6, 7, 1200)
      @test get_boundary_nodes(tri_4) == orig_bn
      @test get_boundary_edge_map(tri_4) == orig_bnn
      @test_throws KeyError DT.get_boundary_edge_map(tri_4, 6, 1200)
      @test_throws KeyError DT.get_boundary_edge_map(tri_4, 1200, 7)
end

@testset "Random sampling of vertices" begin
      tri = triangulate(rand(2, 500); delete_ghosts=false)
      solid = Int64[]
      ghost = Int64[]
      solid_itr = each_solid_vertex(tri)
      ghost_itr = each_ghost_vertex(tri)
      for _ in 1:100_000
            push!(solid, rand(solid_itr))
            push!(ghost, rand(ghost_itr))
      end
      @test all(!DT.is_boundary_index, solid)
      @test all(DT.is_boundary_index, ghost)
end
