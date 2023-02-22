using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StructEquality
using Setfield

@struct_equal Triangulation
@struct_equal DT.Adjacent
@struct_equal DT.Adjacent2Vertex
@struct_equal DT.Graph
@struct_equal ConvexHull
include("../helper_functions.jl")

## Initialising 
pts = rand(2, 500)
tri = Triangulation(pts; IntegerType=Int32)
@test tri == Triangulation(pts,
                           Set{NTuple{3,Int32}}(),
                           DT.Adjacent{Int32,NTuple{2,Int32}}(),
                           DT.Adjacent2Vertex{Int32,Set{NTuple{2,Int32}},NTuple{2,Int32}}(),
                           DT.Graph{Int32}(),
                           Int32[],
                           OrderedDict{Int32,Vector{Int32}}(Int32(DT.BoundaryIndex) => Int32[]),
                           OrderedDict{Int32,UnitRange{Int32}}(-1 => -1:-1),
                           Set{NTuple{2,Int32}}(),
                           ConvexHull(pts, Int32[]))

## Gmsh setup
include("../helper_functions.jl")
x, y = complicated_geometry()
tri = generate_mesh(x, y, 0.1; convert_result=true, add_ghost_triangles=true)
tri_2 = generate_mesh(x[1], y[1], 0.1; convert_result=true, add_ghost_triangles=true)
tri_3 = generate_mesh([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0], 0.1;
                      convert_result=true, add_ghost_triangles=true)
tri_4 = generate_mesh(x[2], y[2], 0.1; convert_result=true, add_ghost_triangles=true)

## Getters 
@test DT.get_points(tri) == tri.points
@test DT.get_triangles(tri) == tri.triangles
@test DT.get_adjacent(tri) == tri.adjacent
@test DT.get_adjacent2vertex(tri) == tri.adjacent2vertex
@test DT.get_graph(tri) == tri.graph
@test DT.get_boundary_nodes(tri) == tri.boundary_nodes
@test DT.get_boundary_map(tri) == tri.boundary_map
@test DT.get_constrained_edges(tri) == tri.constrained_edges
@test DT.get_convex_hull(tri) == tri.convex_hull
@inferred DT.get_points(tri)
@inferred DT.get_triangles(tri)
@inferred DT.get_adjacent(tri)
@inferred DT.get_adjacent2vertex(tri)
@inferred DT.get_graph(tri)
@inferred DT.get_boundary_nodes(tri)
@inferred DT.get_boundary_map(tri)
@inferred DT.get_constrained_edges(tri)
@inferred DT.get_convex_hull(tri)

## Extensions 
# Adjacent
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
@test DT.get_adjacent(tri, 129, DT.BoundaryIndex) == 130
@test DT.get_adjacent(tri, 129, DT.BoundaryIndex; check_existence=Val(true)) == 130
@test DT.get_adjacent(tri, 304, DT.BoundaryIndex - 6; check_existence=Val(true)) == 114
@test DT.get_adjacent(tri, 122, DT.BoundaryIndex - 8; check_existence=Val(true)) == 337
@test DT.get_adjacent(tri, 122, DT.BoundaryIndex - 7) == 337

# Adjacent2Vertex 
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

# Graph 
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

# Convex Hull 
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
circshift!(ch.indices, 1 - shift)
@test __indices == ch.indices
convex_hull!(_tri)
@test length(__indices) == length(_indices)
unique!(__indices)
shift = findfirst(ch.indices .== first(__indices))
circshift!(ch.indices, 1 - shift)
@test __indices == ch.indices
delete_ghost_triangles!(_tri)
convex_hull!(_tri; reconstruct=false)
@test length(__indices) == length(_indices)
unique!(__indices)
unique!(ch.indices)
shift = findfirst(ch.indices .== first(__indices))
circshift!(ch.indices, 1 - shift)
@test __indices == ch.indices
convex_hull!(_tri)
@test length(__indices) == length(_indices)
unique!(__indices)
shift = findfirst(ch.indices .== first(__indices))
circshift!(ch.indices, 1 - shift)
@test __indices == ch.indices
@test !DT.has_ghost_triangles(_tri)

# Boundary Nodes 
@test DT.has_multiple_curves(tri)
@inferred DT.has_multiple_curves(tri)
@test !DT.has_multiple_curves(tri_2)
@test DT.has_multiple_segments(tri)
@test DT.has_multiple_segments(tri_2)
@test !DT.has_multiple_segments(tri_3)
@inferred DT.has_multiple_segments(tri_2)
@test DT.num_curves(tri) == 5
@inferred DT.num_curves(tri)
@test_throws "The" DT.num_segments(tri)
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

# Triangles 
@test DT.triangle_type(tri) == NTuple{3,Int64}
@inferred DT.triangle_type(tri)
@test DT.num_triangles(tri) == length(tri.triangles)
@test DT.each_triangle(tri) == tri.triangles
@inferred DT.num_triangles(tri)
@test DT.contains_triangle(tri, (68, 67, -6))[2]
@inferred DT.contains_triangle(tri, (68, 67, -6))
@test DT.contains_triangle(tri, (900, 832, 1602))[2]
@test !DT.contains_triangle(tri, (1, 1, 5))[2]
@test !DT.contains_triangle(tri, 1, 5, 5)[2]
@test DT.contains_triangle(tri, 3, 140, 1126)[2]
@inferred DT.contains_triangle(tri, 3, 140, 1126)
@test DT.construct_positively_oriented_triangle(tri, 3, 140, 1126) == (3, 140, 1126)
@test DT.construct_positively_oriented_triangle(tri, 140, 3, 1126) == (3, 140, 1126)
@inferred DT.construct_positively_oriented_triangle(tri, 3, 140, 1126)
_solid_itr = each_solid_triangle(tri)
@test Base.IteratorSize(_solid_itr) == Base.SizeUnknown()
@test Base.IteratorEltype(_solid_itr) == Base.HasEltype()
@test Base.eltype(_solid_itr) == NTuple{3,Int64}
@test each_solid_triangle(tri) isa DT.EachSolidTriangle
_solid_tri = collect(_solid_itr)
@inferred collect(_solid_itr)
@test all(!DT.is_ghost_triangle, _solid_tri)
_ghost_itr = each_ghost_triangle(tri)
@test Base.IteratorSize(_ghost_itr) == Base.SizeUnknown()
@test Base.IteratorEltype(_ghost_itr) == Base.HasEltype()
@test Base.eltype(_ghost_itr) == NTuple{3,Int64}
@test each_ghost_triangle(tri) isa DT.EachGhostTriangle
_ghost_tri = collect(_ghost_itr)
@inferred collect(_ghost_itr)
@test all(DT.is_ghost_triangle, _ghost_tri)
@test length(_ghost_tri) + length(_solid_tri) == num_triangles(tri)
___tri = generate_mesh(x, y, 0.1; convert_result=true, add_ghost_triangles=true)
DT.delete_ghost_triangles!(___tri)
@test collect(each_triangle(___tri)) == collect(each_solid_triangle(___tri))
@test length(collect(each_ghost_triangle(___tri))) == 0

# Edges 
@test DT.edge_type(tri) == NTuple{2,Int64}
@inferred DT.edge_type(tri)
@test DT.num_edges(tri) == length(tri.graph.graph.E)
@inferred DT.num_edges(tri)
@test DT.each_edge(tri) == tri.graph.graph.E
@inferred DT.each_edge(tri)

# Points
@test DT.get_point(tri, 2) == Tuple(tri.points[:, 2])
@inferred DT.get_point(tri, 2)
@test DT.get_point(tri, 17) == Tuple(tri.points[:, 17])
@test DT.each_point_index(tri) == axes(tri.points, 2)
@inferred DT.each_point_index(tri)
@test DT.each_point(tri) == eachcol(tri.points)
@inferred DT.each_point(tri)
@test DT.num_points(tri) == size(tri.points, 2)
@inferred DT.num_points(tri)
@test DT.get_point(tri, 2, 5, 6, 10) ==
      ntuple(i -> Tuple(tri.points[:, (2, 5, 6, 10)[i]]), 4)
@inferred DT.get_point(tri, 2, 5, 6, 10)
DT.RepresentativePointList[1] = DT.RepresentativeCoordinates(0.5, 0.3, 2)
DT.RepresentativePointList[2] = DT.RepresentativeCoordinates(2.5, 7.3, 7)
DT.RepresentativePointList[3] = DT.RepresentativeCoordinates(6.5, -0.6, 13)
DT.RepresentativePointList[4] = DT.RepresentativeCoordinates(6.5, -0.66234, 13)
DT.RepresentativePointList[5] = DT.RepresentativeCoordinates(6.534234, -0.6, 13)
@test DT.get_point(tri, -1) == DT.getxy(DT.RepresentativePointList[1])
@test DT.get_point(tri, -2) == DT.getxy(DT.RepresentativePointList[1])
@test DT.get_point(tri, -3) == DT.getxy(DT.RepresentativePointList[1])
@test DT.get_point(tri, -4) == DT.getxy(DT.RepresentativePointList[1])
@test DT.get_point(tri, -5) == DT.getxy(DT.RepresentativePointList[2])
@test DT.get_point(tri, -6) == DT.getxy(DT.RepresentativePointList[3])
@test DT.get_point(tri, -7) == DT.getxy(DT.RepresentativePointList[4])
@test DT.get_point(tri, -8) == DT.getxy(DT.RepresentativePointList[4])
@test DT.get_point(tri, -9) == DT.getxy(DT.RepresentativePointList[4])
@test DT.get_point(tri, -10) == DT.getxy(DT.RepresentativePointList[4])
@test DT.get_point(tri, -11) == DT.getxy(DT.RepresentativePointList[5])
@test DT.get_point(tri_2, -1) == DT.getxy(DT.RepresentativePointList[1])
@inferred DT.get_point(tri, -9)
@test DT.get_point(tri_2, -2) == DT.getxy(DT.RepresentativePointList[1])
@test DT.get_point(tri_2, -3) == DT.getxy(DT.RepresentativePointList[1])
@test DT.get_point(tri_2, -4) == DT.getxy(DT.RepresentativePointList[1])
@test DT.get_point(tri_3, -1) == DT.getxy(DT.RepresentativePointList[1])
@test_throws KeyError DT.get_point(tri, -12)
@test_throws KeyError DT.get_point(tri_2, -5)
@test_throws KeyError DT.get_point(tri_3, -2)
@test get_point(tri, tri.points[:, 2]) == (tri.points[1, 2], tri.points[2, 2])
@inferred get_point(tri, tri.points[:, 2]) == (tri.points[1, 2], tri.points[2, 2])
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

# Miscellaneous
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
@test DT.is_constrained(_tri)
empty!(_tri.boundary_nodes)
_tri = @set _tri.boundary_nodes = Int64[]
@test !DT.is_constrained(_tri)
_tri = @set _tri.boundary_nodes = [2, 3, 4, 5, 6]
@test DT.is_constrained(_tri)
@test collect(DT.all_boundary_indices(_tri)) ==
      [DT.BoundaryIndex, DT.BoundaryIndex - 1, DT.BoundaryIndex - 2, DT.BoundaryIndex - 3]
