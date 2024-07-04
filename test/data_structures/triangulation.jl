using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using StructEquality
using Setfield
using Random
using StatsBase
using StableRNGs
using ..DelaunayTriangulation: add_weight!, get_weight, get_weights

@struct_equal DT.TriangulationCache
include("../helper_functions.jl")

using ..DelaunayTriangulation: Triangulation

global pts = rand(2, 500)
global tri = Triangulation(pts; IntegerType=Int32)

@testset "Initialising a triangulation" begin
      @test tri ⊢ Triangulation(pts, # points
            Set{NTuple{3,Int32}}(), # triangles
            Int32[], # boundary_nodes
            Set{NTuple{2,Int32}}(), # interior_segments
            Set{NTuple{2,Int32}}(), # all_segments
            DT.ZeroWeight(), # weights
            DT.Adjacent{Int32,NTuple{2,Int32}}(), # adjacent
            DT.Adjacent2Vertex{Int32,Set{NTuple{2,Int32}}}(), # adjacent2vertex
            DT.Graph{Int32}(), # graph
            (), # boundary_curves
            Dict{NTuple{2,Int32},Tuple{Vector{Int32},Int32}}(), # boundary_edge_map
            Dict{Int32,Vector{Int32}}(DT.𝒢 => Int32[]), # ghost_vertex_map
            Dict{Int32,UnitRange{Int32}}(-1 => -1:-1), # ghost_vertex_ranges
            DT.ConvexHull(pts, Int32[]), # convex_hull
            Dict{Int32,DT.RepresentativeCoordinates{Int32,Float64}}(), # representative_point_coordinate
            DT.construct_polygon_hierarchy(pts; IntegerType=Int32),
            nothing,
            DT.TriangulationCache(
                  DT.Triangulation(
                        pts,
                        Set{NTuple{3,Int32}}(), # triangles
                        Int32[], # boundary_nodes
                        Set{NTuple{2,Int32}}(), # interior_segments
                        Set{NTuple{2,Int32}}(), # all_segments
                        DT.ZeroWeight(), # weights
                        DT.Adjacent{Int32,NTuple{2,Int32}}(), # adjacent
                        DT.Adjacent2Vertex{Int32,Set{NTuple{2,Int32}}}(), # adjacent2vertex
                        DT.Graph{Int32}(), # graph
                        (), # boundary_curves
                        Dict{NTuple{2,Int32},Tuple{Vector{Int32},Int32}}(), # boundary_edge_map
                        Dict{Int32,Vector{Int32}}(DT.𝒢 => Int32[]), # ghost_vertex_map
                        Dict{Int32,UnitRange{Int32}}(-1 => -1:-1), # ghost_vertex_ranges
                        DT.ConvexHull(pts, Int32[]), # convex_hull
                        Dict{Int32,DT.RepresentativeCoordinates{Int32,Float64}}(), # representative_point_coordinate
                        DT.construct_polygon_hierarchy(pts; IntegerType=Int32),
                        nothing, # boundary_enricher
                        DT.TriangulationCache(nothing, nothing, nothing, nothing, nothing, nothing)
                  ),
                  DT.Triangulation(
                        pts,
                        Set{NTuple{3,Int32}}(), # triangles
                        Int32[], # boundary_nodes
                        Set{NTuple{2,Int32}}(), # interior_segments
                        Set{NTuple{2,Int32}}(), # all_segments
                        DT.ZeroWeight(), # weights
                        DT.Adjacent{Int32,NTuple{2,Int32}}(), # adjacent
                        DT.Adjacent2Vertex{Int32,Set{NTuple{2,Int32}}}(), # adjacent2vertex
                        DT.Graph{Int32}(), # graph
                        (), # boundary_curves
                        Dict{NTuple{2,Int32},Tuple{Vector{Int32},Int32}}(), # boundary_edge_map
                        Dict{Int32,Vector{Int32}}(DT.𝒢 => Int32[]), # ghost_vertex_map
                        Dict{Int32,UnitRange{Int32}}(-1 => -1:-1), # ghost_vertex_ranges
                        DT.ConvexHull(pts, Int32[]), # convex_hull
                        Dict{Int32,DT.RepresentativeCoordinates{Int32,Float64}}(), # representative_point_coordinate
                        DT.construct_polygon_hierarchy(pts; IntegerType=Int32),
                        nothing, # boundary_enricher
                        DT.TriangulationCache(nothing, nothing, nothing, nothing, nothing, nothing)
                  ),
                  Int32[],
                  Set{NTuple{2,Int32}}(),
                  Vector{Int32}(),
                  Set{NTuple{3,Int32}}())
      )
end

_x, _y = complicated_geometry()
global x = _x
global y = _y
boundary_nodes, points = convert_boundary_points_to_indices(x, y)
global tri = triangulate(points; boundary_nodes, delete_ghosts=false)
A = get_area(tri)
refine!(tri; max_area=1e-1A, use_circumcenter=true, use_lens=false)
boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
global tri_2 = triangulate(points; boundary_nodes, delete_ghosts=false)
A = get_area(tri_2)
refine!(tri_2; max_area=1e-1A, use_circumcenter=true, use_lens=false)
boundary_nodes, points = convert_boundary_points_to_indices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
global tri_3 = triangulate(points; boundary_nodes, delete_ghosts=false)
A = get_area(tri_3)
refine!(tri_3; max_area=1e-1A, use_circumcenter=true, use_lens=false)
boundary_nodes, points = convert_boundary_points_to_indices(reverse(reverse.(x[2])), reverse(reverse.(y[2])))
global tri_4 = triangulate(points; boundary_nodes, delete_ghosts=false)
A = get_area(tri_4)
refine!(tri_4; max_area=1e-1A, use_circumcenter=true, use_lens=false)

@testset "Triangulation getters" begin
      @test DT.get_points(tri) == tri.points
      @test DT.get_triangles(tri) == tri.triangles
      @test DT.get_boundary_nodes(tri) == tri.boundary_nodes
      @test DT.get_interior_segments(tri) == tri.interior_segments
      @test DT.get_all_segments(tri) == tri.all_segments
      @test DT.get_weights(tri) == tri.weights
      @test DT.get_adjacent(tri) == tri.adjacent
      @test DT.get_adjacent2vertex(tri) == tri.adjacent2vertex
      @test DT.get_graph(tri) == tri.graph
      @test DT.get_boundary_curves(tri) == tri.boundary_curves
      @test DT.get_boundary_edge_map(tri) == tri.boundary_edge_map
      @test DT.get_ghost_vertex_map(tri) == tri.ghost_vertex_map
      @test DT.get_ghost_vertex_ranges(tri) == tri.ghost_vertex_ranges
      @test DT.get_convex_hull(tri) == tri.convex_hull
      @test DT.get_representative_point_list(tri) == tri.representative_point_list
      @test DT.get_exterior_curve_indices(tri) == keys(tri.polygon_hierarchy.trees)
      @test DT.get_boundary_enricher(tri) === tri.boundary_enricher
      @test compare_trees(DT.get_polygon_hierarchy(tri), DT.construct_polygon_hierarchy(tri.points, tri.boundary_nodes))
      @test compare_edge_vectors(collect(DT.get_all_segments(tri)), collect(tri.all_segments))
      @test compare_edge_vectors(collect(DT.get_all_segments(tri)), collect(DT.merge_segments(get_ghost_vertex_map(tri), get_boundary_nodes(tri), get_interior_segments(tri))))
      @test DT.get_ghost_vertex_map(tri) == tri.ghost_vertex_map == DT.construct_ghost_vertex_map(get_boundary_nodes(tri))
      @test DT.get_ghost_vertex_ranges(tri) == tri.ghost_vertex_ranges == DT.construct_ghost_vertex_ranges(get_boundary_nodes(tri))
      @test DT.get_boundary_edge_map(tri) == tri.boundary_edge_map == DT.construct_boundary_edge_map(get_boundary_nodes(tri))
      @test DT.get_cache(tri) == tri.cache
      @inferred DT.get_points(tri)
      @inferred DT.get_triangles(tri)
      @inferred DT.get_boundary_nodes(tri)
      @inferred DT.get_interior_segments(tri)
      @inferred DT.get_all_segments(tri)
      @inferred DT.get_weights(tri)
      @inferred DT.get_adjacent(tri)
      @inferred DT.get_adjacent2vertex(tri)
      @inferred DT.get_graph(tri)
      @inferred DT.get_boundary_curves(tri)
      @inferred DT.get_boundary_edge_map(tri)
      @inferred DT.get_ghost_vertex_map(tri)
      @inferred DT.get_ghost_vertex_ranges(tri)
      @inferred DT.get_convex_hull(tri)
      @inferred DT.get_representative_point_list(tri)
      @inferred DT.get_exterior_curve_indices(tri)
end

@testset "Forwarded methods" begin
      @testset "Adjacent" begin
            @test DT.get_adjacent(tri, 92, -6) == tri.adjacent.adjacent[(92, -6)]
            @test DT.get_adjacent(tri, (723, 1356)) == get(tri.adjacent.adjacent, (723, 1356), DT.∅)
            @inferred DT.get_adjacent(tri, (723, 1356))
            DT.add_adjacent!(tri, 101117, 20311, 5)
            DT.add_adjacent!(tri, (27311, 50511), 10)
            @test DT.get_adjacent(tri, 101117, 20311) == 5
            @inferred DT.get_adjacent(tri, 101117, 20311)
            @test DT.get_adjacent(tri, 27311, 50511) == 10
            DT.delete_adjacent!(tri, 101117, 20311)
            DT.delete_adjacent!(tri, (27311, 50511))
            @test DT.get_adjacent(tri, 101117, 20311) == DT.∅
            @test DT.get_adjacent(tri, 27311, 50511) == DT.∅
            @test DT.get_adjacent(tri, 711, DT.𝒢) == DT.∅
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
            @test DT.get_edges(tri) == tri.graph.edges == DT.get_edges(tri.graph)
            @inferred DT.get_edges(tri)
            @test DT.get_neighbours(tri, 3) == tri.graph.neighbours[3]
            @inferred DT.get_neighbours(tri, 3)
            DT.add_vertex!(tri, 19998, 23721)
            DT.add_vertex!(tri, 28371)
            @test all(∈(tri.graph.vertices), (19998, 23721, 28371))
            DT.add_neighbour!(tri, 28371, 50912)
            DT.add_neighbour!(tri, 28371, 271, 501)
            @test all(∈(tri.graph.neighbours[28371]), (50912, 271, 501))
            DT.delete_neighbour!(tri, 28371, 50912)
            DT.delete_neighbour!(tri, 28371, 271, 501)
            @test all(∉(tri.graph.neighbours[28371]), (50912, 271, 501))
            DT.delete_vertex!(tri, 19998)
            DT.delete_vertex!(tri, 28371, 3)
            @test all(∉(tri.graph.vertices), (19998, 28371, 3))
            DT.delete_ghost_vertices_from_graph!(tri)
            @test all(∉(tri.graph.vertices), -11:-1)
            @test DT.get_neighbours(tri) == tri.graph.neighbours
            @test DT.get_vertices(tri) == tri.graph.vertices
            @test DT.get_vertices(tri) == each_vertex(tri)
            @test num_vertices(tri) == length(DT.get_vertices(tri))
            @test all((DT.num_neighbours(tri, u) == DT.num_neighbours(DT.get_graph(tri), u) for u in DT.get_vertices(tri)))
      end

      @testset "ConvexHull" begin
            @test DT.get_convex_hull(tri) == tri.convex_hull
            @test DT.get_convex_hull_vertices(tri) == tri.convex_hull.vertices
            _tri = triangulate_rectangle(0.0, 2.0, 5.0, 7.3, 5, 15; single_boundary=true)
            points = get_points(_tri)
            ch = convex_hull(points)
            @test ch.vertices == _tri.boundary_nodes
            _indices = deepcopy(ch.vertices)
            __indices = DT.get_convex_hull_vertices(_tri)
            empty!(__indices)
            DT.convex_hull!(_tri; reconstruct=false)
            @test length(__indices) == length(_indices)
            unique!(__indices)
            unique!(ch.vertices)
            shift = findfirst(ch.vertices .== first(__indices))
            ch.vertices .= circshift(ch.vertices, 1 - shift)
            @test __indices == ch.vertices
            DT.convex_hull!(_tri)
            @test length(__indices) == length(_indices)
            unique!(__indices)
            shift = findfirst(ch.vertices .== first(__indices))
            ch.vertices .= circshift(ch.vertices, 1 - shift)
            @test __indices == ch.vertices
            delete_ghost_triangles!(_tri)
            DT.convex_hull!(_tri; reconstruct=false)
            @test length(__indices) == length(_indices)
            unique!(__indices)
            unique!(ch.vertices)
            shift = findfirst(ch.vertices .== first(__indices))
            ch.vertices .= circshift(ch.vertices, 1 - shift)
            @test __indices == ch.vertices
            DT.convex_hull!(_tri)
            @test length(__indices) == length(_indices)
            unique!(__indices)
            shift = findfirst(ch.vertices .== first(__indices))
            ch.vertices .= circshift(ch.vertices, 1 - shift)
            @test __indices == ch.vertices
            @test !DT.has_ghost_triangles(_tri)
      end

      @testset "Boundary nodes" begin
            @test DT.has_multiple_curves(tri)
            @inferred DT.has_multiple_curves(tri)
            @test !DT.has_multiple_curves(tri_2)
            @test DT.has_multiple_sections(tri)
            @test DT.has_multiple_sections(tri_2)
            @test !DT.has_multiple_sections(tri_3)
            @inferred DT.has_multiple_sections(tri_2)
            @test DT.num_curves(tri) == 5
            @inferred DT.num_curves(tri)
            a, b = 0.0, 5.0
            c, d = 3.0, 7.0
            nx = 12
            ny = 15
            @test DT.num_curves(triangulate_rectangle(0.0, 1.0, 0.0, 1.0, 10, 10; delete_ghosts=false, single_boundary=false)) == 1
            @test DT.num_sections(tri_2) == 4
            @inferred DT.num_sections(tri_2)
            @test DT.get_boundary_nodes(tri, 1) == tri.boundary_nodes[1]
            @test DT.get_boundary_nodes(tri, 1, 3) == tri.boundary_nodes[1][3]
            @test DT.get_boundary_nodes(tri, (5, 1)) == tri.boundary_nodes[5][1]
            @test DT.get_boundary_nodes(tri_2, 1) == tri_2.boundary_nodes[1]
            @test DT.get_boundary_nodes(tri_2, 3) == tri_2.boundary_nodes[3]
            @test DT.get_boundary_nodes(tri_3) == tri_3.boundary_nodes
            @test DT.get_boundary_nodes(tri_2, 3) == tri_2.boundary_nodes[3]
            @test DT.get_curve_index(tri, -1) == 1
            @test DT.get_curve_index(tri, -3) == 1
            @test DT.get_curve_index(tri, -5) == 2
            @test DT.get_curve_index(tri, -11) == 5
            @test DT.get_curve_index(tri_2, -3) == 1
            @test DT.get_curve_index(tri_2, -2) == 1
            @test DT.get_curve_index(tri_3, -1) == 1
            @inferred DT.get_curve_index(tri_3, -1)
      end

      @testset "Triangles" begin
            rng = StableRNG(9882881)
            boundary_nodes, points = convert_boundary_points_to_indices(x, y)
            tri = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
            A = get_area(tri)
            refine!(tri; max_area=1e-1A, rng, use_circumcenter=true)
            @test DT.triangle_type(tri) == NTuple{3,Int}
            @inferred DT.triangle_type(tri)
            @test DT.num_triangles(tri) == length(tri.triangles)
            @test DT.each_triangle(tri) == tri.triangles
            @inferred DT.num_triangles(tri)
            @test DT.contains_triangle(tri, (68, 67, -6))[2]
            @inferred DT.contains_triangle(tri, (68, 67, -6))
            @test !DT.contains_triangle(tri, (1, 1, 5))[2]
            @test !DT.contains_triangle(tri, 1, 5, 5)[2]
            @inferred DT.contains_triangle(tri, 3, 140, 1126)
            @test DT.construct_positively_oriented_triangle(tri, 1, 10, 20) == (10, 1, 20)
            @inferred DT.construct_positively_oriented_triangle(tri, 188, 394, 426)
            _solid_itr = each_solid_triangle(tri)
            @test DelaunayTriangulation.each_triangle(_solid_itr) == _solid_itr
            @test Base.IteratorSize(_solid_itr) == Base.HasLength()
            @test Base.IteratorEltype(_solid_itr) == Base.HasEltype()
            @test Base.eltype(_solid_itr) == NTuple{3,Int}
            @test each_solid_triangle(tri) isa DT.EachSolidTriangle
            _solid_tri = collect(_solid_itr)
            @inferred collect(_solid_itr)
            @test all(!DT.is_ghost_triangle, _solid_tri)
            _ghost_itr = each_ghost_triangle(tri)
            @test Base.IteratorSize(_ghost_itr) == Base.HasLength()
            @test Base.IteratorEltype(_ghost_itr) == Base.HasEltype()
            @test Base.eltype(_ghost_itr) == NTuple{3,Int}
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
            A = get_area(___tri)
            refine!(___tri; max_area=1e-1A, rng, use_circumcenter=true)
            DT.delete_ghost_triangles!(___tri)
            @test collect(each_triangle(___tri)) == collect(each_solid_triangle(___tri))
            @test length(collect(each_ghost_triangle(___tri))) == 0
            @test sort(collect(filter(!DT.is_ghost_triangle, each_triangle(___tri)))) == sort(collect(each_solid_triangle(___tri)))
            @test sort(collect(filter(DT.is_ghost_triangle, each_triangle(___tri)))) == sort(collect(each_ghost_triangle(___tri)))
            @test DelaunayTriangulation.triangle_type(_ghost_itr) == DelaunayTriangulation.triangle_type(_ghost_itr)
            @test DelaunayTriangulation.triangle_type(_solid_itr) == DelaunayTriangulation.triangle_type(_solid_itr)
      end

      @testset "Edges" begin
            rng = StableRNG(998871)
            boundary_nodes, points = convert_boundary_points_to_indices(x, y)
            tri = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
            A = get_area(tri)
            refine!(tri; max_area=1e-1A, rng, use_circumcenter=true)
            @test DT.edge_type(tri) == NTuple{2,Int}
            @inferred DT.edge_type(tri)
            @test DT.num_edges(tri) == length(tri.graph.edges)
            @inferred DT.num_edges(tri)
            @test DT.each_edge(tri) == tri.graph.edges
            @inferred DT.each_edge(tri)
            _solid_itr = each_solid_edge(tri)
            @test DelaunayTriangulation.each_edge(_solid_itr) == _solid_itr
            @test Base.IteratorSize(_solid_itr) == Base.HasLength()
            @test Base.IteratorEltype(_solid_itr) == Base.HasEltype()
            @test Base.eltype(_solid_itr) == NTuple{2,Int}
            @test each_solid_edge(tri) isa DT.EachSolidEdge
            _solid_tri = collect(_solid_itr)
            @inferred collect(_solid_itr)
            @test all(!DT.is_ghost_edge, _solid_tri)
            _ghost_itr = each_ghost_edge(tri)
            @test Base.IteratorSize(_ghost_itr) == Base.HasLength()
            @test Base.IteratorEltype(_ghost_itr) == Base.HasEltype()
            @test Base.eltype(_ghost_itr) == NTuple{2,Int}
            @test each_ghost_edge(tri) isa DT.EachGhostEdge
            _ghost_tri = collect(_ghost_itr)
            @inferred collect(_ghost_itr)
            @test DelaunayTriangulation.each_edge(_ghost_itr) == _ghost_itr
            @test all(DT.is_ghost_edge, _ghost_tri)
            @test length(_ghost_tri) + length(_solid_tri) == num_edges(tri)
            rng = StableRNG(998871)
            boundary_nodes, points = convert_boundary_points_to_indices(x, y)
            ___tri = triangulate(points; rng, boundary_nodes, delete_ghosts=false)
            A = get_area(___tri)
            refine!(___tri; max_area=1e-1A, rng, use_circumcenter=true)
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
            @test reverse(sort(collect(DT.all_ghost_vertices(tri)))) == [DT.𝒢,
                  DT.𝒢 - 1,
                  DT.𝒢 - 2,
                  DT.𝒢 - 3,
                  DT.𝒢 - 4,
                  DT.𝒢 - 5,
                  DT.𝒢 - 6,
                  DT.𝒢 - 7,
                  DT.𝒢 - 8,
                  DT.𝒢 - 9,
                  DT.𝒢 - 10]#(new_point, c′) = (436, (1.8363241117152609, 1.0864161502095007))
            rng = StableRNG(887271)
            boundary_nodes, points = convert_boundary_points_to_indices(x, y)
            _triq = triangulate(points; boundary_nodes, rng)
            A = get_area(_triq)
            refine!(_triq; max_area=1e-1A, rng, use_circumcenter=true, use_lens=false)
            _solid_itr = each_solid_vertex(_triq)
            @test DelaunayTriangulation.each_vertex(_solid_itr) == _solid_itr
            @test Base.IteratorSize(_solid_itr) == Base.HasLength()
            @test Base.IteratorEltype(_solid_itr) == Base.HasEltype()
            @test Base.eltype(_solid_itr) == Int
            @test each_solid_vertex(_triq) isa DT.EachSolidVertex
            _solid_tri = collect(_solid_itr)
            @inferred collect(_solid_itr)
            @test sort(_solid_tri) == sort(filter(i -> DT.has_vertex(_triq, i) && i ≥ 1, 1:DT.num_points(tri)))
            @test all(!DT.is_ghost_vertex, _solid_tri)
            _ghost_itr = each_ghost_vertex(_triq)
            @test Base.IteratorSize(_ghost_itr) == Base.HasLength()
            @test Base.IteratorEltype(_ghost_itr) == Base.HasEltype()
            @test Base.eltype(_ghost_itr) == Int
            @test each_ghost_vertex(_triq) isa DT.EachGhostVertex
            _ghost_tri = collect(_ghost_itr)
            @inferred collect(_ghost_itr)
            @test DelaunayTriangulation.each_vertex(_ghost_itr) == _ghost_itr
            @test all(DT.is_ghost_vertex, _ghost_tri)
            @test reverse(sort(_ghost_tri)) == [DT.𝒢,
                  DT.𝒢 - 1,
                  DT.𝒢 - 2,
                  DT.𝒢 - 3,
                  DT.𝒢 - 4,
                  DT.𝒢 - 5,
                  DT.𝒢 - 6,
                  DT.𝒢 - 7,
                  DT.𝒢 - 8,
                  DT.𝒢 - 9,
                  DT.𝒢 - 10]
            @test length(_ghost_tri) == length(_ghost_itr) == sum(<(0), keys(_triq.adjacent2vertex.adjacent2vertex))
            @test length(_solid_tri) == length(_solid_itr) == sort(filter(i -> DT.has_vertex(_triq, i) && i ≥ 1, 1:DT.num_points(tri))) |> length
            rng = StableRNG(887271)
            ___tri = triangulate(points; boundary_nodes, rng)
            A = get_area(___tri)
            refine!(___tri; max_area=1e-1A, rng, use_circumcenter=true)
            DT.delete_ghost_triangles!(___tri)
            @test sort(collect(filter(!DT.is_ghost_vertex, each_vertex(___tri)))) == sort(collect(each_solid_vertex(___tri)))
            @test sort(collect(filter(DT.is_ghost_vertex, each_vertex(___tri)))) == sort(collect(each_ghost_vertex(___tri)))
            @test length(collect(each_ghost_vertex(___tri))) == num_vertices(___tri) .- length(sort(collect(filter(!DT.is_ghost_vertex, each_vertex(___tri)))))
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
            @test DT.integer_type(tri) == Int
            @test DT.number_type(tri) == Float64
            @test DT.edge_type(tri) == NTuple{2,Int}
            @test DT.edges_type(tri) == Set{NTuple{2,Int}}
            @test DT.triangles_type(tri) == Set{NTuple{3,Int}}
            @test DT.triangle_type(tri) == NTuple{3,Int}
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
            push!(get_all_segments(_tri), (1, 2))
            @test DT.is_constrained(_tri)
            empty!(get_all_segments(_tri))
            @test !DT.is_constrained(_tri)
            @test reverse(sort(collect(DT.all_ghost_vertices(_tri)))) ==
                  [DT.𝒢, DT.𝒢 - 1, DT.𝒢 - 2, DT.𝒢 - 3]
      end
end

@testset "merge_segments" begin
      all_bn = get_boundary_nodes(tri)
      i = rand(1:100000, 50)
      j = rand(1:100000, 50)
      all_ce = Set(((i, j) for (i, j) in zip(i, j)))
      bn_map = get_ghost_vertex_map(tri)
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
      ace = Set{NTuple{2,Int}}()
      for es in (e11, e12, e13, e14, e2, e3, e41, e42, e43, e44, e5)
            for e in es
                  push!(ace, e)
            end
      end
      for e in all_ce
            push!(ace, e)
      end
      @test ace == DT.merge_segments(bn_map, all_bn, all_ce)

      all_bn = get_boundary_nodes(tri_2)
      i = rand(1:100000, 50)
      j = rand(1:100000, 50)
      all_ce = Set(((i, j) for (i, j) in zip(i, j)))
      bn_map = get_ghost_vertex_map(tri_2)
      bn1 = all_bn[1]
      bn2 = all_bn[2]
      bn3 = all_bn[3]
      bn4 = all_bn[4]
      e1 = Set(((bn1[i], bn1[i+1]) for i in 1:(length(bn1)-1)))
      e2 = Set(((bn2[i], bn2[i+1]) for i in 1:(length(bn2)-1)))
      e3 = Set(((bn3[i], bn3[i+1]) for i in 1:(length(bn3)-1)))
      e4 = Set(((bn4[i], bn4[i+1]) for i in 1:(length(bn4)-1)))
      ace = Set{NTuple{2,Int}}()
      for es in (e1, e2, e3, e4)
            for e in es
                  push!(ace, e)
            end
      end
      for e in all_ce
            push!(ace, e)
      end
      @test ace == DT.merge_segments(bn_map, all_bn, all_ce)

      all_bn = get_boundary_nodes(tri_3)
      i = rand(1:100000, 50)
      j = rand(1:100000, 50)
      all_ce = Set(((i, j) for (i, j) in zip(i, j)))
      bn_map = get_ghost_vertex_map(tri_3)
      e = Set(((all_bn[i], all_bn[i+1]) for i in 1:(length(all_bn)-1)))
      ace = Set{NTuple{2,Int}}()
      for e in e
            push!(ace, e)
      end
      for e in all_ce
            push!(ace, e)
      end
      @test ace == DT.merge_segments(bn_map, all_bn, all_ce)

      all_bn = get_boundary_nodes(tri_4)[1]
      i = rand(1:100000, 50)
      j = rand(1:100000, 50)
      all_ce = Set(((i, j) for (i, j) in zip(i, j)))
      bn_map = get_ghost_vertex_map(tri_4)
      e = Set(((all_bn[i], all_bn[i+1]) for i in 1:(length(all_bn)-1)))
      ace = Set{NTuple{2,Int}}()
      for e in e
            push!(ace, e)
      end
      for e in all_ce
            push!(ace, e)
      end
      @test ace == DT.merge_segments(bn_map, [all_bn], all_ce)
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
            tri1 = triangulate(12randn(2, n), delete_ghosts=false)
            tri2 = triangulate(12randn(2, n), delete_ghosts=true)
            for tri in (tri1, tri2)
                  for qi in each_solid_vertex(tri)
                        for k in each_solid_vertex(tri)
                              q = get_point(tri, qi)
                              history = DT.PointLocationHistory{NTuple{3,Int},NTuple{2,Int},Int}()
                              jump_and_march(tri, q;
                                    k,
                                    store_history=true,
                                    history)
                              visited_triangles = history.triangles
                              collinear_segments = history.collinear_segments
                              @test all(T -> DT.is_positively_oriented(DT.triangle_orientation(tri, T...)), visited_triangles)
                              @test all(!DT.is_none, [DT.triangle_line_segment_intersection(tri, T..., (qi, k)...) for T in visited_triangles])
                              @test allunique(visited_triangles)
                              if !isempty(collinear_segments)
                                    @test all(E -> DT.is_collinear(DT.point_position_relative_to_line(tri, qi, k, E[1])), collinear_segments)
                                    @test all(E -> DT.is_collinear(DT.point_position_relative_to_line(tri, qi, k, E[2])), collinear_segments)
                              end
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

@testset "each_segment" begin
      @test each_segment(tric) == each_edge(get_all_segments(tric))
end

@testset "contains_segment" begin
      @test !DT.contains_segment(tric, 12, 17)
      @test DT.contains_segment(tric, 69, 70)
      @test DT.contains_segment(tric, 70, 69)
      @test !DT.contains_segment(tric, 32, 41)
      @test !DT.contains_segment(tric, 45, 38)
      @test DT.contains_segment(tric, 63, 64)
      @test !DT.contains_segment(tric, (45, 38))
      @test !DT.contains_segment(tric, 26, 22)
      @test DT.contains_segment(tric, 64, 65)
      @test DT.contains_segment(tric, 55, 54)
      @test DT.contains_segment(tric, 58, 57)
      @test DT.contains_segment(tric, 59, 60)
      @test !DT.contains_segment(tric, 30, 70)
      @test !DT.contains_segment(tric, 56, 37)
      @test DT.contains_segment(tric, 73, 74)
end

@testset "get_all_boundary_nodes" begin
      x, y = complicated_geometry()
      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x, y)
      tri = triangulate(points; rng, boundary_nodes)
      A = get_area(tri)
      refine!(tri, rng=rng, max_area=1e-1A, use_circumcenter=true, use_lens=false)
      all_bn = DT.get_all_boundary_nodes(tri)
      @test all_bn == Set(reduce(vcat, reduce(vcat, get_boundary_nodes(tri))))
      tri2, label_map, index_map = simple_geometry()
      all_bn = DT.get_all_boundary_nodes(tri2)
      @test all_bn == Set(reduce(vcat, reduce(vcat, get_boundary_nodes(tri2))))
      tri3 = triangulate_rectangle(0, 1, 0, 1, 50, 50; delete_ghosts=false, single_boundary=false)
      all_bn = DT.get_all_boundary_nodes(tri3)
      @test all_bn == Set(reduce(vcat, reduce(vcat, get_boundary_nodes(tri3))))
      tri4 = triangulate_rectangle(0, 1, 0, 1, 50, 50; delete_ghosts=false, single_boundary=true)
      all_bn = DT.get_all_boundary_nodes(tri4)
      @test all_bn == Set(reduce(vcat, reduce(vcat, get_boundary_nodes(tri4))))
      tri = triangulate(rand(2, 50))
      @test isempty(DT.get_all_boundary_nodes(tri))
end

@testset "get_boundary_edge_map" begin
      x, y = complicated_geometry()
      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x, y)
      tri = triangulate(points; rng, boundary_nodes)
      A = get_area(tri)
      refine!(tri, rng=rng, max_area=1e-2A, use_circumcenter=true)
      for (e, (s, i)) in tri.boundary_edge_map
            @test DT.get_boundary_edge_map(tri, e) == (s, i)
            @test DT.get_boundary_edge_map(tri, e...) == (s, i)
            @test get_boundary_nodes(DT.get_boundary_nodes(tri, s), i) == e[1]
      end
end

@testset "split_boundary_edge!" begin
      x, y = complicated_geometry()
      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x, y)
      tri_1 = triangulate(points; rng, boundary_nodes)
      A = get_area(tri_1)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
      tri_2 = triangulate(points; rng, boundary_nodes)
      A = get_area(tri_2)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
      tri_3 = triangulate(points; rng, boundary_nodes)
      A = get_area(tri_3)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(reverse(reverse.(x[2])), reverse(reverse.(y[2])))
      tri_4 = triangulate(points; rng, boundary_nodes)
      A = get_area(tri_4)

      DT.split_boundary_edge!(tri_1, (21, 22), 500)
      @test tri_1.boundary_nodes[2][1][1:11] == [13, 14, 15, 16, 17, 18, 19, 20, 21, 500, 22]
      @test DT.get_boundary_edge_map(tri_1, 21, 500) == ((2, 1), 9)
      @test DT.get_boundary_edge_map(tri_1, 500, 22) == ((2, 1), 10)
      @test_throws KeyError DT.get_boundary_edge_map(tri_1, 21, 22)
      @test !DT.contains_unoriented_edge((21, 22), DT.get_all_segments(tri_1))
      @test DT.contains_unoriented_edge((21, 500), DT.get_all_segments(tri_1))
      @test DT.contains_unoriented_edge((500, 22), DT.get_all_segments(tri_1))

      DT.split_boundary_edge!(tri_1, (7, 8), 5000)
      @test tri_1.boundary_nodes[1][3] == [7, 5000, 8, 9, 10]
      @test DT.get_boundary_edge_map(tri_1, 7, 5000) == ((1, 3), 1)
      @test DT.get_boundary_edge_map(tri_1, 5000, 8) == ((1, 3), 2)
      @test_throws KeyError DT.get_boundary_edge_map(tri_1, 7, 8)
      @test !DT.contains_unoriented_edge((7, 8), DT.get_all_segments(tri_1))
      @test DT.contains_unoriented_edge((7, 5000), DT.get_all_segments(tri_1))
      @test DT.contains_unoriented_edge((5000, 8), DT.get_all_segments(tri_1))

      DT.split_boundary_edge!(tri_2, 8, 9, 300)
      @test tri_2.boundary_nodes[3] == [7, 8, 300, 9, 10]
      @test DT.get_boundary_edge_map(tri_2, 8, 300) == (3, 2)
      @test DT.get_boundary_edge_map(tri_2, 300, 9) == (3, 3)
      @test_throws KeyError DT.get_boundary_edge_map(tri_2, 8, 9)
      @test !DT.contains_unoriented_edge((8, 9), DT.get_all_segments(tri_2))
      @test DT.contains_unoriented_edge((8, 300), DT.get_all_segments(tri_2))
      @test DT.contains_unoriented_edge((300, 9), DT.get_all_segments(tri_2))

      DT.split_boundary_edge!(tri_3, 3, 4, 5000)
      @test tri_3.boundary_nodes == [1, 2, 3, 5000, 4, 1]
      @test DT.get_boundary_edge_map(tri_3, 3, 5000) == (tri_3.boundary_nodes, 3)
      @test DT.get_boundary_edge_map(tri_3, 5000, 4) == (tri_3.boundary_nodes, 4)
      @test_throws KeyError DT.get_boundary_edge_map(tri_3, 3, 4)
      @test !DT.contains_unoriented_edge((3, 4), DT.get_all_segments(tri_3))
      @test DT.contains_unoriented_edge((3, 5000), DT.get_all_segments(tri_3))
      @test DT.contains_unoriented_edge((5000, 4), DT.get_all_segments(tri_3))

      DT.split_boundary_edge!(tri_4, 6, 7, 1200)
      @test DT.get_boundary_edge_map(tri_4, 6, 1200) == (1, 6)
      @test DT.get_boundary_edge_map(tri_4, 1200, 7) == (1, 7)
      @test_throws KeyError DT.get_boundary_edge_map(tri_4, 6, 7)
      @test !DT.contains_unoriented_edge((6, 7), DT.get_all_segments(tri_4))
      @test DT.contains_unoriented_edge((6, 1200), DT.get_all_segments(tri_4))
      @test DT.contains_unoriented_edge((1200, 7), DT.get_all_segments(tri_4))
end

@testset "merge_boundary_edge!" begin
      x, y = complicated_geometry()
      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x, y)
      tri_1 = triangulate(points; rng, boundary_nodes)
      A = get_area(tri_1)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
      tri_2 = triangulate(points; rng, boundary_nodes)
      A = get_area(tri_2)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
      tri_3 = triangulate(points; rng, boundary_nodes)
      A = get_area(tri_3)

      rng = StableRNG(91818)
      boundary_nodes, points = convert_boundary_points_to_indices(reverse(reverse.(x[2])), reverse(reverse.(y[2])))
      tri_4 = triangulate(points; rng, boundary_nodes)
      A = get_area(tri_4)

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
      @test DT.contains_unoriented_edge((21, 22), DT.get_all_segments(tri_1))
      @test !DT.contains_unoriented_edge((21, 170), DT.get_all_segments(tri_1))
      @test !DT.contains_unoriented_edge((170, 22), DT.get_all_segments(tri_1))

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
      @test DT.contains_unoriented_edge((7, 8), DT.get_all_segments(tri_1))
      @test !DT.contains_unoriented_edge((7, 5000), DT.get_all_segments(tri_1))
      @test !DT.contains_unoriented_edge((5000, 8), DT.get_all_segments(tri_1))

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
      @test DT.contains_unoriented_edge((8, 9), DT.get_all_segments(tri_2))
      @test !DT.contains_unoriented_edge((8, 8182), DT.get_all_segments(tri_2))
      @test !DT.contains_unoriented_edge((8182, 9), DT.get_all_segments(tri_2))

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
      @test DT.contains_unoriented_edge((3, 4), DT.get_all_segments(tri_3))
      @test !DT.contains_unoriented_edge((3, 18289), DT.get_all_segments(tri_3))
      @test !DT.contains_unoriented_edge((18289, 4), DT.get_all_segments(tri_3))

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
      @test DT.contains_unoriented_edge((6, 7), DT.get_all_segments(tri_4))
      @test !DT.contains_unoriented_edge((6, 1200), DT.get_all_segments(tri_4))
      @test !DT.contains_unoriented_edge((1200, 7), DT.get_all_segments(tri_4))
end

@testset "get_adjacent concurrency" begin # Shouldn't be an issue anymore since we removed DefaultDict, but let's keep this here anyway. The test here is simply that it doesn't error.
      tri = triangulate(rand(2, 50), delete_ghosts=false)
      Base.Threads.@threads for _ in 1:5000
            get_adjacent(tri, -5, rand(1:1000))
      end
end

@testset "has_vertex and has_ghost_vertices" begin
      tri = triangulate(rand(2, 50), delete_ghosts=false)
      @test DT.has_vertex(tri, 1)
      @test !DT.has_vertex(tri, 57)
      @test DT.has_ghost_vertices(tri)
      @test DT.has_vertex(tri, -1)
      DT.delete_ghost_vertices_from_graph!(tri)
      @test !DT.has_vertex(tri, -1)
      @test !DT.has_ghost_vertices(tri)
end

@testset "Issue #70" begin
      points = [(-1.0, -1.0), (1.0, -1.0), (0.0, 1.0)]
      tri = triangulate(points)
      delete_ghost_triangles!(tri)
      DelaunayTriangulation.delete_ghost_vertices_from_graph!(tri)
      @test collect(each_solid_vertex(tri)) == collect(each_vertex(tri))
      @test !DelaunayTriangulation.has_ghost_vertices(tri)
      @test DelaunayTriangulation.num_ghost_vertices(tri) == 0
      @test DelaunayTriangulation.num_solid_vertices(tri) == 3
      @test isempty(collect(each_ghost_vertex(tri)))
end

@testset "Boundary curve orientation" begin
      tri = triangulate(rand(2, 500))
      @test DT.is_positively_oriented(tri, 1)
      lock_convex_hull!(tri)
      @test DT.is_positively_oriented(tri, 1)

      pts = [
            (-7.36, 12.55), (-9.32, 8.59), (-9.0, 3.0), (-6.32, -0.27),
            (-4.78, -1.53), (2.78, -1.41), (-5.42, 1.45), (7.86, 0.67),
            (10.92, 0.23), (9.9, 7.39), (8.14, 4.77), (13.4, 8.61),
            (7.4, 12.27), (2.2, 13.85), (-3.48, 10.21), (-4.56, 7.35),
            (3.44, 8.99), (3.74, 5.87), (-2.0, 8.0), (-2.52, 4.81),
            (1.34, 6.77), (1.24, 4.15)
      ]
      boundary_points = [
            (0.0, 0.0), (2.0, 1.0), (3.98, 2.85), (6.0, 5.0),
            (7.0, 7.0), (7.0, 9.0), (6.0, 11.0), (4.0, 12.0),
            (2.0, 12.0), (1.0, 11.0), (0.0, 9.13), (-1.0, 11.0),
            (-2.0, 12.0), (-4.0, 12.0), (-6.0, 11.0), (-7.0, 9.0),
            (-6.94, 7.13), (-6.0, 5.0), (-4.0, 3.0), (-2.0, 1.0), (0.0, 0.0)
      ]
      boundary_nodes, pts = convert_boundary_points_to_indices(boundary_points; existing_points=pts)
      tri = triangulate(pts; boundary_nodes, delete_ghosts=false)
      @test DT.is_positively_oriented(tri, 1)

      points = [
            (2.0, 8.0), (6.0, 4.0), (2.0, 6.0),
            (2.0, 4.0), (8.0, 2.0)
      ]
      segment_1 = [(0.0, 0.0), (14.0, 0.0)]
      segment_2 = [(14.0, 0.0), (10.0, 4.0), (4.0, 6.0), (2.0, 12.0), (0.0, 14.0)]
      segment_3 = [(0.0, 14.0), (0.0, 0.0)]
      boundary_points = [segment_1, segment_2, segment_3]
      boundary_nodes, points = convert_boundary_points_to_indices(boundary_points; existing_points=points)
      tri = triangulate(points; boundary_nodes)
      @test DT.is_positively_oriented(tri, 1)

      curve_1 = [[
            (0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (12.0, 0.0), (12.0, 4.0),
            (12.0, 8.0), (14.0, 10.0), (16.0, 12.0), (16.0, 16.0),
            (14.0, 18.0), (12.0, 20.0), (12.0, 24.0), (12.0, 28.0),
            (8.0, 28.0), (4.0, 28.0), (0.0, 28.0), (-2.0, 26.0), (0.0, 22.0),
            (0.0, 18.0), (0.0, 10.0), (0.0, 8.0), (0.0, 4.0), (-4.0, 4.0),
            (-4.0, 0.0), (0.0, 0.0),
      ]]
      curve_2 = [[
            (4.0, 26.0), (8.0, 26.0), (10.0, 26.0), (10.0, 24.0),
            (10.0, 22.0), (10.0, 20.0), (8.0, 20.0), (6.0, 20.0),
            (4.0, 20.0), (4.0, 22.0), (4.0, 24.0), (4.0, 26.0)
      ]]
      curve_3 = [[(4.0, 16.0), (12.0, 16.0), (12.0, 14.0), (4.0, 14.0), (4.0, 16.0)]]
      curve_4 = [[(4.0, 8.0), (10.0, 8.0), (8.0, 6.0), (6.0, 6.0), (4.0, 8.0)]]
      curves = [curve_1, curve_2, curve_3, curve_4]
      points = [
            (2.0, 26.0), (2.0, 24.0), (6.0, 24.0), (6.0, 22.0), (8.0, 24.0), (8.0, 22.0),
            (2.0, 22.0), (0.0, 26.0), (10.0, 18.0), (8.0, 18.0), (4.0, 18.0), (2.0, 16.0),
            (2.0, 12.0), (6.0, 12.0), (2.0, 8.0), (2.0, 4.0), (4.0, 2.0),
            (-2.0, 2.0), (4.0, 6.0), (10.0, 2.0), (10.0, 6.0), (8.0, 10.0), (4.0, 10.0),
            (10.0, 12.0), (12.0, 12.0), (14.0, 26.0), (16.0, 24.0), (18.0, 28.0),
            (16.0, 20.0), (18.0, 12.0), (16.0, 8.0), (14.0, 4.0), (14.0, -2.0),
            (6.0, -2.0), (2.0, -4.0), (-4.0, -2.0), (-2.0, 8.0), (-2.0, 16.0),
            (-4.0, 22.0), (-4.0, 26.0), (-2.0, 28.0), (6.0, 15.0), (7.0, 15.0),
            (8.0, 15.0), (9.0, 15.0), (10.0, 15.0), (6.2, 7.8),
            (5.6, 7.8), (5.6, 7.6), (5.6, 7.4), (6.2, 7.4), (6.0, 7.6),
            (7.0, 7.8), (7.0, 7.4)]
      boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points)
      tri = triangulate(points; boundary_nodes=boundary_nodes)
      @test DT.is_positively_oriented(tri, 1)
      @test !DT.is_positively_oriented(tri, 2)
      @test !DT.is_positively_oriented(tri, 3)
      @test !DT.is_positively_oriented(tri, 4)

      curve_1 = [
            [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0), (15.0, 0.0), (20.0, 0.0), (25.0, 0.0)],
            [(25.0, 0.0), (25.0, 5.0), (25.0, 10.0), (25.0, 15.0), (25.0, 20.0), (25.0, 25.0)],
            [(25.0, 25.0), (20.0, 25.0), (15.0, 25.0), (10.0, 25.0), (5.0, 25.0), (0.0, 25.0)],
            [(0.0, 25.0), (0.0, 20.0), (0.0, 15.0), (0.0, 10.0), (0.0, 5.0), (0.0, 0.0)]
      ] # outer-most boundary: counter-clockwise  
      curve_2 = [
            [(4.0, 6.0), (4.0, 14.0), (4.0, 20.0), (18.0, 20.0), (20.0, 20.0)],
            [(20.0, 20.0), (20.0, 16.0), (20.0, 12.0), (20.0, 8.0), (20.0, 4.0)],
            [(20.0, 4.0), (16.0, 4.0), (12.0, 4.0), (8.0, 4.0), (4.0, 4.0), (4.0, 6.0)]
      ] # inner boundary: clockwise 
      curve_3 = [
            [(12.906, 10.912), (16.0, 12.0), (16.16, 14.46), (16.29, 17.06),
            (13.13, 16.86), (8.92, 16.4), (8.8, 10.9), (12.906, 10.912)]
      ] # this is inside curve_2, so it's counter-clockwise 
      curves = [curve_1, curve_2, curve_3]
      points = [
            (3.0, 23.0), (9.0, 24.0), (9.2, 22.0), (14.8, 22.8), (16.0, 22.0),
            (23.0, 23.0), (22.6, 19.0), (23.8, 17.8), (22.0, 14.0), (22.0, 11.0),
            (24.0, 6.0), (23.0, 2.0), (19.0, 1.0), (16.0, 3.0), (10.0, 1.0), (11.0, 3.0),
            (6.0, 2.0), (6.2, 3.0), (2.0, 3.0), (2.6, 6.2), (2.0, 8.0), (2.0, 11.0),
            (5.0, 12.0), (2.0, 17.0), (3.0, 19.0), (6.0, 18.0), (6.5, 14.5),
            (13.0, 19.0), (13.0, 12.0), (16.0, 8.0), (9.8, 8.0), (7.5, 6.0),
            (12.0, 13.0), (19.0, 15.0)
      ]
      boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points)
      tri = triangulate(points; boundary_nodes=boundary_nodes)
      @test DT.is_positively_oriented(tri, 1)
      @test !DT.is_positively_oriented(tri, 2)
      @test DT.is_positively_oriented(tri, 3)

      θ = LinRange(0, 2π, 20) |> collect
      θ[end] = 0 # need to make sure that 2π gives the exact same coordinates as 0
      xy = Vector{Vector{Vector{NTuple{2,Float64}}}}()
      cx = 0.0
      for i in 1:2
            # Make the exterior circle
            push!(xy, [[(cx + cos(θ), sin(θ)) for θ in θ]])
            # Now the interior circle - clockwise
            push!(xy, [[(cx + 0.5cos(θ), 0.5sin(θ)) for θ in reverse(θ)]])
            cx += 3.0
      end
      boundary_nodes, points = convert_boundary_points_to_indices(xy)
      tri = triangulate(points; boundary_nodes=boundary_nodes)
      @test DT.is_positively_oriented(tri, 1)
      @test !DT.is_positively_oriented(tri, 2)
      @test DT.is_positively_oriented(tri, 3)
      @test !DT.is_positively_oriented(tri, 4)

      C = (15.7109521325776, 33.244486807457)
      D = (14.2705719699703, 32.8530791545746)
      E = (14.3, 27.2)
      F = (14.1, 27.0)
      G = (13.7, 27.2)
      H = (13.4, 27.5)
      I = (13.1, 27.6)
      J = (12.7, 27.4)
      K = (12.5, 27.1)
      L = (12.7, 26.7)
      M = (13.1, 26.5)
      N = (13.6, 26.4)
      O = (14.0, 26.4)
      P = (14.6, 26.5)
      Q = (15.1983491346581, 26.8128534095401)
      R = (15.6, 27.6)
      S = (15.6952958264624, 28.2344688505621)
      T = (17.8088971520274, 33.1192363585346)
      U = (16.3058917649589, 33.0722674401887)
      V = (16.3215480710742, 29.7374742376305)
      W = (16.3841732955354, 29.393035503094)
      Z = (16.6190178872649, 28.9233463196351)
      A1 = (17.0417381523779, 28.5319386667527)
      B1 = (17.5114273358368, 28.3753756055997)
      C1 = (18.1376795804487, 28.3597192994844)
      D1 = (18.7169629067146, 28.5632512789833)
      E1 = (19.2805899268653, 28.8920337074045)
      F1 = (19.26493362075, 28.4536571361762)
      G1 = (20.6426885588962, 28.4223445239456)
      H1 = (20.689657477242, 33.1035800524193)
      I1 = (19.2805899268653, 33.0722674401887)
      J1 = (19.2962462329806, 29.7531305437458)
      K1 = (19.0614016412512, 29.393035503094)
      L1 = (18.7482755189452, 29.236472441941)
      M1 = (18.4508057027546, 29.1425346052493)
      N1 = (18.1689921926793, 29.3147539725175)
      O1 = (17.7932408459121, 29.6278800948235)
      P1 = (22.6466957416542, 35.4207133574833)
      Q1 = (21.2219718851621, 34.9979930923702)
      R1 = (21.2376281912774, 28.4693134422915)
      S1 = (22.6780083538847, 28.4380008300609)
      T1 = (24.5724213938357, 33.1975178891111)
      U1 = (23.3512295168425, 32.8530791545746)
      V1 = (23.3199169046119, 28.4380008300609)
      W1 = (24.6663592305274, 28.3753756055997)
      Z1 = (15.1942940307729, 35.4363696635986)
      A2 = (14.7246048473139, 35.3737444391374)
      B2 = (14.3645098066621, 35.1858687657538)
      C2 = (14.1766341332786, 34.8570863373326)
      D2 = (14.1140089088174, 34.3247719294125)
      E2 = (14.2705719699703, 33.8394264398383)
      F2 = (14.7246048473139, 33.6202381542241)
      G2 = (15.4604512347329, 33.6045818481088)
      H2 = (16.0, 34.0)
      I2 = (15.9771093365377, 34.6848669700643)
      J2 = (15.6170142958859, 35.2328376840997)
      K2 = (24.1653574348379, 35.4520259697138)
      L2 = (23.7739497819555, 35.4363696635986)
      M2 = (23.4608236596496, 35.2641502963303)
      N2 = (23.272947986266, 34.9040552556785)
      O2 = (23.1320412312284, 34.5909291333725)
      P2 = (23.1163849251131, 34.2151777866054)
      Q2 = (23.2886042923813, 33.8081138276077)
      R2 = (23.8209187003014, 33.6045818481088)
      S2 = (24.3062641898756, 33.5576129297629)
      T2 = (24.7602970672192, 33.8550827459536)
      U2 = (25.010797965064, 34.4656786844502)
      V2 = (24.8385785977957, 34.9666804801397)
      W2 = (24.5254524754898, 35.2641502963303)
      Z2 = (25.3708930057158, 37.4716894585871)
      A3 = (24.7916096794498, 37.3464390096648)
      B3 = (24.4471709449133, 36.9550313567823)
      C3 = (24.3062641898756, 36.5636237038999)
      D3 = (24.4941398632592, 35.9999966837492)
      E3 = (25.0264542711793, 35.5929327247515)
      F3 = (25.5587686790994, 35.5929327247515)
      F3 = (25.5587686790994, 35.5929327247515)
      G3 = (26.0, 36.0)
      H3 = (26.1380520053653, 36.5792800100152)
      I3 = (26.0, 37.0)
      J3 = (25.7466443524829, 37.2838137852036)
      K3 = (26.3885529032101, 35.4676822758291)
      L3 = (25.9814889442124, 35.3580881330221)
      M3 = (25.6840191280217, 35.1858687657538)
      N3 = (25.5274560668688, 34.9040552556785)
      O3 = (25.4961434546382, 34.5596165211419)
      P3 = (25.5274560668688, 34.246490398836)
      Q3 = (25.6683628219064, 33.8394264398383)
      R3 = (26.0284578625583, 33.6358944603394)
      S3 = (26.5451159643631, 33.6202381542241)
      T3 = (27.0, 34.0)
      U3 = (27.280962351782, 34.5596165211419)
      V3 = (27.0304614539373, 35.2171813779844)
      W3 = (26.1693646175959, 33.087923746304)
      Z3 = (26.0, 33.0)
      A4 = (25.5274560668688, 32.7278287056522)
      B4 = (25.2612988629087, 32.4147025833463)
      C4 = (25.1830173323322, 32.0702638488098)
      D4 = (25.2299862506781, 31.7727940326191)
      E4 = (25.6527065157911, 31.5222931347744)
      F4 = (26.2946150665183, 31.7258251142732)
      G4 = (26.5607722704784, 32.5086404200381)
      H4 = (27.1557119028596, 32.7434850117675)
      I4 = (27.6097447802033, 32.4929841139228)
      J4 = (27.6410573924338, 32.1015764610403)
      K4 = (27.7193389230103, 31.6005746653509)
      L4 = (27.437525412935, 31.4283552980826)
      M4 = (26.9834925355914, 31.2561359308143)
      N4 = (26.5764285765937, 31.0995728696614)
      O4 = (26.0441141686736, 30.7864467473554)
      P4 = (25.6527065157911, 30.5672584617413)
      Q4 = (25.3239240873699, 30.1915071149741)
      R4 = (25.1673610262169, 29.8783809926682)
      S4 = (25.1047358017558, 29.6122237887082)
      T4 = (25.0890794956405, 29.1895035235952)
      U4 = (25.2926114751393, 28.8294084829433)
      V4 = (25.6840191280217, 28.5632512789833)
      W4 = (26.1537083114806, 28.3753756055997)
      Z4 = (26.8269294744384, 28.391031911715)
      A5 = (27.4844943312809, 28.6102201973292)
      B5 = (27.7342002330051, 28.7239579596219)
      C5 = (27.7264126450755, 28.4202565942047)
      D5 = (29.1825559185446, 28.3922538389457)
      E5 = (29.1545531632856, 32.2146299318021)
      F5 = (29.000538009361, 32.5786657501693)
      G5 = (28.6785063238822, 32.9006974356481)
      H5 = (28.3144705055149, 33.0827153448317)
      I5 = (27.9084305542591, 33.2367304987563)
      J5 = (27.3343740714492, 33.3207387645334)
      K5 = (26.8303244767868, 33.2367304987563)
      L5 = (27.6564057569279, 30.786489413592)
      M5 = (27.6984098898165, 30.3944508399657)
      N5 = (27.6984098898165, 29.7363860913787)
      O5 = (27.5863988687804, 29.4143544059)
      P5 = (27.2643671833016, 29.2043337414573)
      Q5 = (26.9843396307114, 29.1763309861983)
      R5 = (26.6903107004917, 29.3163447624934)
      S5 = (26.5782996794556, 29.7503874690082)
      T5 = (26.7603175886393, 30.3384453294476)
      U5 = (27.3203726938197, 30.7024811478149)
      J_curve = [[C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, C]]
      U_curve = [[T, U, V, W, Z, A1, B1, C1, D1, E1, F1, G1, H1, I1, J1, K1, L1, M1, N1, O1, T]]
      L_curve = [[P1, Q1, R1, S1, P1]]
      I_curve = [[T1, U1, V1, W1, T1]]
      A_curve_outline = [[
            K5, W3, Z3, A4, B4, C4, D4, E4, F4, G4, H4, I4, J4, K4, L4, M4, N4,
            O4, P4, Q4, R4, S4, T4, U4, V4, W4, Z4, A5, B5, C5, D5, E5, F5, G5,
            H5, I5, J5, K5]]
      A_curve_hole = [[L5, M5, N5, O5, P5, Q5, R5, S5, T5, U5, L5]]
      dot_1 = [[Z1, A2, B2, C2, D2, E2, F2, G2, H2, I2, J2, Z1]]
      dot_2 = [[Z2, A3, B3, C3, D3, E3, F3, G3, H3, I3, J3, Z2]]
      dot_3 = [[K2, L2, M2, N2, O2, P2, Q2, R2, S2, T2, U2, V2, W2, K2]]
      dot_4 = [[K3, L3, M3, N3, O3, P3, Q3, R3, S3, T3, U3, V3, K3]]
      curves = [J_curve, U_curve, L_curve, I_curve, A_curve_outline, A_curve_hole, dot_1, dot_2, dot_3, dot_4]
      nodes, points = convert_boundary_points_to_indices(curves)
      tri = triangulate(points; boundary_nodes=nodes)
      @test DT.is_positively_oriented(tri, 1)
      @test DT.is_positively_oriented(tri, 2)
      @test DT.is_positively_oriented(tri, 3)
      @test DT.is_positively_oriented(tri, 4)
      @test DT.is_positively_oriented(tri, 5)
      @test !DT.is_positively_oriented(tri, 6)
      @test DT.is_positively_oriented(tri, 7)
      @test DT.is_positively_oriented(tri, 8)
      @test DT.is_positively_oriented(tri, 9)
      @test DT.is_positively_oriented(tri, 10)
      @test DT.num_curves(tri) == 10
end

@testset "Full constructor" begin
      # Simple example
      tri1 = triangulate(rand(2, 50))
      tri2 = Triangulation(get_points(tri1), each_solid_triangle(tri1), get_convex_hull_vertices(tri1))
      unlock_convex_hull!(tri2)
      @test tri1 == tri2

      # With boundaries 
      points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.2, 0.2), (0.8, 0.2), (0.8, 0.8), (0.2, 0.8)]
      boundary_nodes = [[[1, 2], [2, 3], [3, 4], [4, 1]], [[8, 7], [7, 6], [6, 5], [5, 8]]]
      tri1 = triangulate(points; boundary_nodes)
      tri2 = Triangulation(get_points(tri1), each_solid_triangle(tri1), get_boundary_nodes(tri1))
      @test tri1 == tri2

      # Custom types 
      tri1 = triangulate(rand(2, 50); EdgeType=NTuple{2,Int32}, TrianglesType=Set{NTuple{3,Int}})
      tri2 = Triangulation(get_points(tri1), each_solid_triangle(tri1), get_convex_hull_vertices(tri1); EdgeType=NTuple{2,Int32}, TrianglesType=Set{NTuple{3,Int}})
      unlock_convex_hull!(tri2)
      @test tri1 == tri2
      @test tri1 ⊢ tri2

      # Delete ghosts works properly 
      tri1 = triangulate(rand(2, 50), delete_ghosts=true)
      tri2 = Triangulation(get_points(tri1), each_solid_triangle(tri1), get_convex_hull_vertices(tri1), delete_ghosts=true)
      unlock_convex_hull!(tri2)
      @test tri1 == tri2

      # Weights work properly 
      weights = rand(50)
      tri2 = Triangulation(get_points(tri1), each_solid_triangle(tri1), get_convex_hull_vertices(tri1), weights=weights)
      @test DT.is_weighted(tri2) && get_weights(tri2) == weights
end

@testset "Random sampling of triangles, edges, and vertices" begin
      # Get the triangulation
      A = (0.0, 0.0)
      B = (0.0, 25.0)
      C = (5.0, 25.0)
      D = (5.0, 5.0)
      E = (10.0, 5.0)
      F = (10.0, 10.0)
      G = (25.0, 10.0)
      H = (25.0, 15.0)
      I = (10.0, 15.0)
      J = (10.0, 25.0)
      K = (45.0, 25.0)
      L = (45.0, 20.0)
      M = (40.0, 20.0)
      N = (40.0, 5.0)
      O = (45.0, 5.0)
      P = (45.0, 0.0)
      Q = (10.0, 0.0)
      R = (10.0, -5.0)
      S = (15.0, -5.0)
      T = (15.0, -10.0)
      U = (10.0, -10.0)
      V = (5.0, -10.0)
      W = (5.0, -5.0)
      Z = (5.0, 0.0)
      A1 = (5.0, 2.5)
      B1 = (10.0, 2.5)
      C1 = (38.0, 2.5)
      D1 = (38.0, 20.0)
      E1 = (27.0, 20.0)
      F1 = (27.0, 11.0)
      G1 = (27.0, 4.0)
      H1 = (2.0, 4.0)
      I1 = (2.0, 0.0)
      pts = [A, I1, H1, G1, F1, E1, D1, C1, B1, A1, Z, W, V, U, T, S, R, Q, P, O, N, M, L, K, J, I, H, G, F, E, D, C, B, A]
      J1 = (17.0603265896789, 7.623652007194)
      K1 = (14.8552854162067, 6.5423337394336)
      L1 = (16.6998871670921, 6.9875824379232)
      M1 = (16.0, 6.0)
      N1 = (16.9755173137761, 6.6483453343121)
      O1 = (17.0391242707032, 4.8885528593294)
      P1 = (17.4207660122657, 6.4575244635308)
      Q1 = (17.6327892020226, 4.9945644542079)
      R1 = (22.6789411182379, 6.1818943168468)
      S1 = (21.8096460402344, 6.4787267825065)
      T1 = (26.0, 8.0)
      U1 = (15.0673086059636, 9.086612016517)
      W1 = (15.0, 8.5)
      Z1 = (17.7913089332764, 8.3603005983396)
      inner_pts = [Z1, W1, U1, T1, S1, R1, Q1, P1, O1, N1, M1, L1, K1, J1, Z1]
      boundary_pts = [[pts], [inner_pts]]
      nodes, points = convert_boundary_points_to_indices(boundary_pts)
      push!(points, (20.0, 20.0))
      rng = StableRNG(19191919)
      C = Set{NTuple{2,Int}}()
      for i in 1:50
            θ = 2π * rand(rng)
            r = 4sqrt(rand(rng))
            x = 20 + r * cos(θ)
            y = 20 + r * sin(θ)
            push!(points, (x, y))
            push!(C, (48, 48 + i))
      end
      tri = triangulate(points; boundary_nodes=nodes, segments=C, rng)

      # Vertices 
      ns = 5_000_000
      function test_vertex_sampler(tri, ns)
            solid_vertices = DefaultDict{Int,Float64,Float64}(0.0)
            ghost_vertices = DefaultDict{Int,Float64,Float64}(0.0)
            for _ in 1:ns
                  sv = rand(each_solid_vertex(tri))
                  gv = rand(each_ghost_vertex(tri))
                  solid_vertices[sv] += 1 / ns
                  ghost_vertices[gv] += 1 / ns
            end
            return solid_vertices, ghost_vertices
      end
      solid_vertices, ghost_vertices = test_vertex_sampler(tri, ns)
      @inferred rand(each_solid_vertex(tri))
      @inferred rand(each_ghost_vertex(tri))
      @test all(!DT.is_ghost_vertex, keys(solid_vertices))
      @test all(DT.is_ghost_vertex, keys(ghost_vertices))
      sm = mean(values(solid_vertices)) # sm for solid mean
      gm = mean(values(ghost_vertices))
      @test sm ≈ 1 / DT.num_solid_vertices(tri)
      @test gm ≈ 1 / 2
      @test all(y -> isapprox(y, sm, rtol=1e-1), values(solid_vertices))
      @test all(y -> isapprox(y, gm, rtol=1e-1), values(ghost_vertices))
      for i in 1:250
            for f in (each_vertex, each_solid_vertex, each_ghost_vertex)
                  rng = StableRNG(i)
                  V1 = rand(rng, f(tri))
                  V2 = rand(rng, f(tri))
                  V3 = rand(rng, f(tri))
                  rng = StableRNG(i)
                  @test rand(rng, f(tri)) == V1
                  @test rand(rng, f(tri)) == V2
                  @test rand(rng, f(tri)) == V3
            end
      end

      # Triangles 
      ns = 5_000_000
      function test_triangle_sampler(tri, ns)
            solid_triangles = DefaultDict{NTuple{3,Int},Float64,Float64}(0.0)
            ghost_triangles = DefaultDict{NTuple{3,Int},Float64,Float64}(0.0)
            for _ in 1:ns
                  st = rand(each_solid_triangle(tri))
                  gt = rand(each_ghost_triangle(tri))
                  solid_triangles[st] += 1 / ns
                  ghost_triangles[gt] += 1 / ns
            end
            return solid_triangles, ghost_triangles
      end
      solid_triangles, ghost_triangles = test_triangle_sampler(tri, ns)
      @inferred rand(each_solid_triangle(tri))
      @inferred rand(each_ghost_triangle(tri))
      @test all(!DT.is_ghost_triangle, keys(solid_triangles))
      @test all(DT.is_ghost_triangle, keys(ghost_triangles))
      sm = mean(values(solid_triangles)) # sm for solid mean
      gm = mean(values(ghost_triangles))
      @test sm ≈ 1 / DT.num_solid_triangles(tri)
      @test gm ≈ 1 / DT.num_ghost_triangles(tri)
      @test all(y -> isapprox(y, sm, rtol=1e-1), values(solid_triangles))
      @test all(y -> isapprox(y, gm, rtol=1e-1), values(ghost_triangles))
      for i in 1:250
            for f in (each_triangle, each_solid_triangle, each_ghost_triangle)
                  rng = StableRNG(i)
                  V1 = rand(rng, f(tri))
                  V2 = rand(rng, f(tri))
                  V3 = rand(rng, f(tri))
                  rng = StableRNG(i)
                  @test rand(rng, f(tri)) == V1
                  @test rand(rng, f(tri)) == V2
                  @test rand(rng, f(tri)) == V3
            end
      end

      # Edges 
      ns = 5_000_000
      function test_edge_sampler(tri, ns)
            solid_edges = DefaultDict{NTuple{2,Int},Float64,Float64}(0.0)
            ghost_edges = DefaultDict{NTuple{2,Int},Float64,Float64}(0.0)
            for _ in 1:ns
                  se = rand(each_solid_edge(tri))
                  ge = rand(each_ghost_edge(tri))
                  solid_edges[se] += 1 / ns
                  ghost_edges[ge] += 1 / ns
            end
            return solid_edges, ghost_edges
      end
      solid_edges, ghost_edges = test_edge_sampler(tri, ns)
      @inferred rand(each_solid_edge(tri))
      @inferred rand(each_ghost_edge(tri))
      @test all(!DT.is_ghost_edge, keys(solid_edges))
      @test all(DT.is_ghost_edge, keys(ghost_edges))
      sm = mean(values(solid_edges)) # sm for solid mean
      gm = mean(values(ghost_edges))
      @test sm ≈ 1 / DT.num_solid_edges(tri)
      @test gm ≈ 1 / DT.num_ghost_edges(tri)
      @test all(y -> isapprox(y, sm, rtol=1e-1), values(solid_edges))
      @test all(y -> isapprox(y, gm, rtol=1e-1), values(ghost_edges))
      for i in 1:250
            for f in (each_edge, each_solid_edge, each_ghost_edge)
                  rng = StableRNG(i)
                  V1 = rand(rng, f(tri))
                  V2 = rand(rng, f(tri))
                  rng = StableRNG(i)
                  @test rand(rng, f(tri)) == V1
                  @test rand(rng, f(tri)) == V2
            end
      end
end