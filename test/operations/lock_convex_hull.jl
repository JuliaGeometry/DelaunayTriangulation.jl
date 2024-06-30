using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using DataStructures

Base.include(@__MODULE__, "../helper_functions.jl")

@testset "lock_convex_hull" begin
    for _ in 1:150
        pts = rand(2, 500)
        tri = triangulate(pts)
        tri2 = deepcopy(tri)
        idx = get_convex_hull_vertices(tri)
        bn = get_boundary_nodes(tri)
        resize!(bn, length(idx))
        copyto!(get_boundary_nodes(tri), idx)
        bn_map = get_ghost_vertex_map(tri)
        bnn_map = get_boundary_edge_map(tri)
        IT = DT.integer_type(tri)
        bn_map[IT(DT.𝒢)] = bn
        ne = num_boundary_edges(bn)
        E = DT.edge_type(tri)
        for i in 1:ne
            u = get_boundary_nodes(bn, i)
            v = get_boundary_nodes(bn, i + 1)
            e = DT.construct_edge(E, u, v)
            bnn_map[e] = (bn, i)
        end
        for e in keys(bnn_map)
            add_segment!(tri, e)
        end
        lock_convex_hull!(tri2)
        tri3 = triangulate(pts; boundary_nodes=bn)
        @test tri2.boundary_nodes == tri3.boundary_nodes == bn
        @test tri2.ghost_vertex_map == tri3.ghost_vertex_map == bn_map
        @test tri2.boundary_edge_map == tri3.boundary_edge_map == bnn_map
        @test tri2.convex_hull == tri3.convex_hull == tri.convex_hull
        @test tri2.interior_segments == tri3.interior_segments == tri.interior_segments
        @test tri2.all_segments == tri3.all_segments == tri.all_segments
    end
    p1 = (0.0, 0.0)
    p2 = (1.0, 0.0)
    p3 = (1.0, 1.0)
    p4 = (0.0, 1.0)
    pts = [p1, p2, p3, p4]
    tri = triangulate(pts; boundary_nodes=[1, 2, 3, 4, 1])
    @test_throws ArgumentError("Cannot lock the convex hull of a triangulation with boundary nodes.") lock_convex_hull!(tri)
end

@testset "unlock_convex_hull" begin
    p1 = (0.0, 0.0)
    p2 = (1.0, 0.0)
    p3 = (1.0, 1.0)
    p4 = (0.0, 1.0)
    pts = [p1, p2, p3, p4]
    tri = triangulate(pts)
    add_segment!(tri, 1, 3)
    add_segment!(tri, 3, 1)
    @test tri.interior_segments ≠ Set(((1, 3), (3, 1)))
    @test tri.all_segments ≠ Set(((1, 3), (3, 1)))
    tri2 = deepcopy(tri)
    lock_convex_hull!(tri)
    @test tri.interior_segments == Set(((3, 1),)) || tri.interior_segments == Set(((1, 3),))
    unlock_convex_hull!(tri)
    @test tri.interior_segments == tri2.interior_segments
    @test tri.all_segments == tri2.all_segments
    @test tri.boundary_nodes == tri2.boundary_nodes
    @test tri.ghost_vertex_map == tri2.ghost_vertex_map
    @test tri.boundary_edge_map == tri2.boundary_edge_map
    @test tri.convex_hull == tri2.convex_hull
    @test_throws ArgumentError("Cannot unlock the convex hull of a triangulation without boundary nodes.") unlock_convex_hull!(tri2)
    lock_convex_hull!(tri)
    push!(tri.convex_hull.vertices, 17)
    tri.convex_hull.vertices[1] = 17 # need to be circular
    @test_throws ArgumentError("Cannot unlock the convex hull of a triangulation with boundary nodes that are not the convex hull. If the boundary nodes have been split, consider setting reconstruct=true.") unlock_convex_hull!(tri)
end

@testset "Fixing interior segments that happen to be on the convex hull" begin
    for _ in 1:10
        points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.25, 0.5), (0.75, 0.5)]
        tri = triangulate(points; segments=Set([(1, 2), (2, 3), (5, 6)]))
        lock_convex_hull!(tri)
        @test DT.compare_unoriented_edge_collections(get_interior_segments(tri), Set([(5, 6)]))
        @test DT.compare_unoriented_edge_collections(tri.cache.interior_segments_on_hull, Set([(1, 2), (2, 3)]))
        push!(points, (0.25, 0.0))
        push!(points, (0.75, 0.0))
        push!(points, (0.5, 0.0))
        DT.complete_split_edge_and_legalise!(tri, 1, 2, 7)
        DT.complete_split_edge_and_legalise!(tri, 7, 2, 8)
        DT.complete_split_edge_and_legalise!(tri, 7, 8, 9)
        @test validate_triangulation(tri)
        unlock_convex_hull!(tri; reconstruct=true)
        @test validate_triangulation(tri)
        @test DT.compare_unoriented_edge_collections(get_interior_segments(tri), Set([(5, 6), (1, 7), (7, 9), (9, 8), (8, 2), (2, 3)]))
        @test isempty(tri.cache.interior_segments_on_hull)
        @test !DT.has_boundary_nodes(tri)
        @test DT.circular_equality(get_convex_hull_vertices(tri), [1, 7, 9, 8, 2, 3, 4, 1])
        validate_statistics(tri)

        points = [(0.0, 0.0), (9.0, 0.0), (9.0, 7.0)]
        tri = triangulate(points; segments=Set([(1, 2), (1, 3)]))
        lock_convex_hull!(tri)
        orig_tri = deepcopy(tri)
        orig_points = copy(orig_tri.points)
        push!(points, (4.5, 0.0))
        DT.complete_split_edge_and_legalise!(tri, 1, 2, 4)
        @test validate_triangulation(tri)
        @test collect(get_point(tri, 4)) ≈ [4.5, 0.0]
        unlock_convex_hull!(tri; reconstruct=true)
        @test tri == triangulate([orig_points; get_point(tri, 4)]; segments=Set([(1, 4), (4, 2), (1, 3)]))
    end
end