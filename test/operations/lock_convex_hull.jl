using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using DataStructures

include("../test_setup.jl")
save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

include("../helper_functions.jl")

@testset "lock_convex_hull" begin
    for _ in 1:150
        pts = rand(2, 500)
        tri = triangulate(pts)
        tri2 = deepcopy(tri)
        idx = get_convex_hull_indices(tri)
        bn = get_boundary_nodes(tri)
        resize!(bn, length(idx))
        copyto!(get_boundary_nodes(tri), idx)
        bn_map = get_boundary_map(tri)
        bnn_map = get_boundary_edge_map(tri)
        IT = DT.integer_type(tri)
        bn_map[IT(DT.BoundaryIndex)] = bn
        ne = num_boundary_edges(bn)
        E = DT.edge_type(tri)
        for i in 1:ne
            u = get_boundary_nodes(bn, i)
            v = get_boundary_nodes(bn, i + 1)
            e = DT.construct_edge(E, u, v)
            bnn_map[e] = (bn, i)
        end
        for e in keys(bnn_map)
            add_edge!(tri, e)
        end
        lock_convex_hull!(tri2)
        tri3 = triangulate(pts; boundary_nodes=bn)
        @test tri2.boundary_nodes == tri3.boundary_nodes == bn
        @test tri2.boundary_map == tri3.boundary_map == bn_map
        @test tri2.boundary_edge_map == tri3.boundary_edge_map == bnn_map
        @test tri2.convex_hull == tri3.convex_hull == tri.convex_hull
        @test tri2.constrained_edges == tri3.constrained_edges == tri.constrained_edges
        @test tri2.all_constrained_edges == tri3.all_constrained_edges == tri.all_constrained_edges
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
    add_edge!(tri, 1, 3)
    add_edge!(tri, 3, 1)
    @test tri.constrained_edges ≠ Set(((1, 3), (3, 1)))
    @test tri.all_constrained_edges ≠ Set(((1, 3), (3, 1)))
    tri2 = deepcopy(tri)
    lock_convex_hull!(tri)
    @test tri.constrained_edges == Set(((3, 1),)) || tri.constrained_edges == Set(((1, 3),))
    unlock_convex_hull!(tri)
    @test tri.constrained_edges == tri2.constrained_edges
    @test tri.all_constrained_edges == tri2.all_constrained_edges
    @test tri.boundary_nodes == tri2.boundary_nodes
    @test tri.boundary_map == tri2.boundary_map
    @test tri.boundary_edge_map == tri2.boundary_edge_map
    @test tri.convex_hull == tri2.convex_hull
    @test_throws ArgumentError("Cannot unlock the convex hull of a triangulation without boundary nodes.") unlock_convex_hull!(tri2)
    lock_convex_hull!(tri)
    push!(tri.convex_hull.indices, 17)
    tri.convex_hull.indices[1] = 17 # need to be circular
    @test_throws ArgumentError("Cannot unlock the convex hull of a triangulation with boundary nodes that are not the convex hull.") unlock_convex_hull!(tri)
end