using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using StableRNGs
using ElasticArrays

include("../helper_functions.jl")

@testset "Simple refinements" begin
    for i in 1:10
        @show i
        p1 = (0.0, 0.0)
        p2 = (1.0, 0.0)
        p3 = (0.0, 1.0)
        p4 = (1.0, 1.0)
        pts = [p1, p2, p3, p4]
        tri = triangulate(pts; delete_ghosts=false)
        refine!(tri; max_area=0.0001, maxiters=25_000)
        stats = statistics(tri)
        @test DT.get_smallest_angle(stats) ≥ deg2rad(30.0)
        @test DT.get_largest_area(stats) ≤ 0.0001
        @test !DT.is_constrained(tri)
        @test DT.convex_hull(tri).indices == DT.convex_hull(tri.points).indices
        @test validate_triangulation(tri)
        validate_statistics(tri)

        p1 = (0.0, 0.0)
        p2 = (1.0, 0.0)
        p3 = (0.0, 1.0)
        p4 = (1.0, 1.0)
        pts = [p1, p2, p3, p4]
        tri = triangulate(pts; delete_ghosts=false)
        stats = refine!(tri; max_area=0.0001, max_points=5_000, maxiters=25_000)
        @test num_points(tri) == 5001
        @test DT.get_smallest_angle(stats) ≥ deg2rad(30.0)
        @test DT.get_largest_area(stats) ≤ 0.001
        @test !DT.is_constrained(tri)
        @test DT.convex_hull(tri).indices == DT.convex_hull(tri.points).indices
        @test validate_triangulation(tri)
        validate_statistics(tri, stats)

        tri = triangulate(ElasticMatrix(rand(2, 500)))
        stats = refine!(tri, min_angle=25.0, max_area=0.0001)
        @test DT.get_smallest_angle(stats) ≥ deg2rad(25.0)
        @test DT.get_largest_area(stats) ≤ 0.0001
        @test !DT.is_constrained(tri)
        @test DT.convex_hull(tri).indices == DT.convex_hull(tri.points).indices
        validate_statistics(tri)
    end
end


@testset "Avoiding infinite bouncing with concentric circular shells" begin
    for i in 1:20
        @show i
        p1 = (0.0, 0.0)
        p2 = (1.0, 0.0)
        p3 = (0.0, 1.0)
        p4 = (1.0, 1.0)
        pts = [p1, p2, p3, p4]
        tri = triangulate(pts; delete_ghosts=false)
        add_edge!(tri, 1, 4)
        stats = refine!(tri; max_area=0.001, min_angle=27.3)
        @test DT.get_smallest_angle(stats) ≥ deg2rad(27.3)
        @test DT.get_largest_area(stats) ≤ 0.001
        @test DT.is_constrained(tri)
        @test DT.convex_hull(tri).indices == DT.convex_hull(tri.points).indices
        validate_statistics(tri, stats)
        @test validate_triangulation(tri)
    end
end

@testset "Triangulating with an interior hole" begin
    p1 = (0.0, 0.0)
    p2 = (1.0, 0.0)
    p3 = (1.0, 1.0)
    p4 = (0.0, 1.0)
    p5 = (0.4, 0.4)
    p6 = (0.6, 0.4)
    p7 = (0.6, 0.6)
    p8 = (0.4, 0.6)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8]
    boundary_nodes = [[[1, 2, 3, 4, 1]], [[8, 7, 6, 5, 8]]]
    tri = triangulate(pts; boundary_nodes=boundary_nodes, delete_ghosts=false)
    add_point!(tri, 0.1, 0.8)
    add_point!(tri, 0.3, 0.2)
    add_point!(tri, 0.7, 0.2)
    add_point!(tri, 0.9, 0.8)
    add_edge!(tri, 9, 10)
    add_edge!(tri, 11, 12)
    add_edge!(tri, 9, 12)
    add_edge!(tri, 10, 11)
    stats = refine!(tri; max_area=0.001, max_points=5000)
    @test DT.get_smallest_angle(stats) ≥ deg2rad(30)
    @test DT.get_largest_area(stats) ≤ 0.001
    @test DT.is_constrained(tri)
    @test DT.convex_hull(tri).indices == DT.convex_hull(tri.points).indices
    validate_statistics(tri, stats)
    @test validate_triangulation(tri)
end

@testset "Refining disjoint sets" begin
    θ = LinRange(0, 2π, 30) |> collect
    θ[end] = 0
    xy = Vector{Vector{Vector{NTuple{2,Float64}}}}()
    cx = 0.0
    for i in 1:2
        push!(xy, [[(cx + cos(θ), sin(θ)) for θ in θ]])
        push!(xy, [[(cx + 0.5cos(θ), 0.5sin(θ)) for θ in reverse(θ)]])
        cx += 3.0
    end
    boundary_nodes, points = convert_boundary_points_to_indices(xy)
    tri = triangulate(points; boundary_nodes=boundary_nodes, check_arguments=false)
    A = DT.get_total_area(tri)
    max_area = 0.00001A
    min_area = 1e-9A
    stats = refine!(tri; min_area, max_area, exterior_curve_index=Set((1, 3)))
    validate_statistics(tri, stats)
    @test validate_triangulation(tri; check_planarity=false, check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
end

REFS = [0]

p1 = (0.0, 0.0)
p2 = (1.0, 0.0)
p3 = (0.0, 1.0)
p4 = (1.0, 1.0)
p5 = (0.5, 0.5)
pts = [p1, p2, p3, p4, p5]
C = Set{NTuple{2,Int64}}()
for i in 1:20
    θ = 2π * rand()
    r = 0.5sqrt(rand())
    x = 0.5 + r * cos(θ)
    y = 0.5 + r * sin(θ)
    push!(pts, (x, y))
    push!(C, (5, 5 + i))
end
tri = triangulate(pts; delete_ghosts=false, boundary_nodes=[1, 2, 4, 3, 1], edges=C)
stats = refine!(tri; min_angle=27.3, min_area=0.0, maxiters=15000)

p1 = (0.0, 0.0)
p2 = (2.0, 0.0)
p3 = (2.0, 0.2)
p4 = (0.6464432231964, 0.0646443223196)
p5 = (0.6496673989241, 0.0)
pts = [p1, p2, p3, p4, p5]
A = DT.get_total_area(tri)
tri = triangulate(pts; boundary_nodes=[1, 5, 2, 3, 4, 1])
refine!(tri; max_area=0.001A, min_area=1e-9A)

tri, queue, events, targets, has_ghosts, subsegment_list = DT.initialise_refine(tri)

e = DT.encroachment_dequeue!(queue)
DT.split_subsegment!(tri, queue, events, targets, subsegment_list, e)
DT.split_all_encroached_segments!(tri, queue, events, targets, subsegment_list)


u, v, w = 1, 5, 4
smallest_idx = 2
ratio_flag = true
if smallest_idx == 1
    i, j = u, v
elseif smallest_idx == 2
    i, j = v, w
else
    i, j = w, u
end
shared_flag, common_vertex = DT.edge_lies_on_two_distinct_segments(tri, i, j)
p, q, r = get_point(tri, i, j, common_vertex)
ρ = DT.IndividualTriangleStatistics(p, q, r).radius_edge_ratio
px, py = getxy(p)
qx, qy = getxy(q)
rx, ry = getxy(r)
ℓrp² = (px - rx)^2 + (py - ry)^2
ℓrq² = (qx - rx)^2 + (qy - ry)^2

DT.edge_is_seditious(tri, u, v, w, smallest_idx, ratio_flag, ρ)
DT.assess_triangle_quality(tri, (1, 5, 4), targets)
