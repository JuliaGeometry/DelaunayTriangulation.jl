using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using StableRNGs
using ElasticArrays
using DelimitedFiles

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

        tri = triangulate(pts)
        stats = refine!(tri, min_angle=25.0, max_area=0.0001)
        @test count(rad2deg.(DT.get_minimum_angle.(Ref(stats), each_solid_triangle(tri))) .≥ 25.0) / length(each_solid_triangle(tri)) ≥ 0.99
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
    max_area = 0.001A
    min_area = 1e-9A
    stats = refine!(tri; min_area, max_area, exterior_curve_index=Set((1, 3)))
    validate_statistics(tri, stats)
    @test validate_triangulation(tri; check_planarity=false, check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
end

@testset "Small angles" begin
    ps = 0
    for i in 1:50
        @show i
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
        stats = refine!(tri; min_angle=27.3, min_area=0.0)
        ps += DT.get_largest_angle(stats) ≤ max(π - 2 * deg2rad(17.0), 2asin((sqrt(3) - 1) / 2)) # Corollary 8 of "When and Why Ruppert's Algorithm Works. In Twelfth International Meshing Roundtable, pp. 91–102, Santa Fe, NM, Sept 2003."
        validate_statistics(tri)
        @test validate_triangulation(tri)
    end
    @test ps / 50 ≥ 0.88
end

@testset "Another multiply-connected domain" begin
    t = (collect ∘ LinRange)(0, 2π, 50)
    t[end] = 0
    r = 2 .* (1 .+ cos.(t))
    x = r .* cos.(t)
    y = r .* sin.(t)
    points = [(x, y) for (x, y) in zip(x, y)]
    x = -(x .+ 3)
    inverse_points = [(x, y) for (x, y) in zip(x, y)]
    reverse!(inverse_points)
    curve_points = [[points], [inverse_points]]
    boundary_nodes, points = convert_boundary_points_to_indices(curve_points)
    C = Set([[(50, i) for i in 98:-1:81]..., [(50, i) for i in 69:-1:51]...,
        [(10, i) for i in 36:49]...])
    tri = triangulate(points; boundary_nodes, check_arguments=false, edges=C)
    A = DT.get_total_area(tri)
    refine!(tri, max_area=1e-3A, min_area=0.0, exterior_curve_index=Set((1, 2)))
    @test validate_triangulation(tri, check_planarity=false, check_ghost_triangle_delaunay=false, check_ghost_triangle_orientation=false)
    validate_statistics(tri)
end

@testset "Tight example with a single triangle boundary interior" begin
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
    J1 = (14.0, 6.0)
    K1 = (14.0, 8.0)
    L1 = (22.0, 8.0)
    inner_pts = [J1, K1, L1, J1]
    boundary_pts = [[pts], [inner_pts]]
    nodes, points = convert_boundary_points_to_indices(boundary_pts)
    tri = triangulate(points; boundary_nodes=nodes)
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    A = get_total_area(tri)
    refine!(tri; max_area=1e-3A)
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    validate_statistics(tri)
end

@testset "A complicated example" begin
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
    C = Set{NTuple{2,Int64}}()
    for i in 1:50
        θ = 2π * rand()
        r = 4sqrt(rand())
        x = 20 + r * cos(θ)
        y = 20 + r * sin(θ)
        push!(points, (x, y))
        push!(C, (48, 48 + i))
    end
    tri = triangulate(points; boundary_nodes=nodes, edges=C)
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    A = get_total_area(tri)
    refine!(tri; max_area=1.0)
    @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
    validate_statistics(tri)
end

if !get(ENV, "CI", false)
    @testset "Tasmania" begin
        tassy = readdlm("./test/tassy.txt")
        ymax = @views maximum(tassy[:, 2])
        tassy = [(x, ymax - y) for (x, y) in eachrow(tassy)]
        reverse!(tassy)
        unique!(tassy)
        push!(tassy, tassy[begin])
        boundary_nodes, points = convert_boundary_points_to_indices(tassy)
        tri = triangulate(points; boundary_nodes=boundary_nodes)
        A = get_total_area(tri)
        stats = refine!(tri; max_area=1e-3A)
        @test validate_triangulation(tri; check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
        validate_statistics(tri)
    end
end