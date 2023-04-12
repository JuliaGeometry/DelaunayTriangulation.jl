using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra

include("../helper_functions.jl")

@testset "Defining the targets" begin
    targets1 = DT.RefinementTargets(;
        max_area=10.0,
        max_radius_edge_ratio=0.9,
        max_points=10_000)
    @test targets1.max_area == 10.0
    @test targets1.max_radius_edge_ratio == 0.9
    @test targets1.max_points == 10_000

    targets2 = DT.RefinementTargets(;
        max_area=10.0)
    @test targets2.max_area == 10.0
    @test targets2.max_radius_edge_ratio == 1.0
    @test targets2.max_points == typemax(Int64)

    targets3 = DT.RefinementTargets(;
        max_radius_edge_ratio=0.5)
    @test targets3.max_area == typemax(Float64)
    @test targets3.max_radius_edge_ratio == sqrt(3) / 3
    @test targets3.max_points == typemax(Int64)

    targets4 = DT.RefinementTargets(;
        max_points=10_000,
        min_angle=50.0)
    @test targets4.max_area == typemax(Float64)
    @test targets4.max_radius_edge_ratio == 1 / (2sind(50))
    @test 50.0 ≈ asind(1 / (2targets4.max_radius_edge_ratio))
    @test targets4.max_points == 10_000
end

@testset "Quality comparisons" begin
    targets1 = DT.RefinementTargets(;
        max_area=10.0,
        min_angle=23.7,
        max_points=10_000)
    ρ = 1.0 / (2sind(23.7))
    a = 10.0
    targets2 = DT.RefinementTargets(;
        max_area=(T, p, q, r, A) -> A > 10.0,
        max_radius_edge_ratio=(T, p, q, r, ρ) -> ρ > 1.24,
        max_points=10_000)
    !DT.compare_area(targets1, (1, 2, 3), 0.5, 0.0, 0.0, 0.0)
    DT.compare_area(targets1, (1, 2, 3), 12.0, 0.0, 0.0, 0.0)
    DT.compare_area(targets1, (1, 2, 3), 10.5, 0.0, 0.0, 0.0)
    !DT.compare_area(targets2, (1, 2, 3), 0.5, 0.0, 0.0, 0.0)
    !DT.compare_area(targets2, (1, 2, 3), 7.8, 0.0, 0.0, 0.0)
    DT.compare_area(targets2, (1, 2, 3), 10.5, 0.0, 0.0, 0.0)
    DT.compare_ratio(targets1, (1, 2, 3), 0.5, 0.0, 0.0, 0.0)
    DT.compare_ratio(targets1, (1, 2, 3), 1.5, 0.0, 0.0, 0.0)
    !DT.compare_ratio(targets1, (1, 2, 3), 2.5, 0.0, 0.0, 0.0)
    !DT.compare_ratio(targets1, (1, 2, 3), 3.5, 0.0, 0.0, 0.0)
    DT.compare_ratio(targets2, (1, 2, 3), 0.5, 0.0, 0.0, 0.0)
    DT.compare_ratio(targets2, (1, 2, 3), 1.5, 0.0, 0.0, 0.0)
    !DT.compare_ratio(targets2, (1, 2, 3), 2.5, 0.0, 0.0, 0.0)
    !DT.compare_ratio(targets2, (1, 2, 3), 3.5, 0.0, 0.0, 0.0)
    @test !DT.compare_points(targets1, 5000)
    @test DT.compare_points(targets1, 10_001)
    @test !DT.compare_points(targets2, 5000)
    @test DT.compare_points(targets2, 10_001)
end

@testset "Assessing the quality of a triangle" begin
    tri = poor_triangulation_example()
    max_area = 50.0
    min_angle = 20.0
    targets = DT.RefinementTargets(;
        max_area,
        min_angle
    )
    stats = statistics(tri)
    for T in each_triangle(tri)
        A = DT.get_area(stats, T)
        ρ = DT.get_radius_edge_ratio(stats, T)
        r, flag = DT.assess_triangle_quality(tri, T, targets)
        @test r ≈ ρ
        @test flag == (A > max_area || ρ > 1 / (2sind(min_angle)))
    end
end

@testset "Building the priority queues" begin
    @testset "Unconstrained triangulation" begin
        tri = poor_triangulation_example()
        max_area = 50.0
        min_angle = 20.0
        targets = DT.RefinementTargets(;
            max_area,
            min_angle
        )
        stats = statistics(tri)
        queue = DT.initialise_refinement_queue(tri, targets)
        @test DT.encroachment_queue_is_empty(queue)
        all_encroachment = Any[]
        all_triangle = Any[]
        while !DT.encroachment_queue_is_empty(queue)
            e = DT.encroachment_dequeue!(queue)
            u, v = DT.edge_indices(e)
            p, q = get_point(tri, u, v)
            ℓ² = norm(p .- q)^2
            push!(all_encroachment, (e, ℓ²))
        end
        while !DT.triangle_queue_is_empty(queue)
            T = DT.triangle_dequeue!(queue)
            ρ = DT.get_radius_edge_ratio(stats, T)
            push!(all_triangle, (T, ρ))
        end
        @test issorted(last.(all_encroachment), rev=true)
        @test issorted(last.(all_triangle), rev=true)
        for (e, ℓ²) in all_encroachment
            u, v = DT.edge_indices(e)
            p, q = get_point(tri, u, v)
            ℓ²′ = norm(p .- q)^2
            @test ℓ² ≈ ℓ²′
        end
        for (T, ρ) in all_triangle
            ρ′ = DT.get_radius_edge_ratio(stats, T)
            @test ρ ≈ ρ′
        end
        @test DT.isempty(queue.encroachment_queue)
        @test DT.isempty(queue.triangle_queue)
        @test DT.isempty(queue)
    end

    @testset "Constrained triangulation" begin
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
        uncons_tri = triangulate(points)
        cons_tri = triangulate(points; boundary_nodes=boundary_nodes)
        tri = deepcopy(cons_tri)
        max_area = 50.0
        min_angle = 20.0
        targets = DT.RefinementTargets(;
            max_area,
            min_angle
        )
        stats = statistics(tri)
        queue = DT.initialise_refinement_queue(tri, targets)
        @test !DT.encroachment_queue_is_empty(queue)
        all_encroachment = Any[]
        all_triangle = Any[]
        while !DT.encroachment_queue_is_empty(queue)
            ℓ² = peek(queue.encroachment_queue)[2]
            e = DT.encroachment_dequeue!(queue)
            @test DT.is_encroached(tri, e)
            push!(all_encroachment, (e, ℓ²))
        end
        while !DT.triangle_queue_is_empty(queue)
            ρ = peek(queue.triangle_queue)[2]
            T = DT.triangle_dequeue!(queue)
            push!(all_triangle, (T, ρ))
        end
        @test issorted(last.(all_encroachment), rev=true)
        @test issorted(last.(all_triangle), rev=true)
        for (e, ℓ²) in all_encroachment
            u, v = DT.edge_indices(e)
            p, q = get_point(tri, u, v)
            ℓ²′ = norm(p .- q)^2
            @test ℓ² ≈ ℓ²′
        end
        for (T, ρ) in all_triangle
            ρ′ = DT.get_radius_edge_ratio(stats, T)
            @test ρ ≈ ρ′
            A = DT.get_area(stats, T)
            u, v, w = indices(T)
            p, q, r = get_point(tri, u, v, w)
            A′ = DT.triangle_area(p, q, r)
            ρ2 = DT.triangle_radius_edge_ratio(p, q, r)
            @test A ≈ A′
            @test ρ ≈ ρ2
            @test A ≥ targets.max_area || ρ ≥ targets.max_radius_edge_ratio
        end
        @test DT.isempty(queue.encroachment_queue)
        @test DT.isempty(queue.triangle_queue)
    end
end
