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
    tri = poor_triangulation_example()
    max_area = 50.0
    min_angle = 20.0
    targets = DT.RefinementTargets(;
        max_area,
        min_angle
    )
    stats = statistics(tri)
    queue = DT.initialise_refinement_queue(tri, targets)
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
end