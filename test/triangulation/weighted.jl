using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using StableRNGs
using LinearAlgebra
using ..DelaunayTriangulation: add_weight!, get_weight, get_weights

@testset "ZeroWeight" begin
    zw = DT.ZeroWeight()
    @inferred DT.ZeroWeight()
    @test get_weight(zw, 1) ⊢ 0.0
end

@testset "get_weight" begin
    weights = rand(10)
    @test get_weight(weights, 1) == weights[1]
    @test get_weight(weights, 5) == weights[5]
    tri = Triangulation(rand(2, 10); weights)
    @test get_weight(tri, 1) == weights[1]
    @test get_weight(tri, 5) == weights[5]
    @test DT.get_weights(tri) == weights
    weights = rand(Float32, 10)
    @test get_weight(weights, 2) ⊢ Float64(weights[2])
end

@testset "add_weight!" begin
    weights = rand(10)
    add_weight!(weights, 0.5)
    @test weights[11] == 1 / 2 && length(weights) == 11
    tri = Triangulation(rand(2, 10); weights)
    add_weight!(tri, 27.5)
    @test get_weight(tri, 12) == 27.5 && length(DT.get_weights(tri)) == 12
end

@testset "is_weighted" begin
    tri = Triangulation(rand(2, 10))
    @test !DT.is_weighted(tri)
    tri = Triangulation(rand(2, 10); weights = rand(10))
    @test DT.is_weighted(tri)
    tri = Triangulation(rand(2, 10); weights = DT.ZeroWeight())
    @test !DT.is_weighted(tri)
    tri = Triangulation(rand(2, 10); weights = zeros(10))
    @test DT.is_weighted(tri)
end

@testset "get_lifted_point" begin
    for _ in 1:100
        p, w = rand(2), rand()
        @test DT.get_lifted_point(p, w) == (p..., sum(p .^ 2) - w)
    end
    weights = randn(1000)
    points = randn(2, 1000)
    tri = Triangulation(points; weights)
    for i in 1:1000
        @test DT.get_lifted_point(tri, i) == (points[:, i]..., sum(points[:, i] .^ 2) .- weights[i])
        @inferred DT.get_lifted_point(tri, i)
    end
    points = randn(Float32, 2, 1000)
    tri = Triangulation(points; weights)
    for i in 1:1000
        @test DT.get_lifted_point(tri, i) == (Float32.(points[:, i])..., sum(Float32.(points[:, i]) .^ 2) .- weights[i])
        @inferred DT.get_lifted_point(tri, i)
    end
end

@testset "get_power_distance" begin
    for fi in 1:10
        tri, submerged, nonsubmerged, weights = get_weighted_example(fi)
        for i in eachindex(weights)
            for j in eachindex(weights)
                pᵢ = get_point(tri, i)
                pⱼ = get_point(tri, j)
                wᵢ = weights[i]
                wⱼ = DT.get_weight(tri, j)
                dist = norm(pᵢ .- pⱼ) .^ 2 - wᵢ - wⱼ
                @test dist ≈ DT.get_power_distance(tri, i, j)
            end
        end
    end
end

@testset "get_weighted_nearest_neighbour" begin
    tri = triangulate(rand(2, 50))
    for i in 1:50
        for _ in 1:50
            j = DT.get_weighted_nearest_neighbour(tri, i)
            @test j == i == get_nearest_power_point(tri, i)
        end
    end

    for fi in 1:NUM_WEGT
        tri, submerged, nonsubmerged, weights = get_weighted_example(fi)
        for i in eachindex(weights)
            j = DT.get_weighted_nearest_neighbour(tri, i)
            @inferred DT.get_weighted_nearest_neighbour(tri, i)
            if i ∈ submerged
                @test j == get_nearest_power_point(tri, i)
            else
                @test j == i
            end
        end
    end
end

@testset "get_distance_to_witness_plane" begin
    # solid triangles
    points = [(-2.6, -2.82), (0.0, 0.0), (2.11, 4.42), (1.0, 1.0)]
    _points = [Tuple(rand(2)) for _ in 1:500]
    _points[137] = points[1]
    _points[58] = points[2]
    _points[498] = points[3]
    _points[5] = points[4]
    weights = [3.7, -20.0, 9.2, 2.5]
    _weights = rand(500)
    _weights[137] = weights[1]
    _weights[58] = weights[2]
    _weights[498] = weights[3]
    _weights[5] = weights[4]
    tri = Triangulation(_points; weights = _weights)
    d = DT.get_distance_to_witness_plane(tri, 5, (137, 58, 498))
    @test d ≈ -2.129523129725314
    @test DT.get_distance_to_witness_plane(tri, 5, (137, 58, 498)) ≈
        DT.get_distance_to_witness_plane(tri, 5, (58, 137, 498)) ≈
        DT.get_distance_to_witness_plane(tri, 5, (498, 58, 137)) ≈
        DT.get_distance_to_witness_plane(tri, 5, (137, 498, 58)) ≈
        DT.get_distance_to_witness_plane(tri, 5, (58, 498, 137)) ≈
        DT.get_distance_to_witness_plane(tri, 5, (498, 137, 58))
    _points[5] = (-3.9, 2.01)
    _weights[5] = -15.0
    d = DT.get_distance_to_witness_plane(tri, 5, (137, 58, 498))
    @test d ≈ 5.603226942686893

    # ghost triangles: exterior points
    a, b, c = (0.0, -6.0), (6.0, 2.0), (-7.0, 7.0)
    points = [a, b, c]
    weights = zeros(4)
    tri = triangulate(points; weights)
    d = (3.18, 7.62)
    push!(points, d)
    points[4] = (7.21, 0.5)
    @test DT.get_distance_to_witness_plane(tri, 4, (2, 1, -1)) == -Inf
    points[1] = (0.0, -4.0)
    points[2] = (6.0, 3.0)
    points[3] = (-7.0, 3.0)
    points[4] = (6.19, -1.52)
    @test DT.get_distance_to_witness_plane(tri, 4, (2, 1, -1)) == -Inf
    @test DT.get_distance_to_witness_plane(tri, 4, (1, 3, -1)) == Inf
    @test DT.get_distance_to_witness_plane(tri, 4, (3, 2, -1)) == Inf
    points[4] = (-6.0, 0.0)
    @test DT.get_distance_to_witness_plane(tri, 4, (2, 1, -1)) == Inf
    @test DT.get_distance_to_witness_plane(tri, 4, (1, 3, -1)) == -Inf
    @test DT.get_distance_to_witness_plane(tri, 4, (3, 2, -1)) == Inf
    points[4] = (0.0, 6.0)
    @test DT.get_distance_to_witness_plane(tri, 4, (2, 1, -1)) == Inf
    @test DT.get_distance_to_witness_plane(tri, 4, (1, 3, -1)) == Inf
    @test DT.get_distance_to_witness_plane(tri, 4, (3, 2, -1)) == -Inf

    # ghost triangles: points on the solid unweighted edge
    points[4] = (0.0, 3.0)
    @test DT.get_distance_to_witness_plane(tri, 4, (3, 2, -1)) == DT.get_distance_to_witness_plane(tri, 4, (3, 1, 2))
    points[1] = (2.0, -4.0)
    points[2] = (2.0, 3.0)
    points[3] = (-5.0, 5.0)
    points[4] = (2.0, 0.0)
    @test DT.get_distance_to_witness_plane(tri, 4, (2, 1, -1)) == DT.get_distance_to_witness_plane(tri, 4, (3, 1, 2))
    @test DT.get_distance_to_witness_plane(tri, 4, (1, 3, -1)) == Inf

    # ghost triangles: points on the solid weighted edge
    points[3] = (-2.0, 0.0)
    points[4] = (0.0, -2.0)
    cert = DT.point_position_relative_to_circumcircle(tri, 1, 3, -1, 4)
    @test DT.is_on(cert)
    @test DT.get_distance_to_witness_plane(tri, 4, (1, 3, -1)) == DT.get_distance_to_witness_plane(tri, 4, (3, 1, 2))
    weights .= [2.9, 3.7, -4.4, 0.0]
    cert = DT.point_position_relative_to_circumcircle(tri, 1, 3, -1, 4)
    @test DT.is_on(cert)
    @test DT.get_distance_to_witness_plane(tri, 4, (1, 3, -1)) == DT.get_distance_to_witness_plane(tri, 4, (3, 1, 2))
    weights[4] = -16.0
    cert = DT.point_position_relative_to_circumcircle(tri, 1, 3, -1, 4)
    @test DT.get_distance_to_witness_plane(tri, 4, (1, 3, -1)) == Inf
end

@testset "is_submerged" begin
    for fi in 1:NUM_WEGT
        tri, submerged, nonsubmerged, weights = get_weighted_example(fi)
        for i in submerged
            @test DT.is_submerged(tri, i)
        end
        for i in nonsubmerged
            @test !DT.is_submerged(tri, i)
        end
        for i in 1:DT.num_points(tri)
            @test (DT.is_submerged(tri, i, find_triangle(tri, get_point(tri, i)))) == (i ∈ submerged)
        end
    end
end

@testset "Weighted triangulations with identical/no weights and random sets/collinear sets" begin
    @testset "Triangulation with zero weights" begin
        for n in 3:500
            points = rand(2, n)
            weights = zeros(size(points, 2))
            tri1 = triangulate(points)
            tri2 = triangulate(points; weights)
            @test tri1 == tri2
            @test validate_triangulation(tri2)
        end
        for i in 3:10
            for j in 3:10
                tri = triangulate_rectangle(0, 10, 0, 10, i, j)
                tri = triangulate(get_points(tri); weights = zeros(i * j))
                # @test validate_triangulation(tri) # Why is this failing sometimes? Is validate not branching at weighted triangulations?
            end
        end
    end

    @testset "Triangulation with identical weights" begin
        for n in 3:500
            points = rand(2, n)
            w = randn()
            weights = w * ones(size(points, 2))
            tri1 = triangulate(points; weights)
            tri2 = triangulate(points; weights)
            @test tri1 == tri2
        end
        for i in 3:10
            for j in 3:10
                tri = triangulate_rectangle(0, 10, 0, 10, i, j)
                tri = triangulate(get_points(tri); weights = 10randn() * ones(i * j))
                # @test validate_triangulation(tri) # Why is this failing sometimes? Is validate not branching at weighted triangulations?
            end
        end
    end
end

@testset "Brute force search" begin
    tri = triangulate(rand(2, 50))
    for i in 1:50
        V = DT.brute_force_search_enclosing_circumcircle(tri, i)
        @test !DT.is_outside(DT.point_position_relative_to_circumcircle(tri, V, i))
    end

    for fi in 1:NUM_WEGT
        tri, submerged, nonsubmerged, weights = get_weighted_example(fi)
        for i in eachindex(weights)
            V = DT.brute_force_search_enclosing_circumcircle(tri, i)
            if i ∈ submerged
                @test V == (0, 0, 0)
            else
                @test !DT.is_outside(DT.point_position_relative_to_circumcircle(tri, V, i))
            end
        end
    end
end

@testset "Convex polygons" begin
    for fi in 1:NUM_CWEGT
        @info "Testing triangulation of a weighted convex polygon: $fi"
        tri, submerged, nonsubmerged, weights, S = get_convex_polygon_weighted_example(fi)
        ctri = triangulate_convex(get_points(tri), S; weights)
        @test tri == ctri
        fi < 77 && @test validate_triangulation(ctri) # THIS GETS STUCK IN A LOOP FOR FI == 77 ?!
    end
end

# While we have tested the convex polygons and everything else properly,
# normal triangulations still need to be correctly updated. To test this,
# we will need to test the ones provided in helper_functions.jl properly.
# See the get_weighted_example function in that file and the comments surrounding it.
