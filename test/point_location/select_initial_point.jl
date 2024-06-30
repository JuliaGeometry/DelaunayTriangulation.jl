using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Random
using StableRNGs

include("../helper_functions.jl")

global tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)
DT.compute_representative_points!(tri)
global pts = get_points(tri)

@testset "Getting the default number of samples" begin
    @test DT.default_num_samples(37) == 4
    @test DT.default_num_samples(582) == 9
    @inferred DT.default_num_samples(582)
end

@testset "Getting the closest point" begin
    i1 = 13
    i2 = 24
    qx, qy = pts[5]
    current_dist = norm((qx, qy) .- pts[i1])^2
    current_idx = i1
    next_dist = norm((qx, qy) .- pts[i2])^2
    cd, ci = DT.compare_distance(current_dist, current_idx, pts, i2, qx, qy)
    @test cd ≈ next_dist
    @test ci == i2
    cd, ci = DT.compare_distance(cd, ci, pts, i1, qx, qy)
    @inferred DT.compare_distance(cd, ci, pts, i1, qx, qy)
    @test cd ≈ next_dist
    @test ci == i2
    current_dist = Inf
    current_idx = 1
    for i in DT.each_point_index(pts)
        if i ≠ 5
            current_dist, current_idx = DT.compare_distance(current_dist, current_idx, tri, i,
                qx, qy)
            @inferred DT.compare_distance(current_dist, current_idx, tri, i,
                qx, qy)
            @inferred DT.compare_distance(current_dist, current_idx, pts, i,
                qx, qy)
        end
    end
    cd, ci = findmin([[norm(p .- (qx, qy))^2 + (p == pts[5]) * Inf] for p in pts])
    @test cd[1] ≈ current_dist
    @test ci == current_idx
end

@testset "Selecting random initial points" begin
    @testset "Specific example" begin
        rng = StableRNG(92884881)
        m = DT.default_num_samples(length(each_solid_vertex(tri)))
        tried_idx = [rand(rng, each_solid_vertex(tri)) for _ in 1:m]
        q = (0.59183, -6.2731)
        ci = tried_idx[argmin([norm(pts[i] .- q) for i in tried_idx])]
        rng = StableRNG(92884881)
        ci2 = DT.select_initial_point(tri, q; rng)
        @test ci == ci2
    end

    @testset "Random points" begin
        seeds = rand(1:50000000, 50)
        for s in seeds
            local m, tried_idx, q, ci, ci2
            rng = StableRNG(s)
            m = 9
            tried_idx = [rand(rng, each_solid_vertex(tri)) for _ in 1:m]
            q = (20rand(rng), 20rand(rng))
            ci = tried_idx[argmin([norm(pts[i] .- q) for i in tried_idx])]
            rng = StableRNG(s)
            ci2 = DT.select_initial_point(tri, q; rng, m)
            @inferred DT.select_initial_point(tri, q; rng, m)
            @test ci == ci2
        end

        seeds = rand(1:527700, 50)
        for s in seeds
            rng = StableRNG(s)
            m = 9
            tried_idx = [rand(rng, [1, 5, 8, 19, 20]) for _ in 1:m]
            q = (20rand(rng), 20rand(rng))
            ci = tried_idx[argmin([norm(pts[i] .- q) for i in tried_idx])]
            rng = StableRNG(s)
            ci2 = DT.select_initial_point(tri, q; m, rng, point_indices=(1, 5, 8, 19, 20))
            @inferred DT.select_initial_point(tri, q; m, rng, point_indices=(1, 5, 8, 19, 20))
            @test ci == ci2
            @test ci2 ∈ (1, 5, 8, 19, 20)
        end

        for k in DT.each_point_index(pts)
            @test DT.select_initial_point(tri, k; try_points=k) == k
        end
    end
end
