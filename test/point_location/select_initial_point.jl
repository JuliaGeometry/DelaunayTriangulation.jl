using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Random

include("../helper_functions.jl")

tri, label_map, index_map = simple_geometry()
add_ghost_triangles!(tri)
DT.compute_representative_points!(tri)
pts = get_points(tri)

@test DT.default_num_samples(37) == 4
@test DT.default_num_samples(582) == 9
@inferred DT.default_num_samples(582)

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
for i in each_point_index(pts)
    global current_dist, current_idx
    if i ≠ 5
        current_dist, current_idx = DT.compare_distance(current_dist, current_idx, pts, i,
                                                        qx, qy)
    end
end
cd, ci = findmin([norm(p .- (qx, qy))^2 + (p == pts[5]) * Inf] for p in pts)
@test cd[1] ≈ current_dist
@test ci == current_idx

Random.seed!(92884881)
m = DT.default_num_samples(length(pts))
tried_idx = [rand(eachindex(pts)) for _ in 1:m]
q = (0.59183, -6.2731)
ci = tried_idx[argmin([norm(pts[i] .- q) for i in tried_idx])]
Random.seed!(92884881)
ci2 = DT.select_initial_point(pts, q)
@test ci == ci2

seeds = rand(1:50000000, 50)
for s in seeds
    local m, tried_idx, q, ci, ci2
    Random.seed!(s)
    m = 9
    tried_idx = [rand(eachindex(pts)) for _ in 1:m]
    q = (20rand(), 20rand())
    ci = tried_idx[argmin([norm(pts[i] .- q) for i in tried_idx])]
    Random.seed!(s)
    ci2 = DT.select_initial_point(pts, q; m)
    @inferred DT.select_initial_point(pts, q; m)
    @test ci == ci2
end

seeds = rand(1:500, 50)
for s in seeds
    local m, tried_idx, q, ci, ci2
    Random.seed!(s)
    m = 3
    tried_idx = [rand(eachindex(pts)) for _ in 1:m]
    try_points = (1, 2, 3, 10, 11, 12, 13, 14, 15, 16)
    q = (20rand(), 20rand())
    ci = tried_idx[argmin([norm(pts[i] .- q) for i in tried_idx])]
    ci2 = try_points[argmin([norm(pts[i] .- q) for i in try_points])]
    d1 = norm(pts[ci] .- q)
    d2 = norm(pts[ci2] .- q)
    _ci = d1 < d2 ? ci : ci2
    Random.seed!(s)
    ci2 = DT.select_initial_point(pts, q; m, try_points)
    @test _ci == ci2
end

seeds = rand(1:527700, 50)
for s in seeds
    local m, tried_idx, q, ci, ci2
    Random.seed!(s)
    m = 9
    tried_idx = [rand([1, 5, 8, 19, 20]) for _ in 1:m]
    q = (20rand(), 20rand())
    ci = tried_idx[argmin([norm(pts[i] .- q) for i in tried_idx])]
    Random.seed!(s)
    ci2 = DT.select_initial_point(pts, q; m, point_indices=(1, 5, 8, 19, 20))
    @inferred DT.select_initial_point(pts, q; m, point_indices=(1, 5, 8, 19, 20))
    @test ci == ci2
    @test ci2 ∈ (1, 5, 8, 19, 20)
end

for k in each_point_index(pts)
    @test DT.select_initial_point(pts, k; try_points=k) == k
end
