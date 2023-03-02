using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using Test
using CairoMakie
using StableRNGs

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

include("../test_setup.jl")
include("../helper_functions.jl")

rng = StableRNG(9992881)
pts = rand(rng, 2, 50)
tri = triangulate(pts; rng)
_S = get_convex_hull_indices(tri)
points = get_points(tri)
rng2 = StableRNG(92871)

next, prev, k, S, shuffled_indices = DT.prepare_convex_triangulation_vectors(_S)
@test next == zeros(Int64, length(_S) - 1)
@test prev == zeros(Int64, length(_S) - 1)
@test k == length(_S) - 1
@test shuffled_indices == collect(1:(length(_S)-1))
@test S == _S[begin:(end-1)]

SS = DT.prepare_convex_triangulation_vectors(_S[begin:(end-1)])[4]
@test SS == _S[begin:(end-1)]

DT.reset_convex_triangulation_vectors!(next, prev, shuffled_indices, k, rng2)
@test next == [(2:length(S))..., 1]
@test prev == [length(S), 1:(length(S)-1)...]
@test shuffled_indices == shuffle(StableRNG(92871), collect(1:(length(_S)-1)))

u, v, w = DT.index_shuffled_linked_list(S, next, prev, shuffled_indices, 2)
@test u == S[shuffled_indices[2]]
@test v == S[next[shuffled_indices[2]]]
@test w == S[prev[shuffled_indices[2]]]

for _ in 1:100
    local pts, S, points
    pts = rand(2, 2500)
    tri_orig = triangulate(pts)
    S = get_convex_hull_indices(tri_orig)
    points = get_points(tri_orig)
    tri_bw = @views triangulate(points; skip_points=setdiff(each_point_index(tri_orig), S))
    tri_ch = triangulate_convex(points, S; add_ghost_triangles=false, compute_centers=false, add_convex_hull=false, delete_empty_features=false)
    DT.compare_triangle_collections(get_triangles(tri_bw), get_triangles(tri_ch))

    for (ij, k) in get_adjacent(tri_ch)
        if DT.edge_exists(k)
            @test get_adjacent(tri_bw, ij) == k
        end
    end
    for (w, E) in get_adjacent2vertex(tri_ch)
        for ij in each_edge(E)
            @test DT.contains_edge(ij, get_adjacent2vertex(tri_bw, w))
        end
    end
    for i in get_vertices(tri_ch)
        @test get_neighbours(tri_ch, i) == setdiff(get_neighbours(tri_bw, i), DT.BoundaryIndex)
    end
    @test isempty(get_convex_hull_indices(tri_ch))

    DT.convex_triangulation_post_processing!(tri_ch, S[begin:end-1], false, true, false, false)
    @test get_convex_hull(tri_ch) == get_convex_hull(tri_bw)
    empty!(get_convex_hull_indices(tri_ch))
    _pt = deepcopy(DT.RepresentativePointList)
    DT.convex_triangulation_post_processing!(tri_ch, S[begin:end-1], false, false, true, false)
    @test DT.RepresentativePointList[1].x ≈ _pt[1].x rtol = 1e-1
    @test DT.RepresentativePointList[1].y ≈ _pt[1].y rtol = 1e-1
    @test DT.RepresentativePointList[1].n ≈ _pt[1].n
    @test any(!DT.edge_exists, values(get_adjacent(get_adjacent(tri_ch))))
    DT.convex_triangulation_post_processing!(tri_ch, S[begin:end-1], false, false, false, true)
    @test all(DT.edge_exists, values(get_adjacent(get_adjacent(tri_ch))))
    DT.convex_triangulation_post_processing!(tri_ch, S[begin:end-1], true, false, false, false)
    DT.add_ghost_triangles!(tri_bw)
    @test DT.compare_triangle_collections(get_triangles(tri_ch), get_triangles(tri_bw))
    @test (get_adjacent ∘ get_adjacent)(tri_ch) == (get_adjacent ∘ get_adjacent)(tri_bw)
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_ch) == (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_bw)
    @test (get_graph ∘ get_graph)(tri_ch) == (get_graph ∘ get_graph)(tri_bw)

    tri_ch = triangulate_convex(points, S; add_ghost_triangles=true, delete_empty_features=true, add_convex_hull=false, compute_centers=false)
    DT.add_ghost_triangles!(tri_bw)
    @test DT.compare_triangle_collections(get_triangles(tri_ch), get_triangles(tri_bw))
    @test (get_adjacent ∘ get_adjacent)(tri_ch) == (get_adjacent ∘ get_adjacent)(tri_bw)
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_ch) == (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_bw)
    @test (get_graph ∘ get_graph)(tri_ch) == (get_graph ∘ get_graph)(tri_bw)
    @test !(get_convex_hull(tri_ch) == get_convex_hull(tri_bw))

    tri_ch = triangulate_convex(points, S)
    @test get_convex_hull(tri_ch) == get_convex_hull(tri_bw)
    @test DT.compare_triangle_collections(get_triangles(tri_ch), get_triangles(tri_bw))
    @test (get_adjacent ∘ get_adjacent)(tri_ch) == (get_adjacent ∘ get_adjacent)(tri_bw)
    @test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_ch) == (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_bw)
    @test (get_graph ∘ get_graph)(tri_ch) == (get_graph ∘ get_graph)(tri_bw)
end

p1 = [10.0, 12.0]
p2 = [7.0, 11.0]
p3 = [8.0, 6.0]
p4 = [10.0, 3.0]
p5 = [14.0, 5.0]
p6 = [15.0, 10.0]
p7 = [13.0, 12.0]
pts = [p1, p2, p3, p4, p5, p6, p7]
S = collect(1:7)
tri = triangulate_convex(pts, S)
fig, ax, sc = triplot(tri; plot_convex_hull=false)
SAVE_FIGURE && save("$save_path/convex_triangulation_example.png", fig)

p1 = [8.0, 4.0]
p2 = [10.0, 4.0]
p3 = [12.0, 4.0]
p4 = [14.0, 4.0]
p5 = [14.0, 6.0]
p6 = [14.0, 8.0]
p7 = [14.0, 10.0]
p8 = [12.0, 10.0]
p9 = [10.0, 10.0]
p10 = [8.0, 10.0]
p11 = [8.0, 8.0]
p12 = [8.0, 6.0]
pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12]
for _ in 1:10000
    tri_ch = triangulate_convex(pts, 1:12)
    validate_triangulation(tri_ch)
end

pts = rand(2, 50)
p1 = [0.0, 0.0]
p2 = [1.0, 0.0]
p3 = [0.0, 1.0]
pts[:, 11] .= p1
pts[:, 27] .= p2
pts[:, 5] .= p3
S = [11, 27, 5, 11]
tri_ch = triangulate_convex(pts, S)
tri_bw = triangulate(pts; skip_points=setdiff(1:50, [11, 27, 5]), delete_ghosts=false)
@test get_convex_hull(tri_ch) == get_convex_hull(tri_bw)
@test DT.compare_triangle_collections(get_triangles(tri_ch), get_triangles(tri_bw))
@test (get_adjacent ∘ get_adjacent)(tri_ch) == (get_adjacent ∘ get_adjacent)(tri_bw)
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_ch) == (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_bw)
@test (get_graph ∘ get_graph)(tri_ch) == (get_graph ∘ get_graph)(tri_bw)

pts[:, 28] .= [1.0, 1.0]
S = [11, 27, 28, 5, 11]
tri_ch = triangulate_convex(pts, S)
tri_bw = triangulate(pts; skip_points=setdiff(1:50, [11, 27, 5, 28]), delete_ghosts=false)
@test get_convex_hull(tri_ch) == get_convex_hull(tri_bw)
@test DT.compare_triangle_collections(get_triangles(tri_ch), get_triangles(tri_bw))
@test (get_adjacent ∘ get_adjacent)(tri_ch) == (get_adjacent ∘ get_adjacent)(tri_bw)
@test (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_ch) == (get_adjacent2vertex ∘ get_adjacent2vertex)(tri_bw)
@test (get_graph ∘ get_graph)(tri_ch) == (get_graph ∘ get_graph)(tri_bw)