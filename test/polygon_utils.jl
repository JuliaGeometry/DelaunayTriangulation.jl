using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Random
using CairoMakie

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

include("./helper_functions.jl")

tri, label_map, index_map = simple_geometry()
pts = get_points(tri)

boundary_nodes = [index_map["a"],
      index_map["b"],
      index_map["c"],
      index_map["d"],
      index_map["e"],
      index_map["f"],
      index_map["g"],
      index_map["h"],
      index_map["a"]]
a, (cx, cy) = polygon_features(pts, boundary_nodes)
@inferred polygon_features(pts, boundary_nodes)
@test a ≈ 400.0 && cx ≈ 10.0 && cy ≈ 10.0
boundary_nodes = [index_map["j"],
      index_map["k"],
      index_map["ℓ"],
      index_map["i"],
      index_map["j"]]
a, (cx, cy) = polygon_features(pts, boundary_nodes)
@test a ≈ 40.0 && cx ≈ 6.0 && cy ≈ 11.0
boundary_nodes = [[index_map["a"],
            index_map["b"],
            index_map["c"],
            index_map["d"]],
      [index_map["d"],
            index_map["e"],
            index_map["f"],
            index_map["g"]],
      [index_map["g"],
            index_map["h"],
            index_map["a"]]]
a, (cx, cy) = polygon_features(pts, boundary_nodes)
@inferred polygon_features(pts, boundary_nodes)
@test a ≈ 400.0 && cx ≈ 10.0 && cy ≈ 10.0
boundary_nodes = [[index_map["j"],
            index_map["k"]],
      [index_map["k"],
            index_map["ℓ"],
            index_map["i"]],
      [index_map["i"],
            index_map["j"]]]
a, (cx, cy) = polygon_features(pts, boundary_nodes)
@test a ≈ 40.0 && cx ≈ 6.0 && cy ≈ 11.0
a, (cx, cy) = polygon_features(pts, tri.boundary_nodes)
@inferred polygon_features(pts, tri.boundary_nodes)
a1, a2, a3 = 400.0, 32.0, 40.0
c1, c2, c3 = (10.0, 10.0), (15.58333333333, 7.0), (6.0, 11.0)
@test a ≈ a1 - a2 - a3
@test cx ≈ (c1[1] * a1 - c2[1] * a2 - c3[1] * a3) / (a1 - a2 - a3)
@test cy ≈ (c1[2] * a1 - c2[2] * a2 - c3[2] * a3) / (a1 - a2 - a3)

p1 = [0.0, 0.0]
p2 = [10.0, 0.0]
q = [0.0, 5.0]
d = DT.squared_distance_to_segment(p1..., p2..., q...)
@test sqrt(d) == q[2]
@inferred DT.squared_distance_to_segment(p1..., p2..., q...)
for _ in 1:10000
      local q
      pᵢ = 10randn(2)
      pⱼ = 10randn(2)
      q = 10randn(2)
      pᵢx, pᵢy = pᵢ
      pⱼx, pⱼy = pⱼ
      qx, qy = q
      t = ((pᵢx - pⱼx) * (pᵢx - qx) + (pᵢy - pⱼy) * (pᵢy - qy)) /
          ((pᵢx - pⱼx)^2 + (pᵢy - pⱼy)^2) # solve (d/dt)||q - (pᵢ + t(pⱼ - pᵢ))|| = 0 for t
      if t < 0
            @test DT.squared_distance_to_segment(pᵢ..., pⱼ..., q...) ≈ norm(q - pᵢ)^2
      elseif t > 1
            @test DT.squared_distance_to_segment(pᵢ..., pⱼ..., q...) ≈ norm(q - pⱼ)^2
      else
            @test DT.squared_distance_to_segment(pᵢ..., pⱼ..., q...) ≈
                  norm(q - (pᵢ + t * (pⱼ - pᵢ)))^2
      end
end

boundary_nodes = [index_map["a"],
      index_map["b"],
      index_map["c"],
      index_map["d"],
      index_map["e"],
      index_map["f"],
      index_map["g"],
      index_map["h"],
      index_map["a"]]
q = (28.0, 18.0)
dist = DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@test dist ≈ -8.0
@inferred DT.distance_to_polygon(q, pts, boundary_nodes)
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (26.9116654358588, 23.0120339025522)
dist = DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@test dist ≈ -7.5394606788132
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (9.5687897994641, 11.7840765682329)
dist = DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@test dist ≈ 20.0 - q[2]
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes)
@test DT.distance_to_polygon(q, pts, tri.boundary_nodes) ≈ q[1] - 8.0
q = (2.0, 2.0)
dist = DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@test dist ≈ 2.0
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (-2.0, 0.0)
dist = DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@test dist ≈ -2.0
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (0.0, 0.0)
dist = DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@test dist ≈ 0.0
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (10.0, 0.0)
dist = DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@inferred DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@test dist ≈ 0.0
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (4.6998638334488, 13.8273575177129)
dist = DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@test dist ≈ q[1]
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes)
@test DT.distance_to_polygon(q, pts, tri.boundary_nodes) ≈ -(q[1] - 4.0)
q = [-0.181375963606, 9.9696497047896]
dist = DT.distance_to_polygon_single_segment(q, pts, boundary_nodes)
@test dist ≈ q[1]
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
boundary_nodes = [[index_map["a"],
            index_map["b"],
            index_map["c"],
            index_map["d"]],
      [index_map["d"],
            index_map["e"],
            index_map["f"],
            index_map["g"]],
      [index_map["g"],
            index_map["h"],
            index_map["a"]]]
q = (28.0, 18.0)
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ -8.0
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (26.9116654358588, 23.0120339025522)
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ -7.5394606788132
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (9.5687897994641, 11.7840765682329)
@inferred DT.distance_to_polygon(q, pts, boundary_nodes)
@inferred DT.distance_to_polygon(q, pts, tri.boundary_nodes) ==
          DT.distance_to_polygon(q, pts, tri.boundary_nodes)
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ 20.0 - q[2]
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes)
q = (2.0, 2.0)
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ 2.0
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (-2.0, 0.0)
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ -2.0
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (0.0, 0.0)
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ 0.0
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (10.0, 0.0)
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@inferred DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ 0.0
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes) ==
      DT.distance_to_polygon(q, pts, tri.boundary_nodes)
q = (4.6998638334488, 13.8273575177129)
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ q[1]
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes)
q = [-0.181375963606, 9.9696497047896]
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ q[1]
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes)
q = [14.2840577995064, 9.7320720329939]
dist = DT.distance_to_polygon_multiple_segments(q, pts, boundary_nodes)
@test dist ≈ 20.0 - q[1]
@test dist == DT.distance_to_polygon(q, pts, boundary_nodes)
@test DT.distance_to_polygon(q, pts, tri.boundary_nodes) ≈ -(q[1] - 14.0)

boundary_nodes = [index_map["a"],
      index_map["b"],
      index_map["c"],
      index_map["d"],
      index_map["e"],
      index_map["f"],
      index_map["g"],
      index_map["h"],
      index_map["a"]]
@test DT.polygon_bounds(pts, boundary_nodes) == (0.0, 20.0, 0.0, 20.0)
@inferred DT.polygon_bounds(pts, boundary_nodes)
boundary_nodes = [[index_map["a"],
            index_map["b"],
            index_map["c"],
            index_map["d"]],
      [index_map["d"],
            index_map["e"],
            index_map["f"],
            index_map["g"]],
      [index_map["g"],
            index_map["h"],
            index_map["a"]]]
@test DT.polygon_bounds(pts, boundary_nodes) == (0.0, 20.0, 0.0, 20.0)
@inferred DT.polygon_bounds(pts, boundary_nodes)
@test DT.polygon_bounds(pts, tri.boundary_nodes) == (0.0, 20.0, 0.0, 20.0)
@inferred DT.polygon_bounds(pts, tri.boundary_nodes)

p = DT.Cell(10.0, 10.0, 10.0, pts, tri.boundary_nodes)
@test p.dist == DT.distance_to_polygon((10.0, 10.0), pts, tri.boundary_nodes)
@test p.max_dist ≈ p.dist + 10.0 * sqrt(2)
q = DT.Cell(7.0, 17.0, 0.3, pts, tri.boundary_nodes)
@test (p < q) == (p.max_dist < q.max_dist)
@test (p > q) == (p.max_dist > q.max_dist)
@test (p == q) == (p.max_dist == q.max_dist)
@test (p ≤ q) == (p.max_dist ≤ q.max_dist)
@test (p ≥ q) == (p.max_dist ≥ q.max_dist)
@test p.x == 10.0
@test p.y == 10.0
@test p.half_width == 10.0
@test q.x == 7.0
@test q.y == 17.0
@test q.half_width == 0.3
@inferred p < q

q = DT.CellQueue{Float64}()
@test isempty(q)
@test q.queue.o == Base.Order.Reverse
c = DT.Cell(polygon_features(pts, boundary_nodes)[2]..., 10.0, pts, boundary_nodes)
DT.insert_cell!(q, c)
@test q.queue[c] == c.max_dist
@inferred DT.get_next_cell!(q)

pts = [(0.0, 10.0), (0.0, 8.0), (0.0, 6.0), (0.0, 4.0), (0.0, 2.0),
      (0.0, 0.0), (2.0, 0.0), (4.0, 0.0), (6.0, 0.0), (8.0, 0.0), (10.0, 0.0),
      (10.0, 2.0), (10.0, 4.0), (8.0, 4.0), (6.0, 4.0),
      (4.0, 4.0), (4.0, 6.0), (4.0, 8.0), (4.0, 10.0), (2.0, 10.0),
      (0.0, 10.0)]
boundary_nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1]
x, y = DT.pole_of_inaccessibility(pts, boundary_nodes)
@inferred DT.pole_of_inaccessibility(pts, boundary_nodes)
@test x == 2.5 && y == 2.5

push!(pts,
      (2.4988557436664, 2.4749992804628),
      (1.258244761794, 3.494679539536),
      (1.2242554198249, 1.4553190213896),
      (3.3825786348632, 1.2683776405595),
      (3.3825786348632, 3.4097061846133),
      (2.0, 4.0))
boundary_nodes = [[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            1]],
      [reverse([21, 22, 23, 24, 25, 21])]]
x, y = DT.pole_of_inaccessibility(pts, boundary_nodes)
@inferred DT.pole_of_inaccessibility(pts, boundary_nodes)
@test x ≈ 8.125 && y ≈ 1.875

push!(pts,
      (7.4103156582024, 2.4749992804628),
      (7.0, 1.0),
      (8.8548626918894, 0.8775002079148),
      (9.2797294665033, 2.2370738866791))
boundary_nodes = [[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            1]],
      [reverse([21, 22, 23, 24, 25, 21])],
      [reverse([26, 27, 28, 29, 26])]]
x, y = DT.pole_of_inaccessibility(pts, boundary_nodes)
@test x ≈ 2.5 && y ≈ 7.5

pts = [0.0 8.0
      2.0 5.0
      3.0 7.0
      1.8190689042843 8.1342247183191
      3.2296265960022 8.864995570655
      4.2493068550754 7.7433472856744
      4.5042269198437 5.8739334773735
      3.6714880416006 4.3784024307328
      2.7367811374502 2.6279513193238
      5.5069125079324 1.3873403374514
      8.4299959172756 2.7469140162157
      9.7045962411171 5.5340400576824
      8.565953285152 7.7943312986281
      6.7135341478357 9.0349422805005
      4.1303441581835 9.663745106929
      2.7537758084347 10.3775212882802
      1.0882980519485 10.4964839851721
      -1.1380038470281 9.8336918167745
      -2.2596521320086 8.4571234670257
      -2.7864869325298 5.9419121613117
      -1.3929239117964 3.647631578397
      0.3235378576436 4.9732159151922
      -0.9000784532443 6.6216990006939
      0.9863300260411 9.6807397779135
      0.153591147798 9.5447824100371
      0.2725538446899 8.6610595188403
      2.9067278472957 8.1852087312728
      2.1249729820062 9.4258197131452]'
A = Dict('a':'z' .=> 1:26)
boundary_nodes = [[[A['a'], A['d'], A['c'], A['b']],
            [A['b'], A['i'], A['j'], A['k'], A['h'], A['g'], A['l']],
            [A['l'], A['f'], A['m'], A['e'], A['n'], A['o'], A['p'], A['q'], A['p']],
            [A['p'], A['q'], A['r'], A['s'], A['t'], A['u'], A['v'], A['w'], A['a']]],
      [[28 - 2, 27 - 2, 26 - 2],
            [26 - 2, 30 - 2, 29 - 2, 28 - 2]]]
x, y = DT.pole_of_inaccessibility(pts, boundary_nodes)

@test x ≈ 4.6146922812432685 && y ≈ 3.0953047713990314
len = size(pts, 2)
pts = hcat(pts, [7.274358290326 3.0 5.3369657980869; 2.7978980291693 4.0 1.8801857960034])
boundary_nodes = [[[A['a'], A['d'], A['c'], A['b']],
            [A['b'], A['i'], A['j'], A['k'], A['h'], A['g'], A['l']],
            [A['l'], A['f'], A['m'], A['e'], A['n'], A['o'], A['p'], A['q'], A['p']],
            [A['p'], A['q'], A['r'], A['s'], A['t'], A['u'], A['v'], A['w'], A['a']]],
      [[28 - 2, 27 - 2, 26 - 2],
            [26 - 2, 30 - 2, 29 - 2, 28 - 2]],
      [[len + 1, len + 2, len + 3, len + 1]]]
x, y = DT.pole_of_inaccessibility(pts, boundary_nodes)
@test x ≈ -1.0785224985822 && y ≈ 5.3725906833292

# Making a figure 
pts = [0.0 8.0
      2.0 5.0
      3.0 7.0
      1.81907 8.13422
      3.22963 8.865
      4.24931 7.74335
      4.50423 5.87393
      3.67149 4.3784
      2.73678 2.62795
      5.50691 1.38734
      8.43 2.74691
      9.7046 5.53404
      8.56595 7.79433
      6.71353 9.03494
      4.13034 9.66375
      2.75378 10.3775
      1.0883 10.4965
      -1.138 9.83369
      -2.25965 8.45712
      -2.78649 5.94191
      -1.39292 3.64763
      0.323538 4.97322
      -0.900078 6.6217
      0.98633 9.68074
      0.153591 9.54478
      0.272554 8.66106
      2.90673 8.18521
      2.12497 9.42582
      7.27436 2.7979
      3.0 4.0
      5.33697 1.88019]'
boundary_nodes = [
      [[1, 4, 3, 2], [2, 9, 10, 11, 8, 7, 12], [12, 6, 13, 5, 14, 15, 16, 17, 16], [16, 17, 18, 19, 20, 21, 22, 23, 1]],
      [[26, 25, 24], [24, 28, 27, 26]],
      [[29, 30, 31, 29]]
]
x, y = DT.pole_of_inaccessibility(pts, boundary_nodes)

fig = Figure()
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y")
bn1 = pts[:, unique(reduce(vcat, boundary_nodes[1]))] |> x -> hcat(x, x[:, begin])
bn2 = pts[:, unique(reduce(vcat, boundary_nodes[2]))] |> x -> hcat(x, x[:, begin])
bn3 = pts[:, unique(reduce(vcat, boundary_nodes[3]))] |> x -> hcat(x, x[:, begin])
lines!(ax, bn1, color=:red, linewidth=4)
lines!(ax, bn2, color=:red, linewidth=4)
lines!(ax, bn3, color=:red, linewidth=4)
scatter!(ax, [x], [y], color=:blue, markersize=23)
save("$save_path/pole_of_inaccessibility.png", fig)