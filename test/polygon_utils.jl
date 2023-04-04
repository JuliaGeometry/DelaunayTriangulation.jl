using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using Random
using CairoMakie

save_path = basename(pwd()) == "test" ? "figures" : "test/figures"

include("./test_setup.jl")
include("./helper_functions.jl")

@testset "Getting polygon features" begin
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
end

@testset "Distance to a segment" begin
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
end

@testset "Distance to a polygon" begin
      tri, label_map, index_map = simple_geometry()
      pts = get_points(tri)
      @testset "Single boundary" begin
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
      end

      @testset "Multiple boundary segments" begin
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
      end
end

@testset "Bounding box of a polygon" begin
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

      # this example used to give ymin = ∞!
      pts =
            [
                  (0.42128834958962136, 0.33007028464908217)
                  (0.11356454007618466, 0.7448537954874419)
                  (0.7546603355923669, 0.9543777463196534)
                  (0.4891168285858787, 0.633382367024488)
                  (0.6747495735823583, 0.45029401396930835)
                  (0.4974345692650808, 0.5149317175333161)
                  (0.5484553916212294, 0.5711900118327666)
                  (0.11175023541896634, 0.990159314424705)
                  (0.4879170832093027, 0.08984797306499748)
                  (0.1335657114656048, 0.35758096957091445)
                  (0.7400877461824955, 0.8325280694798072)
                  (0.3299481327824305, 0.3440909795000966)
                  (0.1962438207259194, 0.6775012296614791)
                  (0.3403201981957973, 0.012234115125469014)
                  (0.39090662279892596, 0.6232084829209825)
                  (0.05180909728733263, 0.008306644625064141)
                  (0.4469104766158789, 0.5039047194497466)
                  (0.33193129503638996, 0.1768246437543215)
                  (0.24763476605581736, 0.9547830758766014)
                  (0.8626957317918005, 0.8901670309728742)
                  (0.16962017427458131, 0.8788693051101659)
                  (0.6974737865767218, 0.3655018057608477)
                  (0.5781761692908192, 0.49368701064930676)
                  (0.802284945950765, 0.6391848231098498)
                  (0.24014031334952324, 0.03642844544263135)
                  (0.29635836817046646, 0.49234998547822206)
                  (0.6537526603197776, 0.9534202877086324)
                  (0.22033649109831877, 0.6097719755673441)
                  (0.5794841917252405, 0.6525875695433809)
                  (0.48634161888118377, 0.7185107690604786)
                  (0.5345141678719951, 0.5951779828559485)
                  (0.07485448974139897, 0.3652052168490376)
                  (0.9456233879280223, 0.20388899534798632)
                  (0.27834285268176084, 0.8083123815440214)
                  (0.6267933326859505, 0.39246432872096704)
                  (0.7616653549409313, 0.6567908542485912)
                  (0.7064053508954178, 0.5295025690789412)
                  (0.6402160832134494, 0.7577312997966936)
                  (0.3919353829681529, 0.8457590619098538)
                  (0.9716293296512977, 0.5682387373301687)
            ]
      boundary_nodes = [16, 14, 33, 40, 20, 3, 8, 16]
      _pts = pts[boundary_nodes]
      xmin = minimum(getindex.(_pts, 1))
      xmax = maximum(getindex.(_pts, 1))
      ymin = minimum(getindex.(_pts, 2))
      ymax = maximum(getindex.(_pts, 2))
      _xmin, _xmax, _ymin, _ymax = DT.polygon_bounds(pts, boundary_nodes)
      @test xmin == _xmin 
      @test ymin == _ymin 
      @test xmax == _xmax 
      @test ymax == _ymax
end

@testset "Cell data structure" begin
      tri, label_map, index_map = simple_geometry()
      pts = get_points(tri)
      @testset "Constructor and operations" begin
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
      end

      @testset "CellQueue" begin
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
            q = DT.CellQueue{Float64}()
            @test isempty(q)
            @test q.queue.o == Base.Order.Reverse
            c = DT.Cell(polygon_features(pts, boundary_nodes)[2]..., 10.0, pts, boundary_nodes)
            DT.insert_cell!(q, c)
            @test q.queue[c] == c.max_dist
            @inferred DT.get_next_cell!(q)
      end
end

@testset "Pole of inaccessibility" begin
      pts = [(0.0, 10.0), (0.0, 8.0), (0.0, 6.0), (0.0, 4.0), (0.0, 2.0),
            (0.0, 0.0), (2.0, 0.0), (4.0, 0.0), (6.0, 0.0), (8.0, 0.0), (10.0, 0.0),
            (10.0, 2.0), (10.0, 4.0), (8.0, 4.0), (6.0, 4.0),
            (4.0, 4.0), (4.0, 6.0), (4.0, 8.0), (4.0, 10.0), (2.0, 10.0),
            (0.0, 10.0)]
      @testset "Single segment" begin
            boundary_nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1]
            x, y = DT.pole_of_inaccessibility(pts, boundary_nodes)
            @inferred DT.pole_of_inaccessibility(pts, boundary_nodes)
            @test x == 2.5 && y == 2.5
      end

      @testset "Multiple segments" begin
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
            @test (x ≈ 8.125 || x ≈ 5.625) && y ≈ 1.875
      end

      @testset "Multiply connected" begin
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
      end

      @testset "More complicated examples" begin
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
      end
end

@testset "Making a figure" begin
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
      SAVE_FIGURE && save("$save_path/pole_of_inaccessibility.png", fig)
end

@testset "A previously broken example" begin
      PT = [-3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402,
            3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912,
            2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224,
            2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273,
            1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557,
            0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467,
            -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364,
            -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026,
            -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916,
            -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402,
            3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912,
            2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224,
            2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273,
            1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557,
            0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467,
            -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364,
            -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026,
            -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916,
            -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402,
            3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912,
            2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224,
            2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273,
            1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467, -0.569282706070557,
            0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364, -1.4572528434350467,
            -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026, -2.3452229807995364,
            -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916, -3.233193118164026,
            -2.3452229807995364, -1.4572528434350467, -0.569282706070557, 0.31868743129393273, 1.2066575686584224, 2.094627706022912, 2.982597843387402, 3.8705679807518916,
            -4.927587132008278, -4.927587132008278, -4.927587132008278, -4.927587132008278, -4.927587132008278, -4.927587132008278, -4.927587132008278, -4.927587132008278,
            -4.927587132008278, -3.304869996725775, -3.304869996725775, -3.304869996725775, -3.304869996725775, -3.304869996725775, -3.304869996725775, -3.304869996725775,
            -3.304869996725775, -3.304869996725775, -1.6821528614432721, -1.6821528614432721, -1.6821528614432721, -1.6821528614432721, -1.6821528614432721, -1.6821528614432721,
            -1.6821528614432721, -1.6821528614432721, -1.6821528614432721, -0.05943572616076942, -0.05943572616076942, -0.05943572616076942, -0.05943572616076942, -0.05943572616076942,
            -0.05943572616076942, -0.05943572616076942, -0.05943572616076942, -0.05943572616076942, 1.5632814091217337, 1.5632814091217337, 1.5632814091217337, 1.5632814091217337,
            1.5632814091217337, 1.5632814091217337, 1.5632814091217337, 1.5632814091217337, 1.5632814091217337, 3.185998544404236, 3.185998544404236, 3.185998544404236, 3.185998544404236,
            3.185998544404236, 3.185998544404236, 3.185998544404236, 3.185998544404236, 3.185998544404236, 4.808715679686739, 4.808715679686739, 4.808715679686739, 4.808715679686739,
            4.808715679686739, 4.808715679686739, 4.808715679686739, 4.808715679686739, 4.808715679686739, 6.431432814969242, 6.431432814969242, 6.431432814969242, 6.431432814969242,
            6.431432814969242, 6.431432814969242, 6.431432814969242, 6.431432814969242, 6.431432814969242, 8.054149950251745, 8.054149950251745, 8.054149950251745, 8.054149950251745,
            8.054149950251745, 8.054149950251745, 8.054149950251745, 8.054149950251745, 8.054149950251745, 9.676867085534248, 9.676867085534248, 9.676867085534248, 9.676867085534248,
            9.676867085534248, 9.676867085534248, 9.676867085534248, 9.676867085534248, 9.676867085534248, 11.299584220816751, 11.299584220816751, 11.299584220816751, 11.299584220816751,
            11.299584220816751, 11.299584220816751, 11.299584220816751, 11.299584220816751, 11.299584220816751, 12.922301356099254, 12.922301356099254, 12.922301356099254, 12.922301356099254,
            12.922301356099254, 12.922301356099254, 12.922301356099254, 12.922301356099254, 12.922301356099254, 14.545018491381757, 14.545018491381757, 14.545018491381757, 14.545018491381757,
            14.545018491381757, 14.545018491381757, 14.545018491381757, 14.545018491381757, 14.545018491381757, 16.16773562666426, 16.16773562666426, 16.16773562666426, 16.16773562666426,
            16.16773562666426, 16.16773562666426, 16.16773562666426, 16.16773562666426, 16.16773562666426, 17.790452761946764, 17.790452761946764, 17.790452761946764, 17.790452761946764,
            17.790452761946764, 17.790452761946764, 17.790452761946764, 17.790452761946764, 17.790452761946764, 19.413169897229267, 19.413169897229267, 19.413169897229267, 19.413169897229267,
            19.413169897229267, 19.413169897229267, 19.413169897229267, 19.413169897229267, 19.413169897229267, 21.03588703251177, 21.03588703251177, 21.03588703251177, 21.03588703251177,
            21.03588703251177, 21.03588703251177, 21.03588703251177, 21.03588703251177, 21.03588703251177, 22.658604167794273, 22.658604167794273, 22.658604167794273, 22.658604167794273,
            22.658604167794273, 22.658604167794273, 22.658604167794273, 22.658604167794273, 22.658604167794273, 24.281321303076776, 24.281321303076776, 24.281321303076776, 24.281321303076776,
            24.281321303076776, 24.281321303076776, 24.281321303076776, 24.281321303076776, 24.281321303076776, 25.90403843835928, 25.90403843835928, 25.90403843835928, 25.90403843835928,
            25.90403843835928, 25.90403843835928, 25.90403843835928, 25.90403843835928, 25.90403843835928, 27.52675557364178, 27.52675557364178, 27.52675557364178, 27.52675557364178,
            27.52675557364178, 27.52675557364178, 27.52675557364178, 27.52675557364178, 27.52675557364178, 29.149472708924286, 29.149472708924286, 29.149472708924286, 29.149472708924286,
            29.149472708924286, 29.149472708924286, 29.149472708924286, 29.149472708924286, 29.149472708924286, 30.772189844206785, 30.772189844206785, 30.772189844206785, 30.772189844206785,
            30.772189844206785, 30.772189844206785, 30.772189844206785, 30.772189844206785, 30.772189844206785, 32.39490697948929, 32.39490697948929, 32.39490697948929, 32.39490697948929,
            32.39490697948929, 32.39490697948929, 32.39490697948929, 32.39490697948929, 32.39490697948929]
      PT = [PT[1:216]'; PT[217:end]']
      BN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 18, 27, 36, 45, 54, 63, 72,
            81, 90, 99, 108, 117, 126, 135, 144, 153, 162, 171, 180, 189,
            198, 207, 216, 215, 214, 213, 212, 211, 210, 209, 208, 199,
            190, 181, 172, 163, 154, 145, 136, 127, 118, 109, 100, 91, 82,
            73, 64, 55, 46, 37, 28, 19, 10, 1]
      pc = polylabel(PT, BN)
      @test collect(pc) ≈ collect(DT.polygon_features(PT, BN)[2])
end