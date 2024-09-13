using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using ColorSchemes
using DataStructures
using StableRNGs
import GeometryBasics: Point2f
using Random
using StaticArrays
using LinearAlgebra
using StructEquality
using GeometryBasics
@struct_equal DT.Queue

@testset "Simple example" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for _ in 1:10
            A = (-1.0, 7.0)
            B = (4.0, 4.0)
            C = (-2.0, -1.0)
            D = (-1.0, 3.0)
            E = (3.0, -1.0)
            F = (1.0, 4.0)
            G = (-3.0, 5.0)
            points = [A, B, C, D, E, F, G]
            weights = [0.0, 0.1, 0.2, -2.0, 17.0, 0.9, 1.5]
            tri = triangulate(points; weights, randomise = false, predicates = PT())
            vorn = voronoi(tri)
            @test DT.is_weighted(vorn)
            for (i, p) in DT.get_generators(vorn)
                @test get_point(tri, i) == get_generator(vorn, i) == p
            end
            @test validate_tessellation(vorn; predicates = PT())
            @test DT.get_triangulation(vorn) == tri
            circumcenter_to_triangle = DT.get_circumcenter_to_triangle(vorn)
            triangle_to_circumcenter = DT.get_triangle_to_circumcenter(vorn)
            for V in DT.each_solid_triangle(DT.get_triangulation(vorn))
                V = DT.sort_triangle(V)
                c = DT.get_triangle_to_circumcenter(vorn, V)
                c = get_polygon_point(vorn, c)
                i, j, k = triangle_vertices(V)
                p, q, r = get_point(DT.get_triangulation(vorn), i, j, k)
                cx, cy = DT.triangle_orthocenter(p, q, r, get_weight(tri, i), get_weight(tri, j), get_weight(tri, k))
                @test cx == c[1] && cy == c[2]
            end
            for (c, V) in circumcenter_to_triangle
                @test DT.get_circumcenter_to_triangle(vorn, c) == V
                @test DT.get_triangle_to_circumcenter(vorn, V) == c
            end
            for (V, c) in triangle_to_circumcenter
                @test DT.get_circumcenter_to_triangle(vorn, c) == V
                @test DT.get_triangle_to_circumcenter(vorn, V) == c
            end
            @test isempty(DT.get_boundary_polygons(vorn))
            poly1 = [7, 3, -3, -5, 7]
            poly2 = [3, 6, -1, -3, 3]
            poly3 = [1, 5, -4, -2, 1]
            poly4 = [4, 2, 5, 1, 4]
            poly5 = [1, -2, -1, 6, 4, 1]
            poly6 = [6, 3, 7, 2, 4, 6]
            poly7 = [2, 7, -5, -4, 5, 2]
            polys = [poly1, poly2, poly3, poly4, poly5, poly6, poly7]
            for (i, poly) in enumerate(polys)
                @test DT.circular_equality(get_polygon(vorn, i), poly)
            end
            polypoints = sort(reinterpret(reshape, Float64, DT.get_polygon_points(vorn))', dims = 1)
            _polypoints = sort(
                [
                    2.6333333333333 3.3633333333333
                    -0.765 5.14
                    2.6333333333333 7.4055555555556
                    -3.38 1.745
                    -1.025 4.1
                    -1.18 1.195
                    -0.10833333333333 2.2666666666667
                ], dims = 1,
            )
            @test polypoints ≈ _polypoints
        end
    end
end

@testset "get_nearest_neighbour" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        A = (-1.0, 7.0)
        B = (4.0, 4.0)
        C = (-2.0, -1.0)
        D = (-1.0, 3.0)
        E = (3.0, -1.0)
        F = (1.0, 4.0)
        G = (-3.0, 5.0)
        points = [A, B, C, D, E, F, G]
        weights = [0.0, 0.1, 0.2, -2.0, 17.0, 0.9, 1.5]
        tri = triangulate(points; weights, randomise = false, predicates = PT())
        vorn = voronoi(tri)
        bbox = (-20.0, 20.0, -20.0, 20.0)
        randx = rand(100000) * 40 .- 20
        randy = rand(100000) * 40 .- 20
        randpts = tuple.(randx, randy)
        neighs = [get_nearest_neighbour(vorn, q) for q in randpts]
        for i in each_polygon_index(vorn)
            C = get_polygon_coordinates(vorn, i, bbox)
            pop!(C)
            bnd = [1:length(C); 1]
            all_neighs = findall(==(i), neighs)
            qs = randpts[all_neighs]
            for q in qs
                d = DT.distance_to_polygon(q, C, bnd)
                @test d > 0
                all_dists = [DT.get_power_distance(tri, j, q) for j in DT.each_point_index(tri)]
                @test all_dists[i] == minimum(all_dists)
            end
        end
    end
end

@testset "Clip" begin
    A = (-1.0, 7.0)
    B = (4.0, 4.0)
    C = (-2.0, -1.0)
    D = (-1.0, 3.0)
    E = (3.0, -1.0)
    F = (1.0, 4.0)
    G = (-3.0, 5.0)
    points = [A, B, C, D, E, F, G]
    weights = [0.0, 0.1, 0.2, -12.0, 17.0, 6.9, 1.5]
    tri = triangulate(points; weights, randomise = false)
    vorn = voronoi(tri, clip = true)
    Cs = [get_polygon_coordinates(vorn, i) for i in 1:7]
    @test !DT.has_polygon(vorn, 4)
    Cs1 = [
        (-1.8125, 6.1875)
        (-1.365, 5.74)
        (-0.19736842105263186, 6.518421052631579)
        (-1.0, 7.0)
        (-1.8125, 6.1875)
    ]
    Cs2 = [
        (3.6333333333333333, 3.163333333333333)
        (3.825, 3.125)
        (4.0, 4.0)
        (3.633333333333333, 4.220000000000001)
        (3.6333333333333333, 3.163333333333333)
    ]
    Cs3 = [
        (-2.0, -1.0)
        (-1.1799999999999997, -1.0)
        (-1.1800000000000002, 1.2380000000000009)
        (-2.319565217391304, 1.9217391304347826)
        (-2.4824324324324323, 1.8945945945945946)
        (-2.0, -1.0)
    ]
    Cs4 = NTuple{2, Float64}[]
    Cs5 = [
        (-1.1799999999999997, -1.0)
        (3.0, -1.0)
        (3.825, 3.125)
        (3.6333333333333333, 3.163333333333333)
        (-1.1800000000000002, 1.2380000000000009)
        (-1.1799999999999997, -1.0)
    ]
    Cs6 = [
        (-2.319565217391304, 1.9217391304347826)
        (-1.1800000000000002, 1.2380000000000009)
        (3.6333333333333333, 3.163333333333333)
        (3.633333333333333, 4.220000000000001)
        (-0.19736842105263186, 6.518421052631579)
        (-1.365, 5.74)
        (-2.319565217391304, 1.9217391304347826)
    ]
    Cs7 = [
        (-2.4824324324324323, 1.8945945945945946)
        (-2.319565217391304, 1.9217391304347826)
        (-1.365, 5.74)
        (-1.8125, 6.1875)
        (-3.0, 5.0)
        (-2.4824324324324323, 1.8945945945945946)
    ]
    _Cs = (Cs1, Cs2, Cs3, Cs4, Cs5, Cs6, Cs7)
    for (C, _C) in zip(Cs, _Cs)
        @test DT.circular_equality(C, _C, ⪧)
    end
    @test validate_tessellation(vorn)
    @test isempty(DT.get_unbounded_polygons(vorn))
    @test all([0 ≤ get_area(vorn, i) < Inf for i in each_polygon_index(vorn)])
    @test allunique(DT.get_polygon_points(vorn))
    for i in DT.each_polygon_index(vorn)
        @test DT.has_polygon(vorn, i)
        C = get_polygon(vorn, i)
        for (j, v) in pairs(C)
            δ = DT.distance_to_polygon(get_polygon_point(vorn, v), get_points(tri), get_convex_hull_vertices(tri))
            @test δ ≥ -1.0e-15
        end
    end
end

@testset "A case where no segments of any tile intersect the convex hull" begin
    a = (3.0, 3.0)
    b = (0.0, 3.0)
    c = (0.0, 0.0)
    d = (4.0, 0.0)
    e = (1.0, 1.5)
    pts = [a, b, c, d, e]
    weights = [0.2, -3.0, 0.0, 1.7, 17.0]
    tri = triangulate(pts; weights, randomise = false)
    vorn = voronoi(tri, clip = true)
    @test validate_tessellation(vorn)
    @test isempty(DT.get_unbounded_polygons(vorn))
    @test get_adjacent(vorn, 7, 8) ==
        get_adjacent(vorn, 5, 6) ==
        get_adjacent(vorn, 8, 5) ==
        get_adjacent(vorn, 6, 7) == 5
    @test length(get_adjacent(get_adjacent(vorn))) == 4
    @test DT.get_polygons(vorn) == Dict(5 => [5, 6, 7, 8, 5])
    @test sort(DT.get_polygon_points(vorn)[5:8]) == sort([(4.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)])
    @test get_area(vorn, 5) == 10.5
    @test DT.get_boundary_polygons(vorn) == Set(5)
end

@testset "Single triangle" begin
    points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    weights = [-1.0, 0.0, 1.0]
    tri = triangulate(points; weights, randomise = false)
    vorn = voronoi(tri, clip = true)
    @test validate_tessellation(vorn)
    @test sort(collect.(vorn.polygon_points)) ≈ sort([[0.0, -0.5], [1.0, 0.0], [0.75, 0.25], [0.5, 0.0], [0.0, 1.0], [0.0, 0.0]])
end

@testset "A previously broken example with a non-boundary generator's intersecting edges not being previously detected" begin
    for _ in 1:1000
        a = (0.0, 0.0)
        b = (6.0, 0.0)
        c = (6.0, 6.0)
        d = (0.0, 6.0)
        e = (5.0, 3.0)
        f = (0.2, 5.8)
        g = (0.1, 0.2)
        h = (0.2, 0.1)
        pts = [a, b, c, d, e, f, g, h]
        weights = [
            3.462951732737026
            3.5064686402354086
            6.4774858066700505
            -0.4526740957045017
            6.493383218280129
            -7.892597046583529
            3.4256751890686794
            5.716293846261364
        ]
        PT = AdaptiveKernel
        tri = triangulate(pts; weights)
        vorn = voronoi(tri, predicates = PT(), clip = true)
        @test validate_tessellation(vorn)
        orig_pt = [
            (3.2942423167709074, 0.2669283425828489)
            (2.996373591041801, -17.00945774970529)
            (-7.171506359995909, 3.326302152370127)
            (2.4224866748021205, 5.52848734366764)
            (1.2729267417541226, 3.612554121920976)
            (10.7507029974137, 2.7524152361304464)
            (0.0, 3.569404062878463)
            (3.2896401039677543, 0.0)
            (0.0, 0.0)
            (2.4224866748021205, 6.0)
            (0.0, 6.0)
            (6.0, 1.1688475703258812)
            (6.0, 0.0)
            (6.0, 4.335982901935013)
            (6.0, 6.0)
        ]
        @test !DT.has_polygon(vorn, 1)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 2))), collect.(getindex.(Ref(orig_pt), [1, 8, 13, 12, 1])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 3))), collect.(getindex.(Ref(orig_pt), [14, 15, 10, 4, 14])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 4))), collect.(getindex.(Ref(orig_pt), [7, 5, 4, 10, 11, 7])), ≈)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 5))), collect.(getindex.(Ref(orig_pt), [1, 12, 14, 4, 5, 1])), ≈)
        @test !DT.has_polygon(vorn, 6)
        @test !DT.has_polygon(vorn, 7)
        @test DT.circular_equality(collect.(get_polygon_point.(Ref(vorn), get_polygon(vorn, 8))), collect.(getindex.(Ref(orig_pt), [9, 8, 1, 5, 7, 9])), ≈)
        @test isempty(DT.get_unbounded_polygons(vorn))
        @test all([0 ≤ get_area(vorn, i) < Inf for i in each_polygon_index(vorn)])
        @test allunique(DT.get_polygon_points(vorn))
        for i in each_polygon_index(vorn)
            C = get_polygon(vorn, i)
            for (j, v) in pairs(C)
                δ = DT.distance_to_polygon(get_polygon_point(vorn, v), get_points(tri), get_convex_hull_vertices(tri))
                @test δ ≥ -1.0e-14
            end
        end
        @test DT.get_boundary_polygons(vorn) == Set((4, 3, 5, 2, 8))
    end
end

@testset "Another single triangle example" begin
    (j, n) = (1, 3)
    rng = StableRNG(2^(isqrt(n)) * 3^(isqrt(n)))
    pts = rand(rng, 2, n)
    weights = 10randn(rng, n)
    tri = triangulate(pts; weights, rng)
    PT = AdaptiveKernel
    vorn = voronoi(tri, clip = true, predicates = PT())
    @test sort(vorn.polygon_points) ⪧ sort(
        [
            (-29.461464714041586, -32.54436775973904)
            (0.5358434495871474, 0.3351310497994868)
            (0.4967363663631179, 0.02543632261830825)
            (0.7720557069021632, 0.05488055609536713)
        ],
    )
    @test DT.circular_equality(vorn.polygons[1], [2, 3, 4, 2], ⪧)
    @test vorn.boundary_polygons == Set(1)
    @test validate_tessellation(vorn, predicates = PT())
end

@testset "An example with no intersections but two interior generators" begin
    (j, n) = (1, 7)
    rng = StableRNG(2^(isqrt(n)) * 3^(isqrt(n)))
    pts = rand(rng, 2, n)
    weights = 10randn(rng, n)
    tri = triangulate(pts; weights, rng)
    PT = AdaptiveKernel
    vorn = voronoi(tri, clip = true, predicates = PT())
    @test validate_tessellation(vorn, predicates = PT())
    pt = [
        (0.10462050798907063, 0.4620215842278311)
        (0.5828660494682045, 0.5300729131822177)
        (0.9643636416220359, 0.844658139116969)
        (0.23973629383552075, 0.8818897411817337)
        (0.10462050798907063, 0.4620215842278311)
    ]
    @test DT.circular_equality(vorn.polygon_points[vorn.polygons[3]], pt, ⪧)
    @test isempty(DT.get_unbounded_polygons(vorn))
    @test DT.get_boundary_polygons(vorn) == Set(3)
    @test length(vorn.polygons) == 1
end

@testset "Previously broken example with accidentally overwriting clip_vertices" begin
    (j, n) = (1, 11)
    rng = StableRNG(2^(isqrt(n)) * 3^(isqrt(n)))
    pts = rand(rng, 2, n)
    weights = 10randn(rng, n)
    tri = triangulate(pts; weights, rng)
    PT = AdaptiveKernel
    vorn = voronoi(tri, predicates = PT())
    flag1 = validate_tessellation(vorn, predicates = PT())
    vorn = voronoi(tri, clip = true, predicates = PT())
    flag2 = validate_tessellation(vorn, predicates = PT())
    @test all((flag1, flag2))
    pt = [
        (0.9297574334792393, 0.15491751490212047)
        (0.9042082877857884, 0.7717010771932071)
        (0.22439725232510632, 0.9386516631219246)
        (0.10000931451051631, 0.9048677530226168)
        (0.04890304624527486, 0.017339438428236598)
        (0.9297574334792393, 0.15491751490212047)
    ]
    @test DT.circular_equality(vorn.polygon_points[vorn.polygons[5]], pt, ⪧)
    @test isempty(DT.get_unbounded_polygons(vorn))
    @test DT.get_boundary_polygons(vorn) == Set(5)
    @test length(vorn.polygons) == 1
end

@testset "A case where there are no intersections but the surrounding tile comes from an unbounded polygon" begin
    (j, n) = (1, 15)
    rng = StableRNG(2^(isqrt(n)) * 3^(isqrt(n)))
    pts = rand(rng, 2, n)
    weights = 10randn(rng, n)
    tri = triangulate(pts; weights, rng)
    PT = AdaptiveKernel
    vorn = voronoi(tri, clip = true, predicates = PT())
    @test validate_tessellation(vorn, predicates = PT())
    pt = [
        (0.9297574334792393, 0.15491751490212047)
        (0.9812650498729347, 0.7130412548863212)
        (0.9042082877857883, 0.7717010771932071)
        (0.22439725232510627, 0.9386516631219246)
        (0.10000931451051631, 0.9048677530226168)
        (0.04890304624527487, 0.01733943842823671)
        (0.9297574334792393, 0.15491751490212047)
    ]
    @test DT.circular_equality(vorn.polygon_points[vorn.polygons[3]], pt, ⪧)
    @test isempty(DT.get_unbounded_polygons(vorn))
    @test DT.get_boundary_polygons(vorn) == Set(3)

    (j, n) = (1, 27)
    rng = StableRNG(2^(isqrt(n)) * 3^(isqrt(n)))
    pts = rand(rng, 2, n)
    weights = 10randn(rng, n)
    tri = triangulate(pts; weights, rng)
    PT = AdaptiveKernel
    vorn = voronoi(tri, clip = true, predicates = PT())
    @test validate_tessellation(vorn, predicates = PT())
    pt = [
        (0.04229809654442551, 0.18153429731036597)
        (0.16870777231496592, 0.09895916058959298)
        (0.7634316145681845, 0.15435495983060665)
        (0.8765020075792944, 0.22584838606743632)
        (0.8813286851587872, 0.23540562874585402)
        (0.9531849471555502, 0.41542093777623634)
        (0.9435552587190128, 0.7353245575729825)
        (0.8440901099158531, 0.9456321178014652)
        (0.31996403434288845, 0.9140172639926483)
        (0.033472676603685914, 0.4512476863895618)
        (0.04229809654442551, 0.18153429731036597)
    ]
    @test DT.circular_equality(vorn.polygon_points[vorn.polygons[23]], pt, ⪧)
    @test isempty(DT.get_unbounded_polygons(vorn))
    @test DT.get_boundary_polygons(vorn) == Set(23)
end

@testset "Varying size" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        for n in 3:4:200
            for j in 1:15:250
                rng = StableRNG(2^(isqrt(n)) * 3^(isqrt(n)))
                pts = rand(rng, 2, n)
                weights = 10randn(rng, n)
                tri = triangulate(pts; weights, rng)
                vorn = voronoi(tri, predicates = PT())
                flag1 = validate_tessellation(vorn, predicates = PT())
                vorn = voronoi(tri, clip = true, predicates = PT())
                flag2 = validate_tessellation(vorn, predicates = PT())
                @test flag1
                @test flag2
                flag2 || throw("...")
            end
        end
    end
end

@testset "Generic convex clipping" begin
    pts = [
        0.508812 0.656662 0.785124 0.63427 0.444969 0.609253 0.0826304 0.265388 0.830807 0.658346
        0.647732 0.482994 0.809909 0.482046 0.0170022 0.821742 0.835057 0.591724 0.881006 0.97652
    ]
    circ = CircularArc((0.75, 0.5), (0.75, 0.5), (0.5, 0.5))
    t = LinRange(0, 1, 20)
    clip_points = circ.(t)[begin:(end - 1)]
    clip_vertices = [1:(length(clip_points)); 1]
    clip_polygon = (clip_points, clip_vertices)
    weights = [6.12365, -1.50081, 0.853172, 6.57105, 0.863264, 10.1487, -5.45746, 5.78284, 4.70237, 2.69588]
    tri = triangulate(pts; weights)
    vorn = voronoi(tri; clip = true, clip_polygon = clip_polygon)
    @test length(vorn.polygons) == 1
    pt = [
        (0.28013156219837776, 0.38101315174073164)
        (0.33067960709356464, 0.3160690223317172)
        (0.39957614383675755, 0.27105666833623565)
        (0.4793551636319168, 0.25085387674833254)
        (0.5613713717851997, 0.2576499335151674)
        (0.6367370395306067, 0.2907083804343678)
        (0.6972851273490983, 0.3464468218275831)
        (0.7364543104251585, 0.41882513269882893)
        (0.75, 0.5)
        (0.7364543104251586, 0.5811748673011708)
        (0.6972851273490984, 0.653553178172417)
        (0.6367370395306067, 0.7092916195656321)
        (0.5613713717851998, 0.7423500664848326)
        (0.47935516363191694, 0.7491461232516674)
        (0.39957614383675766, 0.7289433316637643)
        (0.3306796070935648, 0.6839309776682829)
        (0.28013156219837776, 0.6189868482592684)
        (0.2534096741493194, 0.5411486475701835)
        (0.25340967414931936, 0.45885135242981656)
        (0.28013156219837776, 0.38101315174073164)
    ]
    @test DT.circular_equality(vorn.polygon_points[vorn.polygons[6]], pt, ⪧)
    @test vorn.boundary_polygons == Set(6)
    @test isempty(vorn.unbounded_polygons)
    @test validate_tessellation(vorn; check_area = false)

    weights = [
        -1.3249051535428942
        0.624963632666202
        0.4387318402335466
        0.7664689421091694
        -1.260416278903278
        0.6471611481089081
        0.9636964032629605
        0.2885455886737393
        1.191620421858312
        -1.7857456681807637
    ]
    tri = triangulate(pts; weights)
    vorn = voronoi(tri; clip = true, clip_polygon = clip_polygon)
    @test length(vorn.polygons) == 3
    pt4 = [
        (0.33067960709356464, 0.3160690223317172)
        (0.39957614383675755, 0.27105666833623565)
        (0.4793551636319168, 0.25085387674833254)
        (0.5185767219005764, 0.25410386772547405)
        (0.3359642248673101, 0.3440630424032038)
        (0.3237485736949492, 0.32497401422430655)
        (0.33067960709356464, 0.3160690223317172)
    ]
    pt7 = [
        (0.28013156219837776, 0.6189868482592684)
        (0.2534096741493194, 0.5411486475701835)
        (0.25340967414931936, 0.45885135242981656)
        (0.2801315621983778, 0.3810131517407316)
        (0.3237485736949492, 0.32497401422430655)
        (0.3359642248673101, 0.3440630424032038)
        (0.31623139994979865, 0.6653679207504803)
        (0.28013156219837776, 0.6189868482592684)
    ]
    pt9 = [
        (0.5185767219005764, 0.25410386772547405)
        (0.5613713717851996, 0.2576499335151674)
        (0.6367370395306068, 0.29070838043436786)
        (0.6972851273490983, 0.34644682182758296)
        (0.7364543104251587, 0.41882513269882926)
        (0.75, 0.4999999999999999)
        (0.7364543104251586, 0.581174867301171)
        (0.6972851273490985, 0.6535531781724169)
        (0.6367370395306067, 0.7092916195656321)
        (0.5613713717851998, 0.7423500664848326)
        (0.479355163631917, 0.7491461232516674)
        (0.39957614383675766, 0.7289433316637643)
        (0.3306796070935648, 0.6839309776682829)
        (0.31623139994979865, 0.6653679207504803)
        (0.3359642248673101, 0.3440630424032038)
        (0.5185767219005764, 0.25410386772547405)
    ]
    @test DT.circular_equality(vorn.polygon_points[vorn.polygons[4]], pt4, ⪧)
    @test DT.circular_equality(vorn.polygon_points[vorn.polygons[7]], pt7, ⪧)
    @test DT.circular_equality(vorn.polygon_points[vorn.polygons[9]], pt9, ⪧)
    @test vorn.boundary_polygons == Set((4, 7, 9))
    @test isempty(vorn.unbounded_polygons)
    @test validate_tessellation(vorn; check_area = false)
end

@testset "Smoothing" begin
    circ = CircularArc((1 / 2, 0.0), (1 / 2, 0.0), (0.0, 0.0))
    t = LinRange(0, 1, 20)
    clip_points = circ.(t)[begin:(end - 1)]
    clip_vertices = [1:(length(clip_points)); 1]
    clip_polygon = (clip_points, clip_vertices)
    flag = 0
    tot = 0
    for _ in 1:100
        points = 2randn(2, 50)
        weights = 2randn(50)
        tri = triangulate(points; weights)
        vorn = voronoi(tri; clip = true, clip_polygon = clip_polygon, smooth = true)
        @test validate_tessellation(vorn; check_area = false)
        for index in each_polygon_index(vorn)
            c = DT.get_centroid(vorn, index)
            p = get_generator(vorn, index)
            px, py = getxy(p)
            cx, cy = getxy(c)
            dx, dy = px - cx, py - cy
            ℓ = sqrt(dx^2 + dy^2)
            _flag = ℓ ≤ 1.0e-1
            flag += _flag
            tot += 1
        end
    end
    @test flag / tot > 0.9
end
