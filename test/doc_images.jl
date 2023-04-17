using ..DelaunayTriangulation
using ReferenceTests
using CairoMakie
using StableRNGs
using DelimitedFiles
const DT = DelaunayTriangulation

include("./helper_functions.jl")

@testset "Unconstrained Triangulations" begin
    a = [1.5, 4.0]
    b = [0.0, 3.5]
    c = [2.0, 1.5]
    d = [3.0, 2.5]
    e = [2.5, 3.5]
    f = [0.5, 3.0]
    g = [2.5, -2.0]
    h = [0.5, 1.5]
    i = [0.0, 0.5]
    j = [1.5, 3.0]
    pts = [a, b, c, d, e, f, g, h, i, j]
    tri = triangulate(pts)
    fig, ax, sc = triplot(tri)
    @test_reference "../docs/src/triangulations/figs/small_example.png" fig
end

@testset "Constrained Triangulation" begin
    @testset "Constrained edges only" begin
        a = (0.0, 0.0)
        b = (0.0, 1.0)
        c = (0.0, 2.5)
        d = (2.0, 0.0)
        e = (6.0, 0.0)
        f = (8.0, 0.0)
        g = (8.0, 0.5)
        h = (7.5, 1.0)
        i = (4.0, 1.0)
        j = (4.0, 2.5)
        k = (8.0, 2.5)
        pts = [a, b, c, d, e, f, g, h, i, j, k]
        C = Set([(2, 1), (2, 11), (2, 7), (2, 5)])
        uncons_tri = triangulate(pts)
        cons_tri = triangulate(pts; edges=C)
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(a):$ $ Unconstrained", titlealign=:left)
        triplot!(ax, uncons_tri)
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(b):$ $  Constrained", titlealign=:left)
        triplot!(ax, cons_tri)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/simple_constrained.png" fig
    end

    @testset "With an outer boundary" begin
        pts = [
            (-7.36, 12.55), (-9.32, 8.59), (-9.0, 3.0), (-6.32, -0.27),
            (-4.78, -1.53), (2.78, -1.41), (-5.42, 1.45), (7.86, 0.67),
            (10.92, 0.23), (9.9, 7.39), (8.14, 4.77), (13.4, 8.61),
            (7.4, 12.27), (2.2, 13.85), (-3.48, 10.21), (-4.56, 7.35),
            (3.44, 8.99), (3.74, 5.87), (-2.0, 8.0), (-2.52, 4.81),
            (1.34, 6.77), (1.24, 4.15)
        ]
        boundary_points = [
            (0.0, 0.0), (2.0, 1.0), (3.98, 2.85), (6.0, 5.0),
            (7.0, 7.0), (7.0, 9.0), (6.0, 11.0), (4.0, 12.0),
            (2.0, 12.0), (1.0, 11.0), (0.0, 9.13), (-1.0, 11.0),
            (-2.0, 12.0), (-4.0, 12.0), (-6.0, 11.0), (-7.0, 9.0),
            (-6.94, 7.13), (-6.0, 5.0), (-4.0, 3.0), (-2.0, 1.0), (0.0, 0.0)
        ]
        boundary_nodes, pts = convert_boundary_points_to_indices(boundary_points; existing_points=pts)
        uncons_tri = triangulate(pts, delete_ghosts=false)
        cons_tri = triangulate(pts; boundary_nodes, delete_ghosts=false)
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(a):$ $ Unconstrained", titlealign=:left)
        triplot!(ax, uncons_tri)
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(b):$ $  Constrained", titlealign=:left)
        triplot!(ax, cons_tri)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/heart_constrained.png" fig

        add_point!(cons_tri, 0.0, 5.0)
        add_edge!(cons_tri, 40, 26)
        add_edge!(cons_tri, 39, 27)
        add_edge!(cons_tri, 38, 28)
        add_point!(cons_tri, -3.0, 12.0)

        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(a):$ $ Unconstrained", titlealign=:left)
        triplot!(ax, uncons_tri)
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(b):$ $  Constrained", titlealign=:left)
        triplot!(ax, cons_tri)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/heart_add_constrained.png" fig
    end

    @testset "Segmented outer boundary" begin
        points = [
            (2.0, 8.0), (6.0, 4.0), (2.0, 6.0),
            (2.0, 4.0), (8.0, 2.0)
        ]
        segment_1 = [(0.0, 0.0), (14.0, 0.0)]
        segment_2 = [(14.0, 0.0), (10.0, 4.0), (4.0, 6.0), (2.0, 12.0), (0.0, 14.0)]
        segment_3 = [(0.0, 14.0), (0.0, 0.0)]
        boundary_points = [segment_1, segment_2, segment_3]
        boundary_nodes, points = convert_boundary_points_to_indices(boundary_points; existing_points=points)
        uncons_tri = triangulate(points)
        cons_tri = triangulate(points; boundary_nodes)

        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(a):$ $ Unconstrained", titlealign=:left)
        triplot!(ax, uncons_tri)
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(b):$ $  Constrained", titlealign=:left)
        triplot!(ax, cons_tri)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/triangle_triangulation.png" fig

        add_ghost_triangles!(cons_tri)
        add_ghost_triangles!(uncons_tri)
        for x in LinRange(0.1, 13.9, 10) # bottom side
            add_point!(uncons_tri, x, 0.0)
            add_point!(cons_tri, x, 0.0)
        end
        for x in LinRange(13.9, 10.1, 10) # first part of diagonal
            add_point!(uncons_tri, x, 14 - x)
            add_point!(cons_tri, x, 14 - x)
        end
        for x in LinRange(9.9, 4.1, 10) # second part of diagonal 
            add_point!(uncons_tri, x, 22 // 3 - x / 3)
            add_point!(cons_tri, x, 22 // 3 - x / 3)
        end
        for x in LinRange(3.9, 2.1, 10) # third part of diagonal 
            add_point!(uncons_tri, x, 18 - 3x)
            add_point!(cons_tri, x, 18 - 3x)
        end
        for x in LinRange(1.9, 0.1, 10) # last part of diagonal
            add_point!(uncons_tri, x, 14 - x)
            add_point!(cons_tri, x, 14 - x)
        end
        for y in LinRange(13.9, 0.1, 10) # left 
            add_point!(uncons_tri, 0.0, y)
            add_point!(cons_tri, 0.0, y)
        end
        delete_ghost_triangles!(cons_tri)
        delete_ghost_triangles!(uncons_tri)

        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(a):$ $ Unconstrained", titlealign=:left)
        triplot!(ax, uncons_tri)
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(b):$ $  Constrained", titlealign=:left)
        triplot!(ax, cons_tri)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/triangle_triangulation_refined.png" fig
    end

    @testset "Domain with interior holes" begin
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
        rng = StableRNG(1992881)
        uncons_tri = triangulate(points; rng)
        cons_tri = triangulate(points; boundary_nodes=boundary_nodes, rng)

        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(a):$ $ Unconstrained", titlealign=:left)
        triplot!(ax, uncons_tri)
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(b):$ $  Constrained", titlealign=:left)
        triplot!(ax, cons_tri)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/multiply_connected.png" fig
    end

    @testset "Interior holes within interiors" begin
        curve_1 = [
            [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0), (15.0, 0.0), (20.0, 0.0), (25.0, 0.0)],
            [(25.0, 0.0), (25.0, 5.0), (25.0, 10.0), (25.0, 15.0), (25.0, 20.0), (25.0, 25.0)],
            [(25.0, 25.0), (20.0, 25.0), (15.0, 25.0), (10.0, 25.0), (5.0, 25.0), (0.0, 25.0)],
            [(0.0, 25.0), (0.0, 20.0), (0.0, 15.0), (0.0, 10.0), (0.0, 5.0), (0.0, 0.0)]
        ] # outer-most boundary: counter-clockwise  
        curve_2 = [
            [(4.0, 6.0), (4.0, 14.0), (4.0, 20.0), (18.0, 20.0), (20.0, 20.0)],
            [(20.0, 20.0), (20.0, 16.0), (20.0, 12.0), (20.0, 8.0), (20.0, 4.0)],
            [(20.0, 4.0), (16.0, 4.0), (12.0, 4.0), (8.0, 4.0), (4.0, 4.0), (4.0, 6.0)]
        ] # inner boundary: clockwise 
        curve_3 = [
            [(12.906, 10.912), (16.0, 12.0), (16.16, 14.46), (16.29, 17.06),
            (13.13, 16.86), (8.92, 16.4), (8.8, 10.9), (12.906, 10.912)]
        ] # this is inside curve_2, so it's counter-clockwise 
        curves = [curve_1, curve_2, curve_3]
        points = [
            (3.0, 23.0), (9.0, 24.0), (9.2, 22.0), (14.8, 22.8), (16.0, 22.0),
            (23.0, 23.0), (22.6, 19.0), (23.8, 17.8), (22.0, 14.0), (22.0, 11.0),
            (24.0, 6.0), (23.0, 2.0), (19.0, 1.0), (16.0, 3.0), (10.0, 1.0), (11.0, 3.0),
            (6.0, 2.0), (6.2, 3.0), (2.0, 3.0), (2.6, 6.2), (2.0, 8.0), (2.0, 11.0),
            (5.0, 12.0), (2.0, 17.0), (3.0, 19.0), (6.0, 18.0), (6.5, 14.5),
            (13.0, 19.0), (13.0, 12.0), (16.0, 8.0), (9.8, 8.0), (7.5, 6.0),
            (12.0, 13.0), (19.0, 15.0)
        ]
        rng = StableRNG(123)
        boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points)
        uncons_tri = triangulate(points; rng)
        cons_tri = triangulate(points; boundary_nodes=boundary_nodes, rng, check_arguments=false)

        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(a):$ $ Unconstrained", titlealign=:left)
        triplot!(ax, uncons_tri)
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(b):$ $  Constrained", titlealign=:left)
        triplot!(ax, cons_tri)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/multiply_connected_interior_interior.png" fig
    end

    @testset "Disjoint domains" begin
        θ = LinRange(0, 2π, 20) |> collect
        θ[end] = 0 # need to make sure that 2π gives the exact same coordinates as 0
        xy = Vector{Vector{Vector{NTuple{2,Float64}}}}()
        cx = 0.0
        for i in 1:2
            # Make the exterior circle
            push!(xy, [[(cx + cos(θ), sin(θ)) for θ in θ]])
            # Now the interior circle - clockwise
            push!(xy, [[(cx + 0.5cos(θ), 0.5sin(θ)) for θ in reverse(θ)]])
            cx += 3.0
        end
        boundary_nodes, points = convert_boundary_points_to_indices(xy)
        uncons_tri = triangulate(points)
        cons_tri = triangulate(points; boundary_nodes=boundary_nodes, check_arguments=false)

        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(a):$ $ Unconstrained", titlealign=:left)
        triplot!(ax, uncons_tri)
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(b):$ $  Constrained", titlealign=:left)
        triplot!(ax, cons_tri)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/simple_disjoint.png" fig
    end

    @testset "Julia" begin
        C = (15.7109521325776, 33.244486807457)
        D = (14.2705719699703, 32.8530791545746)
        E = (14.3, 27.2)
        F = (14.1, 27.0)
        G = (13.7, 27.2)
        H = (13.4, 27.5)
        I = (13.1, 27.6)
        J = (12.7, 27.4)
        K = (12.5, 27.1)
        L = (12.7, 26.7)
        M = (13.1, 26.5)
        N = (13.6, 26.4)
        O = (14.0, 26.4)
        P = (14.6, 26.5)
        Q = (15.1983491346581, 26.8128534095401)
        R = (15.6, 27.6)
        S = (15.6952958264624, 28.2344688505621)
        T = (17.8088971520274, 33.1192363585346)
        U = (16.3058917649589, 33.0722674401887)
        V = (16.3215480710742, 29.7374742376305)
        W = (16.3841732955354, 29.393035503094)
        Z = (16.6190178872649, 28.9233463196351)
        A1 = (17.0417381523779, 28.5319386667527)
        B1 = (17.5114273358368, 28.3753756055997)
        C1 = (18.1376795804487, 28.3597192994844)
        D1 = (18.7169629067146, 28.5632512789833)
        E1 = (19.2805899268653, 28.8920337074045)
        F1 = (19.26493362075, 28.4536571361762)
        G1 = (20.6426885588962, 28.4223445239456)
        H1 = (20.689657477242, 33.1035800524193)
        I1 = (19.2805899268653, 33.0722674401887)
        J1 = (19.2962462329806, 29.7531305437458)
        K1 = (19.0614016412512, 29.393035503094)
        L1 = (18.7482755189452, 29.236472441941)
        M1 = (18.4508057027546, 29.1425346052493)
        N1 = (18.1689921926793, 29.3147539725175)
        O1 = (17.7932408459121, 29.6278800948235)
        P1 = (22.6466957416542, 35.4207133574833)
        Q1 = (21.2219718851621, 34.9979930923702)
        R1 = (21.2376281912774, 28.4693134422915)
        S1 = (22.6780083538847, 28.4380008300609)
        T1 = (24.5724213938357, 33.1975178891111)
        U1 = (23.3512295168425, 32.8530791545746)
        V1 = (23.3199169046119, 28.4380008300609)
        W1 = (24.6663592305274, 28.3753756055997)
        Z1 = (15.1942940307729, 35.4363696635986)
        A2 = (14.7246048473139, 35.3737444391374)
        B2 = (14.3645098066621, 35.1858687657538)
        C2 = (14.1766341332786, 34.8570863373326)
        D2 = (14.1140089088174, 34.3247719294125)
        E2 = (14.2705719699703, 33.8394264398383)
        F2 = (14.7246048473139, 33.6202381542241)
        G2 = (15.4604512347329, 33.6045818481088)
        H2 = (16.0, 34.0)
        I2 = (15.9771093365377, 34.6848669700643)
        J2 = (15.6170142958859, 35.2328376840997)
        K2 = (24.1653574348379, 35.4520259697138)
        L2 = (23.7739497819555, 35.4363696635986)
        M2 = (23.4608236596496, 35.2641502963303)
        N2 = (23.272947986266, 34.9040552556785)
        O2 = (23.1320412312284, 34.5909291333725)
        P2 = (23.1163849251131, 34.2151777866054)
        Q2 = (23.2886042923813, 33.8081138276077)
        R2 = (23.8209187003014, 33.6045818481088)
        S2 = (24.3062641898756, 33.5576129297629)
        T2 = (24.7602970672192, 33.8550827459536)
        U2 = (25.010797965064, 34.4656786844502)
        V2 = (24.8385785977957, 34.9666804801397)
        W2 = (24.5254524754898, 35.2641502963303)
        Z2 = (25.3708930057158, 37.4716894585871)
        A3 = (24.7916096794498, 37.3464390096648)
        B3 = (24.4471709449133, 36.9550313567823)
        C3 = (24.3062641898756, 36.5636237038999)
        D3 = (24.4941398632592, 35.9999966837492)
        E3 = (25.0264542711793, 35.5929327247515)
        F3 = (25.5587686790994, 35.5929327247515)
        F3 = (25.5587686790994, 35.5929327247515)
        G3 = (26.0, 36.0)
        H3 = (26.1380520053653, 36.5792800100152)
        I3 = (26.0, 37.0)
        J3 = (25.7466443524829, 37.2838137852036)
        K3 = (26.3885529032101, 35.4676822758291)
        L3 = (25.9814889442124, 35.3580881330221)
        M3 = (25.6840191280217, 35.1858687657538)
        N3 = (25.5274560668688, 34.9040552556785)
        O3 = (25.4961434546382, 34.5596165211419)
        P3 = (25.5274560668688, 34.246490398836)
        Q3 = (25.6683628219064, 33.8394264398383)
        R3 = (26.0284578625583, 33.6358944603394)
        S3 = (26.5451159643631, 33.6202381542241)
        T3 = (27.0, 34.0)
        U3 = (27.280962351782, 34.5596165211419)
        V3 = (27.0304614539373, 35.2171813779844)
        W3 = (26.1693646175959, 33.087923746304)
        Z3 = (26.0, 33.0)
        A4 = (25.5274560668688, 32.7278287056522)
        B4 = (25.2612988629087, 32.4147025833463)
        C4 = (25.1830173323322, 32.0702638488098)
        D4 = (25.2299862506781, 31.7727940326191)
        E4 = (25.6527065157911, 31.5222931347744)
        F4 = (26.2946150665183, 31.7258251142732)
        G4 = (26.5607722704784, 32.5086404200381)
        H4 = (27.1557119028596, 32.7434850117675)
        I4 = (27.6097447802033, 32.4929841139228)
        J4 = (27.6410573924338, 32.1015764610403)
        K4 = (27.7193389230103, 31.6005746653509)
        L4 = (27.437525412935, 31.4283552980826)
        M4 = (26.9834925355914, 31.2561359308143)
        N4 = (26.5764285765937, 31.0995728696614)
        O4 = (26.0441141686736, 30.7864467473554)
        P4 = (25.6527065157911, 30.5672584617413)
        Q4 = (25.3239240873699, 30.1915071149741)
        R4 = (25.1673610262169, 29.8783809926682)
        S4 = (25.1047358017558, 29.6122237887082)
        T4 = (25.0890794956405, 29.1895035235952)
        U4 = (25.2926114751393, 28.8294084829433)
        V4 = (25.6840191280217, 28.5632512789833)
        W4 = (26.1537083114806, 28.3753756055997)
        Z4 = (26.8269294744384, 28.391031911715)
        A5 = (27.4844943312809, 28.6102201973292)
        B5 = (27.7342002330051, 28.7239579596219)
        C5 = (27.7264126450755, 28.4202565942047)
        D5 = (29.1825559185446, 28.3922538389457)
        E5 = (29.1545531632856, 32.2146299318021)
        F5 = (29.000538009361, 32.5786657501693)
        G5 = (28.6785063238822, 32.9006974356481)
        H5 = (28.3144705055149, 33.0827153448317)
        I5 = (27.9084305542591, 33.2367304987563)
        J5 = (27.3343740714492, 33.3207387645334)
        K5 = (26.8303244767868, 33.2367304987563)
        L5 = (27.6564057569279, 30.786489413592)
        M5 = (27.6984098898165, 30.3944508399657)
        N5 = (27.6984098898165, 29.7363860913787)
        O5 = (27.5863988687804, 29.4143544059)
        P5 = (27.2643671833016, 29.2043337414573)
        Q5 = (26.9843396307114, 29.1763309861983)
        R5 = (26.6903107004917, 29.3163447624934)
        S5 = (26.5782996794556, 29.7503874690082)
        T5 = (26.7603175886393, 30.3384453294476)
        U5 = (27.3203726938197, 30.7024811478149)

        J_curve = [[C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, C]]
        U_curve = [[T, U, V, W, Z, A1, B1, C1, D1, E1, F1, G1, H1, I1, J1, K1, L1, M1, N1, O1, T]]
        L_curve = [[P1, Q1, R1, S1, P1]]
        I_curve = [[T1, U1, V1, W1, T1]]
        A_curve_outline = [[
            K5, W3, Z3, A4, B4, C4, D4, E4, F4, G4, H4, I4, J4, K4, L4, M4, N4,
            O4, P4, Q4, R4, S4, T4, U4, V4, W4, Z4, A5, B5, C5, D5, E5, F5, G5,
            H5, I5, J5, K5]]
        A_curve_hole = [[L5, M5, N5, O5, P5, Q5, R5, S5, T5, U5, L5]]
        dot_1 = [[Z1, A2, B2, C2, D2, E2, F2, G2, H2, I2, J2, Z1]]
        dot_2 = [[Z2, A3, B3, C3, D3, E3, F3, G3, H3, I3, J3, Z2]]
        dot_3 = [[K2, L2, M2, N2, O2, P2, Q2, R2, S2, T2, U2, V2, W2, K2]]
        dot_4 = [[K3, L3, M3, N3, O3, P3, Q3, R3, S3, T3, U3, V3, K3]]
        curves = [J_curve, U_curve, L_curve, I_curve, A_curve_outline, A_curve_hole, dot_1, dot_2, dot_3, dot_4]
        nodes, points = convert_boundary_points_to_indices(curves)
        uncons_tri = triangulate(points)
        cons_tri = triangulate(points; boundary_nodes=nodes, check_arguments=false)

        fig = Figure()
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(a):$ $ Unconstrained", titlealign=:left)
        triplot!(ax, uncons_tri)
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
            title=L"(b):$ $  Constrained", titlealign=:left)
        triplot!(ax, cons_tri)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/julia.png" fig
    end
end

@testset "Lattice triangulations" begin
    @testset "Four boundaries" begin
        a, b, c, d = 2.0, 10.0, -5.0, 7.5
        nx = 20
        ny = 10
        tri = DT.triangulate_rectangle(a, b, c, d, nx, ny)
        fig, ax, sc = triplot(tri; show_ghost_edges=true)
        xlims!(ax, a - 0.5, b + 0.5)
        ylims!(ax, c - 0.5, d + 0.5)
        lines!(ax, tri.points[get_boundary_nodes(tri, 1)]; linewidth=4)
        lines!(ax, tri.points[get_boundary_nodes(tri, 2)]; linewidth=4)
        lines!(ax, tri.points[get_boundary_nodes(tri, 3)]; linewidth=4)
        lines!(ax, tri.points[get_boundary_nodes(tri, 4)]; linewidth=4)
        @test_reference "../docs/src/triangulations/figs/rectangular_triangulation_1.png" fig
    end

    @testset "Single boundary" begin
        a, b, c, d = 2.0, 10.0, -5.0, 7.5
        nx = 20
        ny = 10
        tri = DT.triangulate_rectangle(a, b, c, d, nx, ny; single_boundary=true)
        fig, ax, sc = triplot(tri; show_ghost_edges=true)
        xlims!(ax, a - 0.5, b + 0.5)
        ylims!(ax, c - 0.5, d + 0.5)
        @test_reference "../docs/src/triangulations/figs/rectangular_triangulation_2.png" fig
    end
end

if !(get(ENV, "CI", "false") == "true")
    @testset "Gmsh" begin
        @testset "Contiguous boundary" begin
            a = 4 / 5
            t = LinRange(0, 2π, 100)
            x = @. a * (2cos(t) + cos(2t))
            y = @. a * (2sin(t) - sin(2t))
            tri = generate_mesh(x, y, 0.1)
            tri2 = generate_mesh(x, y, 1.0)
            fig = Figure()
            ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=300, height=300,
                title=L"(a):$ $ Dense mesh", titlealign=:left)
            triplot!(ax, tri)
            ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", width=300, height=300,
                title=L"(b):$ $  Coarse mesh", titlealign=:left)
            triplot!(ax, tri2)
            resize_to_layout!(fig)
            @test_reference "../docs/src/triangulations/figs/gmsh_example_1.png" fig
        end

        @testset "Single boundary curve with multiple segments" begin
            # The first segment 
            t = LinRange(0, 1 / 4, 25)
            x1 = cos.(2π * t)
            y1 = sin.(2π * t)
            # The second segment 
            t = LinRange(0, -3, 25)
            x2 = collect(t)
            y2 = repeat([1.0], length(t))
            # The third segment 
            t = LinRange(1, 0, 25)
            x3 = -3.0 .+ (1 .- t) .* sin.(t)
            y3 = collect(t)
            # The fourth segment 
            t = LinRange(0, 1, 25)
            x4 = collect(-3.0(1 .- t))
            y4 = collect(0.98t)
            # The fifth segment 
            x5 = [0.073914, 0.0797, 0.1522, 0.1522, 0.2, 0.28128, 0.3659, 0.4127, 0.3922, 0.4068, 0.497, 0.631, 0.728, 0.804, 0.888, 1.0]
            y5 = [0.8815, 0.8056, 0.80268, 0.73258, 0.6, 0.598, 0.5777, 0.525, 0.4346, 0.3645, 0.3032, 0.2886, 0.2623, 0.1367, 0.08127, 0.0]
            # Now combine the vectors 
            x = [x1, x2, x3, x4, x5]
            y = [y1, y2, y3, y4, y5]
            # Mesh 
            tri = generate_mesh(x, y, 0.05)
            fig = Figure()
            ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", width=600, height=300)
            triplot!(ax, tri)
            colors = [:red, :blue, :orange, :purple, :darkgreen]
            bn_map = get_boundary_map(tri)
            for (i, segment_index) in enumerate(values(bn_map))
                bn_nodes = get_boundary_nodes(tri, segment_index)
                lines!(ax, get_points(tri)[:, bn_nodes], color=colors[i], linewidth=4)
            end
            resize_to_layout!(fig)
            @test_reference "../docs/src/triangulations/figs/gmsh_example_2.png" fig
        end

        @testset "Multiple boundaries" begin
            x1 = [collect(LinRange(0, 2, 4)),
                collect(LinRange(2, 2, 4)),
                collect(LinRange(2, 0, 4)),
                collect(LinRange(0, 0, 4))]
            y1 = [collect(LinRange(0, 0, 4)),
                collect(LinRange(0, 6, 4)),
                collect(LinRange(6, 6, 4)),
                collect(LinRange(6, 0, 4))]
            r = 0.5
            h = k = 0.6
            θ = LinRange(2π, 0, 50)
            x2 = [h .+ r .* cos.(θ)]
            y2 = [k .+ r .* sin.(θ)]
            r = 0.2
            h = 1.5
            k = 0.5
            x3 = [h .+ r .* cos.(θ)]
            y3 = [k .+ r .* sin.(θ)]
            x4 = reverse(reverse.([collect(LinRange(1, 1.5, 4)),
                collect(LinRange(1.5, 1.5, 4)),
                collect(LinRange(1.5, 1, 4)),
                collect(LinRange(1, 1, 4))]))
            y4 = reverse(reverse.([collect(LinRange(2, 2, 4)),
                collect(LinRange(2, 5, 4)),
                collect(LinRange(5, 5, 4)),
                collect(LinRange(5, 2, 4))]))
            x5 = [reverse([0.2, 0.5, 0.75, 0.75, 0.2, 0.2])]
            y5 = [reverse([2.0, 2.0, 3.0, 4.0, 5.0, 2.0])]
            x = [x1, x2, x3, x4, x5]
            y = [y1, y2, y3, y4, y5]
            tri = generate_mesh(x, y, 0.2)
            fig, ax, sc = triplot(tri; show_ghost_edges=true, convex_hull_linestyle=:solid, convex_hull_linewidth=4)
            xlims!(ax, -0.5, 2.5)
            ylims!(ax, -0.5, 6.5)
            @test_reference "../docs/src/triangulations/figs/gmsh_example_3.png" fig
        end
    end
end

@testset "Convex polygons" begin
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
    @test_reference "../docs/src/triangulations/figs/convex_triangulation_example.png" fig
end

@testset "Explaining how we compute constrained triangulations" begin
    ## Segment location 
    # Example triangulation
    a = (0.0, 0.0)
    b = (0.0, 1.0)
    c = (0.0, 2.5)
    d = (2.0, 0.0)
    e = (6.0, 0.0)
    f = (8.0, 0.0)
    g = (8.0, 0.5)
    h = (7.5, 1.0)
    i = (4.0, 1.0)
    j = (4.0, 2.5)
    k = (8.0, 2.5)
    pts = [a, b, c, d, e, f, g, h, i, j, k]
    tri = triangulate(pts; delete_ghosts=false, randomise=false)
    fig, ax, sc = triplot(tri)
    let vert = each_solid_vertex(tri)
        text!(ax, collect(get_point(tri, vert...)); text=string.(vert), fontsize=18)
    end
    lines!(ax, [get_point(tri, 2, 7)...], color=:blue, linewidth=2)
    @test_reference "../docs/src/tri_algs/figs/segment_example.png" fig

    # Deleting a cavity 
    e = (2, 7)
    intersecting_triangles, collinear_segments, left_vertices, right_vertices = DelaunayTriangulation.locate_intersecting_triangles(tri, e)
    DelaunayTriangulation.delete_intersected_triangles!(tri, intersecting_triangles)
    fig, ax, sc = triplot(tri)
    let vert = each_solid_vertex(tri)
        text!(ax, collect(get_point(tri, vert...)); text=string.(vert), fontsize=18)
    end
    lines!(ax, [get_point(tri, 2, 7)...], color=:blue, linewidth=2)
    lines!(ax, [get_point(tri, 10, 9)...], color=:black, linestyle=:dash, linewidth=2)
    @test_reference "../docs/src/tri_algs/figs/segment_example_deleted_triangles.png" fig

    ## Triangulating the polygonal cavities 
    tri = triangulate(pts; delete_ghosts=false, randomise=false)
    add_edge!(tri, 2, 7)
    fig, ax, sc = triplot(tri)
    @test_reference "../docs/src/tri_algs/figs/segment_example_completed.png" fig
end

@testset "Pole of inaccessibility" begin
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
    @test_reference "../docs/src/other_features/figs/pole_of_inaccessibility.png" fig
end

@testset "Querying if points are in a polygon" begin
    C = (15.7109521325776, 33.244486807457)
    D = (14.2705719699703, 32.8530791545746)
    E = (14.3, 27.2)
    F = (14.1, 27.0)
    G = (13.7, 27.2)
    H = (13.4, 27.5)
    I = (13.1, 27.6)
    J = (12.7, 27.4)
    K = (12.5, 27.1)
    L = (12.7, 26.7)
    M = (13.1, 26.5)
    N = (13.6, 26.4)
    O = (14.0, 26.4)
    P = (14.6, 26.5)
    Q = (15.1983491346581, 26.8128534095401)
    R = (15.6, 27.6)
    S = (15.6952958264624, 28.2344688505621)
    T = (17.8088971520274, 33.1192363585346)
    U = (16.3058917649589, 33.0722674401887)
    V = (16.3215480710742, 29.7374742376305)
    W = (16.3841732955354, 29.393035503094)
    Z = (16.6190178872649, 28.9233463196351)
    A1 = (17.0417381523779, 28.5319386667527)
    B1 = (17.5114273358368, 28.3753756055997)
    C1 = (18.1376795804487, 28.3597192994844)
    D1 = (18.7169629067146, 28.5632512789833)
    E1 = (19.2805899268653, 28.8920337074045)
    F1 = (19.26493362075, 28.4536571361762)
    G1 = (20.6426885588962, 28.4223445239456)
    H1 = (20.689657477242, 33.1035800524193)
    I1 = (19.2805899268653, 33.0722674401887)
    J1 = (19.2962462329806, 29.7531305437458)
    K1 = (19.0614016412512, 29.393035503094)
    L1 = (18.7482755189452, 29.236472441941)
    M1 = (18.4508057027546, 29.1425346052493)
    N1 = (18.1689921926793, 29.3147539725175)
    O1 = (17.7932408459121, 29.6278800948235)
    P1 = (22.6466957416542, 35.4207133574833)
    Q1 = (21.2219718851621, 34.9979930923702)
    R1 = (21.2376281912774, 28.4693134422915)
    S1 = (22.6780083538847, 28.4380008300609)
    T1 = (24.5724213938357, 33.1975178891111)
    U1 = (23.3512295168425, 32.8530791545746)
    V1 = (23.3199169046119, 28.4380008300609)
    W1 = (24.6663592305274, 28.3753756055997)
    Z1 = (15.1942940307729, 35.4363696635986)
    A2 = (14.7246048473139, 35.3737444391374)
    B2 = (14.3645098066621, 35.1858687657538)
    C2 = (14.1766341332786, 34.8570863373326)
    D2 = (14.1140089088174, 34.3247719294125)
    E2 = (14.2705719699703, 33.8394264398383)
    F2 = (14.7246048473139, 33.6202381542241)
    G2 = (15.4604512347329, 33.6045818481088)
    H2 = (16.0, 34.0)
    I2 = (15.9771093365377, 34.6848669700643)
    J2 = (15.6170142958859, 35.2328376840997)
    K2 = (24.1653574348379, 35.4520259697138)
    L2 = (23.7739497819555, 35.4363696635986)
    M2 = (23.4608236596496, 35.2641502963303)
    N2 = (23.272947986266, 34.9040552556785)
    O2 = (23.1320412312284, 34.5909291333725)
    P2 = (23.1163849251131, 34.2151777866054)
    Q2 = (23.2886042923813, 33.8081138276077)
    R2 = (23.8209187003014, 33.6045818481088)
    S2 = (24.3062641898756, 33.5576129297629)
    T2 = (24.7602970672192, 33.8550827459536)
    U2 = (25.010797965064, 34.4656786844502)
    V2 = (24.8385785977957, 34.9666804801397)
    W2 = (24.5254524754898, 35.2641502963303)
    Z2 = (25.3708930057158, 37.4716894585871)
    A3 = (24.7916096794498, 37.3464390096648)
    B3 = (24.4471709449133, 36.9550313567823)
    C3 = (24.3062641898756, 36.5636237038999)
    D3 = (24.4941398632592, 35.9999966837492)
    E3 = (25.0264542711793, 35.5929327247515)
    F3 = (25.5587686790994, 35.5929327247515)
    F3 = (25.5587686790994, 35.5929327247515)
    G3 = (26.0, 36.0)
    H3 = (26.1380520053653, 36.5792800100152)
    I3 = (26.0, 37.0)
    J3 = (25.7466443524829, 37.2838137852036)
    K3 = (26.3885529032101, 35.4676822758291)
    L3 = (25.9814889442124, 35.3580881330221)
    M3 = (25.6840191280217, 35.1858687657538)
    N3 = (25.5274560668688, 34.9040552556785)
    O3 = (25.4961434546382, 34.5596165211419)
    P3 = (25.5274560668688, 34.246490398836)
    Q3 = (25.6683628219064, 33.8394264398383)
    R3 = (26.0284578625583, 33.6358944603394)
    S3 = (26.5451159643631, 33.6202381542241)
    T3 = (27.0, 34.0)
    U3 = (27.280962351782, 34.5596165211419)
    V3 = (27.0304614539373, 35.2171813779844)
    W3 = (26.1693646175959, 33.087923746304)
    Z3 = (26.0, 33.0)
    A4 = (25.5274560668688, 32.7278287056522)
    B4 = (25.2612988629087, 32.4147025833463)
    C4 = (25.1830173323322, 32.0702638488098)
    D4 = (25.2299862506781, 31.7727940326191)
    E4 = (25.6527065157911, 31.5222931347744)
    F4 = (26.2946150665183, 31.7258251142732)
    G4 = (26.5607722704784, 32.5086404200381)
    H4 = (27.1557119028596, 32.7434850117675)
    I4 = (27.6097447802033, 32.4929841139228)
    J4 = (27.6410573924338, 32.1015764610403)
    K4 = (27.7193389230103, 31.6005746653509)
    L4 = (27.437525412935, 31.4283552980826)
    M4 = (26.9834925355914, 31.2561359308143)
    N4 = (26.5764285765937, 31.0995728696614)
    O4 = (26.0441141686736, 30.7864467473554)
    P4 = (25.6527065157911, 30.5672584617413)
    Q4 = (25.3239240873699, 30.1915071149741)
    R4 = (25.1673610262169, 29.8783809926682)
    S4 = (25.1047358017558, 29.6122237887082)
    T4 = (25.0890794956405, 29.1895035235952)
    U4 = (25.2926114751393, 28.8294084829433)
    V4 = (25.6840191280217, 28.5632512789833)
    W4 = (26.1537083114806, 28.3753756055997)
    Z4 = (26.8269294744384, 28.391031911715)
    A5 = (27.4844943312809, 28.6102201973292)
    B5 = (27.7342002330051, 28.7239579596219)
    C5 = (27.7264126450755, 28.4202565942047)
    D5 = (29.1825559185446, 28.3922538389457)
    E5 = (29.1545531632856, 32.2146299318021)
    F5 = (29.000538009361, 32.5786657501693)
    G5 = (28.6785063238822, 32.9006974356481)
    H5 = (28.3144705055149, 33.0827153448317)
    I5 = (27.9084305542591, 33.2367304987563)
    J5 = (27.3343740714492, 33.3207387645334)
    K5 = (26.8303244767868, 33.2367304987563)
    L5 = (27.6564057569279, 30.786489413592)
    M5 = (27.6984098898165, 30.3944508399657)
    N5 = (27.6984098898165, 29.7363860913787)
    O5 = (27.5863988687804, 29.4143544059)
    P5 = (27.2643671833016, 29.2043337414573)
    Q5 = (26.9843396307114, 29.1763309861983)
    R5 = (26.6903107004917, 29.3163447624934)
    S5 = (26.5782996794556, 29.7503874690082)
    T5 = (26.7603175886393, 30.3384453294476)
    U5 = (27.3203726938197, 30.7024811478149)
    J_curve = [[C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, C]]
    U_curve = [[T, U, V, W, Z, A1, B1, C1, D1, E1, F1, G1, H1, I1, J1, K1, L1, M1, N1, O1, T]]
    L_curve = [[P1, Q1, R1, S1, P1]]
    I_curve = [[T1, U1, V1, W1, T1]]
    A_curve_outline = [[
        K5, W3, Z3, A4, B4, C4, D4, E4, F4, G4, H4, I4, J4, K4, L4, M4, N4,
        O4, P4, Q4, R4, S4, T4, U4, V4, W4, Z4, A5, B5, C5, D5, E5, F5, G5,
        H5, I5, J5, K5]]
    A_curve_hole = [[L5, M5, N5, O5, P5, Q5, R5, S5, T5, U5, L5]]
    dot_1 = [[Z1, A2, B2, C2, D2, E2, F2, G2, H2, I2, J2, Z1]]
    dot_2 = [[Z2, A3, B3, C3, D3, E3, F3, G3, H3, I3, J3, Z2]]
    dot_3 = [[K2, L2, M2, N2, O2, P2, Q2, R2, S2, T2, U2, V2, W2, K2]]
    dot_4 = [[K3, L3, M3, N3, O3, P3, Q3, R3, S3, T3, U3, V3, K3]]
    curves = [J_curve, U_curve, L_curve, I_curve, A_curve_outline, A_curve_hole, dot_1, dot_2, dot_3, dot_4]
    nodes, points = convert_boundary_points_to_indices(curves)
    xmin, xmax, ymin, ymax = DelaunayTriangulation.polygon_bounds(points, nodes, Val(true)) # Val(true) => check all parts of the polygon
    rng = StableRNG(19299)
    query_points = [((xmax - xmin) * rand(rng) + xmin, (ymax - ymin) * rand(rng) + ymin) for _ in 1:1000]

    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, query_points)
    for nodes in nodes
        lines!(ax, points[reduce(vcat, nodes)], color=:magenta, linewidth=3)
    end

    @test_reference "../docs/src/other_features/figs/scattered_julia.png" fig

    is_inside = [DelaunayTriangulation.distance_to_polygon(q, points, nodes) > 0 for q in query_points]
    scatter!(ax, query_points[is_inside], color=:blue)
    scatter!(ax, query_points[.!is_inside], color=:red)

    @test_reference "../docs/src/other_features/figs/point_in_polygon.png" fig
end

@testset "Convex hull" begin
    rng = StableRNG(1929)
    pts = [Tuple(25randn(rng, 2)) for _ in 1:500]
    ch = convex_hull(pts)
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, pts)
    lines!(ax, pts[ch.indices])
    @test_reference "../docs/src/other_features/figs/convex_hull_1.png" fig
end

@testset "Refinement" begin
    @testset "Unconstrained example" begin
        rng = StableRNG(9282881)
        pts = [(rand(rng), rand(rng)) for _ in 1:50]
        tri = triangulate(pts; rng)
        orig_tri = deepcopy(tri)
        A = get_total_area(tri)
        stats = refine!(tri; rng, min_angle=30.0, max_area=0.01A)
        fig = Figure(fontsize=33)
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", title=L"(a):$ $ Original", width=400, height=400, titlealign=:left)
        triplot!(ax, orig_tri) # orig_tri was deepcopy(triangulate(pts)) from before 
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", title=L"(b):$ $ Refined", width=400, height=400, titlealign=:left)
        triplot!(ax, tri)
        areas = get_all_stat(stats, :area) ./ A
        angles = getindex.(get_all_stat(stats, :angles), 1) # 1 is the smallest
        ax = Axis(fig[2, 1], xlabel=L"A/A(\Omega)", ylabel=L"$ $Count", title=L"(c):$ $ Area histogram", width=400, height=400, titlealign=:left)
        hist!(ax, areas, bins=0:0.001:0.01)
        ax = Axis(fig[2, 2], xlabel=L"\theta_{\min}", ylabel=L"$ $Count", title=L"(d):$ $ Angle histogram", width=400, height=400, titlealign=:left)
        hist!(ax, rad2deg.(angles), bins=20:2:60)
        vlines!(ax, [30.0], color=:red)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/unconstrained_refinement.png" fig
        @test validate_triangulation(tri; check_ghost_triangle_delaunay=false)
        validate_statistics(tri)
    end

    @testset "Square" begin
        p1 = (0.0, 0.0)
        p2 = (1.0, 0.0)
        p3 = (1.0, 0.5)
        p4 = (0.0, 0.5)
        pts = [p1, p2, p3, p4, p1]
        boundary_nodes, points = convert_boundary_points_to_indices(pts)
        C = Set(((2, 4),))
        rng = StableRNG(19281)
        tri = triangulate(points; boundary_nodes, edges=C, rng)
        orig_tri = deepcopy(tri)
        A = get_total_area(tri)
        stats = refine!(tri; rng, min_angle=33.0, max_area=0.001A)
        fig = Figure(fontsize=33)
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", title=L"(a):$ $ Original", width=400, height=400, titlealign=:left)
        triplot!(ax, orig_tri) # orig_tri was deepcopy(triangulate(pts)) from before 
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", title=L"(b):$ $ Refined", width=400, height=400, titlealign=:left)
        triplot!(ax, tri)
        areas = get_all_stat(stats, :area) ./ A
        angles = getindex.(get_all_stat(stats, :angles), 1) # 1 is the smallest
        ax = Axis(fig[2, 1], xlabel=L"A/A(\Omega)", ylabel=L"$ $Count", title=L"(c):$ $ Area histogram", width=400, height=400, titlealign=:left)
        hist!(ax, areas, bins=0:0.0001:0.001)
        ax = Axis(fig[2, 2], xlabel=L"\theta_{\min}", ylabel=L"$ $Count", title=L"(d):$ $ Angle histogram", width=400, height=400, titlealign=:left)
        hist!(ax, rad2deg.(angles), bins=0:2:60)
        vlines!(ax, [33.0], color=:red)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/square_constrained_refinement.png" fig
        @test validate_triangulation(tri; check_ghost_triangle_delaunay=false)
        validate_statistics(tri)
    end

    @testset "Multiply-connected" begin
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
        rng = StableRNG(19191919)
        C = Set{NTuple{2,Int64}}()
        for i in 1:50
            θ = 2π * rand(rng)
            r = 4sqrt(rand(rng))
            x = 20 + r * cos(θ)
            y = 20 + r * sin(θ)
            push!(points, (x, y))
            push!(C, (48, 48 + i))
        end
        tri = triangulate(points; boundary_nodes=nodes, edges=C, rng)
        orig_tri = deepcopy(tri)
        A = get_total_area(tri)
        stats = refine!(tri; max_area=0.001A, min_angle=27.3, rng)

        fig = Figure(fontsize=33)
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", title=L"(a):$ $ Original", width=400, height=400, titlealign=:left)
        triplot!(ax, orig_tri) # orig_tri was deepcopy(triangulate(pts)) from before 
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", title=L"(b):$ $ Refined", width=400, height=400, titlealign=:left)
        triplot!(ax, tri)
        areas = get_all_stat(stats, :area) ./ A
        angles = getindex.(get_all_stat(stats, :angles), 1) # 1 is the smallest
        ax = Axis(fig[2, 1], xlabel=L"A/A(\Omega)", ylabel=L"$ $Count", title=L"(c):$ $ Area histogram", width=400, height=400, titlealign=:left)
        hist!(ax, areas, bins=0:0.0000001:0.000001)
        ax = Axis(fig[2, 2], xlabel=L"\theta_{\min}", ylabel=L"$ $Count", title=L"(d):$ $ Angle histogram", width=400, height=400, titlealign=:left)
        hist!(ax, rad2deg.(angles), bins=0:0.5:40)
        vlines!(ax, [27.3], color=:red)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/mc_constrained_refinement.png" fig
    end

    @testset "Julia logo" begin
        C = (15.7109521325776, 33.244486807457)
        D = (14.2705719699703, 32.8530791545746)
        E = (14.3, 27.2)
        F = (14.1, 27.0)
        G = (13.7, 27.2)
        H = (13.4, 27.5)
        I = (13.1, 27.6)
        J = (12.7, 27.4)
        K = (12.5, 27.1)
        L = (12.7, 26.7)
        M = (13.1, 26.5)
        N = (13.6, 26.4)
        O = (14.0, 26.4)
        P = (14.6, 26.5)
        Q = (15.1983491346581, 26.8128534095401)
        R = (15.6, 27.6)
        S = (15.6952958264624, 28.2344688505621)
        T = (17.8088971520274, 33.1192363585346)
        U = (16.3058917649589, 33.0722674401887)
        V = (16.3215480710742, 29.7374742376305)
        W = (16.3841732955354, 29.393035503094)
        Z = (16.6190178872649, 28.9233463196351)
        A1 = (17.0417381523779, 28.5319386667527)
        B1 = (17.5114273358368, 28.3753756055997)
        C1 = (18.1376795804487, 28.3597192994844)
        D1 = (18.7169629067146, 28.5632512789833)
        E1 = (19.2805899268653, 28.8920337074045)
        F1 = (19.26493362075, 28.4536571361762)
        G1 = (20.6426885588962, 28.4223445239456)
        H1 = (20.689657477242, 33.1035800524193)
        I1 = (19.2805899268653, 33.0722674401887)
        J1 = (19.2962462329806, 29.7531305437458)
        K1 = (19.0614016412512, 29.393035503094)
        L1 = (18.7482755189452, 29.236472441941)
        M1 = (18.4508057027546, 29.1425346052493)
        N1 = (18.1689921926793, 29.3147539725175)
        O1 = (17.7932408459121, 29.6278800948235)
        P1 = (22.6466957416542, 35.4207133574833)
        Q1 = (21.2219718851621, 34.9979930923702)
        R1 = (21.2376281912774, 28.4693134422915)
        S1 = (22.6780083538847, 28.4380008300609)
        T1 = (24.5724213938357, 33.1975178891111)
        U1 = (23.3512295168425, 32.8530791545746)
        V1 = (23.3199169046119, 28.4380008300609)
        W1 = (24.6663592305274, 28.3753756055997)
        Z1 = (15.1942940307729, 35.4363696635986)
        A2 = (14.7246048473139, 35.3737444391374)
        B2 = (14.3645098066621, 35.1858687657538)
        C2 = (14.1766341332786, 34.8570863373326)
        D2 = (14.1140089088174, 34.3247719294125)
        E2 = (14.2705719699703, 33.8394264398383)
        F2 = (14.7246048473139, 33.6202381542241)
        G2 = (15.4604512347329, 33.6045818481088)
        H2 = (16.0, 34.0)
        I2 = (15.9771093365377, 34.6848669700643)
        J2 = (15.6170142958859, 35.2328376840997)
        K2 = (24.1653574348379, 35.4520259697138)
        L2 = (23.7739497819555, 35.4363696635986)
        M2 = (23.4608236596496, 35.2641502963303)
        N2 = (23.272947986266, 34.9040552556785)
        O2 = (23.1320412312284, 34.5909291333725)
        P2 = (23.1163849251131, 34.2151777866054)
        Q2 = (23.2886042923813, 33.8081138276077)
        R2 = (23.8209187003014, 33.6045818481088)
        S2 = (24.3062641898756, 33.5576129297629)
        T2 = (24.7602970672192, 33.8550827459536)
        U2 = (25.010797965064, 34.4656786844502)
        V2 = (24.8385785977957, 34.9666804801397)
        W2 = (24.5254524754898, 35.2641502963303)
        Z2 = (25.3708930057158, 37.4716894585871)
        A3 = (24.7916096794498, 37.3464390096648)
        B3 = (24.4471709449133, 36.9550313567823)
        C3 = (24.3062641898756, 36.5636237038999)
        D3 = (24.4941398632592, 35.9999966837492)
        E3 = (25.0264542711793, 35.5929327247515)
        F3 = (25.5587686790994, 35.5929327247515)
        F3 = (25.5587686790994, 35.5929327247515)
        G3 = (26.0, 36.0)
        H3 = (26.1380520053653, 36.5792800100152)
        I3 = (26.0, 37.0)
        J3 = (25.7466443524829, 37.2838137852036)
        K3 = (26.3885529032101, 35.4676822758291)
        L3 = (25.9814889442124, 35.3580881330221)
        M3 = (25.6840191280217, 35.1858687657538)
        N3 = (25.5274560668688, 34.9040552556785)
        O3 = (25.4961434546382, 34.5596165211419)
        P3 = (25.5274560668688, 34.246490398836)
        Q3 = (25.6683628219064, 33.8394264398383)
        R3 = (26.0284578625583, 33.6358944603394)
        S3 = (26.5451159643631, 33.6202381542241)
        T3 = (27.0, 34.0)
        U3 = (27.280962351782, 34.5596165211419)
        V3 = (27.0304614539373, 35.2171813779844)
        W3 = (26.1693646175959, 33.087923746304)
        Z3 = (26.0, 33.0)
        A4 = (25.5274560668688, 32.7278287056522)
        B4 = (25.2612988629087, 32.4147025833463)
        C4 = (25.1830173323322, 32.0702638488098)
        D4 = (25.2299862506781, 31.7727940326191)
        E4 = (25.6527065157911, 31.5222931347744)
        F4 = (26.2946150665183, 31.7258251142732)
        G4 = (26.5607722704784, 32.5086404200381)
        H4 = (27.1557119028596, 32.7434850117675)
        I4 = (27.6097447802033, 32.4929841139228)
        J4 = (27.6410573924338, 32.1015764610403)
        K4 = (27.7193389230103, 31.6005746653509)
        L4 = (27.437525412935, 31.4283552980826)
        M4 = (26.9834925355914, 31.2561359308143)
        N4 = (26.5764285765937, 31.0995728696614)
        O4 = (26.0441141686736, 30.7864467473554)
        P4 = (25.6527065157911, 30.5672584617413)
        Q4 = (25.3239240873699, 30.1915071149741)
        R4 = (25.1673610262169, 29.8783809926682)
        S4 = (25.1047358017558, 29.6122237887082)
        T4 = (25.0890794956405, 29.1895035235952)
        U4 = (25.2926114751393, 28.8294084829433)
        V4 = (25.6840191280217, 28.5632512789833)
        W4 = (26.1537083114806, 28.3753756055997)
        Z4 = (26.8269294744384, 28.391031911715)
        A5 = (27.4844943312809, 28.6102201973292)
        B5 = (27.7342002330051, 28.7239579596219)
        C5 = (27.7264126450755, 28.4202565942047)
        D5 = (29.1825559185446, 28.3922538389457)
        E5 = (29.1545531632856, 32.2146299318021)
        F5 = (29.000538009361, 32.5786657501693)
        G5 = (28.6785063238822, 32.9006974356481)
        H5 = (28.3144705055149, 33.0827153448317)
        I5 = (27.9084305542591, 33.2367304987563)
        J5 = (27.3343740714492, 33.3207387645334)
        K5 = (26.8303244767868, 33.2367304987563)
        L5 = (27.6564057569279, 30.786489413592)
        M5 = (27.6984098898165, 30.3944508399657)
        N5 = (27.6984098898165, 29.7363860913787)
        O5 = (27.5863988687804, 29.4143544059)
        P5 = (27.2643671833016, 29.2043337414573)
        Q5 = (26.9843396307114, 29.1763309861983)
        R5 = (26.6903107004917, 29.3163447624934)
        S5 = (26.5782996794556, 29.7503874690082)
        T5 = (26.7603175886393, 30.3384453294476)
        U5 = (27.3203726938197, 30.7024811478149)

        J_curve = [[C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, C]]
        U_curve = [[T, U, V, W, Z, A1, B1, C1, D1, E1, F1, G1, H1, I1, J1, K1, L1, M1, N1, O1, T]]
        L_curve = [[P1, Q1, R1, S1, P1]]
        I_curve = [[T1, U1, V1, W1, T1]]
        A_curve_outline = [[
            K5, W3, Z3, A4, B4, C4, D4, E4, F4, G4, H4, I4, J4, K4, L4, M4, N4,
            O4, P4, Q4, R4, S4, T4, U4, V4, W4, Z4, A5, B5, C5, D5, E5, F5, G5,
            H5, I5, J5, K5]]
        A_curve_hole = [[L5, M5, N5, O5, P5, Q5, R5, S5, T5, U5, L5]]
        dot_1 = [[Z1, A2, B2, C2, D2, E2, F2, G2, H2, I2, J2, Z1]]
        dot_2 = [[Z2, A3, B3, C3, D3, E3, F3, G3, H3, I3, J3, Z2]]
        dot_3 = [[K2, L2, M2, N2, O2, P2, Q2, R2, S2, T2, U2, V2, W2, K2]]
        dot_4 = [[K3, L3, M3, N3, O3, P3, Q3, R3, S3, T3, U3, V3, K3]]
        curves = [J_curve, U_curve, L_curve, I_curve, A_curve_outline, A_curve_hole, dot_1, dot_2, dot_3, dot_4]
        nodes, points = convert_boundary_points_to_indices(curves)
        tri = triangulate(points; boundary_nodes=nodes, check_arguments=false)
        orig_tri = deepcopy(tri)
        A = get_total_area(tri)
        stats = refine!(tri; min_angle=26.45, max_area=0.005A / 9)
        @test validate_triangulation(tri; check_planarity=false, check_ghost_triangle_orientation=false, check_ghost_triangle_delaunay=false)
        validate_statistics(tri)

        fig = Figure(fontsize=33)
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", title=L"(a):$ $ Original", width=400, height=400, titlealign=:left)
        triplot!(ax, orig_tri, show_convex_hull=false) # orig_tri was deepcopy(triangulate(pts)) from before 
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", title=L"(b):$ $ Refined", width=400, height=400, titlealign=:left)
        triplot!(ax, tri, show_convex_hull=false)
        areas = get_all_stat(stats, :area) ./ (0.005A)
        angles = getindex.(get_all_stat(stats, :angles), 1) # 1 is the smallest
        ax = Axis(fig[2, 1], xlabel=L"A/A(\Omega)", ylabel=L"$ $Count", title=L"(c):$ $ Area histogram", width=400, height=400, titlealign=:left)
        hist!(ax, areas, bins=0:0.001:0.1)
        ax = Axis(fig[2, 2], xlabel=L"\theta_{\min}", ylabel=L"$ $Count", title=L"(d):$ $ Angle histogram", width=400, height=400, titlealign=:left)
        hist!(ax, rad2deg.(angles), bins=0:0.2:40)
        vlines!(ax, [26.45], color=:red)
        resize_to_layout!(fig)
        @test_reference "../docs/src/triangulations/figs/julia_constrained_refinement.png" fig
    end
end

if !(get(ENV, "CI", "false") == "true")
    @testset "Tasmania" begin
        tassy = readdlm("tassy.txt")
        ymax = @views maximum(tassy[:, 2])
        tassy = [(x, ymax - y) for (x, y) in eachrow(tassy)]
        reverse!(tassy)
        unique!(tassy)
        push!(tassy, tassy[begin])
        boundary_nodes, points = convert_boundary_points_to_indices(tassy)
        rng = StableRNG(18181)
        tri = triangulate(points; boundary_nodes=boundary_nodes, rng)
        orig_tri = deepcopy(tri)
        A = get_total_area(tri)
        stats = refine!(tri; max_area=1e-3A, rng)
        fig = Figure(fontsize=33)
        ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", title=L"(a):$ $ Original", width=400, height=400, titlealign=:left)
        triplot!(ax, orig_tri, show_convex_hull=false) # orig_tri was deepcopy(triangulate(pts)) from before 
        ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"y", title=L"(b):$ $ Refined", width=400, height=400, titlealign=:left)
        triplot!(ax, tri, show_convex_hull=false)
        areas = get_all_stat(stats, :area) ./ (1e-3A)
        angles = getindex.(get_all_stat(stats, :angles), 1) # 1 is the smallest
        ax = Axis(fig[2, 1], xlabel=L"A/A(\Omega)", ylabel=L"$ $Count", title=L"(c):$ $ Area histogram", width=400, height=400, titlealign=:left)
        hist!(ax, areas, bins=0:0.01:1)
        ylims!(ax, 0, 1000)
        ax = Axis(fig[2, 2], xlabel=L"\theta_{\min}", ylabel=L"$ $Count", title=L"(d):$ $ Angle histogram", width=400, height=400, titlealign=:left)
        hist!(ax, rad2deg.(angles), bins=0:0.2:60)
        vlines!(ax, [30.0], color=:red)
        resize_to_layout!(fig)
        ylims!(ax, 0, 1000)
        @test_reference "../docs/src/triangulations/figs/tassy_constrained_refinement.png" fig
    end
end

@testset "Custom Inteface" begin
    ## Struct definitions
    struct CustomPoint
        x::Float64
        y::Float64
    end
    struct CustomPoints
        points::Vector{CustomPoint}
    end
    struct CustomEdge
        i::Int32
        j::Int32
    end
    struct CustomEdges
        edges::Set{CustomEdge}
    end
    struct CustomTriangle
        i::Int32
        j::Int32
        k::Int32
    end
    struct CustomTriangles
        triangles::Vector{CustomTriangle}
    end
    struct CustomPolygonSegment
        edges::Vector{CustomEdge}
    end
    struct CustomPolygon
        segments::Vector{CustomPolygonSegment}
    end
    struct CustomPolygons{N}
        polygons::NTuple{N,CustomPolygon}
    end

    ### Defining the methods
    ## Definitions needed for unconstrained triangulations
    # Point
    DT.getx(p::CustomPoint) = p.x
    DT.gety(p::CustomPoint) = p.y
    DT.number_type(::Type{CustomPoint}) = Float64

    # Points
    DT.each_point_index(pts::CustomPoints) = eachindex(pts.points)
    DT.num_points(pts::CustomPoints) = length(pts.points)
    DT.each_point(pts::CustomPoints) = pts.points
    DT.number_type(::Type{CustomPoints}) = DT.number_type(CustomPoint)
    DT.getpoint(pts::CustomPoints, i::Integer) = pts.points[i]

    # Edge
    DT.construct_edge(::Type{CustomEdge}, i, j) = CustomEdge(i, j)
    DT.initial(e::CustomEdge) = e.i
    DT.terminal(e::CustomEdge) = e.j

    # Edges
    DT.initialise_edges(::Type{CustomEdges}) = CustomEdges(Set{CustomEdge}())
    DT.add_to_edges!(edges::CustomEdges, edge::CustomEdge) = push!(edges.edges, edge)
    Base.iterate(edges::CustomEdges, state...) = Base.iterate(edges.edges, state...)
    DT.each_edge(edges::CustomEdges) = edges.edges
    DT.delete_from_edges!(edges::CustomEdges, edge::CustomEdge) = delete!(edges.edges, edge)
    DT.contains_edge(edge::CustomEdge, edges::CustomEdges) = edge ∈ edges.edges
    DT.num_edges(edges::CustomEdges) = length(edges.edges)

    # Triangle
    DT.construct_triangle(::Type{CustomTriangle}, i, j, k) = CustomTriangle(i, j, k)
    DT.geti(tri::CustomTriangle) = tri.i
    DT.getj(tri::CustomTriangle) = tri.j
    DT.getk(tri::CustomTriangle) = tri.k
    DT.integer_type(::Type{CustomTriangle}) = Int32

    # Triangles
    function DT.delete_from_triangles!(tri::CustomTriangles, triangle::CustomTriangle)
        i = findfirst(==(triangle), tri.triangles)
        deleteat!(tri.triangles, i)
        return nothing
    end
    DT.initialise_triangles(::Type{CustomTriangles}) = CustomTriangles(Vector{CustomTriangle}())
    DT.sizehint!(tri::CustomTriangles, n) = sizehint!(tri.triangles, n)
    DT.triangle_type(::Type{CustomTriangles}) = CustomTriangle
    DT.add_to_triangles!(tri::CustomTriangles, triangle::CustomTriangle) = push!(tri.triangles, triangle)
    DT.num_triangles(tri::CustomTriangles) = length(tri.triangles)
    Base.iterate(tri::CustomTriangles, state...) = Base.iterate(tri.triangles, state...)
    Base.in(T, V::CustomTriangles) = T ∈ V.triangles

    ## Definitions needed for constrained segments 
    # Edge 
    DT.integer_type(::Type{CustomEdge}) = Int32

    # Edges 
    DT.edge_type(::Type{CustomEdges}) = CustomEdge

    ## Definitions needed for boundary nodes 
    # Triangles 
    DT.each_triangle(tri::CustomTriangles) = tri.triangles

    # Boundary Nodes 
    DT.has_multiple_curves(::CustomPolygons{N}) where {N} = N > 1
    DT.has_multiple_curves(::CustomPolygon) = false
    DT.has_multiple_curves(::CustomPolygonSegment) = false
    DT.has_multiple_segments(::CustomPolygons) = true
    DT.has_multiple_segments(poly::CustomPolygon) = true
    DT.has_multiple_segments(seg::CustomPolygonSegment) = false
    DT.num_curves(::CustomPolygons{N}) where {N} = N
    DT.getboundarynodes(poly::CustomPolygons, m) = poly.polygons[m] # go down to the mth polygon 
    DT.getboundarynodes(poly::CustomPolygon, m) = poly.segments[m] # go down to the mth segment
    DT.getboundarynodes(seg::CustomPolygonSegment, m) = m > length(seg.edges) ? DT.terminal(seg.edges[m-1]) : DT.initial(seg.edges[m]) # go down to the mth edge and extract the left node
    DT.getboundarynodes(poly::CustomPolygons, (m, n)::NTuple{2,Int32}) = DT.getboundarynodes(DT.getboundarynodes(poly, m), n)
    DT.num_segments(poly::CustomPolygon) = length(poly.segments)
    DT.num_boundary_edges(seg::CustomPolygonSegment) = length(seg.edges)

    ## Definitions needed for refinement 
    # Points 
    DT.push_point!(pts::CustomPoints, x, y) = push!(pts.points, CustomPoint(x, y))
    DT.pop_point!(pts::CustomPoints) = pop!(pts.points)

    # Edges 
    DT.contains_edge(e::CustomEdge, Es::CustomEdges) = e ∈ Es.edges

    # Boundary Nodes 
    function Base.insert!(seg::CustomPolygonSegment, index, node)
        cur_edge = seg.edges[index-1]
        u, v = DT.initial(cur_edge), DT.terminal(cur_edge)
        seg.edges[index-1] = CustomEdge(u, node)
        insert!(seg.edges, index, CustomEdge(node, v))
        return nothing
    end

    ### Example
    ## Build the points 
    p1 = CustomPoint(0.0, 0.0)
    p2 = CustomPoint(1.0, 0.0)
    p3 = CustomPoint(1.0, 1.0)
    p4 = CustomPoint(0.0, 1.0)
    p5 = CustomPoint(0.5, 0.5)
    p6 = CustomPoint(0.25, 0.25)
    p7 = CustomPoint(0.75, 0.25)
    p8 = CustomPoint(0.75, 0.75)
    p9 = CustomPoint(0.25, 0.75)
    p10 = CustomPoint(2.3, 5.5)
    p11 = CustomPoint(-0.5, 2.3)
    points = CustomPoints([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11])

    ## The edges
    edges = CustomEdges(Set{CustomEdge}((CustomEdge(2, 7), CustomEdge(8, 3))))

    ## The boundary nodes
    outer_polygon = CustomPolygon([
        CustomPolygonSegment([CustomEdge(1, 2), CustomEdge(2, 3)]),
        CustomPolygonSegment([CustomEdge(3, 4), CustomEdge(4, 1)]),
    ])
    inner_polygon = CustomPolygon([
        CustomPolygonSegment([CustomEdge(6, 9), CustomEdge(9, 8), CustomEdge(8, 7), CustomEdge(7, 6)]),
    ])
    polygons = CustomPolygons((outer_polygon, inner_polygon))

    ## Triangulate
    rng = StableRNG(982381238)
    tri = triangulate(points; edges=edges, boundary_nodes=polygons,
        IntegerType=Int32,
        EdgeType=CustomEdge,
        TriangleType=CustomTriangle,
        EdgesType=CustomEdges,
        TrianglesType=CustomTriangles,
        randomise=true,
        delete_ghosts=false,
        rng
    )

    ## Refine 
    A = get_total_area(tri)
    stats = refine!(tri; max_area=1e-3A, min_angle=27.3, rng)
    fig, ax, sc = triplot(tri; show_convex_hull=false)
    xlims!(ax, 0, 1)
    ylims!(ax, 0, 1)
    @test validate_triangulation(tri)
    validate_statistics(tri)
    @test_reference "../docs/src/interface/figs/interface_example.png" fig
end