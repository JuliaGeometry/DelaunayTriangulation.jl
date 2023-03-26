using ExactPredicates
using ..DelaunayTriangulation
using StaticArrays
const DT = DelaunayTriangulation

include("../helper_functions.jl")

@testset "Opposite signs" begin
    @test DT.opposite_signs(1, -1)
    @test DT.opposite_signs(-1, 1)
    @test !DT.opposite_signs(0, 0)
    @test !DT.opposite_signs(0, 1)
    @test !DT.opposite_signs(1, 0)
    @test !DT.opposite_signs(0, -1)
    @test !DT.opposite_signs(-1, 0)
    @test !DT.opposite_signs(1, 1)
    @test !DT.opposite_signs(-1, -1)
end

@testset "Basic predicates" begin
    @testset "Random points" begin
        for _ in 1:1500
            p, q, r = eachcol(rand(2, 3))
            @test DT.orient_predicate(p, q, r) == orient(p, q, r)
            @inferred DT.orient_predicate(p, q, r)

            a, b, c, p = eachcol(rand(2, 4))
            @test DT.incircle_predicate(a, b, c, p) == incircle(a, b, c, p)
            @inferred DT.incircle_predicate(a, b, c, p)

            p, a, b = eachcol(rand(2, 3))
            @test DT.sameside_predicate(a, b, p) == sameside(p, a, b)
            @inferred DT.sameside_predicate(a, b, p)

            p, q, a, b = eachcol(rand(2, 4))
            @test DT.meet_predicate(p, q, a, b) == meet(p, q, a, b)
            @inferred DT.meet_predicate(p, q, a, b)
        end
    end

    @testset "Collinearity" begin
        x = [float(i) for i in 0:5, _ in 0:5]
        y = [float(j) for _ in 0:5, j in 0:5]
        for _ in 1:15000
            p = (rand(x), rand(y))
            q = (rand(x), rand(y))
            r = (rand(x), rand(y))
            s = (rand(x), rand(y))
            a = (rand(x), rand(y))
            b = (rand(x), rand(y))
            c = (rand(x), rand(y))
            @test DT.orient_predicate(p, q, r) == orient(p, q, r)
            @test DT.incircle_predicate(a, b, c, p) == incircle(a, b, c, p)
            @test DT.sameside_predicate(a, b, p) == sameside(p, a, b)
            @test DT.meet_predicate(p, q, a, b) == meet(p, q, a, b)
        end
    end

    @testset "Extra orient_predicate tests" begin
        p₁ = [5.7044025422189, 1.801603986463]
        p₂ = [8.3797127128527, 5.8924221838871]
        p₃ = [2.8875415689061, 6.2038339497809]
        @test DT.orient_predicate(p₁, p₂, p₃) == 1
        @test DT.orient_predicate(p₁, p₂, [10.0, 4.0]) == -1
        p₁ = [5.0, 1.0]
        p₂ = [5.0, 6.0]
        @test DT.orient_predicate(p₁, p₂, [5.0, 5.0]) == 0
        @test DT.orient_predicate(p₁, p₂, [5.0, 2.0]) == 0
        @test DT.orient_predicate(p₂, p₁, [5.0, 2.0]) == 0
    end
end

global p1, q1, r1 = (3.0, 1.5), (4.5, 1.5), (4.0, 2.0) # +
global p2, q2, r2 = (3.0, 1.5), (4.5, 1.5), (3.4, 1.2) # - 
global p3, q3, r3 = (2.8, 1.4), (2.8, 2.4), (4.6, 1.4) # - 
global p4, q4, r4 = (2.8, 1.4), (3.8, 1.4), (4.6, 1.4) # 0 
global p5, q5, r5 = (5.0, 1.4), (3.8, 1.4), (4.6, 1.4) # 0 
global p6, q6, r6 = (5.0, 1.4), (3.8, 1.4), (4.4, 0.8) # +
global pqr = ((p1, q1, r1), (p2, q2, r2), (p3, q3, r3),
    (p4, q4, r4), (p5, q5, r5), (p6, q6, r6),
    (p4, q4, r4), (p5, q5, r5), (p6, q6, r6))

@testset "triangle_orientation" begin
    results = [Certificate.PositivelyOriented,
        Certificate.NegativelyOriented,
        Certificate.NegativelyOriented,
        Certificate.Degenerate,
        Certificate.Degenerate,
        Certificate.PositivelyOriented]
    for ((p, q, r), result) in zip(pqr, results)
        @test DT.triangle_orientation(p, q, r) == result
        @inferred DT.triangle_orientation(p, q, r)
    end
    @test DT.is_degenerate(DT.triangle_orientation(p1, p1, q1))
    @test DT.is_degenerate(DT.triangle_orientation(p1, q1, q1))
    @test DT.is_degenerate(DT.triangle_orientation(q1, p1, q1))
    @inferred DT.is_degenerate(DT.triangle_orientation(q1, p1, q1))
end

@testset "point_position_relative_to_circle" begin
    R = 5
    a, b, c = (0.0, 5.0), (-3.0, -4.0), (3.0, 4.0)
    p1, p2, p3, p4, p5 = [(0.0, -5.0),
        (5.0, 0.0),
        (-5.0, 0.0),
        (-3.0, 4.0),
        (3.0, -4.0)]
    for p in (p1, p2, p3, p4, p5)
        cert = DT.point_position_relative_to_circle(a, b, c, p)
        @test DT.is_on(cert)
        @inferred DT.point_position_relative_to_circle(a, b, c, p)
    end
    a, b, c = (3.0, 2.5), (7.5, 2.5), (3.0, 4.5)
    p1, p2, p3 = (4.6438, 4.921), (4.48162, 4.2045), (6.078, 3.17) # in 
    p4, p5, p6 = (2.0, 5.5), (3.13, 0.455), (7.5, -0.5), (7.85, 6.066) # out 
    p7 = (7.5, 4.5)
    for p in (p1, p2, p3)
        cert = DT.point_position_relative_to_circle(a, b, c, p)
        @test DT.is_inside(cert)
        @inferred DT.point_position_relative_to_circle(a, b, c, p)
    end
    for p in (p4, p5, p6)
        cert = DT.point_position_relative_to_circle(a, b, c, p)
        @test DT.is_outside(cert)
        @inferred DT.point_position_relative_to_circle(a, b, c, p)
        @inferred DT.is_outside(cert)
    end
    cert = DT.point_position_relative_to_circle(a, b, c, p7)
    @test DT.is_on(cert)
    a, b, c = (3.439, 4.187), (7.74, 1.36), (7.0, 5.0)
    cert1 = DT.point_position_relative_to_circle(a, b, c, (4.07, 3.98))
    cert2 = DT.point_position_relative_to_circle(a, b, c, (5.5, 3.5))
    cert3 = DT.point_position_relative_to_circle(a, b, c, (4.39, 5.73))
    cert4 = DT.point_position_relative_to_circle(a, b, c, (6.838, 6.049))
    cert5 = DT.point_position_relative_to_circle(a, b, c, (2.36, 0.9767))
    @test DT.is_inside(cert1)
    @test DT.is_inside(cert2)
    @test DT.is_outside(cert3)
    @test DT.is_outside(cert4)
    @test DT.is_outside(cert5)
    @test DT.is_on(DT.point_position_relative_to_circle(a, b, c, a))
    @test DT.is_on(DT.point_position_relative_to_circle(a, b, c, b))
    @test DT.is_on(DT.point_position_relative_to_circle(a, b, c, c))
    @inferred DT.is_on(DT.point_position_relative_to_circle(a, b, c, c))

    p0 = Float64[5, 5]
    p1 = Float64[4.5, 2.5]
    p2 = Float64[2.5, 1.5]
    p3 = Float64[3, 3.5]
    p4 = Float64[0, 2]
    p5 = Float64[1, 5]
    p6 = Float64[1, 3]
    p7 = Float64[4, -1]
    p8 = Float64[-1, 4]
    pts = [p0, p1, p2, p3, p4, p5, p6, p7, p8]
    @test DT.is_inside(DT.point_position_relative_to_circle(get_point(pts, 5, 7, 6, 9)...))
    @test DT.is_outside(DT.point_position_relative_to_circle(get_point(pts, 5, 7, 6, 3)...))
    @test DT.is_outside(DT.point_position_relative_to_circle(get_point(pts, 5, 7, 6, 3)...))
    @test DT.is_on(DT.point_position_relative_to_circle(get_point(pts, 5, 7, 6, 6)...))
    @test DT.is_inside(DT.point_position_relative_to_circle(get_point(pts, 3, 2, 1, 4)...))
    @test DT.is_inside(DT.point_position_relative_to_circle(get_point(pts, 3, 2, 1, 6)...))
    @test DT.is_inside(DT.point_position_relative_to_circle(get_point(pts, 3, 2, 1, 7)...))
    @test DT.is_outside(DT.point_position_relative_to_circle(get_point(pts, 3, 2, 1, 5)...))
    @test DT.is_outside(DT.point_position_relative_to_circle(get_point(pts, 3, 2, 1, 8)...))
    @test DT.is_inside(DT.point_position_relative_to_circle(p4, p6, p5, p8))
end

@testset "point_position_relative_to_line" begin
    results = [Certificate.Left,
        Certificate.Right,
        Certificate.Right,
        Certificate.Collinear,
        Certificate.Collinear,
        Certificate.Left]
    for ((p, q, r), result) in zip(pqr, results)
        @test DT.point_position_relative_to_line(p, q, r) == result
        @inferred DT.point_position_relative_to_line(p, q, r)
    end
    p, q = (5.0, 3.0), (10.0, 3.0)
    @test DT.is_collinear(DT.point_position_relative_to_line(p, q, p))
    @inferred DT.is_collinear(DT.point_position_relative_to_line(p, q, p))
    @test DT.is_collinear(DT.point_position_relative_to_line(p, q, q))
    r = (0.0, 3.0)
    @test DT.is_collinear(DT.point_position_relative_to_line(p, q, r))
    @test DT.is_collinear(DT.point_position_relative_to_line(q, p, r))
    p, q = (3.0, 5.0), (3.0, 3.5)
    a, b, c = (2.77, 4.53), (3.158, 4.4949), (3.0, 4.0)
    @test DT.is_right(DT.point_position_relative_to_line(p, q, a))
    @test DT.is_left(DT.point_position_relative_to_line(p, q, b))
    @test DT.is_collinear(DT.point_position_relative_to_line(p, q, c))

    A = [4.6, 3.2]
    B = [3.2, 2.2]
    C = [3.4, 3.2]
    @test DT.is_left(DT.point_position_relative_to_line(C, B, A))
    @test DT.is_right(DT.point_position_relative_to_line(C, A, B))
    C = [5.8, 3.6]
    @test DT.is_right(DT.point_position_relative_to_line(C, B, A))
    A = [1.0, 7.0]
    B = [1.0, 1.0]
    C = [1.0, 5.0]
    @test DT.is_collinear(DT.point_position_relative_to_line(C, B, A))
    @test DT.is_collinear(DT.point_position_relative_to_line(C, A, B))
    A = [2.123933267613, 7.1892809338214]
    B = [-1.5542939635314, 3.3935384556756]
    C = [2.8172732249214, 5.085758012496]
    @test DT.is_right(DT.point_position_relative_to_line(C, B, A))
    @test DT.is_left(DT.point_position_relative_to_line(C, A, B))
    C = [-2.8172732249214, 5.085758012496]
    @test DT.is_left(DT.point_position_relative_to_line(C, B, A))
    @test DT.is_right(DT.point_position_relative_to_line(C, A, B))
    D = [-6.3, -2.77]
    E = [-3.46, -2.07]
    F = [-5.22, -1.41]
    @test DT.is_left(DT.point_position_relative_to_line(D, E, F))
    @test DT.is_right(DT.point_position_relative_to_line(E, D, F))
    F = [-4.88, -4.57]
    @test DT.is_right(DT.point_position_relative_to_line(D, E, F))
    @test DT.is_left(DT.point_position_relative_to_line(E, D, F))
end

@testset "point_position_on_line_segment" begin
    p, q = (3.0, 3.0), (5.0, 3.0)
    c1, c2, c3 = (4.0, 3.0), (2.0, 3.0), (6.0, 3.0) # in, out, out 
    @test DT.is_on(DT.point_position_on_line_segment(p, q, c1))
    @test DT.is_left(DT.point_position_on_line_segment(p, q, c2))
    @test DT.is_right(DT.point_position_on_line_segment(p, q, c3))
    @inferred DT.is_right(DT.point_position_on_line_segment(p, q, c3))
    @test DT.is_degenerate(DT.point_position_on_line_segment(p, q, p))
    @test DT.is_degenerate(DT.point_position_on_line_segment(p, q, q))
    p, q = (4.0, 5.0), (4.0, 1.50)
    c1, c2, c3, c4, c5, c6, c7, c8 = (4.0, 3.0), (4.0, 2.5), (4.0, 4.0), (4.0, 4.5), (4.0, 1.0),
    (4.0, 6.0), (4.0, 6.5), (4.0, 0.5)
    (4.0, 6.0), (4.0, 6.5), (4.0, 0.5)
    @test DT.is_on(DT.point_position_on_line_segment(p, q, c1))
    @test DT.is_on(DT.point_position_on_line_segment(p, q, c2))
    @test DT.is_on(DT.point_position_on_line_segment(p, q, c3))
    @test DT.is_on(DT.point_position_on_line_segment(p, q, c4))
    @test DT.is_right(DT.point_position_on_line_segment(p, q, c5))
    @test DT.is_left(DT.point_position_on_line_segment(p, q, c6))
    @test DT.is_left(DT.point_position_on_line_segment(p, q, c7))
    @test DT.is_right(DT.point_position_on_line_segment(p, q, c8))
    @test DT.is_degenerate(DT.point_position_on_line_segment(p, q, p))
    @test DT.is_degenerate(DT.point_position_on_line_segment(p, q, q))
    @inferred DT.is_degenerate(DT.point_position_on_line_segment(p, q, q))
end

@testset "line_segment_intersection_type" begin
    p1, q1, a1, b1 = (3.192, 5.59), (5.2, 3.6), (8.2, 5.39), (2.19, 4.34)   # one point 
    p2, q2, a2, b2 = (3.192, 5.59), (5.2, 3.6), (8.2, 5.39), (3.5, 7.0)     # no point 
    p3, q3, a3, b3 = (8.0, 4.5), (8.0, 2.5), (8.0, 5.0), (8.0, 7.5)         # no point 
    p4, q4, a4, b4 = (8.0, 8.0), (8.0, 2.5), (8.0, 5.0), (8.0, 7.5)         # multiple  point 
    p5, q5, a5, b5 = (9.5, 6.0), (7.0, 6.0), (8.0, 5.0), (8.0, 7.5)         # one point 
    p6, q6, a6, b6 = (7.0, 5.0), (9.5, 5.0), (8.0, 5.0), (8.0, 7.5)         # multiple point 
    p7, q7, a7, b7 = (7.0, 5.0), (9.5, 5.0), (8.0, 6.0), (8.0, 7.5)         # no point 
    p8, q8, a8, b8 = (8.0, 5.0), (9.5, 5.0), (8.0, 5.0), (8.0, 7.5)         # on point 
    p9, q9, a9, b9 = (8.0, 5.0), (8.0, 7.5), (8.0, 5.0), (8.0, 7.5)         # multiple point 
    p10, q10, a10, b10 = (6.0, 6.5), (8.0, 6.5), (8.0, 5.0), (8.0, 7.5)     # on point
    p11, q11, a11, b11 = (6.0, 6.5), (8.0, 5.5), (8.0, 5.0), (8.0, 7.5)     # on point
    p12, q12, a12, b12 = (6.0, 6.5), (9.947, 8.137), (8.0, 5.0), (8.0, 7.5) # one point
    p13, q13, a13, b13 = (6.0, 6.5), (9.947, 8.137), (8.0, 5.0), (6.0, 7.0) # one point
    results = [Certificate.Single,
        Certificate.None,
        Certificate.None,
        Certificate.Multiple,
        Certificate.Single,
        Certificate.Touching,
        Certificate.None,
        Certificate.Touching,
        Certificate.Multiple,
        Certificate.Touching,
        Certificate.Touching,
        Certificate.Single,
        Certificate.Single]
    pqab = ((p1, q1, a1, b1), (p2, q2, a2, b2), (p3, q3, a3, b3), (p4, q4, a4, b4),
        (p5, q5, a5, b5), (p6, q6, a6, b6), (p7, q7, a7, b7), (p8, q8, a8, b8),
        (p9, q9, a9, b9), (p10, q10, a10, b10), (p11, q11, a11, b11), (p12, q12, a12, b12),
        (p13, q13, a13, b13))
    for ((p, q, a, b), result) in zip(pqab, results)
        @test DT.line_segment_intersection_type(p, q, a, b) == result
        @inferred DT.line_segment_intersection_type(p, q, a, b)
    end
end

@testset "point_position_relative_to_triangle" begin
    p, q, r = (2.0, 1.0), (5.0, 1.0), (2.0, 5.0)
    c1, c2, c3, c4 = (2.57, 3.35), (2.45, 2.51), (3.11, 2.43), (3.33, 1.83) # inside 
    d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12 = (1.76, 4.06), (1.24, 4.109),
    (2.589, 0.505), (3.19, -0.186),
    (3.875, 4.53), (4.8, 2.698), (2.0, 6.5),
    (2.0, 0.0),
    (0.5, 1.0), (6.0, 1.0), (6.0, 0.0),
    (1.0, 6.0) # outside
    (2.589, 0.505), (3.19, -0.186),
    (3.875, 4.53), (4.8, 2.698), (2.0, 6.5),
    (2.0, 0.0),
    (0.5, 1.0), (6.0, 1.0), (6.0, 0.0),
    (1.0, 6.0) # outside
    e1, e2, e3, e4, e5, e6, e7, e8, e9, e10 = (2.0, 3.0), (2.0, 2.5), (2.0, 4.5), (2.0, 1.5),
    (2.5, 1.0), (3.5, 1.0), (4.5, 1.0), p, q, r # on 
    (2.5, 1.0), (3.5, 1.0), (4.5, 1.0), p, q, r # on 
    for c in (c1, c2, c3, c4)
        @test DT.point_position_relative_to_triangle(p, q, r, c) == Certificate.Inside
        @inferred DT.point_position_relative_to_triangle(p, q, r, c)
    end
    for d in (d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12)
        @test DT.point_position_relative_to_triangle(p, q, r, d) == Certificate.Outside
        @inferred DT.point_position_relative_to_triangle(p, q, r, d)
    end
    for e in (e1, e2, e3, e4, e5, e6, e7, e8, e9, e10)
        @test DT.point_position_relative_to_triangle(p, q, r, e) == Certificate.On
        @inferred DT.point_position_relative_to_triangle(p, q, r, e)
    end

    p1 = Float64[5, 5]
    p2 = Float64[4.5, 2.5]
    p3 = Float64[2.5, 1.5]
    p4 = Float64[3, 3.5]
    p5 = Float64[0, 2]
    p6 = Float64[1, 5]
    p7 = Float64[1, 3]
    p8 = Float64[4, -1]
    p9 = Float64[-1, 4]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
    T1 = DT.construct_triangle(NTuple{3,Int64}, 4, 1, 6)
    T2 = DT.construct_triangle(NTuple{3,Int64}, 4, 2, 1)
    T3 = DT.construct_triangle(NTuple{3,Int64}, 3, 2, 4)
    T4 = DT.construct_triangle(NTuple{3,Int64}, 8, 1, 2)
    T5 = DT.construct_triangle(NTuple{3,Int64}, 8, 2, 3)
    T6 = DT.construct_triangle(NTuple{3,Int64}, 8, 3, 5)
    T7 = DT.construct_triangle(NTuple{3,Int64}, 5, 3, 7)
    T8 = DT.construct_triangle(NTuple{3,Int64}, 3, 4, 7)
    T9 = DT.construct_triangle(NTuple{3,Int64}, 5, 7, 9)
    T10 = DT.construct_triangle(NTuple{3,Int64}, 7, 6, 9)
    T11 = DT.construct_triangle(NTuple{3,Int64}, 7, 4, 6)
    T = [T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11]
    pt1 = (
        A=(3.9551298987095, 4.7489988803935),
        B=(3.2585303811361, 4.3003415639903),
        C=(2.69180534989, 4.4066025073489),
        D=(3.2231100666833, 4.6781582514877),
        E=(2.0424329182538, 4.8198395092992)
    )
    pt2 = (
        A=(3.5, 3.5),
        B=(3.8384340736083, 3.2777159548861),
        C=(4.173119090818, 3.0606229707501),
        D=(4.0736181397556, 3.8475850382431),
        E=(4.390212074954, 3.6938108411468),
        F=(4.390212074954, 4.2546343834981),
        G=(4.7520337151806, 4.2998620885264),
        H=(4.7520337151806, 4.661683728753)
    )
    pt3 = (
        A=(3.114790793155, 2.9520764786821),
        B=(3.3952025643307, 2.8163933635972),
        C=(3.3137926952797, 2.4545717233705),
        D=(2.951971055053, 2.3098430672799),
        E=(2.9157888910304, 2.4726628053818),
        F=(3.0786086291324, 2.6535736254952),
        G=(3.901752860648, 2.6626191665008),
        H=(3.0, 2.0),
        I=(2.7258325299114, 1.8213838529739)
    )
    pt4 = (
        A=(4.5781396673455, 2.7004728027825),
        B=(4.6264236360138, 2.9231155471972),
        C=(4.6693427192745, 3.1618529478347),
        D=(4.70153203172, 3.3871781349531),
        E=(4.5325381413811, 2.437593417811),
        F=(4.4708419591939, 2.1988560171735),
        G=(4.4520648602673, 2.0352270122422),
        H=(4.3501320375233, 1.3082850395147)
    )
    pt5 = (
        A=(3.433718968984, 1.1294534677817),
        B=(3.7382811110409, 1.4137114670348),
        C=(3.8499538964617, 0.9162599683419),
        D=(3.6875207540314, 0.5304812550699),
        E=(3.0377881843101, 1.2614303960064),
        F=(3.7484331824428, -0.0481868148381),
        G=(4.0, 2.0)
    )
    pt6 = (
        A=(2.8956591846836, 0.4086563982472),
        B=(2.4286639001964, 0.90610789694),
        C=(2.0936455439339, 1.2309741818007),
        D=(1.5149774740259, 1.4035593956329),
        E=(0.824636618697, 1.586296680867),
        F=(1.3322401887918, 1.2715824674082),
        G=(1.6063461166429, 1.0177806823609),
        H=(2.0225810441206, 0.713218540304),
        I=(2.3169911147756, 0.4391126124529),
        J=(2.6317053282343, 0.1954628988074),
        K=(3.1697651125348, -0.1395554574552)
    )
    pt7 = (
        A=(1.0581342609406, 2.2766375361959),
        B=(0.9363094041178, 2.550743464047),
        C=(1.2916319031842, 2.3070937504015),
        D=(1.7281709734657, 1.9923795369428),
        E=(1.0, 2.0),
        F=(0.5809869050515, 2.1142043937655)
    )
    pt8 = (
        A=(1.5454336882315, 2.845153534702),
        B=(1.9921248299149, 2.5405913926451),
        C=(2.1545579723453, 2.936522177319),
        D=(2.1951662579528, 2.4187665358224),
        E=(2.4185118287945, 2.8756097489077),
        F=(2.4185118287945, 2.5101351784394),
        G=(2.4489680430002, 2.0431398939523),
        H=(2.6317053282343, 2.9162180345153)
    )
    pt9 = (
        A=(-0.5458930205588, 3.5557985328346),
        B=(0.0733833349568, 3.2512363907778),
        C=(0.3170330486022, 2.9771304629266),
        D=(0.0835354063587, 2.6421121066641),
        E=(0.0, 2.4187665358224),
        F=(0.3576413342098, 2.4695268928319),
        G=(-0.2210267356982, 2.9872825343285),
        H=(-0.4849805921475, 3.2410843193759),
        I=(0.5099224052383, 2.9162180345153)
    )
    pt10 = (
        A=(0.3576413342098, 4.1649228169483),
        B=(0.0, 4.0),
        C=(0.4794661910326, 3.6573192468536),
        D=(0.7028117618743, 3.3527571047967),
        E=(0.6520514048648, 4.4796370304071),
        F=(-0.3530036639228, 4.1243145313408),
        G=(0.0, 3.7689920322744)
    )
    pt11 = (
        A=(1.3931526172031, 4.0735541743313),
        B=(2.0022769013168, 3.718231675265),
        C=(1.3931526172031, 3.5354943900309),
        D=(2.1444059009434, 3.4339736760119),
        E=(1.3017839745861, 4.3172038879768),
        F=(1.3017839745861, 4.5507015302204),
        G=(1.7992354732789, 3.9923376031161),
        H=(1.6875626878581, 3.5151902472271),
        I=(1.4337609028107, 3.809600317882)
    )
    test_pts = [pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9, pt10, pt11]
    for i in eachindex(T)
        @test DT.is_positively_oriented(DT.triangle_orientation(get_point(pts, T[i]...)...))
        for r in eachindex(test_pts[i])
            @test DT.is_inside(DT.point_position_relative_to_triangle(get_point(pts, T[i]...)..., test_pts[i][r]))
        end
    end
    for i in eachindex(T)
        u, v, w = indices(T[i])
        pu, pv, pw = get_point(pts, u, v, w)
        for r in eachindex(test_pts[i])
            pr = test_pts[i][r]
            @test DT.is_inside(DT.point_position_relative_to_triangle(pu, pv, pw, pr))
        end
    end
    p1 = [8.0, 3.0]
    p2 = [9.0, 2.0]
    p3 = [9.0, 3.0]
    p4 = [9.0, 4.0]
    pts = [p1, p2, p3, p4]
    i, j, k = 1, 2, 3
    ℓ = 4
    @test DT.is_outside(DT.point_position_relative_to_triangle(get_point(pts, i, j, k, ℓ)...))
    p1, p2, p3 = ([2.858866215272096, -2.220975375945989], [0.25559192484080395, -0.37340906332046214], [1.3855904656897817, -2.47947044705479])
    pts = [p1, p2, p3]
    τ1 = DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)
    τ2 = DT.construct_triangle(NTuple{3,Int64}, 2, 3, 1)
    τ3 = DT.construct_triangle(NTuple{3,Int64}, 3, 1, 2)
    e = Vector{Any}(undef, 9)
    e[1] = DT.point_position_relative_to_triangle(get_point(pts, τ1..., 1)...)
    e[2] = DT.point_position_relative_to_triangle(get_point(pts, τ2..., 2)...)
    e[3] = DT.point_position_relative_to_triangle(get_point(pts, τ3..., 3)...)
    e[4] = DT.point_position_relative_to_triangle(get_point(pts, τ1..., 1)...)
    e[5] = DT.point_position_relative_to_triangle(get_point(pts, τ2..., 2)...)
    e[6] = DT.point_position_relative_to_triangle(get_point(pts, τ3..., 3)...)
    e[7] = DT.point_position_relative_to_triangle(get_point(pts, τ1..., 1)...)
    e[8] = DT.point_position_relative_to_triangle(get_point(pts, τ2..., 2)...)
    e[9] = DT.point_position_relative_to_triangle(get_point(pts, τ3..., 3)...)
    @test all(DT.is_on, e)
    for _ in 1:5000
        local pts, i, j, k
        n = rand(3:5000)
        pts = rand(SVector{2,Float64}, n)
        i, j, k = rand(1:n, 3)
        while length(unique((i, j, k))) < 3
            i, j, k = rand(1:n, 3)
        end
        local τ1, τ2, τ3, e
        τ = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, i, j, k, pts)
        i, j, k = indices(τ)
        τ1 = DT.construct_triangle(NTuple{3,Int64}, i, j, k)
        τ2 = DT.construct_triangle(NTuple{3,Int64}, j, k, i)
        τ3 = DT.construct_triangle(NTuple{3,Int64}, k, i, j)
        e = Vector{Any}(undef, 9)
        e[1] = DT.point_position_relative_to_triangle(get_point(pts, τ1..., i)...)
        e[2] = DT.point_position_relative_to_triangle(get_point(pts, τ2..., i)...)
        e[3] = DT.point_position_relative_to_triangle(get_point(pts, τ3..., i)...)
        e[4] = DT.point_position_relative_to_triangle(get_point(pts, τ1..., j)...)
        e[5] = DT.point_position_relative_to_triangle(get_point(pts, τ2..., j)...)
        e[6] = DT.point_position_relative_to_triangle(get_point(pts, τ3..., j)...)
        e[7] = DT.point_position_relative_to_triangle(get_point(pts, τ1..., k)...)
        e[8] = DT.point_position_relative_to_triangle(get_point(pts, τ2..., k)...)
        e[9] = DT.point_position_relative_to_triangle(get_point(pts, τ3..., k)...)
        @test all(DT.is_on, e)
    end
    for _ in 1:5000
        local r, pts
        n = rand(3:5000)
        r = 5sqrt.(rand(n))
        θ = 2π * rand(n)
        pts = [@SVector[r * cos(θ), r * sin(θ)] for (r, θ) in zip(r, θ)]
        for _ in 1:50
            local i, j, k, τ1, τ2, τ3, e
            i, j, k = rand(1:n, 3)
            while length(unique((i, j, k))) < 3
                i, j, k = rand(1:n, 3)
            end
            τ = DT.construct_positively_oriented_triangle(NTuple{3,Int64}, i, j, k, pts)
            i, j, k = indices(τ)
            τ1 = DT.construct_triangle(NTuple{3,Int64}, i, j, k)
            τ2 = DT.construct_triangle(NTuple{3,Int64}, j, k, i)
            τ3 = DT.construct_triangle(NTuple{3,Int64}, k, i, j)
            e = Vector{Any}(undef, 9)
            e[1] = DT.point_position_relative_to_triangle(get_point(pts, τ1..., i)...)
            e[2] = DT.point_position_relative_to_triangle(get_point(pts, τ2..., i)...)
            e[3] = DT.point_position_relative_to_triangle(get_point(pts, τ3..., i)...)
            e[4] = DT.point_position_relative_to_triangle(get_point(pts, τ1..., j)...)
            e[5] = DT.point_position_relative_to_triangle(get_point(pts, τ2..., j)...)
            e[6] = DT.point_position_relative_to_triangle(get_point(pts, τ3..., j)...)
            e[7] = DT.point_position_relative_to_triangle(get_point(pts, τ1..., k)...)
            e[8] = DT.point_position_relative_to_triangle(get_point(pts, τ2..., k)...)
            e[9] = DT.point_position_relative_to_triangle(get_point(pts, τ3..., k)...)
            @test all(DT.is_on, e)
            pu, pv, pw = get_point(pts, i, j, k)
            e[1] = DT.point_position_relative_to_triangle(pu, pv, pw, pu)
            e[2] = DT.point_position_relative_to_triangle(pv, pw, pu, pu)
            e[3] = DT.point_position_relative_to_triangle(pw, pu, pv, pu)
            e[4] = DT.point_position_relative_to_triangle(pu, pv, pw, pv)
            e[5] = DT.point_position_relative_to_triangle(pv, pw, pu, pv)
            e[6] = DT.point_position_relative_to_triangle(pw, pu, pv, pv)
            e[7] = DT.point_position_relative_to_triangle(pu, pv, pw, pw)
            e[8] = DT.point_position_relative_to_triangle(pv, pw, pu, pw)
            e[9] = DT.point_position_relative_to_triangle(pw, pu, pv, pw)
            @test all(DT.is_on, e)
        end
    end
end

@testset "point_position_relative_to_oriented_outer_halfplane" begin
    p, q = (3.0, 4.0), (6.0, 4.0)
    c, d, e, f, g, h, i, j, k, ℓ = (2.0, 4.0), (7.0, 4.0), (5.0, 5.0),
    (4.0, 6.0), (6.0, 7.0), (4.0, 4.0), (5.0, 4.0),
    (3.0, 2.0), (5.0, 3.0), (6.0, 1.0)
    (4.0, 6.0), (6.0, 7.0), (4.0, 4.0), (5.0, 4.0),
    (3.0, 2.0), (5.0, 3.0), (6.0, 1.0)
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, c))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, d))
    @test DT.is_inside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, e))
    @test DT.is_inside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, f))
    @test DT.is_inside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, f))
    @test DT.is_on(DT.point_position_relative_to_oriented_outer_halfplane(p, q, h))
    @test DT.is_on(DT.point_position_relative_to_oriented_outer_halfplane(p, q, i))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, j))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, k))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, ℓ))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, p))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, q))
    @inferred DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, ℓ))
    p, q = (3.0, 4.0), (6.0, 7.0)
    c, d, e, f, g, h, i, j, k, ℓ = (3.0, 5.0), (3.0, 7.0), (4.0, 8.0),
    (2.0, 3.0), (5.787, 3.774), (4.128, 1.626),
    (6.95784, 7.715), (4.0, 5.0), (5.0, 6.0),
    (4.143, 6.57)
    (2.0, 3.0), (5.787, 3.774), (4.128, 1.626),
    (6.95784, 7.715), (4.0, 5.0), (5.0, 6.0),
    (4.143, 6.57)
    @test DT.is_inside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, c))
    @test DT.is_inside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, d))
    @test DT.is_inside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, e))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, f))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, f))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, h))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, i))
    @test DT.is_on(DT.point_position_relative_to_oriented_outer_halfplane(p, q, j))
    @test DT.is_on(DT.point_position_relative_to_oriented_outer_halfplane(p, q, k))
    @test DT.is_inside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, ℓ))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, p))
    @test DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, q))
    @inferred DT.is_outside(DT.point_position_relative_to_oriented_outer_halfplane(p, q, ℓ))
end

@testset "is_legal" begin
    @testset "Solid edges" begin
        p1 = [-1.0, 4.0]
        p2 = [4.0, 6.0]
        p3 = [2.0, 2.0]
        p4 = [6.0, -3.0]
        p5 = [7.0, 3.0]
        p6 = [-2.5519459826976, 13.6106262700637]
        p7 = [-21.0502221073507, -23.3458204355075]
        p8 = [23.7055438911111, 19.1906513812123]
        p9 = [31.7813088302274, -22.9759380718838]
        pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
        tri = triangulate(pts)
        @test all(DT.is_legal(DT.is_legal(tri, i, j)) for (i, j) in
                  ((1, 3), (3, 2), (2, 1), (1, 4), (4, 3), (3, 1), (3, 4), (4, 5), (5, 3), (3, 5), (5, 2), (2, 3))
        )
    end

    @testset "Boundary edges" begin
        tri = example_triangulation()
        DT.add_triangle!(tri, 6, 2, 3)
        DT.split_triangle!(tri, 1, 3, 5, 7)
        DT.flip_edge!(tri, 5, 1, 4, 7)
        for ((i, j), v) in get_adjacent(tri)
            @test DT.is_legal(DT.is_legal(tri, i, j))
        end
    end
end
