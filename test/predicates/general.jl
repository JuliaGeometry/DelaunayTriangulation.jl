using ExactPredicates
using ..DelaunayTriangulation
const DT = DelaunayTriangulation

@test DT.opposite_signs(1, -1)
@test DT.opposite_signs(-1, 1)
@test !DT.opposite_signs(0, 0)
@test !DT.opposite_signs(0, 1)
@test !DT.opposite_signs(1, 0)
@test !DT.opposite_signs(0, -1)
@test !DT.opposite_signs(-1, 0)
@test !DT.opposite_signs(1, 1)
@test !DT.opposite_signs(-1, -1)

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

p1, q1, r1 = (3.0, 1.5), (4.5, 1.5), (4.0, 2.0) # +
p2, q2, r2 = (3.0, 1.5), (4.5, 1.5), (3.4, 1.2) # - 
p3, q3, r3 = (2.8, 1.4), (2.8, 2.4), (4.6, 1.4) # - 
p4, q4, r4 = (2.8, 1.4), (3.8, 1.4), (4.6, 1.4) # 0 
p5, q5, r5 = (5.0, 1.4), (3.8, 1.4), (4.6, 1.4) # 0 
p6, q6, r6 = (5.0, 1.4), (3.8, 1.4), (4.4, 0.8) # +
pqr = ((p1, q1, r1), (p2, q2, r2), (p3, q3, r3),
    (p4, q4, r4), (p5, q5, r5), (p6, q6, r6))
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

p, q, r = (2.0, 1.0), (5.0, 1.0), (2.0, 5.0)
c1, c2, c3, c4 = (2.57, 3.35), (2.45, 2.51), (3.11, 2.43), (3.33, 1.83) # inside 
d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12 = (1.76, 4.06), (1.24, 4.109),
(2.589, 0.505), (3.19, -0.186),
(3.875, 4.53), (4.8, 2.698), (2.0, 6.5),
(2.0, 0.0),
(0.5, 1.0), (6.0, 1.0), (6.0, 0.0),
(1.0, 6.0) # outside
e1, e2, e3, e4, e5, e6, e7, e8, e9, e10 = (2.0, 3.0), (2.0, 2.5), (2.0, 4.5), (2.0, 1.5),
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

p, q = (3.0, 4.0), (6.0, 4.0)
c, d, e, f, g, h, i, j, k, ℓ = (2.0, 4.0), (7.0, 4.0), (5.0, 5.0),
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
