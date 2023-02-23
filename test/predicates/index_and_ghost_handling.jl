using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using CairoMakie
using StatsBase

include("../helper_functions.jl")

x, y = complicated_geometry()
tri = generate_mesh(x, y, 2.0; convert_result=true, add_ghost_triangles=true)
DT.compute_representative_points!(tri)
boundary_map = DT.get_boundary_map(tri)
pts = DT.get_points(tri)

for T in each_triangle(tri)
    i, j, k = indices(T)
    cert1 = DT.triangle_orientation(i, j, k, pts, boundary_map)
    cert2 = DT.triangle_orientation(T, pts, boundary_map)
    cert3 = DT.triangle_orientation(tri, T)
    cert4 = DT.triangle_orientation(tri, i, j, k)
    cert5 = DT.triangle_orientation(tri, i, pts[:, j], k)
    @test all(DT.is_positively_oriented, (cert1, cert2, cert3, cert4, cert5))
    @inferred DT.triangle_orientation(i, j, k, pts, boundary_map)
    @inferred DT.triangle_orientation(T, pts, boundary_map)
    @inferred DT.triangle_orientation(tri, T)
    @inferred DT.triangle_orientation(tri, i, j, k)
    @inferred DT.triangle_orientation(tri, i, pts[:, j], k)
end

for T in each_triangle(tri)
    if !DT.is_ghost_triangle(T)
        i, j, k = indices(T)
        for ℓ in (i, j, k)
            cert1 = DT.point_position_relative_to_circumcircle(i, j, k, ℓ, pts,
                                                               boundary_map)
            cert2 = DT.point_position_relative_to_circumcircle(T, ℓ, pts, boundary_map)
            cert3 = DT.point_position_relative_to_circumcircle(tri, T, ℓ)
            cert4 = DT.point_position_relative_to_circumcircle(tri, i, j, k, ℓ)
            cert5 = DT.point_position_relative_to_circumcircle(tri, pts[:, i], j, k,
                                                               pts[:, ℓ])
            cert6 = DT.point_position_relative_to_circumcircle(tri, i, pts[:, j], k, ℓ)
            @test all(DT.is_on, (cert1, cert2, cert3, cert4, cert5, cert6))
            @inferred DT.point_position_relative_to_circumcircle(i, j, k, ℓ, pts,
                                                                 boundary_map)
            @inferred DT.point_position_relative_to_circumcircle(T, ℓ, pts, boundary_map)
            @inferred DT.point_position_relative_to_circumcircle(tri, T, ℓ)
            @inferred DT.point_position_relative_to_circumcircle(tri, i, j, k, ℓ)
            @inferred DT.point_position_relative_to_circumcircle(tri, pts[:, i], j, k,
                                                                 pts[:, ℓ])
            @inferred DT.point_position_relative_to_circumcircle(tri, i, pts[:, j], k, ℓ)
        end
        q = (pts[:, i] .+ pts[:, j] .+ pts[:, k]) ./ 3
        append!(pts, q)
        ℓ = size(pts, 2)
        cert1 = DT.point_position_relative_to_circumcircle(i, j, k, ℓ, pts, boundary_map)
        cert2 = DT.point_position_relative_to_circumcircle(T, ℓ, pts, boundary_map)
        cert3 = DT.point_position_relative_to_circumcircle(tri, T, ℓ)
        cert4 = DT.point_position_relative_to_circumcircle(tri, i, j, k, ℓ)
        @test all(DT.is_inside, (cert1, cert2, cert3, cert4))
        resize!(pts, (2, ℓ - 1))
    end
end
tri, label_map, index_map = simple_geometry()
pts = get_points(tri)
orig_length = length(pts)
DT.compute_representative_points!(tri)
DT.RepresentativePointList[1] = DT.RepresentativeCoordinates(mean(reinterpret(reshape,
                                                                              Float64,
                                                                              pts[unique(reduce(vcat,
                                                                                                tri.boundary_nodes[1]))]);
                                                                  dims=2)..., 0)
DT.RepresentativePointList[2] = DT.RepresentativeCoordinates(mean(reinterpret(reshape,
                                                                              Float64,
                                                                              pts[unique(reduce(vcat,
                                                                                                tri.boundary_nodes[2]))]);
                                                                  dims=2)..., 0)
DT.RepresentativePointList[3] = DT.RepresentativeCoordinates(mean(reinterpret(reshape,
                                                                              Float64,
                                                                              pts[unique(reduce(vcat,
                                                                                                tri.boundary_nodes[3]))]);
                                                                  dims=2)..., 0)
boundary_map = DT.get_boundary_map(tri)
i, j, k = index_map["w"], index_map["m"], index_map["v"]
ℓ = index_map["z"]
T = (i, j, k)
cert1 = DT.point_position_relative_to_circumcircle(i, j, k, ℓ, pts, boundary_map)
cert2 = DT.point_position_relative_to_circumcircle(T, ℓ, pts, boundary_map)
cert3 = DT.point_position_relative_to_circumcircle(tri, T, ℓ)
cert4 = DT.point_position_relative_to_circumcircle(tri, i, j, k, ℓ)
@test all(DT.is_inside, (cert1, cert2, cert3, cert4))
ℓ = index_map["q"]
cert1 = DT.point_position_relative_to_circumcircle(i, j, k, ℓ, pts, boundary_map)
cert2 = DT.point_position_relative_to_circumcircle(T, ℓ, pts, boundary_map)
cert3 = DT.point_position_relative_to_circumcircle(tri, T, ℓ)
cert4 = DT.point_position_relative_to_circumcircle(tri, i, j, k, ℓ)
@test all(DT.is_outside, (cert1, cert2, cert3, cert4))

f1 = (-5.0, 20.0)
g1 = (2.2179294, 23.7869)
push!(pts, f1, g1)
f1_i = length(pts) - 1
g1_i = length(pts)
h1 = (21.517785930112, 22.1693108196262)
i1 = (22.0683240766855, 21.526199282940)
push!(pts, h1, i1)
h1_i = length(pts) - 1
i1_i = length(pts)
j1 = (6.8116906951308, 14.9338268341201)
k1 = (7.3767511907439, 10.1308126214084)
ℓ1 = (6.0582767009799, 9.3303102526232)
m1 = (4.3630952141405, 10.3662544945806)
push!(pts, j1, k1, ℓ1, m1)
j1_i = length(pts) - 3
k1_i = length(pts) - 2
ℓ1_i = length(pts) - 1
m1_i = length(pts)
n1 = (15.5995300899034, 9.481188316425)
o1 = (16.9311921105779, 7.7056389555257)
p1 = (16.4238922931781, 5.2008461071143)
q1 = (14.0, 4.0)
r1 = (14.8385803638037, 3.425296746215)
s1 = (14.3312805464039, 7.4519890468258)
push!(pts, n1, o1, p1, q1, r1, s1)
n1_i = length(pts) - 5
o1_i = length(pts) - 4
p1_i = length(pts) - 3
q1_i = length(pts) - 2
r1_i = length(pts) - 1
s1_i = length(pts)
certs = [("w", "v", "f") => Certificate.Left,
         ("c", "d", "k") => Certificate.Left,
         ("r", "q", "j") => Certificate.Right,
         ("r", "q", "d") => Certificate.Left,
         ("u", "w", "i") => Certificate.Left,
         ("c", "d", "e") => Certificate.Collinear,
         ("g", DT.BoundaryIndex, f1_i) => Certificate.Left,
         ("g", DT.BoundaryIndex, g1_i) => Certificate.Right,
         (DT.BoundaryIndex, "g", f1_i) => Certificate.Right,
         (DT.BoundaryIndex, "g", g1_i) => Certificate.Left,
         ("a", DT.BoundaryIndex, f1_i) => Certificate.Right,
         ("e", DT.BoundaryIndex, h1_i) => Certificate.Left,
         ("e", DT.BoundaryIndex, i1_i) => Certificate.Right,
         ("ℓ", DT.BoundaryIndex - 1, j1_i) => Certificate.Right,
         (DT.BoundaryIndex - 1, "ℓ", j1_i) => Certificate.Left,
         ("i", DT.BoundaryIndex - 1, j1_i) => Certificate.Left,
         ("i", DT.BoundaryIndex - 1, m1_i) => Certificate.Right,
         (DT.BoundaryIndex - 1, "i", m1_i) => Certificate.Left,
         ("j", DT.BoundaryIndex - 1, "b1") => Certificate.Left,
         ("j", DT.BoundaryIndex - 1, ℓ1_i) => Certificate.Right,
         ("k", DT.BoundaryIndex - 1, ℓ1_i) => Certificate.Left,
         ("k", DT.BoundaryIndex - 1, k1_i) => Certificate.Right,
         (DT.BoundaryIndex - 1, "k", k1_i) => Certificate.Left,
         ("r", DT.BoundaryIndex - 2, n1_i) => Certificate.Right,
         ("r", DT.BoundaryIndex - 2, o1_i) => Certificate.Left,
         (DT.BoundaryIndex - 3, "r", n1_i) => Certificate.Left,
         (DT.BoundaryIndex - 2, "q", o1_i) => Certificate.Left,
         ("q", DT.BoundaryIndex - 2, o1_i) => Certificate.Right,
         ("p", DT.BoundaryIndex - 3, q1_i) => Certificate.Left,
         ("q", DT.BoundaryIndex - 3, p1_i) => Certificate.Left,
         ("o", DT.BoundaryIndex - 2, q1_i) => Certificate.Right,
         ("m", DT.BoundaryIndex - 2, s1_i) => Certificate.Right,
         (DT.BoundaryIndex - 2, "m", s1_i) => Certificate.Left,
         (DT.BoundaryIndex - 2, "p", r1_i) => Certificate.Left,
         ("p", DT.BoundaryIndex - 2, r1_i) => Certificate.Right]
for ((i, j, u), cert) in certs
    i = i isa String ? index_map[i] : i
    j = j isa String ? index_map[j] : j
    u = u isa String ? index_map[u] : u
    @test DT.point_position_relative_to_line(i, j, pts[u], pts, boundary_map) == cert ==
          DT.point_position_relative_to_line(i, j, u, pts, boundary_map)
    @inferred DT.point_position_relative_to_line(i, j, u, pts, boundary_map)
    @inferred DT.point_position_relative_to_line(i, j, pts[u], pts, boundary_map)
end

t1 = (20.0, 25.0)
u1 = (20.0, -5.0)
push!(pts, t1, u1)
t1_i = length(pts) - 1
u1_i = length(pts)
i, j = index_map["c"], index_map["e"]
k = index_map["d"]
cert1 = DT.point_position_on_line_segment(i, j, k, pts)
cert2 = DT.point_position_on_line_segment(i, j, t1_i, pts)
cert3 = DT.point_position_on_line_segment(i, j, u1_i, pts)
cert4 = DT.point_position_on_line_segment(j, i, k, pts)
cert5 = DT.point_position_on_line_segment(j, i, t1_i, pts)
cert6 = DT.point_position_on_line_segment(j, i, u1_i, pts)
cert7 = DT.point_position_on_line_segment(j, i, i, pts)
cert8 = DT.point_position_on_line_segment(j, i, j, pts)
cert9 = DT.point_position_on_line_segment(i, j, i, pts)
cert10 = DT.point_position_on_line_segment(i, j, j, pts)
cert11 = DT.point_position_on_line_segment(i, j, pts[j], pts)
cert12 = DT.point_position_on_line_segment(pts[i], j, i, pts)
@test all(DT.is_on, (cert1, cert4))
@test all(DT.is_left, (cert3, cert5))
@test all(DT.is_right, (cert2, cert6))
@test all(DT.is_degenerate, (cert7, cert8, cert9, cert10, cert11, cert12))
v1 = (-3.0, 0.0)
w1 = (22.0, 0.0)
push!(pts, v1, w1)
v1_i = length(pts) - 1
w1_i = length(pts)
i, j = index_map["a"], index_map["c"]
k = index_map["b"]
cert1 = DT.point_position_on_line_segment(i, j, k, pts)
cert2 = DT.point_position_on_line_segment(i, j, v1_i, pts)
@inferred DT.point_position_on_line_segment(i, j, v1_i, pts)
cert3 = DT.point_position_on_line_segment(i, j, w1_i, pts)
cert4 = DT.point_position_on_line_segment(j, i, k, pts)
cert5 = DT.point_position_on_line_segment(j, i, v1_i, pts)
cert6 = DT.point_position_on_line_segment(j, i, w1_i, pts)
cert7 = DT.point_position_on_line_segment(j, i, i, pts)
cert8 = DT.point_position_on_line_segment(j, i, j, pts)
cert9 = DT.point_position_on_line_segment(i, j, i, pts)
cert10 = DT.point_position_on_line_segment(i, j, j, pts)
@inferred DT.point_position_on_line_segment(i, j, j, pts)
@inferred DT.point_position_on_line_segment(i, j, pts[j], pts)
@inferred DT.point_position_on_line_segment(i, pts[j], j, pts)
@test all(DT.is_on, (cert1, cert4))
@test all(DT.is_left, (cert2, cert6))
@test all(DT.is_right, (cert3, cert5))
@test all(DT.is_degenerate, (cert7, cert8, cert9, cert10))

z1 = (47.93, 17.33)
a2 = (53.135, 10.66)
b2 = (54.8, 18.16)
d2 = (45.7057, 12.1682)
push!(pts, z1, a2, b2, d2)
z1_i = length(pts) - 3
a2_i = length(pts) - 2
b2_i = length(pts) - 1
d2_i = length(pts)
@test DT.is_single(DT.line_segment_intersection_type(z1_i, a2_i, b2_i, d2_i, pts))
@test DT.is_single(DT.line_segment_intersection_type(b2_i, d2_i, z1_i, a2_i, pts))
@test DT.is_single(DT.line_segment_intersection_type(b2_i, d2, z1_i, a2, pts))
@inferred DT.line_segment_intersection_type(z1_i, a2_i, b2_i, d2_i, pts)
@inferred DT.line_segment_intersection_type(z1, a2, b2_i, d2_i, pts)
e2 = (48.0, 6.0)
f2 = (54.0, 6.0)
g2 = (52.0, 8.0)
h2 = (52.0, 6.0)
push!(pts, e2, f2, g2, h2)
e2_i = length(pts) - 3
f2_i = length(pts) - 2
g2_i = length(pts) - 1
h2_i = length(pts)
@test DT.is_touching(DT.line_segment_intersection_type(e2_i, f2_i, g2_i, h2_i, pts))
@test DT.is_touching(DT.line_segment_intersection_type(g2_i, h2_i, e2_i, f2_i, pts))
@test DT.is_touching(DT.line_segment_intersection_type(e2_i, h2_i, h2_i, g2_i, pts))
@test DT.is_touching(DT.line_segment_intersection_type(g2_i, h2_i, h2_i, e2_i, pts))
@inferred DT.line_segment_intersection_type(g2_i, h2_i, e2_i, f2_i, pts)
i2 = (53.58, 35.45)
j2 = (57.11, 27.3)
k2 = (49.61, 37.12)
ℓ2 = (56.39, 40.42)
push!(pts, i2, j2, k2, ℓ2)
i2_i = length(pts) - 3
j2_i = length(pts) - 2
k2_i = length(pts) - 1
ℓ2_i = length(pts)
@test DT.is_none(DT.line_segment_intersection_type(i2_i, j2_i, k2_i, ℓ2_i, pts))
@test DT.is_none(DT.line_segment_intersection_type(k2_i, ℓ2_i, i2_i, j2_i, pts))
m2 = (50.0, 30.0)
n2 = (50.0, 24.0)
o2 = (50.0, 28.0)
p2 = (50.0, 20.0)
push!(pts, m2, n2, o2, p2)
m2_i = length(pts) - 3
n2_i = length(pts) - 2
o2_i = length(pts) - 1
p2_i = length(pts)
@test DT.is_multiple(DT.line_segment_intersection_type(m2_i, n2_i, o2_i, p2_i, pts))
@test DT.is_multiple(DT.line_segment_intersection_type(n2_i, o2_i, m2_i, p2_i, pts))
@test DT.is_multiple(DT.line_segment_intersection_type(o2_i, p2_i, m2_i, n2_i, pts))
@test DT.is_none(DT.line_segment_intersection_type(p2_i, n2_i, o2_i, m2_i, pts))

resize!(pts, orig_length)
f1 = (0.5742824217282, 13.5106352620416)
g1 = (1.3454554411888, 9.8289060078422)
h1 = (3.3355793623777, 12.01804232115)
i1 = (2.5146532448873, 15.4510060852008)
j1 = (2.1415050096644, 6.9183497731035)
k1 = (3.1863200682886, 5.5501395772861)
ℓ1 = (1.0220603039957, 4.0326700873796)
m1 = (4, 17.0679817711668)
n1 = (6.0, 2.0)
o1 = (4.0, 20.0)
p1 = (8.0, 20.0)
q1 = (0.0, 2.0)
r1 = (2.0, 0.0)
push!(pts, f1, g1, h1, i1, j1, k1, ℓ1, m1, n1, o1, p1, q1, r1)
r1_i = length(pts)
q1_i = length(pts) - 1
p1_i = length(pts) - 2
o1_i = length(pts) - 3
n1_i = length(pts) - 4
m1_i = length(pts) - 5
ℓ1_i = length(pts) - 6
k1_i = length(pts) - 7
j1_i = length(pts) - 8
i1_i = length(pts) - 9
h1_i = length(pts) - 10
g1_i = length(pts) - 11
f1_i = length(pts) - 12
certs = [("g", "h", "b1", f1_i) => Certificate.Inside,
         ("b1", "i", "g", f1_i) => Certificate.Outside,
         ("b1", "h", "j", g1_i) => Certificate.Inside,
         ("i", "b1", "j", h1_i) => Certificate.Inside,
         ("i", "ℓ", "f", h1_i) => Certificate.Outside,
         ("h", "a", "s", ℓ1_i) => Certificate.Inside,
         ("h", "s", "j", k1_i) => Certificate.Inside,
         ("s", "j", "h", j1_i) => Certificate.Inside,
         ("i", "g", "b1", i1_i) => Certificate.Inside,
         ("a", "b", "t", r1_i) => Certificate.On,
         ("b", "t", "a", r1_i) => Certificate.On,
         ("a", "t", "s", n1_i) => Certificate.On,
         ("i", "f", "a1", m1_i) => Certificate.On,
         ("g", "i", "a1", m1_i) => Certificate.On,
         ("g", "a1", "f", o1_i) => Certificate.On,
         ("a1", "f", "g", o1_i) => Certificate.On,
         ("h", "a", "s", q1_i) => Certificate.On,
         ("a1", "f", "g", p1_i) => Certificate.On,
         ("b1", "j", "i", "v") => Certificate.Outside]
for ((i, j, k, u), cert) in certs
    local cert1, cert2, cert3, cert4, cert5, cert6
    i = i isa String ? index_map[i] : i
    j = j isa String ? index_map[j] : j
    k = k isa String ? index_map[k] : k
    u = u isa String ? index_map[u] : u
    T1 = (i, j, k)
    T2 = [i, j, k]
    if !DT.is_boundary_index(i)
        cert1 = DT.point_position_relative_to_triangle(i, j, k, u, pts, tri.boundary_map)
    else
        cert1 = DT.point_position_relative_to_triangle(pts[i], j, pts[k], u, pts,
                                                       tri.boundary_map)
        @inferred DT.point_position_relative_to_triangle(pts[i], j, pts[k], u, pts,
                                                         tri.boundary_map)
    end
    cert2 = DT.point_position_relative_to_triangle(T1, u, pts, tri.boundary_map)
    cert3 = DT.point_position_relative_to_triangle(T2, u, pts, tri.boundary_map)
    cert4 = DT.point_position_relative_to_triangle(tri, i, j, k, u)
    cert5 = DT.point_position_relative_to_triangle(tri, T1, u)
    cert6 = DT.point_position_relative_to_triangle(tri, T2, pts[u])
    @test all(==(cert), (cert1, cert2, cert3, cert4, cert5, cert6))
    @inferred DT.point_position_relative_to_triangle(i, j, k, u, pts, tri.boundary_map)
    @inferred DT.point_position_relative_to_triangle(T1, u, pts, tri.boundary_map)
    @inferred DT.point_position_relative_to_triangle(T2, u, pts, tri.boundary_map)
    @inferred DT.point_position_relative_to_triangle(tri, i, j, k, u)
    @inferred DT.point_position_relative_to_triangle(tri, T1, u)
    @inferred DT.point_position_relative_to_triangle(tri, T2, u)
    @inferred DT.point_position_relative_to_triangle(tri, T2, pts[u])
end
resize!(pts, orig_length)
f1 = (-4.0, 20.0)
g1 = (-2.0, 18.0)
h1 = (-6.0, 12.0)
i1 = (-6.0, 22.0)
j1 = (0.0, 18.0)
k1 = (0.0, 12.0)
ℓ1 = (-2.0, 10.0)
m1 = (-4.0, 24.0)
n1 = (-4.0, 8.0)
o1 = (0.0, 4.0)
p1 = (0.0, 22.0)
push!(pts, f1, g1, h1, i1, j1, k1, ℓ1, m1, n1, o1, p1)
p1_i = length(pts)
o1_i = length(pts) - 1
n1_i = length(pts) - 2
m1_i = length(pts) - 3
ℓ1_i = length(pts) - 4
k1_i = length(pts) - 5
j1_i = length(pts) - 6
i1_i = length(pts) - 7
h1_i = length(pts) - 8
g1_i = length(pts) - 9
f1_i = length(pts) - 10
certs = [("h", "g", DT.BoundaryIndex, f1_i) => Certificate.Inside,
         ("h", "g", DT.BoundaryIndex, g1_i) => Certificate.Inside,
         ("h", "g", DT.BoundaryIndex, h1_i) => Certificate.Inside,
         (DT.BoundaryIndex, "h", "g", i1_i) => Certificate.Inside,
         ("g", DT.BoundaryIndex, "h", j1_i) => Certificate.On,
         ("h", "g", DT.BoundaryIndex, k1_i) => Certificate.On,
         ("h", "g", DT.BoundaryIndex, ℓ1_i) => Certificate.On,
         (DT.BoundaryIndex, "h", "g", m1_i) => Certificate.On,
         ("h", "g", DT.BoundaryIndex, n1_i) => Certificate.Outside,
         ("h", "g", DT.BoundaryIndex, o1_i) => Certificate.Outside,
         ("h", "g", DT.BoundaryIndex, p1_i) => Certificate.Outside,
         ("h", "g", DT.BoundaryIndex, "b1") => Certificate.Outside]
DT.RepresentativePointList[1].x = 10.0
DT.RepresentativePointList[1].y = 10.0
for ((i, j, k, u), cert) in certs
    local cert1, cert2, cert3, cert4, cert5, cert6
    i = i isa String ? index_map[i] : i
    j = j isa String ? index_map[j] : j
    k = k isa String ? index_map[k] : k
    u = u isa String ? index_map[u] : u
    T1 = (i, j, k)
    T2 = [i, j, k]
    cert1 = DT.point_position_relative_to_triangle(i, j, k, u, pts, tri.boundary_map)
    cert2 = DT.point_position_relative_to_triangle(T1, u, pts, tri.boundary_map)
    cert3 = DT.point_position_relative_to_triangle(T2, u, pts, tri.boundary_map)
    cert4 = DT.point_position_relative_to_triangle(tri, i, j, k, u)
    cert5 = DT.point_position_relative_to_triangle(tri, T1, u)
    cert6 = DT.point_position_relative_to_triangle(tri, T2, u)
    @test all(==(cert), (cert1, cert2, cert3, cert4, cert5, cert6))
    @inferred DT.point_position_relative_to_triangle(i, j, k, u, pts, tri.boundary_map)
    @inferred DT.point_position_relative_to_triangle(T1, u, pts, tri.boundary_map)
    @inferred DT.point_position_relative_to_triangle(T2, u, pts, tri.boundary_map)
    @inferred DT.point_position_relative_to_triangle(tri, i, j, k, u)
    @inferred DT.point_position_relative_to_triangle(tri, T1, u)
    @inferred DT.point_position_relative_to_triangle(tri, T2, u)
end
resize!(pts, orig_length)
f1 = (2.98004, 22.706)
g1 = (6.449, 23.343)
h1 = (10.0, 24.0)
i1 = (12.4666, 22.777)
j1 = (18.555, 22.600)
k1 = (24.0, 24.0)
push!(pts, f1, g1, h1, i1, j1, k1)
f1_i = length(pts) - 5
g1_i = length(pts) - 4
h1_i = length(pts) - 3
i1_i = length(pts) - 2
j1_i = length(pts) - 1
k1_i = length(pts)
certs = [("g", "f", DT.BoundaryIndex, f1_i) => Certificate.Inside,
         ("f", "e", DT.BoundaryIndex, f1_i) => Certificate.Outside,
         ("g", "f", DT.BoundaryIndex, g1_i) => Certificate.Inside,
         ("e", DT.BoundaryIndex, "f", h1_i) => Certificate.On,
         (DT.BoundaryIndex, "f", "e", h1_i) => Certificate.On,
         ("f", "e", DT.BoundaryIndex, i1_i) => Certificate.Inside,
         ("e", DT.BoundaryIndex, "f", j1_i) => Certificate.Inside,
         (DT.BoundaryIndex, "f", "e", k1_i) => Certificate.On,
         ("e", "d", DT.BoundaryIndex, k1_i) => Certificate.On,
         ("b", "a", DT.BoundaryIndex, "s") => Certificate.Outside,
         (DT.BoundaryIndex, "b", "a", "v") => Certificate.Outside]
for ((i, j, k, u), cert) in certs
    local cert1, cert2, cert3, cert4, cert5, cert6
    i = i isa String ? index_map[i] : i
    j = j isa String ? index_map[j] : j
    k = k isa String ? index_map[k] : k
    u = u isa String ? index_map[u] : u
    T1 = (i, j, k)
    T2 = [i, j, k]
    cert1 = DT.point_position_relative_to_triangle(i, j, k, u, pts, tri.boundary_map)
    cert2 = DT.point_position_relative_to_triangle(T1, u, pts, tri.boundary_map)
    cert3 = DT.point_position_relative_to_triangle(T2, u, pts, tri.boundary_map)
    cert4 = DT.point_position_relative_to_triangle(tri, i, j, k, u)
    cert5 = DT.point_position_relative_to_triangle(tri, T1, u)
    cert6 = DT.point_position_relative_to_triangle(tri, T2, u)
    @test all(==(cert), (cert1, cert2, cert3, cert4, cert5, cert6))
end
resize!(pts, orig_length)
f1 = (5.765, 12.7676)
g1 = (6.0, 15.0)
h1 = (7.0, 12.0)
i1 = (6.0, 9.0)
j1 = (6.0, 8.0)
k1 = (5.0, 12.0)
ℓ1 = (4.0, 13.0)
m1 = (8.0, 13.0)
n1 = (7.0, 10.0)
o1 = (6.0, 10.0)
push!(pts, f1, g1, h1, i1, j1, k1, ℓ1, m1, n1, o1, p1)
p1_i = length(pts)
o1_i = length(pts) - 1
n1_i = length(pts) - 2
m1_i = length(pts) - 3
ℓ1_i = length(pts) - 4
k1_i = length(pts) - 5
j1_i = length(pts) - 6
i1_i = length(pts) - 7
h1_i = length(pts) - 8
g1_i = length(pts) - 9
f1_i = length(pts) - 10
certs = [("ℓ", "i", DT.BoundaryIndex - 1, f1_i) => Certificate.Inside,
         ("i", DT.BoundaryIndex - 1, "ℓ", g1_i) => Certificate.Inside,
         (DT.BoundaryIndex - 1, "ℓ", "i", k1_i) => Certificate.Outside,
         ("i", "j", DT.BoundaryIndex - 1, k1_i) => Certificate.Inside,
         ("i", "j", DT.BoundaryIndex - 1, ℓ1_i) => Certificate.On,
         ("j", "k", DT.BoundaryIndex - 1, o1_i) => Certificate.Inside,
         ("j", "k", DT.BoundaryIndex - 1, h1_i) => Certificate.Outside,
         ("i", DT.BoundaryIndex - 1, "ℓ", i1_i) => Certificate.Outside,
         ("k", DT.BoundaryIndex - 1, "j", j1_i) => Certificate.Inside,
         (DT.BoundaryIndex - 1, "j", "k", i1_i) => Certificate.Inside,
         (DT.BoundaryIndex - 1, "k", "ℓ", k1_i) => Certificate.Outside,
         ("k", "ℓ", DT.BoundaryIndex - 1, n1_i) => Certificate.Inside,
         ("k", "ℓ", DT.BoundaryIndex - 1, h1_i) => Certificate.Inside,
         ("ℓ", DT.BoundaryIndex - 1, "k", m1_i) => Certificate.On,
         ("ℓ", DT.BoundaryIndex - 1, "k", "b1") => Certificate.Outside,
         ("j", "k", DT.BoundaryIndex - 1, "s") => Certificate.Outside]
for ((i, j, k, u), cert) in certs
    local cert1, cert2, cert3, cert4, cert5, cert6
    i = i isa String ? index_map[i] : i
    j = j isa String ? index_map[j] : j
    k = k isa String ? index_map[k] : k
    u = u isa String ? index_map[u] : u
    T1 = (i, j, k)
    T2 = [i, j, k]
    cert1 = DT.point_position_relative_to_triangle(i, j, k, u, pts, tri.boundary_map)
    cert2 = DT.point_position_relative_to_triangle(T1, u, pts, tri.boundary_map)
    cert3 = DT.point_position_relative_to_triangle(T2, u, pts, tri.boundary_map)
    cert4 = DT.point_position_relative_to_triangle(tri, i, j, k, u)
    cert5 = DT.point_position_relative_to_triangle(tri, T1, u)
    cert6 = DT.point_position_relative_to_triangle(tri, T2, u)
    @test all(==(cert), (cert1, cert2, cert3, cert4, cert5, cert6))
end
resize!(pts, orig_length)
f1 = (15.0, 9.0)
g1 = (15.5, 8.5)
h1 = (15.0, 7.5)
i1 = (16.0, 7.0)
j1 = (17.5, 8.0)
k1 = (17.5, 9.5)
ℓ1 = (17.0, 5.0)
m1 = (15.5, 3.5)
n1 = (13.5, 2.5)
o1 = (12.5, 3.5)
p1 = (14.21, 8.53)
q1 = (13.49, 4.22)
r1 = (15.438, 5.66)
s1 = (13.26, 5.439)
t1 = (12.0, 4.5)
push!(pts, f1, g1, h1, i1, j1, k1, ℓ1, m1, n1, o1, p1, q1, r1, s1, t1)
t1_i = length(pts)
s1_i = length(pts) - 1
r1_i = length(pts) - 2
q1_i = length(pts) - 3
p1_i = length(pts) - 4
o1_i = length(pts) - 5
n1_i = length(pts) - 6
m1_i = length(pts) - 7
ℓ1_i = length(pts) - 8
k1_i = length(pts) - 9
j1_i = length(pts) - 10
i1_i = length(pts) - 11
h1_i = length(pts) - 12
g1_i = length(pts) - 13
f1_i = length(pts) - 14
certs = [("o", "p", DT.BoundaryIndex - 2, o1_i) => Certificate.On,
         ("p", DT.BoundaryIndex - 2, "o", n1_i) => Certificate.On,
         (DT.BoundaryIndex - 2, "o", "p", q1_i) => Certificate.Inside,
         (DT.BoundaryIndex - 3, "o", "p", t1_i) => Certificate.Outside,
         ("o", DT.BoundaryIndex - 3, "n", s1_i) => Certificate.Outside,
         ("n", "o", DT.BoundaryIndex - 3, t1_i) => Certificate.Outside,
         ("p", "q", DT.BoundaryIndex - 2, p1_i) => Certificate.Outside,
         ("p", "q", DT.BoundaryIndex - 3, r1_i) => Certificate.Inside,
         ("q", DT.BoundaryIndex - 3, "p", m1_i) => Certificate.On,
         (DT.BoundaryIndex - 2, "p", "q", ℓ1_i) => Certificate.On,
         ("r", DT.BoundaryIndex - 3, "q", k1_i) => Certificate.Inside,
         ("r", DT.BoundaryIndex - 3, "q", "d") => Certificate.Outside,
         ("q", "r", DT.BoundaryIndex - 2, j1_i) => Certificate.Inside,
         ("r", DT.BoundaryIndex - 3, "q", i1_i) => Certificate.Inside,
         ("r", DT.BoundaryIndex - 3, "q", "m") => Certificate.Outside,
         ("r", "m", DT.BoundaryIndex - 2, f1_i) => Certificate.Inside,
         ("m", DT.BoundaryIndex - 3, "r", g1_i) => Certificate.Inside,
         ("m", DT.BoundaryIndex - 3, "r", r1_i) => Certificate.Outside,
         ("m", DT.BoundaryIndex - 3, "r", h1_i) => Certificate.Inside,
         ("m", "n", DT.BoundaryIndex - 2, p1_i) => Certificate.Inside,
         (DT.BoundaryIndex - 2, "m", "n", s1_i) => Certificate.Outside]
for ((i, j, k, u), cert) in certs
    local cert1, cert2, cert3, cert4, cert5, cert6
    i = i isa String ? index_map[i] : i
    j = j isa String ? index_map[j] : j
    k = k isa String ? index_map[k] : k
    u = u isa String ? index_map[u] : u
    T1 = (i, j, k)
    T2 = [i, j, k]
    cert1 = DT.point_position_relative_to_triangle(i, j, k, u, pts, tri.boundary_map)
    cert2 = DT.point_position_relative_to_triangle(T1, u, pts, tri.boundary_map)
    cert3 = DT.point_position_relative_to_triangle(T2, u, pts, tri.boundary_map)
    cert4 = DT.point_position_relative_to_triangle(tri, i, j, k, u)
    cert5 = DT.point_position_relative_to_triangle(tri, T1, u)
    cert6 = DT.point_position_relative_to_triangle(tri, T2, u)
    @test all(==(cert), (cert1, cert2, cert3, cert4, cert5, cert6))
end