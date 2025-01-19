# # Point-in-Polygon Testing 
# 

# This tutorial shows how we can perform point-in-polygon testing, and how we can find a polygon 
# containing a point. First, let us build this ploygon.
using DelaunayTriangulation
using CairoMakie
using StableRNGs
using ReferenceTests #src
using StatsBase #src
fig_path = joinpath(@__DIR__, "../figures") #src

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
A_curve_outline = [
    [
        K5, W3, Z3, A4, B4, C4, D4, E4, F4, G4, H4, I4, J4, K4, L4, M4, N4,
        O4, P4, Q4, R4, S4, T4, U4, V4, W4, Z4, A5, B5, C5, D5, E5, F5, G5,
        H5, I5, J5, K5,
    ],
]
A_curve_hole = [[L5, M5, N5, O5, P5, Q5, R5, S5, T5, U5, L5]]
dot_1 = [[Z1, A2, B2, C2, D2, E2, F2, G2, H2, I2, J2, Z1]]
dot_2 = [[Z2, A3, B3, C3, D3, E3, F3, G3, H3, I3, J3, Z2]]
dot_3 = [[K2, L2, M2, N2, O2, P2, Q2, R2, S2, T2, U2, V2, W2, K2]]
dot_4 = [[K3, L3, M3, N3, O3, P3, Q3, R3, S3, T3, U3, V3, K3]]
curves = [J_curve, U_curve, L_curve, I_curve, A_curve_outline, A_curve_hole, dot_1, dot_2, dot_3, dot_4]
nodes, points = convert_boundary_points_to_indices(curves)
xmin, xmax, ymin, ymax = DelaunayTriangulation.polygon_bounds(points, nodes, Val(true)) # Val(true) => check all parts of the polygon
rng = StableRNG(123)
query_points = [((xmax - xmin) * rand(rng) + xmin, (ymax - ymin) * rand(rng) + ymin) for _ in 1:1000]
fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, query_points)
for nodes in nodes
    lines!(ax, points[reduce(vcat, nodes)], color = :magenta, linewidth = 3)
end
fig

# To now determine which of these query points are inside any of the magenta regions, we have several methods:
# - Using [`DelaunayTriangulation.distance_to_polygon`](@ref), check if the distance to the polygon is positive (inside) or negative (outside).
# - Triangulate the domain and use [`DelaunayTriangulation.dist`](@ref) to get the distance from a point to the triangulation.
# - Triangulate the domain and use [`find_polygon`](@ref) to find a polygon, or identify the lack of a polygon, containing the point. 
# - Similar to the previous method, but constructing the [`DelaunayTriangulation.PolygonHierarchy`](@ref) to avoid the need for triangulating (this is what the third method does).

# The most efficient approach is the first (since all of them essentially just use the first method), but for the sake of demonstration, we will show all four methods here since, 
# for example, you might already have a triangulation and want to use that. We also emphasise that the tests we perform 
# here do not use exact arithmetic, unlike other predicates in this package, so there may be some robustness issues 
# for points very close to the boundary. Here is the first approach:
is_inside = [DelaunayTriangulation.distance_to_polygon(q, points, nodes) > 0 for q in query_points]
scatter!(ax, query_points[is_inside], color = :blue)
scatter!(ax, query_points[.!is_inside], color = :red)
fig
@test_reference joinpath(fig_path, "point_in_polygon_ex_1.png") fig #src

# Here is the second method. 
tri = triangulate(points; boundary_nodes = nodes)
is_inside_2 = [DelaunayTriangulation.dist(tri, q) > 0 for q in query_points];
@test is_inside == is_inside_2 #src

# The third method is to use [`find_polygon`](@ref) to find the polygon containing the point. If no such polygon exists, `find_polygon` returns 
# `0`, so this is what we use to determine if a point is inside or outside the polygon.
is_inside_3 = [find_polygon(tri, q) ≠ 0 for q in query_points];
@test mean(is_inside) ≈ mean(is_inside_3) atol = 1.0e-2 #src

# This test is not exactly the same as the previous one (with a difference of about five points) due to points near the boundary. 
# The fourth method is:
hierarchy = DelaunayTriangulation.construct_polygon_hierarchy(points, nodes)
is_inside_4 = [find_polygon(hierarchy, points, nodes, q) ≠ 0 for q in query_points];
@test is_inside_4 == is_inside_3 #src
