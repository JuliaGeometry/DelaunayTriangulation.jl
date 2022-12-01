include("functions.jl")

p1 = (5.0, 6.0)
p2 = (9.0, 6.0)
p3 = (13.0, 5.0)
p4 = (10.38, 0.0)
p5 = (12.64, -1.69)
p6 = (2.0, -2.0)
p7 = (3.0, 4.0)
p8 = (7.5, 3.53)
p9 = (4.02, 1.85)
p10 = (4.26, 0.0)
pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0))
xlims!(ax, 1.0, 14.0)
ylims!(ax, -3.0, 7.0)
hidedecorations!(ax)
text!(ax,
    [p1, p2, p3, (p4[1] + 0.2, p4[2] - 0.1),
        (p5[1] + 0.1, p5[2] - 0.1),
        (p6[1] - 0.53, p6[2] - 0.2),
        (p7[1] - 0.46, p7[2] - 0.2),
        (p8[1] - 0.15, p8[2] - 1.1),
        (p9[1], p9[2] + 0.1),
        (p10[1], p10[2] - 1.0)];
    text=[L"%$s" for s in 1:10], textsize=55)
p11 = (6.0, 2.5)
scatter!(ax, p11, color=:blue, markersize=17)
push!(pts, p11)
text!(ax, (p11[1] - 0.74, p11[2] - 0.75); text=L"11", textsize=55)
save("writeups/figures/excavated_cavity_triangulation_a.pdf", fig)

bw_in = NTuple{3,Int64}[]
for T in T
    DelaunayTriangulation.isincircle(T, pts, 11) == 1 && push!(bw_in, T)
end

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
for T in bw_in
    i, j, k = T
    pti, ptj, ptk = pts[i], pts[j], pts[k]
    circx, circy = circle_data(pti, ptj, ptk)
    lines!(ax, circx, circy, color=(:darkgreen, 0.9), linewidth=3, linestyle=:dash)
end
poly!(ax, pts, [bw_in[i][j] for i in eachindex(bw_in), j in 1:3],
    strokewidth=2, color=(:lightblue, 0.5))
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0))
xlims!(ax, 1.0, 14.0)
ylims!(ax, -3.0, 7.0)
hidedecorations!(ax)
text!(ax,
    [p1, p2, p3, (p4[1] + 0.2, p4[2] - 0.1),
        (p5[1] + 0.1, p5[2] - 0.1),
        (p6[1] - 0.53, p6[2] - 0.2),
        (p7[1] - 0.46, p7[2] - 0.2),
        (p8[1] - 0.15, p8[2] - 1.1),
        (p9[1], p9[2] + 0.1),
        (p10[1], p10[2] - 1.0)];
    text=[L"%$s" for s in 1:10], textsize=55)
scatter!(ax, p11, color=:blue, markersize=17)
text!(ax, (p11[1] - 0.74, p11[2] - 0.75); text=L"11", textsize=55)
save("writeups/figures/excavated_cavity_triangulation_b.pdf", fig)

Tmat = T
[DT.delete_triangle!(Tmat, T) for T in bw_in]
DT.add_triangle!(Tmat,
    (1, 11, 8),
    (7, 11, 1),
    (7, 9, 11),
    (9, 10, 11),
    (10, 4, 11),
    (4, 8, 11))

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0))
xlims!(ax, 1.0, 14.0)
ylims!(ax, -3.0, 7.0)
hidedecorations!(ax)
text!(ax,
    [p1, p2, p3, (p4[1] + 0.2, p4[2] - 0.1),
        (p5[1] + 0.1, p5[2] - 0.1),
        (p6[1] - 0.53, p6[2] - 0.2),
        (p7[1] - 0.46, p7[2] - 0.2),
        (p8[1] - 0.15, p8[2] - 1.1),
        (p9[1], p9[2] + 0.1),
        (p10[1], p10[2] - 1.0)];
    text=[L"%$s" for s in 1:10], textsize=55)
scatter!(ax, p11, color=:blue, markersize=17)
text!(ax, (p11[1] - 0.24, p11[2] - 1.2); text=L"11", textsize=55)
fig
lines!(ax, [p11, p1, (NaN, NaN),
        p11, p4, (NaN, NaN),
        p11, p10, (NaN, NaN),
        p11, p9, (NaN, NaN),
        p11, p7, (NaN, NaN),
        p11, p8, (NaN, NaN)],
    color=:blue, linewidth=4)
save("writeups/figures/excavated_cavity_triangulation_c.pdf", fig)

## Setup
p1 = (5.0, 6.0)
p2 = (9.0, 6.0)
p3 = (13.0, 5.0)
p4 = (10.38, 0.0)
p5 = (12.64, -1.69)
p6 = (2.0, -2.0)
p7 = (3.0, 4.0)
p8 = (7.5, 3.53)
p9 = (4.02, 1.85)
p10 = (4.26, 0.0)
pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
Random.seed!(928881)
T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
p11 = (6.0, 2.5)
push!(pts, p11)

DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 11)
_T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
DT.clear_empty_keys!(adj, DG)
DT.clear_empty_keys!(_adj, _DG)
@test DT.compare_triangle_sets(T, _T) &&
      adjacent(adj) == adjacent(_adj) &&
      adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
      graph(DG) == graph(_DG)

p12 = (10.3, 2.85)
push!(pts, p12)
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 12)
_T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
DT.clear_empty_keys!(adj, DG)
DT.clear_empty_keys!(_adj, _DG)
@test DT.compare_triangle_sets(T, _T) &&
      adjacent(adj) == adjacent(_adj) &&
      adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
      graph(DG) == graph(_DG)

p13 = (7.5, 3.5)
push!(pts, p13)
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 13)
_T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
DT.clear_empty_keys!(adj, DG)
DT.clear_empty_keys!(_adj, _DG)
@test DT.compare_triangle_sets(T, _T) &&
      adjacent(adj) == adjacent(_adj) &&
      adjacent2vertex(adj2v) == adjacent2vertex(_adj2v) &&
      graph(DG) == graph(_DG)

begin
    fig = Figure(fontsize=55)
    ax = Axis(fig[1, 1])
    triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0))
    xlims!(ax, 1.0, 14.0)
    ylims!(ax, -3.0, 7.0)
    hidedecorations!(ax)
    text!(ax,
        [p1, p2, p3, (p4[1] + 0.2, p4[2] - 0.1),
            (p5[1] + 0.1, p5[2] - 0.1),
            (p6[1] - 0.53, p6[2] - 0.2),
            (p7[1] - 0.46, p7[2] - 0.2),
            (p8[1] - 0.15, p8[2] - 1.1),
            (p9[1], p9[2] + 0.1),
            (p10[1], p10[2] - 1.0)];
        text=[L"%$s" for s in 1:10], textsize=55)
    scatter!(ax, p11, color=:blue, markersize=17)
    text!(ax, (p11[1] - 0.24, p11[2] - 1.2); text=L"11", textsize=55)
    fig
end

R = 10.0
n = 1381
pts = [(rand((-1,1)) * R * rand(), rand((-1,1)) * R * rand()) for _ in 1:n] # rand((-1, 1)) to get random signs
pushfirst!(pts, (-11.0,-11.0),(11.0,-11.0),(11.0,11.0),(-11.0,11.0)) # bounding box to guarantee interior insertions only
T, adj, adj2v, DG, HG = @views DT.triangulate_berg(pts[1:7])
_T, _adj, _adj2v, _DG, _HG = @views DT.triangulate_berg(pts[1:7])
n = length(pts)
for i in 8:n
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, i)
    DT.add_point_berg!(_T, _adj, _adj2v, _DG, _HG, pts, i)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
end


fig = Figure()
ax = Axis(fig[1, 1], title="BW")
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0))
ax = Axis(fig[1, 2], title="BERG")
triplot!(ax, _T, pts; strokewidth=2, color=(:white, 0.0))
fig