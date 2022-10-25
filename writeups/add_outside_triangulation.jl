include("functions.jl")

p1 = @SVector[-3.32, 3.53]
p2 = @SVector[-5.98, 2.17]
p3 = @SVector[-6.36, -1.55]
p4 = @SVector[-2.26, -4.31]
p5 = @SVector[6.34, -3.23]
p6 = @SVector[-3.24, 1.01]
p7 = @SVector[0.14, -1.51]
p8 = @SVector[0.2, 1.25]
p9 = @SVector[1.0, 4.0]
p10 = @SVector[4.74, 2.21]
p11 = @SVector[2.32, -0.27]
pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
DT.compute_centroid!(pts)
T = Set{NTuple{3,Int64}}()
adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
DG = DT.DelaunayGraph{Int64}()
for (i, j, k) in (
    (1, 2, 6),
    (1, 6, 8),
    (9, 1, 8),
    (9, 8, 10),
    (10, 8, 11),
    (8, 7, 11),
    (8, 6, 7),
    (6, 2, 3),
    (6, 3, 4),
    (6, 4, 7),
    (7, 4, 5),
    (11, 7, 5),
    (10, 11, 5)
)
    DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true)
end

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG)
xlims!(ax, -8.0, 8.0)
ylims!(ax, -8.0, 8.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1],
        @SVector[p3[1] - 0.5, p3[2]],
        @SVector[p4[1] + 0.2, p4[2] - 1.35],
        @SVector[p5[1] + 0.1, p5[2]],
        @SVector[p6[1] + 0.23, p6[2]],
        @SVector[p7[1] - 0.1, p7[2] - 1.7],
        @SVector[p8[1] + 0.25, p8[2] + 0.1],
        @SVector[p9[1] + 0.3, p9[2]],
        @SVector[p10[1] + 0.3, p10[2] - 0.6],
        @SVector[p11[1] + 0.45, p11[2] - 0.75]];
    text=[L"%$s" for s in 1:11], textsize=55)
p12 = @SVector[4.382, 3.2599]
push!(pts, p12)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] + 0.2, p12[2])]; text=L"12", textsize=55)
save("writeups/figures/point_to_be_added_outside_a.pdf", fig)

b = NTuple{3,Int64}[]
for _T in T
    if DT.isincircle(_T, pts, 12) == 1
        push!(b, _T)
    end
end

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, b[2:3], pts; markersize=0, color=(:lightblue, 0.5))
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG)
(; x, y) = DT.CentroidCoordinates
pc = @SVector[x, y]
g9 = 10 * pts[9] + (1 - 10) * pc
g10 = 10 * pts[10] + (1 - 10) * pc
poly!(ax, [pts[9], pts[10], g10, g9]; color=(:lightblue, 0.5))
xlims!(ax, -8.0, 8.0)
ylims!(ax, -8.0, 8.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1],
        @SVector[p3[1] - 0.5, p3[2]],
        @SVector[p4[1] + 0.2, p4[2] - 1.35],
        @SVector[p5[1] + 0.1, p5[2]],
        @SVector[p6[1] + 0.23, p6[2]],
        @SVector[p7[1] - 0.1, p7[2] - 1.7],
        @SVector[p8[1] + 0.25, p8[2] + 0.1],
        @SVector[p9[1] + 0.3, p9[2]],
        @SVector[p10[1] + 0.3, p10[2] - 0.6],
        @SVector[p11[1] + 0.45, p11[2] - 0.75]];
    text=[L"%$s" for s in 1:11], textsize=55)
p12 = @SVector[4.382, 3.2599]
push!(pts, p12)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] + 0.2, p12[2])]; text=L"12", textsize=55)
save("writeups/figures/point_to_be_added_outside_b.pdf", fig)

#=
DT.delete_triangle!(T, (9, 8, 10))
DT.delete_triangle!(T, (8, 11, 10))
DT.delete_triangle!(T, (9, 10, DT.BoundaryIndex))
DT.add_triangle!(T, (9, 8, 12))
DT.add_triangle!(T, (8, 11, 12))
DT.add_triangle!(T, (9, 12, DT.BoundaryIndex))
DT.add_triangle!(T, (10, 12, DT.BoundaryIndex))
DT.add_triangle!(T, (11, 10, 12))
=#
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 12)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG)
xlims!(ax, -8.0, 8.0)
ylims!(ax, -8.0, 8.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1],
        @SVector[p3[1] - 0.5, p3[2]],
        @SVector[p4[1] + 0.2, p4[2] - 1.35],
        @SVector[p5[1] + 0.1, p5[2]],
        @SVector[p6[1] + 0.23, p6[2]],
        @SVector[p7[1] - 0.1, p7[2] - 1.7],
        @SVector[p8[1] + 0.25, p8[2] + 0.1],
        @SVector[p9[1] + 0.3, p9[2]],
        @SVector[p10[1] + 0.3, p10[2] - 0.6],
        @SVector[p11[1] + 0.45, p11[2] - 0.75]];
    text=[L"%$s" for s in 1:11], textsize=55)
p12 = @SVector[4.382, 3.2599]
push!(pts, p12)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] + 0.2, p12[2])]; text=L"12", textsize=55)
lines!(ax, [p9, p12, @SVector[NaN, NaN],
        p8, p12, @SVector[NaN, NaN],
        p11, p12, @SVector[NaN, NaN],
        p10, p12, @SVector[NaN, NaN]],
    color=:blue, linewidth=4)
save("writeups/figures/point_to_be_added_outside_c.pdf", fig)
fig


include("functions.jl")

p1 = @SVector[-3.32, 3.53]
p2 = @SVector[-5.98, 2.17]
p3 = @SVector[-6.36, -1.55]
p4 = @SVector[-2.26, -4.31]
p5 = @SVector[6.34, -3.23]
p6 = @SVector[-3.24, 1.01]
p7 = @SVector[0.14, -1.51]
p8 = @SVector[0.2, 1.25]
p9 = @SVector[1.0, 4.0]
p10 = @SVector[4.74, 2.21]
p11 = @SVector[2.32, -0.27]
pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
DT.compute_centroid!(pts)
T = Set{NTuple{3,Int64}}()
adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
DG = DT.DelaunayGraph{Int64}()
for (i, j, k) in (
    (1, 2, 6),
    (1, 6, 8),
    (9, 1, 8),
    (9, 8, 10),
    (10, 8, 11),
    (8, 7, 11),
    (8, 6, 7),
    (6, 2, 3),
    (6, 3, 4),
    (6, 4, 7),
    (7, 4, 5),
    (11, 7, 5),
    (10, 11, 5)
)
    DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true)
end
p12 = @SVector[4.382, 3.2599]
push!(pts, p12)
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 12)
_T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
DT.add_ghost_triangles!(_T, _adj, _adj2v, _DG)
@test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
@test DT.compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)

p13 = @SVector[-5.25303, 4.761]
p14 = @SVector[-9.83801, 0.562]
p15 = @SVector[-7.15986, -5.99]
p16 = @SVector[4.79, 2.74]
p17 = @SVector[3.77, 2.7689]
push!(pts, p13, p14, p15, p16, p17)
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 13)
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 14)
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 15)
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 16)
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 17)
_T, _adj, _adj2v, _DG, _HG = DT.triangulate_berg(pts)
DT.add_ghost_triangles!(_T, _adj, _adj2v, _DG)
@test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
@test DT.compare_deberg_to_bowyerwatson(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)

