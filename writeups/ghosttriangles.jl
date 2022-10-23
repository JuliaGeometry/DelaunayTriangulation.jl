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
centr = mean(pts)
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

fig = Figure(fontsize=24)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0))
xlims!(ax, -8.0, 8.0)
ylims!(ax, -8.0, 8.0)
#hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1],
        @SVector[p3[1] - 0.5, p3[2]],
        @SVector[p4[1] + 0.2, p4[2] - 0.7],
        @SVector[p5[1] + 0.1, p5[2]],
        @SVector[p6[1] + 0.23, p6[2]],
        @SVector[p7[1] - 0.1, p7[2] - 0.75],
        @SVector[p8[1] + 0.25, p8[2] + 0.1],
        @SVector[p9[1] + 0.3, p9[2]],
        @SVector[p10[1] + 0.3, p10[2] - 0.6],
        @SVector[p11[1] + 0.35, p11[2] - 0.25]];
    text=[L"%$s" for s in 1:11], textsize=24)
p12 = @SVector[3.382, 4.3599]
push!(pts, p12)
scatter!(ax, p12, color=:darkgreen, markersize=11)
scatter!(ax, [centr[1]], [centr[2]], color=:blue, markersize=11)
for i in [1, 9, 10, 5, 4, 3, 2]
    p = pts[i]
    pt2 = 10 * p + (1 - 10) * centr
    lines!(ax, [centr, p, pt2], color=:blue, linestyle=:dash, linewidth=2)
end
text!(ax, [(centr[1], centr[2] - 0.9)]; text=L"p_c", textsize=24)
text!(ax, [(p12[1] + 0.2, p12[2])]; text=L"q", textsize=24)
fig
save("writeups/figures/ghost_triangles.pdf", fig)

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
centr = mean(pts)
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
