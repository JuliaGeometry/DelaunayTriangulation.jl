include("functions.jl")

Random.seed!(25500551)
r = 5sqrt.(rand(50))
θ = 2π * rand(50)
pts = [@SVector[r * cos(θ), r * sin(θ)] for (r, θ) in zip(r, θ)]

T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0), markersize=0)
q = @SVector[0.7, 2.0]
text!(ax, [q[1] - 0.4], [q[2] - 0.75]; text=L"q", textsize=38)
scatter!(ax, [q[1]], [q[2]], color=:blue, markersize=14)
xlims!(ax, -5, 5.0)
ylims!(ax, -5.0, 5.0)
hidedecorations!(ax)
save("writeups/figures/point_location_initialisation.pdf", fig)
fig

n = length(pts)
m = ceil(Int64, n^(1 / 3))
Random.seed!(2889682)
random_pts = rand(1:n, m)
s = pts[random_pts]
colors = [:red, :green, :magenta, :purple]
[lines!(ax, [s[i], q], color=colors[i], linewidth=5) for i in 1:4]
[scatter!(ax, s[i], color=colors[i], markersize=14) for i in 1:4]
scatter!(ax, [q[1]], [q[2]], color=:blue, markersize=14)
xlims!(ax, -5, 5.0)
ylims!(ax, -5.0, 5.0)
hidedecorations!(ax)
save("writeups/figures/point_location_selected_vertices.pdf", fig)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0), markersize=0)
q = @SVector[0.7, 2.0]
xlims!(ax, -5, 5.0)
ylims!(ax, -5.0, 5.0)
hidedecorations!(ax)
scatter!(ax, s[1], color=colors[1], markersize=14)
marked = [(46, 38, 41),
    (38, 3, 41),
    (38, 34, 3),
    (34, 27, 3),
    (27, 36, 3),
    (36, 4, 3)]
triplot!(ax, marked, pts; strokewidth=2, color=(:lightblue, 0.5), markersize=0)
lines!(ax, [s[1], q], color=colors[1], linewidth=5)
tq = text!(ax, [q[1] - 0.4], [q[2] - 0.6]; text=L"q", textsize=38)
scatter!(ax, [q[1]], [q[2]], color=:blue, markersize=14)
save("writeups/figures/point_location_marked_simplices.pdf", fig)

xlims!(ax, -2.0, 3.0)
ylims!(ax, 0.8, 5.3)
text!(ax, pts; text=[L"%$s" for s in 1:n], textsize=55, color=:red)
delete!(ax, tq)
tq = text!(ax, [q[1] - 0.2], [q[2] - 0.3]; text=L"q", textsize=38)
save("writeups/figures/point_location_marked_simplices_zoomed.pdf", fig)

## Example initialisation 
pts = [
    (0.0, 3.0),
    (2.264, 1.475),
    (4.285, -1.0482),
    (3.2818, -3.123),
    (-1.127, -4.0599),
    (-3.65, -1.875),
    (-2.9178, 1.759),
    (0.0, 0.0)
]
T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
q = (4.09, 4.1882)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; markersize=11, strokewidth=2, color=(:white, 0))
xlims!(ax, -6, 6.0)
ylims!(ax, -6.0, 6.0)
hidedecorations!(ax)
text!(ax, [(pts[6][1] - 0.7, pts[6][2] - 0.2),
        (pts[7][1] - 0.4, pts[7][2] + 0.1),
        (pts[8][1] - 0.1, pts[8][2] - 1.2)]; text=[L"p_j", L"p_i", L"p_k"], textsize=55)
ablines!(ax, [0.0], [q[2] / q[1]], color=:blue, linewidth=3, linestyle=:dash)
scatter!(ax, q, color=:blue, markersize=11)
lines!(ax, [pts[end], q], color=:blue, linewidth=3)
scatter!(ax, pts[end], color=:blue, markersize=11)
text!(ax, (q[1] + 0.1, q[2] - 0.5); text=[L"q"], textsize=55)
triplot!(ax, [(8, 7, 6)], pts; strokewidth=2, color=(:lightblue, 0.5))
save("writeups/figures/point_location_clockwise_a.pdf", fig)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; markersize=11, strokewidth=2, color=(:white, 0))
xlims!(ax, -6, 6.0)
ylims!(ax, -6.0, 6.0)
hidedecorations!(ax)
text!(ax, [pts[1],
        (pts[7][1] - 0.4, pts[7][2] + 0.1),
        (pts[8][1] - 0.1, pts[8][2] - 1.2)]; text=[L"p_i", L"p_j", L"p_k"], textsize=55)
ablines!(ax, [0.0], [q[2] / q[1]], color=:blue, linewidth=3, linestyle=:dash)
scatter!(ax, q, color=:blue, markersize=11)
lines!(ax, [pts[end], q], color=:blue, linewidth=3)
scatter!(ax, pts[end], color=:blue, markersize=11)
text!(ax, (q[1] + 0.1, q[2] - 0.5); text=[L"q"], textsize=55)
triplot!(ax, [(8, 7, 1)], pts; strokewidth=2, color=(:lightblue, 0.5))
save("writeups/figures/point_location_clockwise_b.pdf", fig)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; markersize=11, strokewidth=2, color=(:white, 0))
xlims!(ax, -6, 6.0)
ylims!(ax, -6.0, 6.0)
hidedecorations!(ax)
text!(ax, [pts[1],
        pts[2],
        (pts[8][1] - 0.1, pts[8][2] - 1.2)]; text=[L"p_j", L"p_i", L"p_k"], textsize=55)
triplot!(ax, [(8, 1, 2)], pts; strokewidth=2, color=(:lightblue, 0.5))
ablines!(ax, [0.0], [q[2] / q[1]], color=:blue, linewidth=3, linestyle=:dash)
scatter!(ax, q, color=:blue, markersize=11)
lines!(ax, [pts[end], q], color=:blue, linewidth=3)
scatter!(ax, pts[end], color=:blue, markersize=11)
text!(ax, (q[1] + 0.1, q[2] - 0.5); text=[L"q"], textsize=55)
save("writeups/figures/point_location_clockwise_c.pdf", fig)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; markersize=11, strokewidth=2, color=(:white, 0))
xlims!(ax, -6, 6.0)
ylims!(ax, -6.0, 6.0)
hidedecorations!(ax)
ablines!(ax, [0.0], [q[2] / q[1]], color=:blue, linewidth=3, linestyle=:dash)
scatter!(ax, q, color=:blue, markersize=11)
lines!(ax, [pts[end], q], color=:blue, linewidth=3)
scatter!(ax, pts[end], color=:blue, markersize=11)
text!(ax, (q[1] + 0.1, q[2] - 0.5); text=[L"q"], textsize=55)
triplot!(ax, [(8, 3, 4)], pts; strokewidth=2, color=(:lightblue, 0.5))
text!(ax, [(pts[4][1] + 0.2, pts[4][2] - 0.5),
        (pts[3]),
        (pts[8][1] - 0.1, pts[8][2] - 1.2)]; text=[L"p_i", L"p_j", L"p_k"], textsize=55)
save("writeups/figures/point_location_counterclockwise_a.pdf", fig)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; markersize=11, strokewidth=2, color=(:white, 0))
xlims!(ax, -6, 6.0)
ylims!(ax, -6.0, 6.0)
hidedecorations!(ax)
ablines!(ax, [0.0], [q[2] / q[1]], color=:blue, linewidth=3, linestyle=:dash)
scatter!(ax, q, color=:blue, markersize=11)
lines!(ax, [pts[end], q], color=:blue, linewidth=3)
scatter!(ax, pts[end], color=:blue, markersize=11)
text!(ax, (q[1] + 0.1, q[2] - 0.5); text=[L"q"], textsize=55)
triplot!(ax, [(8, 3, 4)], pts; strokewidth=2, color=(:lightblue, 0.5))
text!(ax, [(pts[4][1] + 0.2, pts[4][2] - 0.5),
        (pts[3]),
        (pts[8][1] - 0.1, pts[8][2] - 1.2)]; text=[L"p_i", L"p_j", L"p_k"], textsize=55)
save("writeups/figures/point_location_counterclockwise_b.pdf", fig)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; markersize=11, strokewidth=2, color=(:white, 0))
xlims!(ax, -6, 6.0)
ylims!(ax, -6.0, 6.0)
hidedecorations!(ax)
text!(ax, [pts[1],
        pts[2],
        (pts[8][1] - 0.1, pts[8][2] - 1.2)]; text=[L"p_j", L"p_i", L"p_k"], textsize=55)
triplot!(ax, [(8, 1, 2)], pts; strokewidth=2, color=(:lightblue, 0.5))
ablines!(ax, [0.0], [q[2] / q[1]], color=:blue, linewidth=3, linestyle=:dash)
scatter!(ax, q, color=:blue, markersize=11)
lines!(ax, [pts[end], q], color=:blue, linewidth=3)
scatter!(ax, pts[end], color=:blue, markersize=11)
text!(ax, (q[1] + 0.1, q[2] - 0.5); text=[L"q"], textsize=55)
save("writeups/figures/point_location_counterclockwise_c.pdf", fig)

