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

n = length(pts)
m = ceil(Int64, n^(1 / 3))
Random.seed!(292988888881)
random_edges = rand(edges(DG), m)
colors = [:red, :green, :magenta, :purple]
for k in 1:m
    i, j = random_edges[k]
    pi, pj = pts[i], pts[j]
    midpt = (pi + pj) / 2
    push!(dists, (midpt[1] - q[1])^2 + (midpt[2] - q[2])^2)
    lines!(ax, [pi, pj], color=colors[k], linewidth=5)
    scatter!(ax, midpt, color=colors[k], markersize=14)
    lines!(ax, [midpt, q], color=colors[k], linewidth=5)
end
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
i, j = random_edges[3]
pi, pj = pts[i], pts[j]
midpt = (pi + pj) / 2
scatter!(ax, midpt, color=colors[3], markersize=14)
scatter!(ax, [q[1]], [q[2]], color=:blue, markersize=14)
midpts = [
    (pts[24] + pts[1]) / 2,
    (pts[28] + pts[1]) / 2,
    (pts[1] + pts[36]) / 2,
    (pts[36] + pts[48]) / 2,
    (pts[36] + pts[4]) / 2,
    q
]
marked = [
    (24, 1, 28),
    (28, 1, 36),
    (1, 48, 36),
    (36, 48, 4),
    (36, 4, 3)
]
triplot!(ax, marked, pts; strokewidth=2, color=(:lightblue, 0.5), markersize=0)
lines!(ax, [midpts[1], q], color=(colors[1],0.5), linewidth=5)
lines!(ax, midpts, color=colors[3], linewidth=5)
text!(ax, [q[1] - 0.2], [q[2] - 0.6]; text=L"q", textsize=38)
save("writeups/figures/point_location_marked_simplices.pdf", fig)

xlims!(ax, -4.5, 1.0)
ylims!(ax, -2.0, 3.0)
scatter!(ax, pts, color = :red, markersize = 11)
text!(ax, pts; text = [L"%$s" for s in 1:n], textsize=55, color = :blue)
save("writeups/figures/point_location_marked_simplices_zoomed.pdf", fig)

