include("functions.jl")

rmin = 0.01
rmax = 10.0
pts = SVector{2,Float64}[]
t = LinRange(0.0001, 5.0, 500)
θmin = 0.0
θmax = 8π
for (i, τ) in enumerate(t)
    r = rmin * (1 - τ) + rmax * τ
    θ = θmin * (1 - τ) + θmax * τ
    rei = sqrt(r) * exp(2im * θ)
    x = real(rei)
    y = imag(rei)
    push!(pts, @SVector[x, y])
end

T, adj, adj2v, DG = DT.triangulate_bowyer(pts[1:3]; randomise=false, trim=false)

fig = Figure()
ax = Axis(fig[1, 1])
record(fig, "writeups/figures/animation.mp4", 4:length(pts); framerate=27) do r
    empty!(ax)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, r; initial_search_point=r - 1)
    triplot!(ax, T, pts; strokewidth=0.2,
        color=(:white, 0.0), markersize=1,
        plot_ghost_edges=true, DG=DG,
        recompute_centroid=false)
    hidedecorations!(ax)
    xlims!(ax, -12, 12)
    ylims!(ax, -12, 12)
end

pts = 25 .* rand(SVector{2, Float64}, 500)
T, adj, adj2v, DG = DT.triangulate_bowyer(pts[1:3]; randomise=false, trim=false)
fig = Figure()
ax = Axis(fig[1, 1])
record(fig, "writeups/figures/animation_shuffled.mp4", 4:length(pts); framerate=27) do r
    empty!(ax)
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, r; initial_search_point=r - 1)
    triplot!(ax, T, pts; strokewidth=0.2,
        color=(:white, 0.0), markersize=1,
        plot_ghost_edges=true, DG=DG,
        recompute_centroid=false)
    hidedecorations!(ax)
    xlims!(ax, -2, 27)
    ylims!(ax, -2, 27)
end