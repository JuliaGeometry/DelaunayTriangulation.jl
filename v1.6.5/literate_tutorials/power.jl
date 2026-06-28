# # Power Diagrams 
# 

# In this tutorial, we demonstrate how we can construct 
# power diagrams (also called weighted Voronoi tessellations).
# These are dual to the weighted Delaunay triangulation and, 
# instead of being based on the Euclidean metric like Voronoi tessellations,
# the power distance is used to define the Voronoi tiles. The power distance 
# is defined by $\pi(p, q) = d(p, q)^2 - w_p - w_q$, where $d(p, q)$
# is the Euclidean distance between points $p$ and $q$, and $w_p$ and $w_q$
# are _weights_ associated with $p$ and $q$. See [this page](../math/power.md)
# for more details. To start with the tutorial, we load in the packages we'll need. 
using DelaunayTriangulation
using CairoMakie
using StableRNGs
using StatsBase #src
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

# To build a power diagram, you need to start from a weighted 
# Delaunay triangulation as described in [this tutorial](weighted.md).
points = [
    (-3.0, 7.0), (1.0, 6.0), (-1.0, 3.0),
    (-2.0, 4.0), (3.0, -2.0), (5.0, 5.0),
    (-4.0, -3.0), (3.0, 8.0),
]
weights = [0.2, 0.0, -3.0, 2.0, 7.5, 2.3, 0.0, -0.5]
rng = StableRNG(123)
tri = triangulate(points; weights, rng)
pvorn = voronoi(tri; rng)

# Let's compare the power diagram to its unweighted counterpart.
vorn = voronoi(triangulate(points; rng); rng)
fig = Figure()
ax1 = Axis(fig[1, 1], title="Unweighted", width=300, height=400)
ax2 = Axis(fig[1, 2], title="Weighted", width=300, height=400)
voronoiplot!(ax1, vorn, show_generators=true, colormap=:matter, strokewidth=4, clip=(-5, 5, -5, 10))
voronoiplot!(ax2, pvorn, show_generators=true, colormap=:matter, strokewidth=4, clip=(-5, 5, -5, 10))
resize_to_layout!(fig)
fig
@test_reference joinpath(fig_path, "power_diagram_ex_1.png") fig #src

# Notice that, unlike the unweighted tessellation, the generators in the power diagram 
# don't actually have to live in their own tile, and more than one generator
# may inhibit a tile. This is because the power diagram is based on the power distance,
# not the Euclidean distance.

# We can easily look at the affect of the weights on the power diagram.
A, B, C, D, E, F, G, H,
I, J, K, L, M, N = (-1.0, 3.0), (1.0, 3.0),
(2.0, 2.0), (2.0, 0.0), (1.0, -1.0), (-1.0, -1.0),
(-2.0, 0.0), (-2.0, 2.0), (-1.0, 2.0),
(0.0, 1.5), (1.0, 2.5), (-1.0, 0.5),
(1.0, 0.0), (1.5, 1.0)
points = [A, B, C, D, E, F, G, H, I, J, K, L, M, N]
w = Observable(-10.0)
weights = @lift (wts = zeros(length(points)); wts[10] = $w; wts)
tri = @lift tri = triangulate(points; weights=$weights)
vorn = @lift voronoi($tri)
weight_itr_base = LinRange(-10, 10, 30 * 5)
weight_itr = vcat(weight_itr_base, reverse(weight_itr_base))
title_obs = lift(w -> L"w_{10} = %$(round(w, sigdigits = 4))", w)
fig, ax, sc = voronoiplot(vorn,
    axis=(title=title_obs, titlealign=:left),
    figure=(fontsize=24,),
    clip = (-4, 4, -4, 4),
    color = [:red, :blue, :green, :purple, :orange, :yellow, :cyan, :white, :black, :magenta, :gray, :brown, :pink, :lightblue]
)
scatter!(ax, [J], color=:red, markersize=13)
record(fig, "varying_weight_power.mp4", weight_itr; framerate=30) do _w
    w[] = _w
end;

# ![](varying_weight_power.mp4)

# See that, for small weights, the Voronoi tile of the 10th point isn't even present. As the weight increases, the tile 
# grows and eventually dominates the tessellation.

# We also note that, just like standard Voronoi tessellations, we can also 
# apply clipping and smoothing to power diagrams. These can be done in exactly the same manner.