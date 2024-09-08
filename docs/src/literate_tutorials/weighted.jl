# # Weighted Triangulations
#

# In this tutorial, we demonstrate how to construct a weighted Delaunay triangulation.
# For more information about weighted Delaunay triangulations, see the [math details here](../math/weighted.md).
# To summarise, weighted Delaunay triangulations associate with each vertex $p_i$ a scalar 
# weight $w_i$. To pass weights into `triangulate`, use the `weights` keyword argument, 
# which by default is `ZeroWeight()`.

# Let's first consider a simple example. We use a triangulation with random weights and compare 
# it to the standard unweighted Delaunay triangulation.
using DelaunayTriangulation, CairoMakie
A, B, C, D, E, F, G, H,
I, J, K, L, M, N = (-1.0, 3.0), (1.0, 3.0),
(2.0, 2.0), (2.0, 0.0), (1.0, -1.0), (-1.0, -1.0),
(-2.0, 0.0), (-2.0, 2.0), (-1.0, 2.0),
(0.0, 1.5), (1.0, 2.5), (-1.0, 0.5),
(1.0, 0.0), (1.5, 1.0)
points = [A, B, C, D, E, F, G, H, I, J, K, L, M, N]
weights = randn(length(points))
tri = triangulate(points; weights)

#-
tri_unweighted = triangulate(points)

#-
fig = Figure()
ax = Axis(fig[1, 1], width = 400, height = 400, title = "Unweighted", titlealign = :left)
triplot!(ax, tri_unweighted, strokecolor = :black)
ax2 = Axis(fig[1, 2], width = 400, height = 400, title = "Weighted", titlealign = :left)
triplot!(ax2, tri, strokecolor = :red)
resize_to_layout!(fig) 
fig

# The triangles are of course not the same. Also, note that not all 
# vertices appear in the weighted triangulation. To see the effect of the 
# weights, let's reset all the weights to zero and see how the triangulation changes 
# as we vary the weight of the point at $(x, y) = (0, 1.5)$, starting with a 
# weight of $-10$.
w = Observable(-10.0)
weights = @lift (wts = zeros(length(points)); wts[10] = $w; wts)
tri = @lift triangulate(points; weights = $weights)
weight_itr_base = LinRange(-10, 10, 30*5)
weight_itr = vcat(weight_itr_base, reverse(weight_itr_base))
title_obs = lift(w -> L"w_{10} = %$(round(w, sigdigits = 4))", w)
fig, ax, sc = triplot(tri, 
axis = (title = title_obs, titlealign = :left),
figure = (fontsize = 24,))
scatter!(ax, [J], color = :red, markersize = 13)
record(fig, "varying_weight.mp4", weight_itr; framerate = 30) do _w 
    w[] = _w
end;

# ![](varying_weight.mp4)

# See that, once the weight gets so large, it essentially dominates the triangulation.