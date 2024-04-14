# Centroidal Voronoi Tessellations

We now discuss centroidal Voronoi tessellations. These are tessellations in which the generator of each Voronoi polygon is given by the polygon's centroid. The algorithm we use for this is very simple. Suppose we have a point set $\mathcal P$ and we want to compute the centroidal Voronoi tessellation of $\mathcal P$, denoted $\mathcal C\mathcal V(\mathcal P)$. We do the following:

1. Let $\mathcal P' = \mathcal P$.
2. Compute the clippped Voronoi tessellation $\tilde{\mathcal V}(\mathcal P') = \mathcal V(\mathcal P') \cap \mathcal C\mathcal H(\mathcal P')$.

    Then, until terminating, we:

2. For each generator $v_i' \in \mathcal P'$ not on the boundary, move $v_i'$ to the centroid $v_i''$ of $\mathcal V_i(\mathcal P')$.
3. Compute $\delta = \max_i \|p_i' - p_i''\|$, i.e. the maximum displacement of a point.
4. If $\delta < \varepsilon$, terminate. Otherwise, set $\mathcal P' = \{v_1'', \ldots, v_n''\}$ and go to step 2 after computing $\tilde{\mathcal V}(\mathcal P')$ once again.. Here, $\varepsilon$ is some displacement tolerance. In this package, we default to $\varepsilon = 10^{-4}$, where $h$ is the maximum edge length of the bounding box box $\tilde{\mathcal V}(\mathcal P)$.

That is the entire description of the algorithm. We give an example of this below.

```@setup centroidal 
using DelaunayTriangulation
using CairoMakie 
using StableRNGs
rng = StableRNG(12365)
points = 5randn(rng, 2, 125)
tri = triangulate(points; rng)
vorn = voronoi(tri, clip = true, rng = rng)
centroidal_smooth(vorn)
vorn_obs = Observable(vorn)
point_obs = Observable(points)
iter = 0
F = DelaunayTriangulation.number_type(vorn)
max_dist = typemax(F)
tri = DelaunayTriangulation.get_triangulation(vorn)
has_ghost = DelaunayTriangulation.has_ghost_triangles(tri)
!has_ghost && DelaunayTriangulation.add_ghost_triangles!(tri)
has_bnds = DelaunayTriangulation.has_boundary_nodes(tri)
!has_bnds && DelaunayTriangulation.lock_convex_hull!(tri)
set_of_boundary_nodes = DelaunayTriangulation.get_all_boundary_nodes(tri)
points = (deepcopy ∘ get_points)(tri)
fig = Figure()
ax = Axis(fig[1, 1], width = 300, height = 300)
hidedecorations!(ax)
voronoiplot!(ax, vorn_obs, show_generators = false)
scatter!(ax, point_obs, color = :black, markersize = 7)
resize_to_layout!(fig)
record(fig, "centroidal_voronoi.mp4", 1:200; framerate = 24) do _ # actual number of iterations needed is around 168
    global vorn, points
    vorn, _max_dist = DelaunayTriangulation._centroidal_smooth_itr(vorn, set_of_boundary_nodes, points, rng)
    vorn_obs[] = vorn 
    point_obs[] = points
end
```

![](centroidal_voronoi.mp4)
