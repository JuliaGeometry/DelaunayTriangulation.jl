# Mesh Refinement

Now we describe how we refine meshes. We do not discuss curve-bounded domains here, leaving this discussion to the [curve-bounded section](curve_bounded.md). The mesh refinement algorithm we implement is Ruppert's Delaunay refinement algorithm with some modifications from Shewchuk, following the presentation in Chapter 6 of the book [_Delaunay Mesh Generation_ by Cheng, Dey, and Shewchuk (2013)](https://people.eecs.berkeley.edu/~jrs/meshbook.html). 

To introduce the mesh refinement algorithm, we need to discuss a few preliminaries. We let $\mathcal M$ denote a mesh $\mathcal D\mathcal T(\mathcal P, \mathcal S, \mathcal B)$, where $\mathcal P$ is the point set, $\mathcal S$ is the segment set, and $\mathcal B$ is the boundary.

## Triangle Quality 

We need a way to measure the "quality" of a triangle in order to know what triangles need refinement. The quality of a triangle is typically measured using the _radius-edge ratio_, which is the ratio of the radius of the triangle's circumference to its minimmu edge length. Symbolically, $\rho = R/\ell_{min}$, where $R$ is the triangle's circumradius, $\ell_{\min}$ its minimum edge length, and $\rho$ is the radius-edge ratio. The lower the radius-edge ratio, the better we say the triangle's quality is. We can also view reducing $\rho$ as increasing the minimum angle, $\theta_{\min}$, since $\rho$ is related to $\theta_{\min}$ by 

```math 
\rho = \frac{R}{\ell_{\min}} = \frac{1}{2\sin\theta_{\min}}.
```

We say that a triangle is of bad quality if $\rho > \bar\rho$ for some $\bar\rho > 0$ or, equivalently, $\theta_{\min} < \bar\theta$. Typically, $\bar\theta \approx 33.9^{\circ}$ is about the limit of convergence for the algorithm. In addition to controlling $\rho$, we can also control the area $A$ of the triangle, saying a triangle is too large if $A > A_{\max}$ and too small if $A < A_{\min}$. 

Thus, given a triangle $T$ and constants $\bar\rho$, $A_{\min}$, and $A_{\max}$, we say that $T$ is of bad quality if any of the following are true: (1) $\rho > \bar \rho$ or (2) $A > A_{\max}$. For $A_{\min}$, we of course can't split $T$ in order to increase its area. Instead, the $A_{\min}$ threshold will be used to stop the splitting of any bad quality triangle if its area would be reduced to below $A_{\min}$, and any triangle whose area is already below $A_{\min}$ will never be refined further.

## Edge Encroachment 

Mesh refinement algorithms typically refine bad quality triangles by inserting their circumcenter into the triangulation (there are other variants, but we will not discuss them here). There is one problem with this approach: the circumcenter of a triangle might lie outside of $\mathcal M$, and so we would never be able to insert it. To get around this, we introduce the concept of _edge encroachment_. We say that a vertex $v$ _encroaches upon_ a segment $e$ if it is inside the closed diametral circle of $e$ but is not a vertex of $e$, and $v$ is visible from $e$.

Checking for encroachment is relatively straightforward: The only vertices that might encroach upon a segment $e$ are those adjacent to it. In particular, let the segment be $e = e_{uv}$ and find the triangles $T_{uvw}$ and $T_{uvx}$ adjacent to $e$. If either $p_w$ or $p_x$ are in the diametral circle of $e_{uv}$, then $e$ is encroached. Checking whether a point is inside a diametral circle is also simple: By Thales' theorem, the angle $\angle p_upp_v$ is a right angle when $p$ is on the diametral circle, so $p$ is inside the diametral circle if and only if $\angle p_upp_v \leq 90^{\circ}$. To efficiently check if $\theta = \angle p_upp_v \leq 90^{\circ}$, we can use a dot product. Since 

```math 
\cos \theta = \frac{\langle p - p_u, p - p_v \rangle}{\|p - p_u\|\|p - p_v\|}
```

and $\cos 90^{\circ} = 0$, we can check if $\theta \leq 90^{\circ}$ by checking if $\langle p - p_u, p - p_v \rangle \geq 0$. Thus:

1. If $\langle p - p_u, p - p_v\rangle > 0$, then $\theta < 90^{\circ}$.
2. If $\langle p - p_u, p - p_v\rangle = 0$, then $\theta = 90^{\circ}$.
3. If $\langle p - p_u, p - p_v\rangle < 0$, then $\theta > 90^{\circ}$.

The first two cases would imply that $p$ encroaches upon $e_{uv}$.

Once we have identifed an edge as being encroached, we need to split it. We do this by simply inserting the midpoint $p_m = (p_u + p_v) / 2$ into the triangulation, and replacing the segment $e_{uv}$ by the two new subsegments, legalising any new edges as needed to restore the Delaunay property.

One issue with this definition of encroachment is that there may be a large number of vertices that have to be inserted to deal with all encroached edges as the diametral circle is quite large. To avoid this, we generalise the definition of a diametral circle to that of a _diametral lens_ defined by some lens angle $\theta_\ell$. Firstly, let $p$ be a point on the perpendicular bisector $L_{uv}$ of $e_{uv}$, and left of $e_{uv}$, such that $\angle p_upp_v = \theta_\ell$, and similarly for a point $q$ right of $L$. Now draw two circles through $p_u, p_v, p$ and $p_u, p_v, q$, respectively, and compute their intersection. The diametral lens is this intersection. We show an example of a diametral lens below.

```@setup refex 
using DelaunayTriangulation
const DT = DelaunayTriangulation
using LinearAlgebra
using CairoMakie 
using StableRNGs 
# from helper_functions.jl in test
function reflect_point_across_line(p, q, u)
    orig_p = p
    q = collect(q .- p)
    u = collect(u .- p)
    p = zeros(2)
    pq_θ = -atan(q[2], q[1])
    cwrot_θ = [cos(pq_θ) -sin(pq_θ); sin(pq_θ) cos(pq_θ)]
    rotu = cwrot_θ * u
    # pq is now a horizontal line on the x-axis. So, to reflect rotu about this line, just do 
    rotu = [rotu[1], -rotu[2]]
    # Now, rotate back
    rotu = cwrot_θ' * rotu
    u = rotu .+ orig_p
    return Tuple(u)
end
function compute_diametral_circle(p, q)
    m = (p .+ q) ./ 2
    θ = LinRange(0, 2π, 250)
    rad = norm(p .- q) ./ 2
    circle = tuple.(m[1] .+ rad * cos.(θ), m[2] .+ rad * sin.(θ))
    lower_circle = circle[1:125, :] |> vec
    upper_circle = circle[126:end, :] |> vec
    unique!(circle)
    push!(circle, circle[begin])
    return upper_circle, lower_circle, circle
end
function compute_diametral_lens(p, q, lens_angle)
    orig_p = p
    orig_q = q
    q = collect(q .- p)
    p = zeros(2)
    pq_θ = -atan(q[2], q[1])
    cwrot_θ = [cos(pq_θ) -sin(pq_θ); sin(pq_θ) cos(pq_θ)]
    p = cwrot_θ * p
    q = cwrot_θ * q
    m = (p .+ q) ./ 2
    pq_normal = (q[2] - p[2], q[1] - p[1])
    pq_normal = pq_normal ./ norm(pq_normal)
    a = norm(p .- q) / 2
    c = a / cosd(lens_angle)
    b = c * sind(lens_angle)
    u = m .+ b .* pq_normal
    circumcenter = DT.triangle_circumcenter(p, q, u)
    circumradius = DT.triangle_circumradius(p, q, u)
    p_angle = atan(p[2] - circumcenter[2], p[1] - circumcenter[1])
    q_angle = atan(q[2] - circumcenter[2], q[1] - circumcenter[1])
    θ = LinRange(p_angle, q_angle, 250)
    lens = hcat(circumcenter[1] .+ circumradius * cos.(θ), circumcenter[2] .+ circumradius * sin.(θ))
    reflected_lens = reflect_point_across_line.(Ref(p), Ref(q), eachrow(lens))
    reflected_lens = hcat(first.(reflected_lens), last.(reflected_lens))
    combined_lens = vcat(lens, reverse(reflected_lens))
    rerotated_and_shifted_lens = (cwrot_θ' * lens')'
    rerotated_and_shifted_reflected_lens = (cwrot_θ' * reflected_lens')'
    rerotated_and_shifted_combined_lens = (cwrot_θ' * combined_lens')'
    rerotated_and_shifted_lens = tuple.(rerotated_and_shifted_lens[:, 1] .+ orig_p[1], rerotated_and_shifted_lens[:, 2] .+ orig_p[2])
    rerotated_and_shifted_reflected_lens = tuple.(rerotated_and_shifted_reflected_lens[:, 1] .+ orig_p[1], rerotated_and_shifted_reflected_lens[:, 2] .+ orig_p[2])
    rerotated_and_shifted_combined_lens = vcat(rerotated_and_shifted_lens, reverse(rerotated_and_shifted_reflected_lens))
    unique!(rerotated_and_shifted_combined_lens)
    push!(rerotated_and_shifted_combined_lens, rerotated_and_shifted_combined_lens[begin])
    return rerotated_and_shifted_lens, rerotated_and_shifted_reflected_lens, rerotated_and_shifted_combined_lens
end
function circ(p, q, r)
    x1, y1 = getxy(p)
    x2, y2 = getxy(q)
    x3, y3 = getxy(r)
    x12, x13, x31, x21 = x1 - x2, x1 - x3, x3 - x1, x2 - x1
    y12, y13, y31, y21 = y1 - y2, y1 - y3, y3 - y1, y2 - y1
    sx13 = x1^2 - x3^2 
    sy13 = y1^2 - y3^2
    sx21 = x2^2 - x1^2
    sy21 = y2^2 - y1^2
    f = (sx13 * x12 + sy13 * x12 + sx21 * x13 + sy21 * x13) / (2 * (y31 * x12 - y21 * x13))
    g = (sx13 * y12 + sy13 * y12 + sx21 * y13 + sy21 * y13) / (2 * (x31 * y12 - x21 * y13))
    c = -x1^2 - y1^2 - 2g*x1 - 2f*y1 
    h = -g 
    k = -f
    r = sqrt(h^2 + k^2 - c)
    return (h, k, r)
end
p, q = (0.0, 0.0), (2.0, 0.0)
_, _, circle = compute_diametral_circle(p, q)
t = 30.0
_, _, lens = compute_diametral_lens(p, q, t)
fig = Figure()
ax = Axis(fig[1, 1], width = 300, height = 300)
lines!(ax, [p, q], color = :black, linewidth = 3)
lines!(ax, circle, color = :blue)
scatter!(ax, [p, q], color = :black)
xlims!(ax, -0.3, 2.3)
ylims!(ax, -1.3, 1.3)
m = (p .+ q) ./ 2 
vlines!(ax, [m[1]], color = (:grey, 0.4), linestyle = :dash)
text!(ax, [(1.05, -1.3)], text = [L"L"], fontsize = 36)
du = norm(p .- q) * tand(t) / 2
u = (1.0, du)
scatter!(ax, [u, (1.0, -du)], color = :red)
lines!(ax, [p, u], color = :red, linestyle = :dash)
lines!(ax, [p, (1.0, -du)], color = :red, linestyle = :dash)
text!(ax, [(0.5, 0.0), (0.5, -0.35)], text = [L"\theta_\ell", L"\theta_\ell"], fontsize = 36)
H, K, R = circ(p, q, u)
arc!(ax, (H, K), R, 0.0, 2pi, color = (:magenta, 0.7))
H, K, R = circ(p, q, (1.0, -du))
arc!(ax, (H, K), R, 0.0, 2pi, color = (:magenta, 0.7))
lines!(ax, lens, color = :red)
text!(ax, [p .- (0.3, 0.05), q .- (0.0, 0.2)], text = [L"p_u", L"p_v"], fontsize = 36)
text!(ax, [u .+ (0.05, 0.05)], text = [L"p"], fontsize = 36)
text!(ax, [(1.0, -du) .+ (0.05, 0.05)], text = [L"q"], fontsize = 36)
hidedecorations!(ax)
resize_to_layout!(fig)
```

```@example refex
fig #hide
```

In this figure, the blue circle shows the diametral circle, and the red shape shows the diametral lens, the grey line is the perpendicular bisector of $e_{uv}$, and the magenta circles show the two circles whose intersection defines the diametral lens. Checking if a point $p$ is inside the diametral lens can be done using an extended version of Thales' theorem (see, for example, Theorem 9 in [Lisboa's thesis](https://repositorio.ufmg.br/bitstream/1843/RHCT-7GMJR6/1/adriano_chaves_lisboa.pdf)). We define the quantity 

```math
\Delta(p_u, p_v, p) = \langle p_u - p, p_v - p \rangle^2 - \|p_u - p\|^2 \|p_v - p\|^2 \cos^2(\theta_\ell).
```

Then:

- If $\Delta(p_u, p_v, p) > 0$, then $p$ is inside the diametral lens.
- If $\Delta(p_u, p_v, p) = 0$, then $p$ is on the boundary of the diametral lens.
- If $\Delta(p_u, p_v, p) < 0$, then $p$ is outside the diametral lens.

One issue with diametral lens is that the final mesh is not guaranteed to be Delaunay, but many more subsegment splits will be avoided than if we had used diametral circles so that the final mesh has fewer triangles. It is possible that circumcenters are outside of the boundary when using diametral lens. In this case, we insert the triangle's centroid rather than its circumcenter.

## Deleting Free Vertices 

In addition to splitting a subsegment whenever it is encroached, we can also use an idea from Chew's algorithm for mesh refinement to improve our refinement. When we split a subsegment, points inside the original segment's diametral circle might cause unduly short edges to be created, leading to bad quality triangles that will just have to be split once again. To overcome this, we delete all _free vertices_ inside the diametral circle, except those that are not visible to the segment (i.e., a segment in $\mathcal S$ occludes the visibility between the segment and the vertex), before splitting the segment. A free vertex is any vertex not belonging to a segment that was inserted into the triangulation through refinement, i.e. a vertex that was not originally in $\mathcal P$ and is not on the boundary or on an interior segment.

To delete these free vertices, we apply the following routine to each of $e_{uv}$ and $e_{vu}$. We write these details for $e_{uv}$ only.

1. Get the vertex $w$ adjacent to $e_{uv}$ using the adjacent map.
2. If $w$ is a free vertex and is either inside or on the diametral circle of $e_{uv}$, delete $w$ from the triangulation and return to step 1. Otherwise, stop the routine here.

Using this idea, we can reduce the number of triangles present in the final refined mesh while still retaining a high quality output. Moreover, since we use diametral lens instead of diametral circles, we avoid introducing as many vertices that would just be deleted by this procedure anyway.

## Splitting a Triangle

Now let's discuss how we actually insert a circumcenter into the triangulation to improve the quality of a triangle $T_{uvw}$. The procedure is simple:

1. Let $c$ be the circumcenter of a bad quality $T_{uvw}$.
2. If $c$ encroaches upon some subsegment $e \in \mathcal S$, split $e$. Otherwise, insert $c$ into $\mathcal M$.

Inserting $c$ into $\mathcal M$ is simple using the Bowyer-Watson algorithm. Note that, for the point location step, we already know that the triangle $T_{uvw}$ is a triangle containg $c$ in its circumcircle, so we can skip the point location step.

Checking if $c$ encroaches upon a subsegment $e$ can be expensive. The cheapest way to do this is to actually just insert $c$ into $\mathcal M$, and simply check if any of the edges of triangles containing $c$ are encroahced. If they are, we undo the insertion of $c$ and return to the original $\mathcal M$ prior to the insertion of $c$, and split the marked encroached edges. To undo this insertion efficiently, we store a list of all changes to the triangulation made during the insertion of $c$.

## Small Angles

A very crucial issue to notice with our refinement algorithm thus far is that it may fail to handle small angles. The first problem is called _ping-pong encroachment_, encountered when segments share a vertex and meet at an angle less than $45^{\circ}$. Consider the example below; we illustrate this using diametral circles, but the same problem can be encountered when using diametral lenses.

```@setup refex 
p = (0.0, 0.0)
q = (1.0, 0.0)
r = (0.7, 0.4)
fig = Figure()
ax = Axis(fig[1, 1], width = 300, height = 300)
lines!(ax, [p, q], color = :black, linewidth = 3)
lines!(ax, [p, r], color = :black, linewidth = 3)
t = LinRange(0, 1, 250)
arc!(ax, (p .+ q) ./ 2, norm(p .- q) / 2, 0.0, 2pi, color = :red)
ax2 = Axis(fig[1, 2], width = 300, height = 300)
pq1 = (p .+ q) ./ 2 
lines!(ax2, [p, q], color = :black, linewidth = 3)
lines!(ax2, [p, r], color = :black, linewidth = 3)
scatter!(ax2, [pq1], color = :black, markersize = 17)
arc!(ax2, (p .+ r) ./ 2, norm(p .- r) / 2, 0.0, 2pi, color = :red)
ax3 = Axis(fig[1, 3], width = 300, height = 300)
pr1 = (p .+ r) ./ 2
lines!(ax3, [p, q], color = :black, linewidth = 3)
lines!(ax3, [p, r], color = :black, linewidth = 3)
scatter!(ax3, [pq1, pr1], color = :black, markersize = 17)
arc!(ax3, (p .+ pq1) ./ 2, norm(p .- pq1) / 2, 0.0, 2pi, color = :red)
ax4 = Axis(fig[1, 4], width = 300, height = 300)
ppq11 = (p .+ pq1) ./ 2
lines!(ax4, [p, q], color = :black, linewidth = 3)
lines!(ax4, [p, r], color = :black, linewidth = 3)
scatter!(ax4, [pq1, pr1, ppq11], color = :black, markersize = 17)
arc!(ax4, (p .+ pr1) ./ 2, norm(p .- pr1) / 2, 0.0, 2pi, color = :red)
ax5 = Axis(fig[1, 5], width = 300, height = 300)
lines!(ax5, [p, q], color = :black, linewidth = 3)
lines!(ax5, [p, r], color = :black, linewidth = 3)
m1 = (p .+ q) ./ 2
m2 = (p .+ r) ./ 2
points1 = [m1]
points2 = [m2]
for _ in 1:100
    push!(points1, (p .+ points1[end]) ./ 2)
    push!(points2, (p .+ points2[end]) ./ 2)
end
scatter!(ax5, points1, color = :black, markersize = 17)
scatter!(ax5, points2, color = :black, markersize = 17)
hidedecorations!.((ax, ax2, ax3, ax4, ax5))
xlims!.((ax, ax2, ax3, ax4, ax5), -0.1, 1.1)
ylims!.((ax, ax2, ax3, ax4, ax5), -0.3, 0.5)
resize_to_layout!(fig)
```

```@example refex
fig #hide
```

In the first figure, the bottom segment $e_1$ encroaches upon the vertex of the other adjoining segment $e_2$, so we split $e_1$ at its midpoint. Once we insert this new midpoint, the segment $e_2$ is encroached upon, and so we need to split $e_2$. The third figure then shows how the new segment of $e_1$ is encroached upon by the vertex adding onto $e_2$, and so yet again we must split this subsegment. We can continue this process again to obtain the fourth figure. This will repeat indefinitely, leading to many points added as shown in the last figure. This is the ping-pong encroachment problem.

To overcome this problem, we use _concentric circular shells_. The basic idea is to imagine that each input vertex is surrounded by concentric circles whose radii are all the powers of two. We will still always split a segment initially at its midpoint, but for any future subsegments we need to make use of the concentric circles. When an encroached subsegment adjoins another segment at an acute angle, we split it at one of the circular shells centred at the shared vertex, so that one of the new subsegments has a power of two length. The shell we choose to split at is the one that guarantees that the two new subsegments produced by the split are between $1/3$ and $2/3$ the length of the split subsegment. With this approach, we can avoid the ping-pong encroachment. Notice that the choice of having the shells be powers of two implies that, for any future splits, the most balanced split for the power-of-two length subsegment will always be at the midpoint. If both vertices of a segment adjoin other segments, then the segment could be split twice at each end. To deal with this, just chosoe one vertex arbitrarily and split it so that the subsegment adjoining that vertex has a power-of-two length between $1/4$ and $1/2$ the length of the split subsegment. The other subsegment could still undergo another off-center split, but eventually all subsegment splits are bisections. This solves our ping-pong encroachment problem since adjoining subsegments of equal length cannot encroach upon each other. An example of this adaptation is shown below.

```@setup refex
p = (0.0, 0.0)
q = (1.0, 0.0)
r = (0.7, 0.4)
fig = Figure()
ax = Axis(fig[1, 1], width = 300, height = 300)
for i in -7:7
    arc!(ax, p, 2.0^i, 0.0, 2pi, color = (:grey, 0.4))
end
xlims!(ax, -0.5, 1.5)
ylims!(ax, -1, 1)
lines!(ax, [p, q], color = :black, linewidth = 3)
lines!(ax, [p, r], color = :black, linewidth = 3)
t = LinRange(0, 1, 250)
arc!(ax, (p .+ q) ./ 2, norm(p .- q) / 2, 0.0, 2pi, color = :red)
ax2 = Axis(fig[1, 2], width = 300, height = 300)
for i in -7:7
    arc!(ax2, p, 2.0^i, 0.0, 2pi, color = (:grey, 0.4))
end
xlims!(ax2, -0.5, 1.5)
ylims!(ax2, -1, 1)
pq1 = (p .+ q) ./ 2 
lines!(ax2, [p, q], color = :black, linewidth = 3)
lines!(ax2, [p, r], color = :black, linewidth = 3)
scatter!(ax2, [pq1], color = :black, markersize = 17)
arc!(ax2, (p .+ r) ./ 2, norm(p .- r) / 2, 0.0, 2pi, color = :red)
ax3 = Axis(fig[1, 3], width = 300, height = 300)
for i in -7:7
    arc!(ax3, p, 2.0^i, 0.0, 2pi, color = (:grey, 0.4))
end
xlims!(ax3, -0.5, 1.5)
ylims!(ax3, -1, 1)
pr1 = (p .+ r) ./ 2
lines!(ax3, [p, q], color = :black, linewidth = 3)
lines!(ax3, [p, r], color = :black, linewidth = 3)
scatter!(ax3, [pq1, pr1], color = :black, markersize = 17)
arc!(ax3, (p .+ pq1) ./ 2, norm(p .- pq1) / 2, 0.0, 2pi, color = :red)
ax4 = Axis(fig[1, 4], width = 300, height = 300)
for i in -7:7
    arc!(ax4, p, 2.0^i, 0.0, 2pi, color = (:grey, 0.4))
end
xlims!(ax4, -0.5, 1.5)
ylims!(ax4, -1, 1)
ppq11 = p .+ (1 - DelaunayTriangulation.compute_concentric_shell_ternary_split_position(p, pq1)) .* pq1
lines!(ax4, [p, q], color = :black, linewidth = 3)
lines!(ax4, [p, r], color = :black, linewidth = 3)
scatter!(ax4, [pq1, pr1, ppq11], color = :black, markersize = 17)
arc!(ax4, (p .+ pr1) ./ 2, norm(p .- pr1) / 2, 0.0, 2pi, color = :red)
ax5 = Axis(fig[1, 5], width = 300, height = 300)
for i in -7:7
    arc!(ax5, p, 2.0^i, 0.0, 2pi, color = (:grey, 0.4))
end
xlims!(ax5, -0.5, 1.5)
ylims!(ax5, -1, 1)
ppr11 = p .+ DelaunayTriangulation.compute_concentric_shell_ternary_split_position(p, pr1) .* pr1
lines!(ax5, [p, q], color = :black, linewidth = 3)
lines!(ax5, [p, r], color = :black, linewidth = 3)
scatter!(ax5, [pq1, pr1, ppq11, ppr11], color = :black, markersize = 17)
arc!(ax5, (p .+ ppq11) ./ 2, norm(p .- ppq11) / 2, 0.0, 2pi, color = :red)
hidedecorations!.((ax, ax2, ax3, ax4, ax5))
resize_to_layout!(fig)
```

```@example refex
fig #hide
```

In the second and third figures we have two midpoint splits since the segments are the input segments. For the subsegment in the fourth figure, we see that the new point is being put onto the concentric circles surrounding the input vertex, and similarly for the last figure. In the last figure, we finally see that the newly inserted vertx is no longer encroaching upon the other subsegment, and so the ping-pong encroachment stops.

One other improvement made by Ruppert is to avoid splitting triangles that are nestled in the corner of a small input angle. For a triangle $T_{uvw}$, suppose that $e_{uv}$ is its shortest edge so that the smallest angle of $T_{uvw}$ is opposite $e_{uv}$. Then, if $e_{wu}$ and $e_{wv}$ are both segments and the triangle is skinny, it is considered to be a nestled triangle and so the triangle will never be split.

The last improvement we consider involves _seditious triangles_. If two adjoining subsegments meet at a very small angle, then splitting them may lead to a new edge that is shorter than the previously shortest edge in the mesh, leading to a bad quality mesh. Moreover, this short edge will cause more triangles to be refined as the adjoining triangles will necessarily be skinny, leading to more short edges, thus leading to an infinite loop. To avoid this, we need to prevent these short edges from causing more refinement. We say that an edge is _seditious_ if its two vertices lie on two distinct segments that meet each other at an angle less than $60^{\circ}$ (in this package, the default definition for a seditious edge actually uses an angle of $20^{\circ}$), they lie on the same concentric shell (we don't check this requirement in this package), and the two vertices are true midpoints (not off-center splits). With this definition, we say that a triangle is _seditious_ if its shortest edge is seditious, and refuse to split any skinny triangle that is seditious, thus preventing seditious edges from infesting the rest of the mesh.

## The Complete Algorithm

Now that we have an understanding of all the pieces involved in the refinement algorithm, we can list the complete algorithm.

1. Start by identifying all encroached segments and placing them into a priority queue, prioritising the longer segments first.
2. For each encroached subsegment $e$: Delete all free vertices in the diametral circle (or lens) of $e$ and then split $e$ at a position depending on whether $e$ is an input segment or meets another subsegment at a small angle.
3. Next, identify all triangles that need to be refined, i.e. any triangle $T_{uvw}$ with $\rho > \bar\rho$ or $A > A_{\max}$, ignoring any of those which are nestled or seditious or $A < A_{\min}$. Store these triangles in a priority queue, prioriting the triangles with the largest radius-edge ratio first.
4. Next, while there are any bad quality triangles: Attempt to split the bad quality triangle $T$ by inserting its circumcenter $c$ into $\mathcal M$ (or centroid, if $c$ is outside of the domain in case diametral lenses are used). If $c$ encroaches on any new edges, undo the insertion and then split all those encroached segments as in step 2. If the insertion was successful, check all the newly added triangles for bad quality and add them to the priority queue if needed.
5. Once there are no more bad quality triangles to split, the algorithm is complete.