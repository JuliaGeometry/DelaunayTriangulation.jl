include("functions.jl")

qa = (-1.4, 0.61)
qb = (1.0, 1.0)
qc = (0.0, -2.0)
qd = (1.0, 3.0)
qe = (0.0, 0.0)
qf = (-0.56, 1.53)
qg = (2.22, 1.43)
qh = (2.46, -0.57)
qi = (-0.68, -1.07)
pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]
T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
cents, tri_to_idx, idx_to_tri = DT.circumcenters(T, pts)


all_cent_x = getindex.(collect(values(cents)), 1)
all_cent_y = getindex.(collect(values(cents)), 2)

fig = Figure()
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0))
scatter!(ax, all_cent_x, all_cent_y, color=:purple)
fig

a1
b2
c3
d4
e5
f6
g7
h8
i9

cell_order = Int64[]
i = 5
nghs = DT.get_edge(adj2v, i)
j, k = iterate(nghs)[1]
τ = (i, j, k)
cent_idx = tri_to_idx[τ]
push!(cell_order, j, cents[cent_idx])
j = DT.get_edge(adj, i, k)
push!(cell_order, j)
j = DT.get_edge(adj, i, j)
push!(cell_order, j)
j = DT.get_edge(adj, i, j)
push!(cell_order, j)
j = DT.get_edge(adj, i, j)
push!(cell_order, j)
j = DT.get_edge(adj, i, j)

cell_order = Int64[]
i = 5
j = 6# = rand(DT.get_neighbour(DG, i))
k = DT.get_edge(adj, i, j)
τ = (i, j, k)
cent_idx = tri_to_idx[τ]
push!(cell_order, cent_idx)
j = k
k = DT.get_edge(ajd, i, j)
τ = (i, j, k)

i = 5
num_ngh = DT.num_neighbours(DG, i)
cell_idx = zeros(Int64, num_ngh)
j = rand(DT.get_neighbour(DG, i))
for r in 1:num_ngh
    k = DT.get_edge(adj, i, j)
    τ = (i, j, k)
    if τ ∉ T
        τ = DT.shift_triangle_1(τ)
    end
    if τ ∉ T
        τ = DT.shift_triangle_1(τ)
    end
    cent_idx = tri_to_idx[τ]
    cell_idx[r] = cent_idx
    j = k
end

qa = (-2.84, 1.28)
qb = (1.2, 3.8)
qc = (5.72, 2.74)
qd = (10.22, -3.72)
qe = (6.26, -2.38)
qf = (7.0, -6.0)
qg = (2.54, -3.90)
qh = (-4.48, -7.7)
qi = (-7.9, -4.46)
qj = (-2.42, -2.14)
qk = (-7.92, 0.42)
pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi, qj, qk]
T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
cents = DT.circumcenters(T, pts)
all_cent_x = getindex.(collect(values(cents)), 1)
all_cent_y = getindex.(collect(values(cents)), 2)

fig = Figure()
ax = Axis(fig[1, 1])
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0))
scatter!(ax, all_cent_x, all_cent_y, color=:purple)
fig



Random.seed!(2992991)
qa = (-1.4, 0.61)
qb = (1.0, 1.0)
qc = (0.0, -2.0)
qd = (1.0, 3.0)
qe = (0.0, 0.0)
qf = (-0.56, 1.53)
qg = (2.22, 1.43)
qh = (2.46, -0.57)
qi = (-0.68, -1.07)
pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]


pts = 2rand(SVector{2,Float64}, 22)
T, adj, adj2v, _DG = DT.triangulate_bowyer(pts)
vorn = DT.voronoi(T, adj, adj2v, _DG, pts)

fig = Figure()
ax = Axis(fig[1, 1])
mkp = voronoiplot!(ax, vorn, pts, _DG; markersize=2, strokewidth=2)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=6)
xlims!(ax, 0, 2)
ylims!(ax, 0, 2)
fig

a = 0.0
b = 2.0
c = 0.0
d = 2.0
nx = 5
ny = 5
T, adj, adj2v, DG, pts = triangulate_structured(a, b, c, d, nx, ny)
T, adj, adj2v, DG = triangulate_bowyer(pts)
vorn = voronoi(T, adj, adj2v, DG, pts)
fig = Figure()
ax = Axis(fig[1, 1])
mkp = voronoiplot!(ax, vorn, pts, DG; markersize=2, strokewidth=2)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=6)
xlims!(ax, a - 1 / 2, b + 1 / 2)
ylims!(ax, c - 1 / 2, d + 1 / 2)
fig




Random.seed!(2992991)
qa = (-1.4, 0.61)
qb = (1.0, 1.0)
qc = (0.0, -2.0)
qd = (1.0, 3.0)
qe = (0.0, 0.0)
qf = (-0.56, 1.53)
qg = (2.22, 1.43)
qh = (2.46, -0.57)
qi = (-0.68, -1.07)
pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]

T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=true)

fig = Figure()
ax = Axis(fig[1, 1])
mkp = voronoiplot!(ax, vorn, pts, _DG; markersize=2, strokewidth=2)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=6)
xlims!(ax, -4, 4)
ylims!(ax, -4, 4)
fig

scatter!(ax, vorn.circumcenters, color=:purple)

fig

vorn.polygons[1]

scatter!(ax, pts[1])
fig




Random.seed!(2992991)
qa = (-1.4, 0.61)
qb = (1.0, 1.0)
qc = (-0.5, -1.5)
qd = (1.0, 3.0)
qe = (0.0, 0.0)
qf = (-0.56, 1.23)
qg = (2.22, 1.43)
qh = (3.46, -0.2)
qi = (-0.68, -1.07)
pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]

T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=false)

fig = Figure()
ax = Axis(fig[1, 1])
mkp = voronoiplot!(ax, vorn, pts, _DG; markersize=0, strokewidth=0.5)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=0)
xlims!(ax, -4, 4)
ylims!(ax, -4, 4)
fig

scatter!(ax, vorn.circumcenters, color=:purple)

fig

vorn.polygons[1]

scatter!(ax, pts[1])
fig

j = DT.find_circumcenters_outside_convex_hull(vorn, adj, adj2v)

scatter!(ax, vorn.circumcenters[j], color=:black, markersize=17)
fig

Random.seed!(2992991)
qa = (-1.4, 0.61)
qb = (1.0, 1.0)
qc = (-0.5, -1.5)
qd = (1.0, 3.0)
qe = (0.0, 0.0)
qf = (-0.56, 1.23)
qg = (2.22, 1.43)
qh = (3.46, -0.2)
qi = (-0.68, -1.07)
pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]

T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=true)
_vorn = deepcopy(vorn)
DT.truncate_bounded_polygon_outside_convex_hull(vorn, adj, adj2v, _DG)
DT.trim_voronoi_cell!(vorn, adj, adj2v, _DG)


j = DT.find_circumcenters_outside_convex_hull(vorn, adj, adj2v)
_j = 10

bounded_polygon = DT.find_bounded_polygons(vorn, _j; find_first=true)
p = vorn.polygons[bounded_polygon]
_j_idx = findfirst(p .== _j)
next_idx = _j_idx == lastindex(p) ? first(p) : p[_j_idx+1]
prev_idx = _j_idx == firstindex(p) ? last(p) : p[_j_idx-1]
boundary_τ = vorn.idx_to_triangle[_j]
u, v, w = indices(boundary_τ)
ub, vb, wb = DT.is_boundary_point.((u, v, w), Ref(adj), Ref(_DG))
if ub && vb
    u, v = u, v
elseif ub && wb
    u, v = w, u
elseif vb && wu
    u, v = v, w
elseif ub && vb && wb
    u, v = u, v
end
pᵤ, pᵥ = get_point(pts, u, v)
pⱼ = vorn.circumcenters[_j]
pᵣ = vorn.circumcenters[next_idx]
pₛ = vorn.circumcenters[prev_idx]

coord_1 = DT.intersection_of_two_line_segments(pᵤ, pᵥ, pⱼ, pᵣ)
coord_2 = DT.intersection_of_two_line_segments(pᵤ, pᵥ, pⱼ, pₛ)
coord_1_idx = length(vorn.circumcenters) + 1
coord_2_idx = length(vorn.circumcenters) + 2
push!(vorn.circumcenters, coord_1)
push!(vorn.circumcenters, coord_2)

other_polys = setdiff(vorn.circumcenter_to_polygons[_j], bounded_polygon)
P, Q = other_polys
P_verts = vorn.polygons[P]
Q_verts = vorn.polygons[Q]
delete!(vorn.circumcenter_to_polygons[_j], P)
delete!(vorn.circumcenter_to_polygons[_j], Q)
P_j = findfirst(P_verts .== _j)
Q_j = findfirst(Q_verts .== _j)

if next_idx ∈ P_verts#⟹ prev_idx ∈ Q
    P_verts[P_j] = coord_1_idx
    DT.add_polygon!(vorn.circumcenter_to_polygons, coord_1_idx, P)
    Q_verts[Q_j] = coord_2_idx
    DT.add_polygon!(vorn.circumcenter_to_polygons, coord_2_idx, Q)
else#prev_idx ∈ P_verts ⟹ next_idx ∈ Q
    P_verts[P_j] = coord_2_idx
    DT.add_polygon!(vorn.circumcenter_to_polygons, coord_2_idx, P)
    Q_verts[Q_j] = coord_1_idx
    DT.add_polygon!(vorn.circumcenter_to_polygons, coord_1_idx, Q)
end

p[_j_idx] = coord_2
p[next_idx]

vorn.circumcenter_to_polygons[_j]
P_verts
Q_verts
lines!(ax, vorn.circumcenters[[3, 16, 17]], linewidth=4, color=:white)
lines!(ax, vorn.circumcenters[[18, 4, 2]], linewidth=4, color=:white)


Random.seed!(29291919919)
fig = Figure()
ax = Axis(fig[1, 1])
mkp = voronoiplot!(ax, vorn, pts, _DG; markersize=0, strokewidth=0.5)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=0)
xlims!(ax, -4, 4)
ylims!(ax, -4, 4)
fig
ax = Axis(fig[1, 2])
Random.seed!(29291919919)
mkp = voronoiplot!(ax, _vorn, pts, _DG; markersize=0, strokewidth=0.5)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=0)
xlims!(ax, -4, 4)
ylims!(ax, -4, 4)
fig

vorn.polygons[1]
vorn.polygons[2]
vorn.polygons[3]
vorn.polygons[4]
vorn.polygons[5]
vorn.polygons[6]
vorn.polygons[7]
vorn.polygons[8]
vorn.polygons[9]



lines!(ax, vorn.circumcenters[[20, 12, 1, 25]], color=:white, linewidth=4)
lines!(ax, vorn.circumcenters[[8, 16, 3, 14, 5, 8]], color=:red, linewidth=4)
lines!(ax, vorn.circumcenters[[23, 4, 19]], color=:blue, linewidth=4)
lines!(ax, vorn.circumcenters[[26, 5, 14]], color=:green, linewidth=4)
lines!(ax, vorn.circumcenters[[22, 16, 8, 1, 12, 4, 23, 22]], color=:orange, linewidth=4)
lines!(ax, vorn.circumcenters[[25, 1, 8, 5, 26]], color=:black, linewidth=4)
lines!(ax, vorn.circumcenters[[14, 3]], color=:brown, linewidth=4)
lines!(ax, vorn.circumcenters[[3, 16, 22]], color=:magenta, linewidth=4)
lines!(ax, vorn.circumcenters[[19, 4, 12, 20, 19]], color=:purple, linewidth=4)

lines!(ax, vorn.circumcenters[vorn.polygons[1]], color=:red, linewidth=4)

scatter!(ax, [pⱼ], color=:blue)
fig
scatter!(ax, [pᵤ, pᵥ], color=:green)
scatter!(ax, [pᵣ, pₛ], color=:magenta)
scatter!(ax, [coord_1, coord_2], color=:black)
fig


[lines!(ax, vorn.circumcenters[vorn.polygons[i]], linewidth=6) for i in keys(vorn.polygons)]
fig

Random.seed!(2992991)
qa = (-1.4, 0.61)
qb = (1.0, 1.0)
qc = (-0.5, -1.5)
qd = (1.0, 3.0)
qe = (0.0, 0.0)
qf = (-0.56, 1.23)
qg = (2.22, 1.43)
qh = (3.46, -0.2)
qi = (-0.68, -1.07)
pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]
pts = rand(SVector{2,Float64}, 20)
pts = Tuple.(pts)

T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=false)
j = DT.find_circumcenters_outside_convex_hull(vorn, adj, adj2v)

Random.seed!(29291919919)
fig = Figure()
ax = Axis(fig[1, 1])
mkp = voronoiplot!(ax, vorn, pts, _DG; markersize=0, strokewidth=0.5)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=0, linecolor=:blue)
xlims!(ax, -3, 3)
ylims!(ax, -3, 3)
scatter!(ax, vorn.circumcenters[collect(j)], color=:black, markersize=14)
fig


function triangle_to_mat(T)
    _T = collect(T)
    V = [_T[i][j] for i in eachindex(_T), j in 1:3]
    return V
end
function triangle_area(p, q, r)
    return 0.5 * (p[1] * q[2] + q[1] * r[2] + r[1] * p[2] - p[1] * r[2] - r[1] * q[2] - q[1] * p[2])
end
function triangle_area(T, r)
    i, j, k = indices(T)
    pᵢ, pⱼ, pₖ = get_point(r, i, j, k)
    return triangle_area(pᵢ, pⱼ, pₖ)
end
function compute_cell_densities!(q, T, r::AbstractMatrix{F}) where {F}
    for T in T
        A = triangle_area(T, r)
        i, j, k = indices(T)
        q[i] += A
        q[j] += A
        q[k] += A
    end
    q .= 1.0 ./ q
    return nothing
end
function compute_cell_densities(T, r::AbstractMatrix{F}) where {F}
    q = zeros(F, size(r, 2))
    compute_cell_densities!(q, T, r)
    return q
end

Random.seed!(2993381881)
a, b, c, d, nx, ny = 0.0, 10.0, 0.0, 10.0, 10, 10
T, adj, adj2v, _DG, pts = triangulate_structured(a, b, c, d, nx, ny)
_T, _adj, _adj2v, __DG = triangulate_bowyer(pts)
vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=false)
_vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=true)

q1 = 1.0 ./ compute_cell_densities(T, pts)
q2 = 1.0 ./ compute_cell_densities(_T, pts)

fig = Figure(fontsize=43, resolution=(1600, 800))
ax = Axis(fig[1, 1],
    xlabel=L"x", ylabel=L"y",
    title=L"(a):$ $ Structured", titlealign=:left,
    xticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    yticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    aspect=1)
xlims!(ax, a, b)
ylims!(ax, c, d)
mesh!(ax, pts, triangle_to_mat(T), color=q1, colormap=:matter, colorrange=(0, 10))
triplot!(ax, T, pts; strokewidth=0.6, color=(:white, 0.0), plot_ghost_edges=false, markersize=0, linecolor=:blue)

ax2 = Axis(fig[1, 2],
    xlabel=L"x", ylabel=L"y",
    title=L"(b):$ $ Randomised", titlealign=:left,
    xticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    yticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    aspect=1)
xlims!(ax2, a, b)
ylims!(ax2, c, d)
mesh!(ax2, pts, triangle_to_mat(_T), color=q2, colormap=:matter, colorrange=(0, 10))
triplot!(ax2, _T, pts; strokewidth=0.6, color=(:white, 0.0), plot_ghost_edges=false, markersize=0, linecolor=:blue)

Colorbar(fig[1, 3], limits=(0, 5), colormap=:matter, vertical=true, label=L"$ $Area", ticks=(0:5:10, [L"%$s" for s in 0:5:10]), labelsize=63)

resize_to_layout!(fig)

save("writeups/figures/density_comparisons_areas.pdf", fig)

fig = Figure(fontsize=43, resolution=(1600, 800))
ax = Axis(fig[1, 1],
    xlabel=L"x", ylabel=L"y",
    title=L"(a):$ $ Unbounded", titlealign=:left,
    xticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    yticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    aspect=1)
xlims!(ax, a - 1, b + 1)
ylims!(ax, c - 1, d + 1)
Random.seed!(2993381881)
mkp = voronoiplot!(ax, vorn, pts, _DG; markersize=0, strokewidth=0)
triplot!(ax, T, pts; strokewidth=0.6, color=(:white, 0.0), plot_ghost_edges=false, markersize=0, linecolor=:blue)

ax2 = Axis(fig[1, 2],
    xlabel=L"x", ylabel=L"y",
    title=L"(b):$ $ Bounded", titlealign=:left,
    xticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    yticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    aspect=1)
xlims!(ax, a - 1, b + 1)
ylims!(ax, c - 1, d + 1)
Random.seed!(2993381881)
mkp = voronoiplot!(ax2, _vorn, pts, _DG; markersize=0, strokewidth=0)
triplot!(ax2, T, pts; strokewidth=0.6, color=(:white, 0.0), plot_ghost_edges=false, markersize=0, linecolor=:blue)

save("writeups/figures/dvoronoi_tessellations_areas.pdf", fig)

Random.seed!(2993381881)
a, b, c, d, nx, ny = 0.0, 10.0, 0.0, 10.0, 10, 10
T, adj, adj2v, _DG, pts = triangulate_structured(a, b, c, d, nx, ny)
_T, _adj, _adj2v, __DG = triangulate_bowyer(pts)
vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=true)
_vorn = DT.voronoi(_T, _adj, _adj2v, __DG, pts; trim=true)

q1 = DT.area(vorn)
q2 = DT.area(_vorn)

fig = Figure(fontsize=43, resolution=(1600, 800))
ax = Axis(fig[1, 1],
    xlabel=L"x", ylabel=L"y",
    title=L"(a):$ $ Structured", titlealign=:left,
    xticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    yticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    aspect=1)
xlims!(ax, a, b)
ylims!(ax, c, d)
mesh!(ax, pts, triangle_to_mat(T), color=q1, colormap=:matter, colorrange=(0, 10))
triplot!(ax, T, pts; strokewidth=0.6, color=(:white, 0.0), plot_ghost_edges=false, markersize=0, linecolor=:blue)

ax2 = Axis(fig[1, 2],
    xlabel=L"x", ylabel=L"y",
    title=L"(b):$ $ Randomised", titlealign=:left,
    xticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    yticks=(0:5:10, [L"%$s" for s in 0:5:10]),
    aspect=1)
xlims!(ax2, a, b)
ylims!(ax2, c, d)
mesh!(ax2, pts, triangle_to_mat(_T), color=q2, colormap=:matter, colorrange=(0, 10))
triplot!(ax2, _T, pts; strokewidth=0.6, color=(:white, 0.0), plot_ghost_edges=false, markersize=0, linecolor=:blue)

Colorbar(fig[1, 3], limits=(0, 5), colormap=:matter, vertical=true, label=L"$ $Area", ticks=(0:5:10, [L"%$s" for s in 0:5:10]), labelsize=63)

resize_to_layout!(fig)

save("writeups/figures/density_comparisons_areas_voronoi.pdf", fig)








Random.seed!(2992991)
qa = (-1.4, 0.61)
qb = (1.0, 1.0)
qc = (-0.5, -1.5)
qd = (1.0, 3.0)
qe = (0.0, 0.0)
qf = (-0.56, 1.23)
qg = (2.22, 1.43)
qh = (3.46, -0.2)
qi = (-0.68, -1.07)
pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]
pts = rand(SVector{2,Float64}, 20)
pts = Tuple.(pts)

T, adj, adj2v, _DG = DT.triangulate_bowyer(pts; trim=false)
vorn = DT.voronoi(T, adj, adj2v, _DG, pts; trim=true)
j = DT.find_circumcenters_outside_convex_hull(vorn, T, adj, adj2v, _DG)

Random.seed!(29291919919)
fig = Figure()
ax = Axis(fig[1, 1])
mkp = voronoiplot!(ax, vorn, pts, _DG; markersize=0, strokewidth=0.5)
triplot!(ax, T, pts; strokewidth=0.5, color=(:white, 0.0), plot_ghost_edges=false, markersize=0, linecolor=:blue)
xlims!(ax, -3, 3)
ylims!(ax, -3, 3)
scatter!(ax, vorn.circumcenters[collect(j)], color=:black, markersize=14)
fig
