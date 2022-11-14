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

