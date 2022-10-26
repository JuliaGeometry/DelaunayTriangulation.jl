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
DT.compute_centroid!(pts)
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

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG)
xlims!(ax, -8.0, 8.0)
ylims!(ax, -8.0, 8.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1],
        @SVector[p3[1] - 0.5, p3[2]],
        @SVector[p4[1] + 0.2, p4[2] - 1.35],
        @SVector[p5[1] + 0.1, p5[2]],
        @SVector[p6[1] + 0.23, p6[2]],
        @SVector[p7[1] - 0.1, p7[2] - 1.7],
        @SVector[p8[1] + 0.25, p8[2] + 0.1],
        @SVector[p9[1] + 0.3, p9[2]],
        @SVector[p10[1] + 0.3, p10[2] - 0.6],
        @SVector[p11[1] + 0.45, p11[2] - 0.75]];
    text=[L"%$s" for s in 1:11], textsize=55)
p12 = @SVector[4.382, 3.2599]
push!(pts, p12)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] + 0.2, p12[2])]; text=L"12", textsize=55)
save("writeups/figures/point_to_be_added_outside_a.pdf", fig)

b = NTuple{3,Int64}[]
for _T in T
    if DT.isincircle(_T, pts, 12) == 1
        push!(b, _T)
    end
end

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, b[2:3], pts; markersize=0, color=(:lightblue, 0.5))
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG)
(; x, y) = DT.CentroidCoordinates
pc = @SVector[x, y]
g9 = 10 * pts[9] + (1 - 10) * pc
g10 = 10 * pts[10] + (1 - 10) * pc
poly!(ax, [pts[9], pts[10], g10, g9]; color=(:lightblue, 0.5))
xlims!(ax, -8.0, 8.0)
ylims!(ax, -8.0, 8.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1],
        @SVector[p3[1] - 0.5, p3[2]],
        @SVector[p4[1] + 0.2, p4[2] - 1.35],
        @SVector[p5[1] + 0.1, p5[2]],
        @SVector[p6[1] + 0.23, p6[2]],
        @SVector[p7[1] - 0.1, p7[2] - 1.7],
        @SVector[p8[1] + 0.25, p8[2] + 0.1],
        @SVector[p9[1] + 0.3, p9[2]],
        @SVector[p10[1] + 0.3, p10[2] - 0.6],
        @SVector[p11[1] + 0.45, p11[2] - 0.75]];
    text=[L"%$s" for s in 1:11], textsize=55)
p12 = @SVector[4.382, 3.2599]
push!(pts, p12)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] + 0.2, p12[2])]; text=L"12", textsize=55)
save("writeups/figures/point_to_be_added_outside_b.pdf", fig)

DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 12)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG)
xlims!(ax, -8.0, 8.0)
ylims!(ax, -8.0, 8.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1],
        @SVector[p3[1] - 0.5, p3[2]],
        @SVector[p4[1] + 0.2, p4[2] - 1.35],
        @SVector[p5[1] + 0.1, p5[2]],
        @SVector[p6[1] + 0.23, p6[2]],
        @SVector[p7[1] - 0.1, p7[2] - 1.7],
        @SVector[p8[1] + 0.25, p8[2] + 0.1],
        @SVector[p9[1] + 0.3, p9[2]],
        @SVector[p10[1] + 0.3, p10[2] - 0.6],
        @SVector[p11[1] + 0.45, p11[2] - 0.75]];
    text=[L"%$s" for s in 1:11], textsize=55)
p12 = @SVector[4.382, 3.2599]
push!(pts, p12)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] + 0.2, p12[2])]; text=L"12", textsize=55)
lines!(ax, [p9, p12, @SVector[NaN, NaN],
        p8, p12, @SVector[NaN, NaN],
        p11, p12, @SVector[NaN, NaN],
        p10, p12, @SVector[NaN, NaN]],
    color=:blue, linewidth=4)
save("writeups/figures/point_to_be_added_outside_c.pdf", fig)
fig

## collinear points 
Random.seed!(288282) # there is a set of cocircular points, hence non-unique, so let's set the seed
p1 = (-6.0, 8.0)
p2 = (0.0, 10.0)
p3 = (-6.0, 0.0)
p4 = (-2.0, 0.0)
p5 = (2.0, 8.0)
p6 = (-8.0, 6.0)
p7 = (4.0, 4.0)
p8 = (6.0, 6.0)
p9 = (10.0, -2.0)
p10 = (2.0, 6.0)
p11 = (-4.0, -4.0)
pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
T, adj, adj2v, DG = DT.triangulate_bowyer(pts)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG, recompute_centroid=false)
xlims!(ax, -10.0, 12.0)
ylims!(ax, -5.0, 12.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1].data,
        @SVector[p3[1] - 1.5, p3[2]].data,
        @SVector[p4[1] - 0.4, p4[2] - 1.7].data,
        @SVector[p5[1], p5[2]].data,
        @SVector[p6[1] - 0.65, p6[2]].data,
        @SVector[p7[1] - 0.6, p7[2] - 1.9].data,
        @SVector[p8[1] + 0.25, p8[2] + 0.1].data,
        @SVector[p9[1] + 0.3, p9[2]].data,
        @SVector[p10[1] - 3.1, p10[2] - 1.4].data,
        @SVector[p11[1] - 2.3, p11[2] - 0.75].data];
    text=[L"%$s" for s in 1:11], textsize=55)
p12 = (0.0, 8.0)
push!(pts, p12)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] - 2.0, p12[2])]; text=L"12", textsize=55)
fig
save("writeups/figures/collinear_point_to_be_added_outside_a.pdf", fig)

b = NTuple{3,Int64}[]
for _T in T
    if DT.isincircle(_T, pts, 12) == 1
        push!(b, _T)
    end
end

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, b, pts; markersize=0, color=(:lightblue, 0.5))
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG, recompute_centroid=false)
#=
(; x, y) = DT.CentroidCoordinates
pc = (x, y)
g9 = 10 .* pts[9] + (1 - 10) .* pc
g10 = 10 .* pts[10] + (1 - 10) .* pc
poly!(ax, [pts[9], pts[10], g10, g9]; color=(:lightblue, 0.5))
=#
xlims!(ax, -10.0, 12.0)
ylims!(ax, -5.0, 12.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1].data,
        @SVector[p3[1] - 1.5, p3[2]].data,
        @SVector[p4[1] - 0.4, p4[2] - 1.7].data,
        @SVector[p5[1], p5[2]].data,
        @SVector[p6[1] - 0.65, p6[2]].data,
        @SVector[p7[1] - 0.6, p7[2] - 1.9].data,
        @SVector[p8[1] + 0.25, p8[2] + 0.1].data,
        @SVector[p9[1] + 0.3, p9[2]].data,
        @SVector[p10[1] - 3.1, p10[2] - 1.4].data,
        @SVector[p11[1] - 2.3, p11[2] - 0.75].data];
    text=[L"%$s" for s in 1:11], textsize=55)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] - 2.0, p12[2])]; text=L"12", textsize=55)
save("writeups/figures/collinear_point_to_be_added_outside_b.pdf", fig)

DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 12)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG, recompute_centroid=false)
xlims!(ax, -10.0, 12.0)
ylims!(ax, -5.0, 12.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1].data,
        @SVector[p3[1] - 1.5, p3[2]].data,
        @SVector[p4[1] - 0.4, p4[2] - 1.7].data,
        @SVector[p5[1], p5[2]].data,
        @SVector[p6[1] - 0.65, p6[2]].data,
        @SVector[p7[1] - 0.6, p7[2] - 1.9].data,
        @SVector[p8[1] + 0.25, p8[2] + 0.1].data,
        @SVector[p9[1] + 0.3, p9[2]].data,
        @SVector[p10[1] - 3.1, p10[2] - 1.4].data,
        @SVector[p11[1] - 2.3, p11[2] - 0.75].data];
    text=[L"%$s" for s in 1:11], textsize=55)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] - 2.0, p12[2])]; text=L"12", textsize=55)
lines!(ax, [p2, p12, (NaN, NaN),
        p1, p12, (NaN, NaN),
        p5, p12, (NaN, NaN),
        p4, p12, (NaN, NaN),
        p10, p12, (NaN, NaN)],
    color=:blue, linewidth=4)
save("writeups/figures/collinear_point_to_be_added_outside_c.pdf", fig)
fig

## boundary collinear 
Random.seed!(288282) # there is a set of cocircular points, hence non-unique, so let's set the seed
p1 = (-6.0, 8.0)
p2 = (0.0, 10.0)
p3 = (-6.0, 0.0)
p4 = (-2.0, 0.0)
p5 = (2.0, 8.0)
p6 = (-8.0, 6.0)
p7 = (4.0, 4.0)
p8 = (6.0, 6.0)
p9 = (10.0, -2.0)
p10 = (2.0, 6.0)
p11 = (-4.0, -4.0)
pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12]
T, adj, adj2v, DG = DT.triangulate_bowyer(pts)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG, recompute_centroid=false)
xlims!(ax, -10.0, 12.0)
ylims!(ax, -5.0, 12.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1].data,
        @SVector[p3[1] - 1.5, p3[2]].data,
        @SVector[p4[1] - 0.4, p4[2] - 1.7].data,
        @SVector[p5[1], p5[2]].data,
        @SVector[p6[1] - 0.65, p6[2]].data,
        @SVector[p7[1] - 0.6, p7[2] - 1.9].data,
        @SVector[p8[1] + 0.25, p8[2] + 0.1].data,
        @SVector[p9[1] + 0.3, p9[2]].data,
        @SVector[p10[1] - 3.1, p10[2] - 1.4].data,
        @SVector[p11[1] - 2.3, p11[2] - 0.75].data];
    text=[L"%$s" for s in 1:11], textsize=55)
p12 = (0.0, 8.0)
push!(pts, p12)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] - 2.0, p12[2])]; text=L"12", textsize=55)
fig
save("writeups/figures/collinear_point_to_be_added_outside_d.pdf", fig)

b = NTuple{3,Int64}[]
for _T in T
    if DT.isincircle(_T, pts, 12) == 1
        push!(b, _T)
    end
end

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, b, pts; markersize=0, color=(:lightblue, 0.5))
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG, recompute_centroid=false)
#=
(; x, y) = DT.CentroidCoordinates
pc = (x, y)
g9 = 10 .* pts[9] + (1 - 10) .* pc
g10 = 10 .* pts[10] + (1 - 10) .* pc
poly!(ax, [pts[9], pts[10], g10, g9]; color=(:lightblue, 0.5))
=#
xlims!(ax, -10.0, 12.0)
ylims!(ax, -5.0, 12.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1].data,
        @SVector[p3[1] - 1.5, p3[2]].data,
        @SVector[p4[1] - 0.4, p4[2] - 1.7].data,
        @SVector[p5[1], p5[2]].data,
        @SVector[p6[1] - 0.65, p6[2]].data,
        @SVector[p7[1] - 0.6, p7[2] - 1.9].data,
        @SVector[p8[1] + 0.25, p8[2] + 0.1].data,
        @SVector[p9[1] + 0.3, p9[2]].data,
        @SVector[p10[1] - 3.1, p10[2] - 1.4].data,
        @SVector[p11[1] - 2.3, p11[2] - 0.75].data];
    text=[L"%$s" for s in 1:11], textsize=55)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] - 2.0, p12[2])]; text=L"12", textsize=55)
save("writeups/figures/collinear_point_to_be_added_outside_b.pdf", fig)

DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 12)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1], aspect=1)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), plot_ghost_edges=true, DG=DG, recompute_centroid=false)
xlims!(ax, -10.0, 12.0)
ylims!(ax, -5.0, 12.0)
hidedecorations!(ax)
text!(ax,
    [p1,
        @SVector[p2[1], p2[2] + 0.1].data,
        @SVector[p3[1] - 1.5, p3[2]].data,
        @SVector[p4[1] - 0.4, p4[2] - 1.7].data,
        @SVector[p5[1], p5[2]].data,
        @SVector[p6[1] - 0.65, p6[2]].data,
        @SVector[p7[1] - 0.6, p7[2] - 1.9].data,
        @SVector[p8[1] + 0.25, p8[2] + 0.1].data,
        @SVector[p9[1] + 0.3, p9[2]].data,
        @SVector[p10[1] - 3.1, p10[2] - 1.4].data,
        @SVector[p11[1] - 2.3, p11[2] - 0.75].data];
    text=[L"%$s" for s in 1:11], textsize=55)
scatter!(ax, p12, color=:blue, markersize=17)
text!(ax, [(p12[1] - 2.0, p12[2])]; text=L"12", textsize=55)
lines!(ax, [p2, p12, (NaN, NaN),
        p1, p12, (NaN, NaN),
        p5, p12, (NaN, NaN),
        p4, p12, (NaN, NaN),
        p10, p12, (NaN, NaN)],
    color=:blue, linewidth=4)
save("writeups/figures/collinear_point_to_be_added_outside_c.pdf", fig)
fig

p13 = (9.0,0.0)
push!(pts, p13)
DT.add_point_bowyer!(T, adj, adj2v, DG, pts, 13)

r = 13
pt_idx = graph(DG).V
m = ceil(Int64, length(pt_idx)^(1 / 3))
initial_search_point = DT.select_initial_point(pts, r; pt_idx, m)
V = jump_and_march(r, adj, adj2v, DG, pts; k=initial_search_point)
i, j, k = indices(V)
DT.delete_triangle!(i, j, k, T, adj, adj2v, DG; protect_boundary=true, update_ghost_edges=false)

i, j, k = k, i, j
ℓ = DT.get_edge(adj, j, i) # (j, i, ℓ) is the triangle on the other side of the edge (i, j) from r 




## smaller collinear example

begin
    p1 = (-2.0, 4.0)
    p2 = (-6.0, 2.0)
    p3 = (-6.0, 0.0)
    p4 = (-2.0, 0.0)
    p5 = (-2.0, 2.0)
    pts = [p1, p2, p3, p4, p5][[4, 5, 1, 2, 3]]
    IntegerType, EdgeType, TriangleType,
    EdgesType, TrianglesType,
    try_last_inserted_point = Int64, NTuple{2,Int64},
    NTuple{3,Int64}, Set{NTuple{2,Int64}},
    Set{NTuple{3,Int64}}, true
    I, E, V, Es, Ts = IntegerType, EdgeType,
    TriangleType, EdgesType, TrianglesType
    randomise, trim = false, true
    pt_order = randomise ? shuffle(eachindex(pts)) : eachindex(pts)
    T = Ts()
    adj = DT.Adjacent{I,E}()
    adj2v = DT.Adjacent2Vertex{I,Es,E}()
    DG = DT.DelaunayGraph{I}()
    initial_triangle = DT.construct_positively_oriented_triangle(TriangleType,
        I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]),
        pts
    )
    if DT.isoriented(initial_triangle, pts) == 0 # We will need to randomise the insertion order here - we've started with a degenerate triangle 
        if !randomise # If we didn't randomise, then shuffle! may not work on eachindex(pts)
            pt_order = collect(pt_order)
        end
        while DT.isoriented(initial_triangle, pts) == 0 # We cannot start with a degenerate triangle
            shuffle!(pt_order)
            initial_triangle = DT.construct_positively_oriented_triangle(TriangleType,
                I(pt_order[begin]), I(pt_order[begin+1]), I(pt_order[begin+2]),
                pts
            )
        end
    end
    u, v, w = indices(initial_triangle)
    DT.add_triangle!(u, v, w, T, adj, adj2v, DG; update_ghost_edges=true)
    DT.compute_centroid!(@view pts[pt_order[begin:(begin+2)]])

    num_points = 1
    new_point = pt_order[begin+3]
    last_inserted_point_number = num_points + 3 - 1 # + 3 for the first three points already inserted
    last_inserted_point = pt_order[last_inserted_point_number]
    pt_idx = @view pt_order[begin:last_inserted_point_number] # + 3 for the first three points already inserted
    m = ceil(Int64, length(pt_idx)^(1 / 3))
    initial_search_point = I(DT.select_initial_point(pts, new_point; pt_idx, m, try_points=try_last_inserted_point ? last_inserted_point : ()))
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, new_point;
        pt_idx, m, initial_search_point)
    DT.update_centroid_after_new_point!(pts, new_point)

    num_points = 2
    new_point = pt_order[begin+4]
    last_inserted_point_number = num_points + 3 - 1 # + 3 for the first three points already inserted
    last_inserted_point = pt_order[last_inserted_point_number]
    pt_idx = @view pt_order[begin:last_inserted_point_number] # + 3 for the first three points already inserted
    m = ceil(Int64, length(pt_idx)^(1 / 3))
    initial_search_point = I(DT.select_initial_point(pts, new_point; pt_idx, m, try_points=try_last_inserted_point ? last_inserted_point : ()))
    DT.add_point_bowyer!(T, adj, adj2v, DG, pts, new_point;
        pt_idx, m, initial_search_point)
    DT.update_centroid_after_new_point!(pts, new_point)

    trim && DT.remove_ghost_triangles!(T, adj, adj2v, DG)

    fig = Figure(fontsize=38)
    ax = Axis(fig[1, 1])
    triplot!(ax, T, pts; strokewidth=2, color=(:white, 0.0), markersize=4)
    fig
end



T, adj, adj2v, DG = DT.triangulate_bowyer(pts; trim=true, randomise=false)

