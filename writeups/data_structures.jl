using .DelaunayTriangulation
using CairoMakie
using StaticArrays

function triplot!(ax, tri; kwargs...)
    Tmat = zeros(Int64, num_triangles(tri), 3)
    for (i, T) in enumerate(triangles(tri))
        Tmat[i, :] = [geti(T), getj(T), getk(T)]
    end
    pmat = zeros(num_points(tri), 2)
    for (i, p) in enumerate(points(tri))
        pmat[i, :] = [getx(p), gety(p)]
    end
    poly!(ax, pmat, Tmat; kwargs...)
    scatter!(ax, pmat, color=:red, markersize=11)
end

p1 = @SVector[0.0, 1.0]
p2 = @SVector[3.0, -1.0]
p3 = @SVector[2.0, 0.0]
p4 = @SVector[-1.0, 2.0]
p5 = @SVector[4.0, 2.0]
p6 = @SVector[-2.0, -1.0]
pts = [p1, p2, p3, p4, p5, p6]
tri = triangulate(pts)

fig = Figure(fontsize=28)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 4.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, tri; strokewidth=2, color=(:white, 0))
text!(ax, [(0.0, 1.05), (3.1, -1.05), (2.1, -0.09),
        (-1.0, 2.02), (4.0, 2.02), (-2.2, -1.05)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6"], textsize=28)
hidedecorations!(ax)
save("writeups/figures/example_triangulation.pdf", fig)

p7 = @SVector[2.0, 1.0]
lines!(ax, [p1[1], p3[1], p7[1], p1[1]], [p1[2], p3[2], p7[2], p1[2]], color=:blue, linestyle=:dash, linewidth=3)
scatter!(ax, p7, color=:red, markersize=11)
text!(ax, [(2.0, 1.0)]; text=[L"7"], textsize=28)
save("writeups/figures/addition_triangulation.pdf", fig)

lines!(ax, [p1[1], p5[1], p7[1], p1[1]], [p1[2], p5[2], p7[2], p1[2]], color=:blue, linestyle=:dash, linewidth=3)
lines!(ax, [p3[1], p5[1], p7[1], p3[1]], [p3[2], p5[2], p7[2], p3[2]], color=:blue, linestyle=:dash, linewidth=3)
scatter!(ax, [p1, p5, p7, p3], color=:red, markersize=11)
save("writeups/figures/split_triangulation.pdf", fig)

newtri = deepcopy(tri)
push!(points(newtri), p7)
DelaunayTriangulation.split_triangle!(triangles(newtri),
    pointlocation(newtri), adjacent(newtri),
    adjacent2vertex(newtri), graph(newtri),
    (1, 3, 5), 7)

fig = Figure(fontsize=28)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 4.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, newtri; strokewidth=2, color=(:white, 0))
text!(ax, [(0.0, 1.05), (3.1, -1.05), (2.1, -0.09),
        (-1.0, 2.02), (4.0, 2.02), (-2.2, -1.05)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6"], textsize=28)
hidedecorations!(ax)
lines!(ax, [p1, p3, p7, p1], color=:red, linestyle=:dash, linewidth=3)
scatter!(ax, [p1, p5, p7, p3], color=:red, markersize=11)
text!(ax, [(2.0, 1.0)]; text=[L"7"], textsize=28)
save("writeups/figures/delete_triangulation.pdf", fig)

newtri2 = deepcopy(newtri)
DelaunayTriangulation.delete_triangle!(triangles(newtri2),
    adjacent(newtri2),
    adjacent2vertex(newtri2), graph(newtri2),
    1, 3, 7)
DelaunayTriangulation.delete_triangle!(triangles(newtri2),
    adjacent(newtri2),
    adjacent2vertex(newtri2), graph(newtri2),
    3, 2, 5)

fig = Figure(fontsize=28)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 4.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, newtri2; strokewidth=2, color=(:white, 0))
text!(ax, [(0.0, 1.05), (3.1, -1.05), (2.1, -0.09),
        (-1.0, 2.02), (4.0, 2.02), (-2.2, -1.05)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6"], textsize=28)
hidedecorations!(ax)
text!(ax, [(2.0, 1.0)]; text=[L"7"], textsize=28)
fig
save("writeups/figures/delete_triangulation_2.pdf", fig)


fig = Figure(fontsize=28)
DelaunayTriangulation.flip_edge!(triangles(tri), pointlocation(tri),
    adjacent(tri), adjacent2vertex(tri), graph(tri),
    1, 3, 5, 6)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 4.5)
ylims!(ax, -1.5, 2.5)
text!(ax, [(0.0, 1.05), (3.1, -1.05), (2.1, -0.09),
        (-1.0, 2.02), (4.0, 2.02), (-2.2, -1.05)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6"], textsize=28)
DelaunayTriangulation.delete_triangle!(triangles(tri), (1, 6, 3), (1, 3, 5))
DelaunayTriangulation.add_triangle!(triangles(tri), (6, 5, 1), (6, 3, 5))
lines!(ax, [p1, p3], color=:blue, linestyle=:dash, linewidth=3)
triplot!(ax, tri; strokewidth=2, color=(:white, 0))
hidedecorations!(ax)
save("writeups/figures/flipped_triangulation.pdf", fig)

p1 = @SVector[0.0, 1.0]
p2 = @SVector[3.0, -1.0]
p3 = @SVector[2.0, 0.0]
p4 = @SVector[-1.0, 2.0]
p5 = @SVector[4.0, 2.0]
p6 = @SVector[-2.0, -1.0]
pts = [p1, p2, p3, p4, p5, p6]
tri = triangulate(pts)
newtri = deepcopy(tri)
newtri = deepcopy(tri)
push!(points(newtri), p7)
DelaunayTriangulation.split_triangle!(triangles(newtri),
    pointlocation(newtri), adjacent(newtri),
    adjacent2vertex(newtri), graph(newtri),
    (1, 3, 5), 7)
fig = Figure(fontsize=28)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 4.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, newtri; strokewidth=2, color=(:white, 0))
text!(ax, [(0.0, 1.05), (3.1, -1.05), (2.1, -0.09),
        (-1.0, 2.02), (4.0, 2.02), (-2.2, -1.05)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6"], textsize=28)
hidedecorations!(ax)
text!(ax, [(2.0, 1.0)]; text=[L"7"], textsize=28)
lines!(ax, [p1, p5, p3, p1], color=:blue, linewidth=3)
scatter!(ax, [p1, p5, p3], color=:red, markersize=11)
save("writeups/figures/split_triangulation_2.pdf", fig)

function circle_three_points(a, b, c)
    ax, ay = a
    bx, by = b
    cx, cy = c
    d1x, d1y = [by - ay, ax - bx]
    d2x, d2y = [cy - ay, ax - cx]
    k = d2x * d1y - d2y * d1x
    s1x, s1y = [(ax + bx) / 2, (ay + by) / 2]
    s2x, s2y = [(ax + cx) / 2, (ay + cy) / 2]
    ℓ = d1x * (s2y - s1y) - d1y * (s2x - s1x)
    m = ℓ / k
    centx, centy = [s2x + m * d2x, s2y + m * d2y]
    dx = centx - ax
    dy = centy - ay
    r = sqrt(dx^2 + dy^2)
    return [centx, centy], r
end
function circle_data(a, b, c)
    θ = LinRange(0, 2π, 500)
    (h, k), r = circle_three_points(a, b, c)
    circx = h .+ r * sin.(θ)
    circy = k .+ r * cos.(θ)
    return circx, circy
end


fig = Figure(fontsize=28)
ax = Axis(fig[1, 1])
circ1x, circ1y = circle_data(p1, p3, p7)
circ2x, circ2y = circle_data(p3, p5, p7)
circ3x, circ3y = circle_data(p5, p1, p7)
lines!(ax, circ1x, circ1y, color=(:darkgreen, 0.9), linewidth=3, linestyle=:dash)
lines!(ax, circ2x, circ2y, color=(:darkgreen, 0.9), linewidth=3, linestyle=:dash)
lines!(ax, circ3x, circ3y, color=(:darkgreen, 0.9), linewidth=3, linestyle=:dash)
xlims!(ax, -2.5, 4.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, newtri; strokewidth=2, color=(:white, 0))
text!(ax, [(0.0, 1.05), (3.1, -1.05), (2.1, -0.09),
        (-1.0, 2.02), (4.0, 2.02), (-2.2, -1.05)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6"], textsize=28)
hidedecorations!(ax)
text!(ax, [(2.0, 1.0)]; text=[L"7"], textsize=28)
lines!(ax, [p1, p5, p3, p1], color=:blue, linewidth=3)
scatter!(ax, [p1, p5, p3], color=:red, markersize=11)
save("writeups/figures/split_triangulation_circles.pdf", fig)

p1 = (5.0, 6.0)
p2 = (9.0, 6.0)
p3 = (13.0, 5.0)
p4 = (10.38, 0.0)
p5 = (12.64, -1.69)
p6 = (2.0, -2.0)
p7 = (3.0, 4.0)
p8 = (7.5, 3.53)
p9 = (4.02, 1.85)
p10 = (4.26, 0.0)
pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
tri = triangulate(pts)
fig = Figure(fontsize=28)
ax = Axis(fig[1, 1])
triplot!(ax, tri; strokewidth=2, color=(:white, 0.0))
xlims!(ax, 1.0, 14.0)
ylims!(ax, -3.0, 7.0)
hidedecorations!(ax)
text!(ax,
    [p1, p2, p3, (p4[1] + 0.2, p4[2] - 0.1),
        (p5[1] + 0.1, p5[2] - 0.1),
        (p6[1] - 0.53, p6[2] - 0.2),
        (p7[1] - 0.46, p7[2] - 0.2),
        (p8[1] - 0.15, p8[2] - 0.85),
        (p9[1], p9[2] + 0.1),
        (p10[1], p10[2] - 0.8)];
    text=[L"%$s" for s in 1:10], textsize=43)
p11 = (6.0, 2.5)
scatter!(ax, p11, color=:blue, markersize=17)
push!(pts, p11)
text!(ax, (p11[1] - 0.74, p11[2] - 0.75); text=L"11", textsize=43)
save("writeups/figures/excavated_cavity_triangulation_a.pdf", fig)

bw_in = NTuple{3,Int64}[]
for T in triangles(tri)
    DelaunayTriangulation.isincircle(T, pts, 11) == 1 && push!(bw_in, T)
end

fig = Figure(fontsize=28)
ax = Axis(fig[1, 1])
for T in bw_in
    i, j, k = T
    pti, ptj, ptk = pts[i], pts[j], pts[k]
    circx, circy = circle_data(pti, ptj, ptk)
    lines!(ax, circx, circy, color=(:darkgreen, 0.9), linewidth=3, linestyle=:dash)
end
poly!(ax, pts, [bw_in[i][j] for i in eachindex(bw_in), j in 1:3],
    strokewidth=2, color=(:lightblue, 0.5))
triplot!(ax, tri; strokewidth=2, color=(:white, 0.0))
xlims!(ax, 1.0, 14.0)
ylims!(ax, -3.0, 7.0)
hidedecorations!(ax)
text!(ax,
    [p1, p2, p3, (p4[1] + 0.2, p4[2] - 0.1),
        (p5[1] + 0.1, p5[2] - 0.1),
        (p6[1] - 0.53, p6[2] - 0.2),
        (p7[1] - 0.46, p7[2] - 0.2),
        (p8[1] - 0.15, p8[2] - 0.85),
        (p9[1], p9[2] + 0.1),
        (p10[1], p10[2] - 0.8)];
    text=[L"%$s" for s in 1:10], textsize=43)
scatter!(ax, p11, color=:blue, markersize=17)
text!(ax, (p11[1] - 0.74, p11[2] - 0.75); text=L"11", textsize=43)
save("writeups/figures/excavated_cavity_triangulation_b.pdf", fig)

Tmat = triangles(tri)
[DelaunayTriangulation.delete_triangle!(Tmat, T) for T in bw_in]
DelaunayTriangulation.add_triangle!(Tmat,
    (1, 11, 8),
    (7, 11, 1),
    (7, 9, 11),
    (9, 10, 11),
    (10, 4, 11),
    (4, 8, 11))

fig = Figure(fontsize=28)
ax = Axis(fig[1, 1])
triplot!(ax, tri; strokewidth=2, color=(:white, 0.0))
xlims!(ax, 1.0, 14.0)
ylims!(ax, -3.0, 7.0)
hidedecorations!(ax)
text!(ax,
    [p1, p2, p3, (p4[1] + 0.2, p4[2] - 0.1),
        (p5[1] + 0.1, p5[2] - 0.1),
        (p6[1] - 0.53, p6[2] - 0.2),
        (p7[1] - 0.46, p7[2] - 0.2),
        (p8[1] - 0.15, p8[2] - 0.85),
        (p9[1], p9[2] + 0.1),
        (p10[1], p10[2] - 0.8)];
    text=[L"%$s" for s in 1:10], textsize=43)
scatter!(ax, p11, color=:blue, markersize=17)
text!(ax, (p11[1] - 0.24, p11[2] - 0.93); text=L"11", textsize=43)
fig
lines!(ax, [p11, p1, (NaN, NaN),
        p11, p4, (NaN, NaN),
        p11, p10, (NaN, NaN),
        p11, p9, (NaN, NaN),
        p11, p7, (NaN, NaN),
        p11, p8, (NaN, NaN)],
    color=:blue, linewidth=4)
save("writeups/figures/excavated_cavity_triangulation_c.pdf", fig)

p1 = @SVector[0.0, 1.0]
p2 = @SVector[3.0, -1.0]
p3 = @SVector[2.0, 0.0]
p4 = @SVector[-1.0, 2.0]
p5 = @SVector[4.0, 2.0]
p6 = @SVector[-2.0, -1.0]
pts = [p1, p2, p3, p4, p5, p6]
tri = triangulate(pts)

fig = Figure(fontsize=36)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 4.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, tri; strokewidth=2, color=(:white, 0))
text!(ax, [(0.0, 1.05), (3.1, -1.05), (2.1, -0.09),
        (-1.0, 2.02), (4.0, 2.02), (-2.2, -1.05)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6"], textsize=28)
hidedecorations!(ax)
p7 = @SVector[2.0, 1.0]
lines!(ax, [p1[1], p3[1], p7[1], p1[1]], [p1[2], p3[2], p7[2], p1[2]], color=:blue, linestyle=:dash, linewidth=3)
scatter!(ax, p7, color=:red, markersize=11)
text!(ax, [(2.0, 1.0)]; text=[L"7"], textsize=36)
save("writeups/figures/addition_triangulation.pdf", fig)

p7 = @SVector[5.0, 1.0]
fig = Figure(fontsize=36)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 5.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, tri; strokewidth=2, color=(:white, 0))
text!(ax, [(0.0, 1.05), (3.1, -1.05), (2.1, -0.09),
        (-1.0, 2.02), (4.0, 2.02), (-2.2, -1.05)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6"], textsize=36)
hidedecorations!(ax)
lines!(ax, [p5, p2, p7, p5], color=:blue, linestyle=:dash, linewidth=3)
scatter!(ax, [p5, p2, p7, p5], color=:red, markersize=11)
text!(ax, [(5.0, 0.7)]; text=[L"7"], textsize=36)
save("writeups/figures/triangulation_boundary_addition_a.pdf", fig)

DelaunayTriangulation.delete_triangle!(triangles(tri), (6, 2, 3))
fig = Figure(fontsize=36)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 5.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, tri; strokewidth=2, color=(:white, 0))
text!(ax, [(0.0, 1.05), (3.1, -1.05), (2.1, -0.09),
        (-1.0, 2.02), (4.0, 2.02), (-2.2, -1.05)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6"], textsize=36)
hidedecorations!(ax)
lines!(ax, [p3, p6, p2, p3], color=:blue, linestyle=:dash, linewidth=3)
scatter!(ax, [p3, p6, p2, p3], color=:red, markersize=11)
save("writeups/figures/triangulation_boundary_addition_b.pdf", fig)