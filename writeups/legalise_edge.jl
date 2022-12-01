include("functions.jl")

pts, T, DG, adj, adj2v = example_triangulation()
p1, p2, p3, p4, p5, p6 = pts
p7 = @SVector[2.0, 1.0]
DT.add_point!(pts, p7)
DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)
DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 4.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0))
t1 = text!(ax, [(-0.05, 0.55), (3.1, -1.2), (2.2, -0.15),
        (-1.0, 2.02), (4.0, 2.02), (-2.4, -1.15), (2.0, 1.0)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6", L"7"], textsize=55)
hidedecorations!(ax)
ℓ = lines!(ax, [p1, p3, p5, p1], color=:blue, linewidth=4)
save("writeups/figures/legalise_triangulation_1.pdf", fig)

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

circ1x, circ1y = circle_data(p1, p3, p7)
circ2x, circ2y = circle_data(p3, p5, p7)
circ3x, circ3y = circle_data(p5, p1, p7)
lines!(ax, circ1x, circ1y, color=(:darkgreen, 0.9), linewidth=3, linestyle=:dash)
lines!(ax, circ2x, circ2y, color=(:darkgreen, 0.9), linewidth=3, linestyle=:dash)
lines!(ax, circ3x, circ3y, color=(:darkgreen, 0.9), linewidth=3, linestyle=:dash)
save("writeups/figures/legalise_triangulation_2.pdf", fig)

## Test 
pts, T, DG, adj, adj2v = example_triangulation()
hg = DT.HistoryGraph{NTuple{3,Int64}}()
DT.add_triangle!(hg, T...)
p1, p2, p3, p4, p5, p6 = pts
p7 = @SVector[2.0, 1.0]
DT.add_point!(pts, p7)
DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)
DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG)
@test DT.islegal(1, 3, adj, pts)
@test DT.islegal(3, 5, adj, pts)
i, j, r = 5, 1, 7
e = DT.get_edge(adj, j, i)
DT.legalise_edge!(i, j, r, T, hg, adj, adj2v, DG, pts)
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (1, 3, 7),
    (3, 5, 7),
    (6, 3, 1),
    (4, 6, 1),
    (6, 2, 3),
    (7, 5, 4),
    (7, 4, 1)
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
        (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (7, 5) => 4, (5, 4) => 7, (4, 7) => 5,
        (7, 4) => 1, (4, 1) => 7, (1, 7) => 4,
        (4, 5) => DT.BoundaryIndex,
        (5, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex,
        (1, 5) => DT.DefaultAdjacentValue
    )
)
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(6, 3), (3, 7), (7, 4), (4, 6)]),
    2 => Set{NTuple{2,Int64}}([(5, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(2, 5), (5, 7), (7, 1), (1, 6), (6, 2)]),
    4 => Set{NTuple{2,Int64}}([(6, 1), (1, 7), (7, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 7), (7, 3), (3, 2)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
    7 => Set{NTuple{2,Int64}}([(3, 5), (5, 4), (4, 1), (1, 3)])
)
true_DG = UndirectedGraph([
    0 0 1 1 0 1 1
    0 0 1 0 1 1 0
    1 1 0 0 1 1 1
    1 0 0 0 1 1 1
    0 1 1 1 0 0 1
    1 1 1 1 0 0 0
    1 0 1 1 1 0 0
])
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG
@test all(DT.isoriented(T, pts) == 1 for T in T)
for ((i, j), v) in adjacent(adj)
    if v ≠ DT.DefaultAdjacentValue
        @test DT.islegal(i, j, adj, pts)
    end
end