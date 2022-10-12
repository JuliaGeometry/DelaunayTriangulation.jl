include("functions.jl")

pts, T, DG, adj, adj2v = example_triangulation()
hg = DT.HistoryGraph{NTuple{3,Int64}}()
p1, p2, p3, p4, p5, p6 = pts
p7 = @SVector[2.0, 1.0]
DT.add_point!(pts, p7)
p8 = @SVector[5.0, 1.0]
DT.add_point!(pts, p8)
p9 = @SVector[-1.0, 1.0]
DT.add_point!(pts, p9)
p10 = @SVector[4.0, 1.0]
DT.add_point!(pts, p10)
DT.add_triangle!(5, 2, 8, T, adj, adj2v, DG)
DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 5.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0))
t1 = text!(ax, [(0.0, 1.05), (3.1, -1.2), (2.2, -0.15),
        (-1.0, 2.02), (4.0, 2.02), (-2.4, -1.15), (2.0, 1.0), (5.0, 0.6),
        (-1.0, 0.6), (4.0, 1.0)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6", L"7", L"8", L"9", L"10"], textsize=55)
hidedecorations!(ax)
save("writeups/figures/split_triangulation_a.pdf", fig)

ℓ = lines!(ax, [p1, p3, p7, p1, @SVector[NaN, NaN],
        p3, p7, p5, p3, @SVector[NaN, NaN],
        p5, p7, p1, p5], color=:blue, linewidth=4)
save("writeups/figures/split_triangulation_b.pdf", fig)

ℓ.color = :black
ℓ.linewidth = 2
ℓ2 = lines!(ax, [p4, p9, p1, p4, @SVector[NaN, NaN],
        p6, p9, p1, p6, @SVector[NaN, NaN],
        p6, p9, p4, p6], color=:blue, linewidth=4)
save("writeups/figures/split_triangulation_c.pdf", fig)

ℓ2.color = :black
ℓ2.linewidth = 2
ℓ3 = lines!(ax, [p2, p5, p10, p2, @SVector[NaN, NaN],
        p2, p10, p8, p2, @SVector[NaN, NaN],
        p5, p10, p8, p5], color=:blue, linewidth=4)
save("writeups/figures/split_triangulation_d.pdf", fig)

## Resetup
pts, T, DG, adj, adj2v = example_triangulation()
p1, p2, p3, p4, p5, p6 = pts
p7 = @SVector[2.0, 1.0]
DT.add_point!(pts, p7)
p8 = @SVector[5.0, 1.0]
DT.add_point!(pts, p8)
p9 = @SVector[-1.0, 1.0]
DT.add_point!(pts, p9)
p10 = @SVector[4.0, 1.0]
DT.add_point!(pts, p10)
DT.add_triangle!(5, 2, 8, T, adj, adj2v, DG)
DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)

## Interior splitting 
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (4, 1, 5),
    (5, 2, 8),
    (6, 3, 1),
    (4, 6, 1),
    (1, 3, 7), (3, 5, 7), (5, 1, 7),
    (6, 2, 3),
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
        (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
        (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
        (5, 1) => 7, (1, 7) => 5, (7, 5) => 1,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (4, 5) => DT.BoundaryIndex,
        (5, 8) => DT.BoundaryIndex,
        (8, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex
    )
)
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(5, 4), (4, 6), (6, 3), (3, 7), (7, 5)]),
    2 => Set{NTuple{2,Int64}}([(8, 5), (5, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(5, 7), (7, 1), (1, 6), (6, 2), (2, 5)]),
    4 => Set{NTuple{2,Int64}}([(6, 1), (1, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 1), (1, 7), (7, 3), (3, 2), (2, 8)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
    7 => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 1)]),
    8 => Set{NTuple{2,Int64}}([(5, 2)])
)
true_DG = UndirectedGraph([
    0 0 1 1 1 1 1 0
    0 0 1 0 1 1 0 1
    1 1 0 0 1 1 1 0
    1 0 0 0 1 1 0 0
    1 1 1 1 0 0 1 1
    1 1 1 1 0 0 0 0
    1 0 1 0 1 0 0 0
    0 1 0 0 1 0 0 0
])
true_HG = DirectedGraph{NTuple{3,Int64}}()
[add!(true_HG, (i, j, k)) for (i, j, k) in ((1, 3, 5), (1, 3, 7), (3, 5, 7), (5, 1, 7))]
[add!(true_HG, (1, 3, 5), (i, j, k)) for (i, j, k) in ((1, 3, 7), (3, 5, 7), (5, 1, 7))]
DT.add_triangle!(hg, (1, 3, 5))
DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG, hg)
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG
@test graph(hg) == true_HG

## Splitting a triangle with one boundary edge 
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (4, 1, 5),
    (5, 2, 8),
    (6, 3, 1),
    (4, 6, 9), (6, 1, 9), (1, 4, 9),
    (1, 3, 7), (3, 5, 7), (5, 1, 7),
    (6, 2, 3),
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
        (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
        (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
        (5, 1) => 7, (1, 7) => 5, (7, 5) => 1,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (4, 6) => 9, (6, 9) => 4, (9, 4) => 6,
        (6, 1) => 9, (1, 9) => 6, (9, 6) => 1,
        (9, 1) => 4, (1, 4) => 9, (4, 9) => 1,
        (4, 5) => DT.BoundaryIndex,
        (5, 8) => DT.BoundaryIndex,
        (8, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex
    )
)
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(5, 4), (4, 9), (9, 6), (6, 3), (3, 7), (7, 5)]),
    2 => Set{NTuple{2,Int64}}([(8, 5), (5, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(5, 7), (7, 1), (1, 6), (6, 2), (2, 5)]),
    4 => Set{NTuple{2,Int64}}([(6, 9), (9, 1), (1, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 1), (1, 7), (7, 3), (3, 2), (2, 8)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 9), (9, 4)]),
    7 => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 1)]),
    8 => Set{NTuple{2,Int64}}([(5, 2)]),
    9 => Set{NTuple{2,Int64}}([(6, 1), (1, 4), (4, 6)])
)
true_DG = UndirectedGraph([
    0 0 1 1 1 1 1 0 1
    0 0 1 0 1 1 0 1 0
    1 1 0 0 1 1 1 0 0
    1 0 0 0 1 1 0 0 1
    1 1 1 1 0 0 1 1 0
    1 1 1 1 0 0 0 0 1
    1 0 1 0 1 0 0 0 0
    0 1 0 0 1 0 0 0 0
    1 0 0 1 0 1 0 0 0
])
[add!(true_HG, (i, j, k)) for (i, j, k) in ((4, 6, 1), (4, 6, 9), (6, 1, 9), (1, 4, 9))]
[add!(true_HG, (4, 6, 1), (i, j, k)) for (i, j, k) in ((4, 6, 9), (6, 1, 9), (1, 4, 9))]
DT.add_triangle!(hg, (4, 6, 1))
DT.split_triangle!(4, 6, 1, 9, T, adj, adj2v, DG, hg)
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG
@test graph(hg) == true_HG

## Splitting two boundary edges 
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (4, 1, 5),
    (5, 2, 10), (2, 8, 10), (8, 5, 10),
    (6, 3, 1),
    (4, 6, 9), (6, 1, 9), (1, 4, 9),
    (1, 3, 7), (3, 5, 7), (5, 1, 7),
    (6, 2, 3),
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
        (5, 2) => 10, (2, 10) => 5, (10, 5) => 2,
        (2, 8) => 10, (8, 10) => 2, (10, 2) => 8,
        (8, 5) => 10, (5, 10) => 8, (10, 8) => 5,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (1, 3) => 7, (3, 7) => 1, (7, 1) => 3,
        (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
        (5, 1) => 7, (1, 7) => 5, (7, 5) => 1,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (4, 6) => 9, (6, 9) => 4, (9, 4) => 6,
        (6, 1) => 9, (1, 9) => 6, (9, 6) => 1,
        (9, 1) => 4, (1, 4) => 9, (4, 9) => 1,
        (4, 5) => DT.BoundaryIndex,
        (5, 8) => DT.BoundaryIndex,
        (8, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex
    )
)
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(5, 4), (4, 9), (9, 6), (6, 3), (3, 7), (7, 5)]),
    2 => Set{NTuple{2,Int64}}([(8, 10), (10, 5), (5, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(5, 7), (7, 1), (1, 6), (6, 2), (2, 5)]),
    4 => Set{NTuple{2,Int64}}([(6, 9), (9, 1), (1, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 1), (1, 7), (7, 3), (3, 2), (2, 10), (10, 8)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 9), (9, 4)]),
    7 => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 1)]),
    8 => Set{NTuple{2,Int64}}([(5, 10), (10, 2)]),
    9 => Set{NTuple{2,Int64}}([(6, 1), (1, 4), (4, 6)]),
    10 => Set{NTuple{2,Int64}}([(5, 2), (2, 8), (8, 5)])
)
true_DG = UndirectedGraph([
    0 0 1 1 1 1 1 0 1 0
    0 0 1 0 1 1 0 1 0 1
    1 1 0 0 1 1 1 0 0 0
    1 0 0 0 1 1 0 0 1 0
    1 1 1 1 0 0 1 1 0 1
    1 1 1 1 0 0 0 0 1 0
    1 0 1 0 1 0 0 0 0 0
    0 1 0 0 1 0 0 0 0 1
    1 0 0 1 0 1 0 0 0 0
    0 1 0 0 1 0 0 1 0 0
])
[add!(true_HG, (i, j, k)) for (i, j, k) in ((5, 2, 8), (5, 2, 10), (2, 8, 10), (8, 5, 10))]
[add!(true_HG, (5, 2, 8), (i, j, k)) for (i, j, k) in ((5, 2, 10), (2, 8, 10), (8, 5, 10))]
DT.add_triangle!(hg, (5, 2, 8))
DT.split_triangle!(5, 2, 8, 10, T, adj, adj2v, DG, hg)
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG
@test graph(hg) == true_HG