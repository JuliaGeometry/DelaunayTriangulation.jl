include("functions.jl")

pts, T, DG, adj, adj2v = example_triangulation()
p1, p2, p3, p4, p5, p6 = pts
p7 = @SVector[2.0, 1.0]
DT.add_point!(pts, p7)
p8 = @SVector[5.0, 1.0]
DT.add_point!(pts, p8)
DT.add_triangle!(1, 3, 7, T, adj, adj2v, DG)
DT.add_triangle!(5, 2, 8, T, adj, adj2v, DG)
DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 5.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0))
t1 = text!(ax, [(0.0, 1.05), (3.1, -1.2), (2.2, -0.15),
        (-1.0, 2.02), (4.0, 2.02), (-2.4, -1.15), (2.0, 1.0), (5.0, 0.6)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6", L"7", L"8"], textsize=55)
hidedecorations!(ax)
save("writeups/figures/deletion_triangulation_a.pdf", fig)

alphabet = join('a':'z')
for (r, (i, j, k)) in zip(2:4, ((1, 3, 7), (5, 2, 8), (6, 2, 3)))
    DT.delete_triangle!(T, (i, j, k))
    fig = Figure(fontsize=55)
    ax = Axis(fig[1, 1])
    xlims!(ax, -2.5, 5.5)
    ylims!(ax, -1.5, 2.5)
    triplot!(ax, T, pts; strokewidth=2, color=(:white, 0))
    t1 = text!(ax, [(0.0, 1.05), (3.1, -1.2), (2.2, -0.15),
            (-1.0, 2.02), (4.0, 2.02), (-2.4, -1.15), (2.0, 1.0), (5.0, 0.6)];
        text=[L"1", L"2", L"3", L"4", L"5", L"6", L"7", L"8"], textsize=55)
    hidedecorations!(ax)
    PI, PJ, PK = pts[i], pts[j], pts[k]
    lines!(ax, [PI, PJ, PK, PI], color=:red, linewidth=4, linestyle=:dash)
    save("writeups/figures/deletion_triangulation_$(alphabet[r]).pdf", fig)
end

## Resetup the problem 
pts, T, DG, adj, adj2v = example_triangulation()
p1, p2, p3, p4, p5, p6 = pts
p7 = @SVector[2.0, 1.0]
DT.add_point!(pts, p7)
p8 = @SVector[5.0, 1.0]
DT.add_point!(pts, p8)
DT.add_triangle!(1, 3, 7, T, adj, adj2v, DG)
DT.add_triangle!(5, 2, 8, T, adj, adj2v, DG)
DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)

## Now delete the interior 
i, j, k = 1, 3, 7
DT.delete_triangle!(i, j, k, T, adj, adj2v, DG)
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (4, 1, 5),
    (6, 3, 1),
    (4, 6, 1),
    (5, 1, 3),
    (5, 2, 8),
    (6, 2, 3)
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (5, 1) => 3, (3, 5) => 1,
        (1, 7) => DT.DefaultAdjacentValue,
        (7, 3) => DT.DefaultAdjacentValue,
        (5, 2) => 8, (2, 8) => 5, (8, 5) => 2,
        (2, 3) => 6, (3, 6) => 2, (6, 2) => 3,
        (4, 5) => DT.BoundaryIndex, (5, 8) => DT.BoundaryIndex,
        (8, 2) => DT.BoundaryIndex, (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex
    )
)
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 8), (8, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(3, 5), (6, 3), (5, 4), (4, 6)]),
    2 => Set{NTuple{2,Int64}}([(5, 3), (8, 5), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(2, 5), (5, 1), (1, 6), (6, 2)]),
    4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
    5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2), (2, 8)]),
    6 => Set{NTuple{2,Int64}}([(3, 1), (1, 4), (2, 3)]),
    7 => Set{NTuple{2,Int64}}([]),
    8 => Set{NTuple{2,Int64}}([(5, 2)])
)
true_DG = UndirectedGraph(
    [
        0 0 1 1 1 1 0 0
        0 0 1 0 1 1 0 1
        1 1 0 0 1 1 0 0
        1 0 0 0 1 1 0 0
        1 1 1 1 0 0 0 1
        1 1 1 1 0 0 0 0
        0 0 0 0 0 0 0 0
        0 1 0 0 1 0 0 0
    ]
)
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG

## Now delete the triangle with two boundary edges 
for (i, j, k) in ((5, 2, 8), (2, 8, 5), (8, 5, 2))
    Tc, adjc, adj2vc, DGc = deepcopy(T), deepcopy(adj), deepcopy(adj2v), deepcopy(DG)
    DT.delete_triangle!(i, j, k, Tc, adjc, adj2vc, DGc)
    true_T = Set{NTuple{3,Int64}}([
        (3, 2, 5),
        (4, 1, 5),
        (6, 3, 1),
        (4, 6, 1),
        (5, 1, 3),
        (6, 2, 3)
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (3, 5) => 1,
            (1, 7) => DT.DefaultAdjacentValue,
            (7, 3) => DT.DefaultAdjacentValue,
            (2, 3) => 6, (3, 6) => 2, (6, 2) => 3,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex, (2, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(3, 5), (6, 3), (5, 4), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3), (3, 6)]),
        3 => Set{NTuple{2,Int64}}([(2, 5), (5, 1), (1, 6), (6, 2)]),
        4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(3, 1), (1, 4), (2, 3)]),
        7 => Set{NTuple{2,Int64}}([]),
        8 => Set{NTuple{2,Int64}}([])
    )
    true_DG = UndirectedGraph(
        [
            0 0 1 1 1 1 0 0
            0 0 1 0 1 1 0 0
            1 1 0 0 1 1 0 0
            1 0 0 0 1 1 0 0
            1 1 1 1 0 0 0 0
            1 1 1 1 0 0 0 0
            0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0
        ]
    )
    @test Tc == true_T
    @test adjacent(adjc) == true_adj
    @test adjacent2vertex(adj2vc) == true_adj2v
    @test graph(DGc) == true_DG
end
i, j, k = 5, 2, 8
DT.delete_triangle!(i, j, k, T, adj, adj2v, DG)

## Now delete the triangle with a single boundary edge 
for (i, j, k) in ((6, 2, 3), (2, 3, 6), (3, 6, 2))
    Tc, adjc, adj2vc, DGc = deepcopy(T), deepcopy(adj), deepcopy(adj2v), deepcopy(DG)
    DT.delete_triangle!(i, j, k, Tc, adjc, adj2vc, DGc)
    true_T = Set{NTuple{3,Int64}}([
        (3, 2, 5),
        (4, 1, 5),
        (6, 3, 1),
        (4, 6, 1),
        (5, 1, 3),
    ])
    true_adj = DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (3, 5) => 1,
            (1, 7) => DT.DefaultAdjacentValue,
            (7, 3) => DT.DefaultAdjacentValue,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex,
            (2, 3) => DT.BoundaryIndex, (3, 6) => DT.BoundaryIndex,
        )
    )
    true_adj2v = Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(3, 5), (6, 3), (5, 4), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3)]),
        3 => Set{NTuple{2,Int64}}([(2, 5), (5, 1), (1, 6)]),
        4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(3, 1), (1, 4)]),
        7 => Set{NTuple{2,Int64}}([]),
        8 => Set{NTuple{2,Int64}}([])
    )
    true_DG = UndirectedGraph(
        [
            0 0 1 1 1 1 0 0
            0 0 1 0 1 0 0 0
            1 1 0 0 1 1 0 0
            1 0 0 0 1 1 0 0
            1 1 1 1 0 0 0 0
            1 0 1 1 0 0 0 0
            0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0
        ]
    )
    @test Tc == true_T
    @test adjacent(adjc) == true_adj
    @test adjacent2vertex(adj2vc) == true_adj2v
    @test graph(DGc) == true_DG
end
i, j, k = 6, 2, 3
DT.delete_triangle!(i, j, k, T, adj, adj2v, DG)

## Test that we can add to an empty triangulation correctly 
pts, T, DG, adj, adj2v = example_empty_triangulation()
p1, p2, p3 = pts
true_T = Set{NTuple{3,Int64}}([(1, 2, 3)])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict((1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
        (3, 2) => DT.BoundaryIndex, (2, 1) => DT.BoundaryIndex, (1, 3) => DT.BoundaryIndex))
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(3, 2), (2, 1), (1, 3)]),
    1 => Set{NTuple{2,Int64}}([(2, 3)]),
    2 => Set{NTuple{2,Int64}}([(3, 1)]),
    3 => Set{NTuple{2,Int64}}([(1, 2)])
)
true_DG = UndirectedGraph([1 1 1; 1 1 1; 1 1 1])
DT.add_triangle!(1, 2, 3, T, adj, adj2v, DG)
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG

## Test that we can add delete a single triangle correctly
true_T = Set{NTuple{3,Int64}}([])
true_adj = DefaultDict(DT.DefaultAdjacentValue,Dict())
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}(),
    1 => Set{NTuple{2,Int64}}(),
    2 => Set{NTuple{2,Int64}}(),
    3 => Set{NTuple{2,Int64}}()
)
true_DG = UndirectedGraph([0 0 0; 0 0 0; 0 0 0])
DT.delete_triangle!(1, 2, 3, T, adj, adj2v, DG)
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG