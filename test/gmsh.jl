x = [2.0, -4.0, -5.0, 5.0, 10.0]
y = [3.0, 3.0, -5.0, -5.0, 2.0]
const GMSH_PATH = "./gmsh-4.9.4-Windows64/gmsh.exe"
for r in 0.05:10.0
    (T, adj, adj2v, DG, nodes), BN = DT.generate_mesh(x, y, 0.2; gmsh_path=GMSH_PATH)

    function DT._get_point(pts::AbstractMatrix, i)
        return @view pts[:, i]
    end
    function DT._eachindex(pts::AbstractMatrix)
        return axes(pts, 2)
    end

    _T, _adj, _adj2v, _DG = DT.triangulate_bowyer(nodes; trim=true)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
end

using LinearAlgebra

pts1 = [LinRange(0, 3, 500)'; repeat([0.0], 500)']
pts2 = [repeat([3.0], 500)'; LinRange(0, 2, 500)']
pts3 = [LinRange(3, 0, 500)'; repeat([2.0], 500)']
pts4 = [repeat([0.0], 500)'; LinRange(2, 0, 500)']
x1, y1 = pts1[1, :], pts1[2, :]
x2, y2 = pts2[1, :], pts2[2, :]
x3, y3 = pts3[1, :], pts3[2, :]
x4, y4 = pts4[1, :], pts4[2, :]
r = 0.5
(T, adj, adj2v, DG, nodes), BN = DT.generate_mesh([x1, x2, x3, x4], [y1, y2, y3, y4], r; gmsh_path=GMSH_PATH)
_T, _adj, _adj2v, _DG = DT.triangulate_bowyer(nodes; trim=true)
@test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, _T, _adj, _adj2v, _DG)
@test norm(nodes[2, BN[1]]) ≈ 0.0 atol = 1e-13
@test norm(nodes[1, BN[2]] .- 3.0) ≈ 0.0 atol = 1e-13
@test norm(nodes[2, BN[3]] .- 2.0) ≈ 0.0 atol = 1e-13
@test norm(nodes[1, BN[4]]) ≈ 0.0 atol = 1e-13
@test length(BN) == 4
@test length(BN[1]) > 1
@test nodes ≈ [0.0 0.0
    3.0 0.0
    3.0 2.0
    0.0 2.0
    0.5 0.0
    1.0 0.0
    1.5 0.0
    2.0 0.0
    2.5 0.0
    3.0 0.5
    3.0 1.0
    3.0 1.5
    2.5 2.0
    2.0 2.0
    1.5 2.0
    1.0 2.0
    0.5 2.0
    0.0 1.5
    0.0 1.0
    0.0 0.5
    1.25 1.56699
    1.74591 0.430653
    0.752119 0.423386
    2.20243 1.57738
    2.52036 0.795066
    0.395353 1.28162
    1.24967 0.431015
    1.49926 0.844114
    1.99049 0.829466
    1.73487 1.26758
    0.997105 0.807507
    2.23179 0.403535
    2.58725 1.25461
    0.720683 1.61378
    0.447379 0.785961
    1.25628 1.14716
    1.75 1.6941
    2.63397 0.366025
    0.366025 1.63397
    0.366025 0.366025
    2.63397 1.63397
    0.820982 1.20019
    2.20456 1.16649]' atol = 1e-4

a = 0.0
b = 2.0
c = 0.0
d = 2.5
ref = 1.2
tri, BN = generate_mesh(a, b, c, d, ref)
T, adj, adj2v, DG, pts = tri

@test DT.validate_triangulation(T, adj, adj2v, DG, pts)
@test pts ≈ [0.0 2.0 2.0 0.0 1.0 2.0 2.0 1.0 0.0 0.0 1.25052 0.553123 1.36904 1.36904 0.695448 0.695448
    0.0 0.0 2.5 2.5 0.0 0.833333 1.66667 2.5 1.66667 0.833333 1.25 1.25 0.573816 1.92618 1.84881 0.651192] rtol = 1e-3
@test BN[1] == [1, 5, 2, 6, 7, 3, 8, 4, 9, 10]

ref = 1.1
tri, BN = generate_mesh(a, b, c, d, ref; single_boundary=false)
T, adj, adj2v, DG, pts = tri
@test DT.validate_triangulation(T, adj, adj2v, DG, pts)
@test pts ≈ [0.0 2.0 2.0 0.0 1.0 2.0 2.0 1.0 0.0 0.0 1.25052 0.553123 1.36904 1.36904 0.695448 0.695448
    0.0 0.0 2.5 2.5 0.0 0.833333 1.66667 2.5 1.66667 0.833333 1.25 1.25 0.573816 1.92618 1.84881 0.651192] rtol = 1e-3
@test BN[1] == [1, 5, 2]
@test BN[2] == [2, 6, 7, 3]
@test BN[3] == [3, 8, 4]
@test BN[4] == [4, 9, 10, 1]