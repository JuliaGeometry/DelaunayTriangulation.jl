using DelaunayTriangulation, Test, Accessors, StatsBase, Random
const DT = DelaunayTriangulation

@testset "separate!" begin
    v = rand(1000)
    p = x -> x < 1 / 2
    i = DT.separate!(v, p)
    @test all(p, v[1:i])
    @test all(!p, v[i+1:end])
    @test i == count(p, v)

    v = tuple.(randn(100000), randn(100000))
    vc = copy(v)
    p = q -> getx(q) < 0
    pts = DT.PointsWrapper(v)
    i = DT.separate!(pts, p)
    @test all(p, pts[1:i])
    @test all(!p, pts[i+1:end])
    @test i == count(p, pts)
    @test sort(v) == sort(vc) # make sure no elements were lost

    v = tuple.(randn(100000), randn(100000))
    vc = copy(v)
    indices = collect(eachindex(v))
    p = q -> getx(q) < 0
    pts = DT.PointsWrapper(v)
    i = DT.separate!(pts, p; carry=indices)
    @test vc[indices] == v
end

@testset "triseparate!" begin
    v = rand(-100:100, 250000)
    p = sign
    i, j = DT.triseparate!(v, p)
    @test all(x -> p(x) < 0, v[1:i-1])
    @test all(x -> p(x) == 0, v[i:j-1])
    @test all(x -> p(x) > 0, v[j:end])
    @test i == count(x -> p(x) < 0, v) + 1
    @test j == i + count(x -> p(x) == 0, v)

    v = tuple.(rand(-100:100, 250000), rand(-100:100, 250000))
    vc = copy(v)
    p = q -> sign(getx(q) - gety(q))
    pts = DT.PointsWrapper(v)
    carry = collect(eachindex(v))
    i, j = DT.triseparate!(pts, p; carry=carry)
    @test all(x -> p(x) < 0, pts[1:i-1])
    @test all(x -> p(x) == 0, pts[i:j-1])
    @test all(x -> p(x) > 0, pts[j:end])
    @test i == count(x -> p(x) < 0, pts) + 1
    @test j == i + count(x -> p(x) == 0, pts)
    @test vc[carry] == v
    @test sort(v) == sort(vc)
    vcs = sort(vc, by=p)
    @test p.(vcs) == p.(v)
end

@testset "get_median!" begin
    _op = (p, q) -> p < q ? -1 : p > q ? 1 : 0

    v = rand(1001)
    vc = copy(v)
    m = DT.get_median!(v, _op)
    @test m == median(vc)

    v = rand(1000)
    vc = copy(v)
    m = DT.get_median!(v, _op)
    DT.sort!(vc)
    @test m == vc[(length(vc)+1)÷2]

    for j in 1:1271
        v = randn(j)
        vc = copy(v)
        carry = collect(eachindex(v))
        m = DT.get_median!(v, _op; carry)
        vcs = DT.sort(vc)
        @test m == vcs[(length(vcs)+1)÷2]
        @test vc[carry] == v
    end

    for j in 1:1200
        v = tuple.(rand(j), rand(j))
        vc = copy(v)
        carry = collect(eachindex(v))
        m = DT.get_median!(v, (p, q) -> _op(getx(p), getx(q)); carry)
        vcs = DT.sort(vc, lt=(p, q) -> <(getx(p), getx(q)))
        @test m == vcs[(length(vcs)+1)÷2]
        @test vc[carry] == v
        vcs = DT.sort(getx.(vc))
        @test getx(m) == vcs[(length(vcs)+1)÷2]

        v = tuple.(rand(j), rand(j))
        vc = copy(v)
        carry = collect(eachindex(v))
        m = DT.get_median!(v, (p, q) -> _op(gety(p), gety(q)); carry)
        vcs = DT.sort(vc, lt=(p, q) -> <(gety(p), gety(q)))
        @test m == vcs[(length(vcs)+1)÷2]
        @test vc[carry] == v
        vcs = DT.sort(gety.(vc))
        @test gety(m) == vcs[(length(vcs)+1)÷2]
    end

    for _ in 1:10000
        v = [(22.9722771127131, 22.7197112165985), (22.9722771127131, 0.7755048895852)]
        m = DT.get_median!(v, (p, q) -> _op(gety(p), gety(q)))
        @test m == (22.9722771127131, 0.7755048895852)
        v = [(22.9722771127131, 22.7197112165985), (22.9722771127131, 0.7755048895852)]
        m = DT.get_median!(v, (p, q) -> _op(getx(p), getx(q)))
        @test m == (22.9722771127131, 22.7197112165985)
    end
end

@testset "X/Y" begin
    @test DT.next(DT.X()) == DT.Y()
    @test DT.next(DT.Y()) == DT.X()
    @test DT.get(DT.X()) === DT.getx
    @test DT.get(DT.Y()) === DT.gety
end

@testset "Left/NotLeft" begin
    @test DT.op(DT.Left()) === <
    @test DT.op(DT.NotLeft()) === ≥
    @test DT.reverse(DT.Left()) === DT.NotLeft()
    @test DT.reverse(DT.NotLeft()) === DT.Left()
end

@testset "HilbertOp" begin
    hop = DT.HilbertOp(DT.Left(), DT.X(), 1 / 2)
    @test hop((1 / 4, 0.0))
    @test hop((1 / 4, 2.0))
    @test !hop((3 / 4, 0.0))
    @test !hop((3 / 4, 2.0))
    @test !hop((1 / 2, 0.0))
    hop = DT.HilbertOp(DT.Left(), DT.Y(), 1 / 2)
    @test hop((0.0, 1 / 4))
    @test hop((2.0, 1 / 4))
    @test !hop((0.0, 3 / 4))
    @test !hop((2.0, 3 / 4))
    @test !hop((0.0, 1 / 2))
    hop = DT.HilbertOp(DT.NotLeft(), DT.X(), 1 / 2)
    @test !hop((1 / 4, 0.0))
    @test !hop((1 / 4, 2.0))
    @test hop((3 / 4, 0.0))
    @test hop((3 / 4, 2.0))
    @test hop((1 / 2, 0.0))
    hop = DT.HilbertOp(DT.NotLeft(), DT.Y(), 1 / 2)
    @test !hop((0.0, 1 / 4))
    @test !hop((2.0, 1 / 4))
    @test hop((0.0, 3 / 4))
    @test hop((2.0, 3 / 4))
    @test hop((0.0, 1 / 2))
end

@testset "get_split" begin
    bbox = (1.0, 2.3, 5.3, 10.3)
    @test DT.get_middle_split(bbox) == (DT.midpoint(bbox[1:2]...), DT.midpoint(bbox[3:4]...))

    points = rand(2, 2500)
    sorter = DT.HilbertSorter(points; splitter=DT.Middle())
    @test DT.get_split(sorter) == (DT.midpoint(sorter.bbox[1:2]...), DT.midpoint(sorter.bbox[3:4]...))

    _op = (p, q) -> p < q ? -1 : p > q ? 1 : 0
    sorter = DT.HilbertSorter(points; splitter=DT.Median())
    split = DT.get_split(sorter)
    @test split[1] == DT.get_median!(copy(DT.PointsWrapper(points)), (p, q) -> _op(getx(p), getx(q)))[1]
    @test split[2] == DT.get_median!(copy(DT.PointsWrapper(points)), (p, q) -> _op(gety(p), gety(q)))[2]
    sorter = @set sorter.directions.coord = DT.Y()
    split = DT.get_split(sorter)
    @test split[1] == DT.get_median!(copy(DT.PointsWrapper(points)), (p, q) -> _op(gety(p), gety(q)))[2]
    @test split[2] == DT.get_median!(copy(DT.PointsWrapper(points)), (p, q) -> _op(getx(p), getx(q)))[1]


    rng = Xoshiro(1)
    rng2 = copy(rng)
    sorter = DT.HilbertSorter(points; splitter=DT.Median(), rng)
    @test rng == rng2
    split = DT.get_split(sorter)
    @test rng != rng2
end

@testset "span/check_bbox" begin
    points = rand(2, 25)
    sorter = DT.HilbertSorter(
        DT.PointsWrapper(points),
        Int[],
        (-1.0, 1.5, 2.3, 5.7),
        1,
        25,
        DT.HilbertDirections(DT.X(), DT.Left(), DT.Left()),
        DT.Middle(),
        5,
        Random.default_rng())
    @test DT.span(sorter) == 25
    @test !DT.check_bbox(sorter)
    sorter = @set sorter.bbox = (eps(), 2eps(), eps(), 2eps())
    @test DT.check_bbox(sorter)
    @test DT.stop_sorter(sorter)

    sorter = DT.HilbertSorter(
        DT.PointsWrapper(points),
        Int[],
        (-1.0, 1.5, 2.3, 5.7),
        1,
        25,
        DT.HilbertDirections(DT.X(), DT.Left(), DT.Left()),
        DT.Middle(),
        25,
        Random.default_rng()
    )
    @test DT.stop_sorter(sorter)
end

@testset "quadrant_sort!" begin
    points = randn(2, 250000)
    sorter = DT.HilbertSorter(
        DT.PointsWrapper(points),
        collect(axes(points, 2)),
        (-2.0, 2.0, -2.0, 2.0),
        1, 250000,
        DT.HilbertDirections(DT.X(), DT.Left(), DT.Left()),
        DT.Middle(),
        5,
        Random.default_rng()
    )
    i₁, i₂, i₃, i₄, i₅ = DT.quadrant_sort!(sorter, 0.0, 0.0)
    @test all(p -> getx(p) < 0 && gety(p) < 0, sorter.points[i₁:i₂-1])
    @test all(p -> getx(p) < 0 && gety(p) ≥ 0, sorter.points[i₂:i₃-1])
    @test all(p -> getx(p) ≥ 0 && gety(p) ≥ 0, sorter.points[i₃:i₄-1])
    @test all(p -> getx(p) ≥ 0 && gety(p) < 0, sorter.points[i₄:i₅])

    points = randn(2, 250000)
    sorter = DT.HilbertSorter(
        DT.PointsWrapper(points),
        collect(axes(points, 2)),
        (-2.0, 2.0, -2.0, 2.0),
        581, 672,
        DT.HilbertDirections(DT.X(), DT.Left(), DT.Left()),
        DT.Middle(),
        5,
        Random.default_rng()
    )
    cpoints = copy(points)
    i₁, i₂, i₃, i₄, i₅ = DT.quadrant_sort!(sorter, 0.0, 0.0)
    @test all(p -> getx(p) < 0 && gety(p) < 0, sorter.points[i₁:i₂-1])
    @test all(p -> getx(p) < 0 && gety(p) ≥ 0, sorter.points[i₂:i₃-1])
    @test all(p -> getx(p) ≥ 0 && gety(p) ≥ 0, sorter.points[i₃:i₄-1])
    @test all(p -> getx(p) ≥ 0 && gety(p) < 0, sorter.points[i₄:i₅])
    @test points[:, 1:580] == cpoints[:, 1:580]
    @test points[:, 673:250000] == cpoints[:, 673:250000]
end

@testset "HilbertSorter" begin
    points = rand(2, 250000)
    perm = collect(axes(points, 2))
    bbox = (extrema(points[1, :])..., extrema(points[2, :])...)
    first = 1
    last = 250000
    directions = DT.HilbertDirections(DT.X(), DT.Left(), DT.Left())
    splitter = DT.Middle()
    capacity = 8
    sorter = DT.HilbertSorter(DT.PointsWrapper(points), perm, bbox, first, last, directions, splitter, capacity, Random.default_rng())
    sorter2 = DT.HilbertSorter(points; capacity, splitter)
    @test sorter.points == sorter2.points
    @test !(sorter.points === sorter2.points)
    @test !(sorter.points === points)
end

@testset "hilbert_sort" begin
    E = (12.8728090546081, 1.057828011921)
    F = (2.9914997728536, 4.9162440171775)
    G = (4.9677616292045, 10.8920834399529)
    H = (2.7091766505178, 21.1498235514886)
    I = (11.0847626131478, 18.9853462802471)
    J = (17.0263810513976, 26.8690055903223)
    K = (23.2289344966548, 24.9012989801027)
    L = (22.9722771127131, 22.7197112165985)
    M = (25.0683124149035, 22.9763686005402)
    N = (25.1110886455604, 16.6027102326552)
    O = (24.9399837229326, 15.1055421596621)
    P = (31.014208476219, 9.3307510209744)
    Q = (25.0683124149035, 6.8497296428715)
    R = (28.6938091798674, 0.903833581556)
    S = (22.9722771127131, 0.7755048895852)
    T = (30.6719986309634, 31.1038524253599)
    points = [E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T]
    cpoints = copy(points)
    spoints, perm, order = hilbert_sort(points; capacity=1)

    @test cpoints == points # make sure there was no mutation of the original points 
    @test !(spoints === points)
    @test points[perm] == spoints
    ids = hilbert_ids(spoints, order)
    @test issorted(ids)
    @test ids == [2, 5, 15, 20, 31, 35, 37, 40, 43, 46, 47, 50, 52, 60, 62, 64]

    @test_throws ArgumentError hilbert_sort(points; capacity=1, splitter=DT.Median())

    for i in 1:10
        for j in [1:100; rand(500:1000, 5)...]
            @show i, j
            r = rand() < 1 / 2 ? rand : randn
            pts = r(2, j)
            @test issorted(hilbert_ids(hilbert_sort(pts; capacity=1)[[1, 3]]...))
        end
    end

    for _ in 1:15
        points = randn(2, 25000)
        bboxes, bbox_lines, splits, fas, _order, (_to_square, _to_point) = _hilbert_sort(points; capacity=4)
        spoints, perm, order = hilbert_sort(points; capacity=4)
        @test order == _order
        @test maximum(length.(values(_to_point))) ≤ 4
    end
end

@testset "assign_rounds" begin
    indices = 1:25000
    rounds = DT.assign_rounds(indices)
    distr = check_distr(indices)
    for j in eachindex(distr)
        dist = distr[end-j+1]
        prob = (1 / 2)^j
        @test dist ≈ prob rtol = 1e-2 atol = 1e-3
    end

    indices = 1:25000
    rng = Xoshiro(2)
    rounds = DT.assign_rounds(indices; rng)
    rounds2 = DT.assign_rounds(indices; rng=Xoshiro(2))
    @test rounds == rounds2
    rounds3 = DT.assign_rounds(indices)
    @test rounds != rounds3
end

@testset "hilbert_sort_rounds" begin
    points = rand(2, 2500)
    cpoints = copy(points)
    rounds = hilbert_sort_rounds(points; rng=Xoshiro(1))
    rounds2 = hilbert_sort_rounds(points; rng=Xoshiro(1))
    @test rounds == rounds2
    @test points == cpoints

    rounds = hilbert_sort_rounds(points; rng=Xoshiro(1))
    _rounds = DT.assign_rounds(1:2500; rng=Xoshiro(1))
    spoints = [hilbert_sort(points[:, r]; capacity=2)[1] for r in _rounds]
    _spoints = [points[:, r] for r in rounds]
    @test spoints == _spoints
end