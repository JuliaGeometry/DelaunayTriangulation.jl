using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using Random
using StableRNGs

@testset "hilbert_index" begin
    for order in 1:6
        n = 2^order
        cells = Dict{UInt64, NTuple{2, Int}}()
        for x in 0:(n - 1), y in 0:(n - 1)
            cells[DT.hilbert_index(x, y, order)] = (x, y)
        end
        # Each cell has its own position, the curve goes from (0, 0) to (n - 1, 0), and consecutive cells are adjacent
        @test sort(collect(keys(cells))) == 0:(n^2 - 1)
        @test cells[0] == (0, 0) && cells[n^2 - 1] == (n - 1, 0)
        @test all(d -> sum(abs.(cells[d] .- cells[d + 1])) == 1, 0:(n^2 - 2))
    end
    @test DT.hilbert_index(typemax(UInt32), 0, DT.HILBERT_ORDER) == typemax(UInt64)
end

@testset "hilbert_cell" begin
    scale = Float64(typemax(UInt32))
    @test DT.hilbert_cell(0.0, 0.0, scale) == 0
    @test DT.hilbert_cell(1.0, 0.0, scale) == typemax(UInt32)
    @test DT.hilbert_cell(2.0, 0.0, scale) == typemax(UInt32)
    @test DT.hilbert_cell(-1.0, 0.0, scale) == 0
    @test DT.hilbert_cell(NaN, 0.0, scale) == 0
    @test DT.hilbert_cell(Inf, 0.0, scale) == 0
    @test DT.hilbert_cell(0.5f0, 0.0, scale) == DT.hilbert_cell(0.5, 0.0, scale)
end

@testset "hilbert_bounds" begin
    @test DT.hilbert_bounds([(1.0, 2.0), (3.0, -1.0), (NaN, 5.0), (0.5, Inf)], 1:4) == (0.5, 3.0, -1.0, 5.0)
    @test DT.hilbert_bounds([(1.0, 2.0), (3.0, -1.0)], [2]) == (3.0, 3.0, -1.0, -1.0)
    @test DT.hilbert_bounds([(NaN, NaN)], 1:1) == (0.0, 0.0, 0.0, 0.0)
    @test DT.hilbert_bounds([1.0 3.0; 2.0 -1.0], 1:2) == (1.0, 3.0, -1.0, 2.0)
end

@testset "brio_order" begin
    for n in (0, 1, 2, 3, 4, 5, 17, 1000)
        points = [(rand(rng), rand(rng)) for rng in (StableRNG(n),) for _ in 1:n]
        order = DT.brio_order(points, collect(1:n), StableRNG(1))
        @test sort(order) == 1:n
        @test order == DT.brio_order(points, collect(1:n), StableRNG(1))
    end
    @test DT.brio_order([(0.0, 0.0)], Int32[], StableRNG(1)) == Int32[]
    @test eltype(DT.brio_order([(0.0, 0.0), (1.0, 1.0)], Int32[2, 1], StableRNG(1))) == Int32

    # Degenerate bounding boxes and coordinates that are not finite
    @test sort(DT.brio_order([(1.0, y) for y in 1:10], collect(1:10), StableRNG(1))) == 1:10
    @test sort(DT.brio_order([(1.0, 1.0) for _ in 1:10], collect(1:10), StableRNG(1))) == 1:10
    @test sort(DT.brio_order([(NaN, 1.0), (Inf, 2.0), (0.0, 3.0)], collect(1:3), StableRNG(1))) == 1:3

    # The rounds: each vertex is in the last round with probability 1/2, the round before with probability 1/4, and so on.
    # The rounds are sorted along the Hilbert curve, so most consecutive vertices should be close to each other
    rng = StableRNG(2)
    points = [(rand(rng), rand(rng)) for _ in 1:10_000]
    order = DT.brio_order(points, collect(1:10_000), StableRNG(3))
    random_order = shuffle(StableRNG(3), collect(1:10_000))
    total_length(order) = sum(i -> hypot((points[order[i]] .- points[order[i + 1]])...), 1:(length(order) - 1))
    @test total_length(order) < total_length(random_order) / 10
end

@testset "get_insertion_order with the spatial order" begin
    rng = StableRNG(4)
    points = rand(rng, 2, 500)
    order = DT.get_insertion_order(points, true, (3, 7), Int32, StableRNG(1))
    @test eltype(order) == Int32
    @test sort(order) == setdiff(1:500, (3, 7))
    @test order == DT.get_insertion_order(points, true, (3, 7), Int32, StableRNG(1))
    @test DT.get_insertion_order(points, false, (3, 7), Int32, StableRNG(1)) == setdiff(1:500, (3, 7))
end

@testset "The spatial order gives the same triangulation" begin
    for n in (10, 100, 5000)
        rng = StableRNG(n)
        points = [(rand(rng), rand(rng)) for _ in 1:n]
        tri = triangulate(points; rng = StableRNG(1))
        tri_shuffled = triangulate(points; rng = StableRNG(1), insertion_order = shuffle(StableRNG(1), collect(1:n)), num_sample_rule = DT.default_num_samples)
        @test validate_triangulation(tri)
        @test DT.compare_triangle_collections(get_triangles(tri), get_triangles(tri_shuffled))
    end
    # Float32 coordinates, points to skip, and degenerate inputs
    rng = StableRNG(5)
    points = rand(rng, 2, 500)
    @test validate_triangulation(triangulate(Float32.(points); rng = StableRNG(1)))
    tri = triangulate(points; rng = StableRNG(1), skip_points = (3, 7))
    @test validate_triangulation(tri) && !DT.has_vertex(tri, 3) && !DT.has_vertex(tri, 7)
    @test validate_triangulation(triangulate(get_points(triangulate_rectangle(0, 1, 0, 1, 30, 30)); rng = StableRNG(1)))
    for n in 3:8
        rng = StableRNG(n)
        @test validate_triangulation(triangulate([(rand(rng), rand(rng)) for _ in 1:n]; rng = StableRNG(1)))
    end
    @test_throws DT.NonFiniteCoordinatesError triangulate([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (NaN, 0.5)])
end

@testset "num_sample_rule defaults" begin
    @test DT.spatial_num_samples(10) === 1
    @test DT.spatial_num_samples(Int32(10)) === Int32(1)
end
