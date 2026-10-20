using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using ReferenceTests
using StableRNGs
const ASSETS = joinpath(@__DIR__, "..", "..", "..", "assets", "test", "data_structures")


@testset "Computing statistics" begin
    _x, _y = complicated_geometry()
    x = _x
    y = _y
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_1 = triangulate(points; boundary_nodes, delete_ghosts = false)
    A = get_area(tri_1) # 2356
    refine!(tri_1; max_area = 1.0e-3A, use_circumcenter = true)
    boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
    tri_2 = triangulate(points; boundary_nodes, delete_ghosts = false)
    A = get_area(tri_2)
    refine!(tri_2; max_area = 1.0e-3A, use_circumcenter = true)
    boundary_nodes, points = convert_boundary_points_to_indices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
    tri_3 = triangulate(points; boundary_nodes, delete_ghosts = false)
    A = get_area(tri_3)
    refine!(tri_3; max_area = 1.0e-3A, use_circumcenter = true)
    boundary_nodes, points = convert_boundary_points_to_indices(reverse(reverse.(x[2])), reverse(reverse.(y[2])))
    tri_4 = triangulate(points; boundary_nodes, delete_ghosts = false)
    A = get_area(tri_4)
    refine!(tri_4; max_area = 1.0e-3A, use_circumcenter = true)
    a, b = 0.0, 5.0
    c, d = 3.0, 7.0
    nx = 3
    ny = 3
    tri_5 = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts = false, single_boundary = false)
    tri_6 = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts = false, single_boundary = true)
    for (i, tri) in enumerate((tri_1, tri_2, tri_3, tri_4, tri_5, tri_6))
        validate_statistics(tri)
        @inferred statistics(tri)
        stats = statistics(tri)
        validate_statistics(tri, stats)
        delete_ghost_triangles!(tri)
        validate_statistics(tri)
        @inferred statistics(tri)
        stats = statistics(tri)
        validate_statistics(tri, stats)
    end

    pts = [
        (-7.36, 12.55), (-9.32, 8.59), (-9.0, 3.0), (-6.32, -0.27),
        (-4.78, -1.53), (2.78, -1.41), (-5.42, 1.45), (7.86, 0.67),
        (10.92, 0.23), (9.9, 7.39), (8.14, 4.77), (13.4, 8.61),
        (7.4, 12.27), (2.2, 13.85), (-3.48, 10.21), (-4.56, 7.35),
        (3.44, 8.99), (3.74, 5.87), (-2.0, 8.0), (-2.52, 4.81),
        (1.34, 6.77), (1.24, 4.15),
    ]
    boundary_points = [
        (0.0, 0.0), (2.0, 1.0), (3.98, 2.85), (6.0, 5.0),
        (7.0, 7.0), (7.0, 9.0), (6.0, 11.0), (4.0, 12.0),
        (2.0, 12.0), (1.0, 11.0), (0.0, 9.13), (-1.0, 11.0),
        (-2.0, 12.0), (-4.0, 12.0), (-6.0, 11.0), (-7.0, 9.0),
        (-6.94, 7.13), (-6.0, 5.0), (-4.0, 3.0), (-2.0, 1.0), (0.0, 0.0),
    ]
    boundary_nodes, pts = convert_boundary_points_to_indices(boundary_points; existing_points = pts)
    uncons_tri = triangulate(pts, delete_ghosts = false)
    cons_tri = triangulate(pts; boundary_nodes, delete_ghosts = false)
    add_point!(cons_tri, 0.0, 5.0)
    add_segment!(cons_tri, 40, 26)
    add_segment!(cons_tri, 39, 27)
    add_segment!(cons_tri, 38, 28)
    add_point!(cons_tri, -3.0, 12.0)
    validate_statistics(uncons_tri)
    validate_statistics(cons_tri)
    delete_ghost_triangles!(uncons_tri)
    delete_ghost_triangles!(cons_tri)
    validate_statistics(uncons_tri)
    validate_statistics(cons_tri)
end

@testset "tri method for circumcenter" begin
    tri = triangulate(rand(2, 50))
    T = rand(get_triangles(tri))
    @test DT.triangle_circumcenter(tri, T) == DT.triangle_circumcenter(get_point(tri, T...)...)
end

@testset "Sinks" begin
    # For a lattice triangulation, the circumcenters are all on the the edges, so all triangles should be sinks.
    tri1 = triangulate_rectangle(0, 10, 0, 10, 14, 14)
    tri2 = triangulate_rectangle(0, 1, 0, 6, 250, 23)
    for tri in (tri1, tri2)
        stats = statistics(tri)
        sinks = Dict(DT.sort_triangle.(each_solid_triangle(tri)) .=> DT.triangle_sink.(Ref(tri), each_solid_triangle(tri)))
        for T in each_solid_triangle(tri)
            @test is_sink(tri, T)
            @test collect(sinks[DT.sort_triangle(T)]) ≈ collect(DT.triangle_circumcenter(tri, T)) ≈ collect(DT.get_sink(stats, T)) ≈ collect(DT.triangle_sink(tri, T))
            @inferred DT.triangle_sink(tri, T)
        end
    end

    # Make the number of sinks makes sense 
    for _ in 1:100
        tri = triangulate(rand(2, 50))
        sinks = Dict(DT.sort_triangle.(each_solid_triangle(tri)) .=> DT.triangle_sink.(Ref(tri), each_solid_triangle(tri)))
        stats = statistics(tri)
        num_sinks = count(T -> is_sink(tri, T), each_solid_triangle(tri))
        @test length(unique(values(sinks))) == num_sinks
    end

    # A specific example with holes
    tri, label_map, index_map = simple_geometry()
    sinks = Dict(DT.sort_triangle.(each_solid_triangle(tri)) .=> DT.triangle_sink.(Ref(tri), each_solid_triangle(tri)))
    stats = statistics(tri)
    for T in each_solid_triangle(tri) # (23,11,21) bounces back and forth, hence the need for prev_T in the function definitions
        @test collect(sinks[DT.sort_triangle(T)]) ≈ collect(DT.get_sink(stats, T)) ≈ collect(DT.triangle_sink(tri, T))
        @inferred DT.triangle_sink(tri, T)
    end

    # A simpler example 
    tri = example_with_special_corners()
    stats = statistics(tri)
    @test length(get_all_stat(stats, :sink)) == num_solid_triangles(tri)
    p1p2p3 = Vector{Float64}[]
    sink = Vector{Float64}[]
    push!(p1p2p3, [-5.0, 5.0], [-1.0, 2.0], [-1.0, 6.0])
    push!(sink, [-2.625, 4.0])
    push!(p1p2p3, [-5.0, 9.0], [-5.0, 5.0], [-1.0, 6.0])
    push!(sink, [-3.375, 7.0])
    push!(p1p2p3, [0.0, 9.87], [-4.0, 10.0], [-1.0, 6.0])
    push!(sink, collect(DT.triangle_circumcenter([0.0, 9.87], [-4.0, 10.0], [-1.0, 6.0])))
    push!(p1p2p3, [0.0, 9.87], [-1.0, 6.0], [1.0, 7.0])
    push!(sink, collect(DT.triangle_circumcenter([0.0, 9.87], [-4.0, 10.0], [-1.0, 6.0])))
    push!(p1p2p3, [3.0, 8.0], [2.0, 4.0], [5.0, 5.0])
    push!(sink, [3.045454545454545, 5.86363636363636])
    push!(p1p2p3, [-5.0, 11.0], [0.0, 9.87], [3.0, 11.0])
    push!(sink, collect(DT.triangle_circumcenter([-5.0, 11.0], [0.0, 9.87], [3.0, 11.0])))
    push!(p1p2p3, [-5.0, 11.0], [-4.0, 10.0], [0.0, 9.87])
    push!(sink, collect(DT.triangle_circumcenter([-5.0, 11.0], [0.0, 9.87], [3.0, 11.0])))
    push!(p1p2p3, [-5.0, 11.0], [-5.0, 9.0], [-4.0, 10.0])
    push!(sink, [-5.0, 10.0])
    push!(p1p2p3, [-5.0, 9.0], [-1.0, 6.0], [-4.0, 10.0])
    push!(sink, [-2.7857142857143, 7.7857142857143])
    push!(p1p2p3, [5.0, 5.0], [6.0, 1.0], [8.0, 3.0])
    push!(sink, [5.9, 3.1])
    push!(p1p2p3, [5.0, 5.0], [3.0, 3.0], [6.0, 1.0])
    push!(sink, [5.1, 2.9])
    push!(p1p2p3, [6.0, 1.0], [3.0, 3.0], [-1.0, 2.0])
    push!(sink, collect(DT.triangle_circumcenter([6.0, 1.0], [3.0, 3.0], [-1.0, 2.0])))
    push!(p1p2p3, [3.0, 3.0], [2.0, 4.0], [-1.0, 2.0])
    push!(sink, collect(DT.triangle_circumcenter([6.0, 1.0], [3.0, 3.0], [-1.0, 2.0])))
    push!(p1p2p3, [-1.0, 6.0], [-1.0, 2.0], [2.0, 4.0])
    push!(sink, [-0.16666666666666667, 4.0])
    push!(p1p2p3, [1.0, 7.0], [-1.0, 6.0], [2.0, 4.0])
    push!(sink, [0.6428571428571, 5.2142857142857])
    push!(p1p2p3, [5.0, 5.0], [2.0, 4.0], [3.0, 3.0])
    push!(sink, [3.5, 4.5])
    push!(p1p2p3, [6.0, 1.0], [8.0, 1.0], [8.0, 3.0])
    push!(sink, [7.0, 2.0])
    push!(p1p2p3, [7.0, 5.0], [5.0, 5.0], [8.0, 3.0])
    push!(sink, collect(DT.triangle_circumcenter([6.0, 1.0], [8.0, 3.0], [5.0, 5.0])))
    push!(p1p2p3, [3.0, 8.0], [1.0, 7.0], [2.0, 4.0])
    push!(sink, [3.04545454545454545, 5.86363636363636])
    push!(p1p2p3, [0.0, 9.87], [1.0, 7.0], [3.0, 8.0])
    push!(sink, [1.3793100890208, 8.7413798219585])
    push!(p1p2p3, [0.0, 9.87], [3.0, 8.0], [3.0, 11.0])
    push!(sink, [1.852183333333, 9.5])
    push!(p1p2p3, [3.0, 11.0], [3.0, 8.0], [6.0, 8.0])
    push!(sink, [4.5, 9.5])
    push!(p1p2p3, [3.0, 8.0], [5.0, 5.0], [6.0, 8.0])
    push!(sink, [4.5, 6.833333333333333333])
    push!(p1p2p3, [7.0, 5.0], [6.0, 8.0], [5.0, 5.0])
    push!(sink, [6.0, 6.33333333333333])
    push!(p1p2p3, [7.0, 5.0], [8.0, 3.0], [6.0, 8.0])
    push!(sink, [24.5, 12.5])
    p1p2p3 = reshape(p1p2p3, 3, :)
    for ((p1, p2, p3), sink) in zip(eachcol(p1p2p3), sink)
        i, j, k = get_idx(tri, p1, p2, p3)
        @test collect(DT.get_sink(stats, (i, j, k))) ≈ collect(DT.triangle_sink(tri, (i, j, k))) ≈ collect(sink)
    end

    # Single sink 
    r = 1.0
    gθ = θ -> sin(7θ)
    θ = LinRange(0, 2π, 250) |> collect
    θ[end] = θ[begin]
    x = r .* (1 .+ (1 / 10) .* gθ.(θ)) .* cos.(θ)
    y = r .* (1 .+ (1 / 10) .* gθ.(θ)) .* sin.(θ)
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri = triangulate(points; boundary_nodes)
    stats = statistics(tri)
    interior_sinks = NTuple{2, Float64}[]
    for T in each_solid_triangle(tri)
        is_interior_sink(tri, T) && push!(interior_sinks, DT.get_sink(stats, T))
    end
    @test length(interior_sinks) == 1 == count(T -> is_interior_sink(tri, T), each_solid_triangle(tri))
    fig, ax, sc = scatter(get_all_stat(stats, :sink))
    scatter!(ax, interior_sinks, color = :red)
    triplot!(ax, tri)
    fig
    @test_reference reference_path(ASSETS, "sink_figures.png") fig
end

@testset "triangle_offcenter" begin
    p, q, r = (0.0, 2.0), (0.0, 0.0), (6.0, 7.0)
    _validate_offcenter(p, q, r, 1.2)
    _validate_offcenter(p, q, r, 1.2)
    _validate_offcenter(p, q, r, 1.2)
    _validate_offcenter(p, q, r, sqrt(2))
    _validate_offcenter(p, q, r, sqrt(2))
    # For this problem, the circumcenter lies right on the triangle, which caused a problem in the past.
    p, q, r = (0.0, 1.0), (0.0, 0.0), (1.0, 0.0)
    cx, cy = DT.triangle_circumcenter(p, q, r)
    ox, oy = DT.triangle_offcenter(p, q, r)
    _validate_offcenter(p, q, r, 1.0)
    _validate_offcenter(p, q, r, sqrt(2))
    _validate_offcenter(p, q, r, 1.0)
    _validate_offcenter(p, q, r, 1.0)

    # Some random triangles
    for _ in 1:50
        tri = triangulate(rand(2, 50))
        for T in each_solid_triangle(tri)
            p, q, r = get_point(tri, T...)
            _validate_offcenter(p, q, r, 1.0)
            _validate_offcenter(p, q, r, sqrt(2))
        end
    end

    # Triangle with equal edges 
    p, q, r = (0.0, 0.0), (1.0, 1.0), (0.0, 1.0)
    cx, cy = DT.triangle_circumcenter(p, q, r)
    ℓ² = 1.0
    c₁ = (cx, cy)
    β = 1.0
    ℓ₁², ℓ₂², _, idx = DT.squared_triangle_lengths_and_smallest_index(p, q, r)
    if idx == 2 # off-centers are defined relative to the shortest edge, so make pq the shortest edge
        p, q, r = q, r, p
    elseif idx == 3
        p, q, r = r, p, q
    end
    _p, _q, _r = DT.select_shortest_edge_for_offcenter(p, q, r, c₁, ℓ²)
    @test _p == (0.0, 1.0)
    @test _q == (0.0, 0.0)
    @test _r == (1.0, 1.0)
    _p, _q, _r = DT.select_shortest_edge_for_offcenter(q, r, p, c₁, ℓ²)
    @test _p == (0.0, 1.0)
    @test _q == (0.0, 0.0)
    @test _r == (1.0, 1.0)
    _validate_offcenter(p, q, r, β)

    # An equilateral triangle
    p, q, r = (0.0, 0.0), (1.0, 0.0), (0.5, sqrt(3) / 2)
    cx, cy = DT.triangle_circumcenter(p, q, r)
    ℓ² = 1.0
    c₁ = (cx, cy)
    β = 1.0
    α = 0.999
    _p, _q, _r = DT.select_shortest_edge_for_offcenter(p, q, r, c₁, ℓ²)
    @test _p == p # p is the lexicographically smallest point
    @test _q == q
    @test _r == r
    _p, _q, _r = DT.select_shortest_edge_for_offcenter(q, r, p, c₁, ℓ²)
    @test _p == p
    @test _q == q
    @test _r == r
    _p, _q, _r = DT.select_shortest_edge_for_offcenter(r, p, q, c₁, ℓ²)
    @test _p == p
    @test _q == q
    @test _r == r

    # Some lattice triangles 
    tri = triangulate_rectangle(0, 10, 0, 10, 11, 11)
    for T in each_solid_triangle(tri)
        p, q, r = get_point(tri, T...)
        _validate_offcenter(p, q, r, 1.0)
        _validate_offcenter(p, q, r, sqrt(2))
        _validate_offcenter(p, q, r, 1.0)
        _validate_offcenter(p, q, r, 1.0)
        _validate_offcenter(p, q, r, 1.0)
    end
end

@testset "Previously broken centroid" begin
    p = (1.25, -0.8947368421052633)
    q = (1.6666666666666667, 1.2105263157894735)
    r = (2.0833333333333335, 3.3157894736842106)
    c = DT.triangle_centroid(p, q, r)
    @test DT.is_on(DT.point_position_relative_to_triangle(p, q, r, c))
end

@testset "Statistics of a single triangle" begin
    p, q, r = (0.0, 0.0), (3.0, 0.0), (0.0, 4.0)
    @test collect(DT.triangle_lengths(p, q, r)) ≈ [3.0, 4.0, 5.0]
    @test collect(DT.triangle_angles(p, q, r)) ≈ [atan(3, 4), atan(4, 3), π / 2]
    tri = triangulate([p, q, r])
    stats = statistics(tri)
    T = first(each_solid_triangle(tri))
    @test DT.get_minimum_angle(stats, T) ≈ atan(3, 4)
    @test DT.get_median_angle(stats, T) ≈ atan(4, 3)
    @test DT.get_maximum_angle(stats, T) ≈ π / 2
    @test collect(DT.get_lengths(stats, T)) ≈ [3.0, 4.0, 5.0]
    @test DT.get_area(stats, T) ≈ 6.0
    @test DT.get_inradius(stats, T) ≈ 1.0
    @test DT.get_circumradius(stats, T) ≈ 2.5
    @test collect(DT.get_circumcenter(stats, T)) ≈ [1.5, 2.0]
end

@testset "Statistics of a triangulation without any triangles" begin
    tri = Triangulation([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
    stats = statistics(tri)
    @test DT.num_solid_triangles(stats) == 0
    @test isempty(DT.get_all_stat(stats, :area))
    @test DT.get_all_stat(stats, :area) isa Vector{Float64}
    @test DT.get_all_stat(stats, :angles) isa Vector{NTuple{3, Float64}}
end

@testset "Getters for a triangle that is not in the triangulation" begin
    tri = triangulate([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]; rng = StableRNG(1))
    stats = statistics(tri)
    @test_throws KeyError DT.get_area(stats, (1, 2, 100))
    @test_throws KeyError DT.get_minimum_angle(stats, (1, 2, 100))
    T = first(each_solid_triangle(tri))
    @test DT.get_area(stats, T) ≈ 0.5
    @test DT.get_area(stats, (T[2], T[3], T[1])) ≈ 0.5 # rotations of the triangle are found too
end

@testset "distance_to_offcenter" begin
    # The triangle with the shortest edge (p, q) and the offcenter as its third vertex has radius-edge ratio β
    for β in (0.5, 1 / sqrt(3), 1.0, 2.5), ℓ in (1.0e-6, 1.0, 1.0e6)
        h = DT.distance_to_offcenter(β, ℓ)
        p, q, c = (0.0, 0.0), (ℓ, 0.0), (ℓ / 2, h)
        @test DT.triangle_circumradius(p, q, c) / ℓ ≈ β
    end
    @test DT.distance_to_offcenter(1.0, 1.0) ≈ (2 + sqrt(3)) / 2
    @test_throws DomainError DT.distance_to_offcenter(0.4, 1.0)
end
