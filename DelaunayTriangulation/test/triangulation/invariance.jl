using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StableRNGs
using Test

# Multiplying every coordinate by a power of two is exact in floating point arithmetic, so it does not change the
# result of any predicate, and every computed coordinate is scaled by the same factor. An algorithm that only uses
# predicates and relative tolerances should therefore make exactly the same decisions after such a scaling.
const SCALES = (2.0^-20, 2.0^-10, 1.0, 2.0^10, 2.0^20)
scale_points(points, s) = [s .* p for p in points]

@testset "Unconstrained triangulations" begin
    for seed in 1:10
        rng = StableRNG(seed)
        points = [(rand(rng), rand(rng)) for _ in 1:100]
        tri = triangulate(points; rng = StableRNG(seed))
        for s in SCALES
            tri_s = triangulate(scale_points(points, s); rng = StableRNG(seed))
            @test validate_triangulation(tri_s; print_result = false)
            @test DT.compare_triangle_collections(get_triangles(tri_s), get_triangles(tri))
        end
        # Translating the points loses some precision, so the triangulation itself can change
        tri_t = triangulate([p .+ (1.0e8, 1.0e8) for p in points]; rng = StableRNG(seed))
        @test validate_triangulation(tri_t; print_result = false)
    end
end

@testset "Constrained triangulations" begin
    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.3, 0.3), (0.7, 0.3), (0.7, 0.7), (0.3, 0.7)]
    rng = StableRNG(123)
    append!(points, [(rand(rng), rand(rng)) for _ in 1:50])
    filter!(p -> !(0.3 < p[1] < 0.7 && 0.3 < p[2] < 0.7), points) # remove the points inside the hole
    boundary_nodes = [[[1, 2, 3, 4, 1]], [[5, 8, 7, 6, 5]]]
    tri = triangulate(points; boundary_nodes, rng = StableRNG(1))
    for s in SCALES
        tri_s = triangulate(scale_points(points, s); boundary_nodes, rng = StableRNG(1))
        @test validate_triangulation(tri_s; print_result = false)
        @test DT.compare_triangle_collections(get_triangles(tri_s), get_triangles(tri))
        @test get_area(tri_s) ≈ s^2 * get_area(tri)
    end
end

@testset "Clipped Voronoi tessellations" begin
    rng = StableRNG(123)
    points = [(rand(rng), rand(rng)) for _ in 1:50]
    for s in SCALES
        tri = triangulate(scale_points(points, s); rng = StableRNG(1))
        vorn = voronoi(tri; clip = true, rng = StableRNG(1))
        @test validate_tessellation(vorn)
        @test num_polygons(vorn) == 50
        @test sum(get_area(vorn, i) for i in each_polygon_index(vorn)) ≈ get_area(tri)
    end
end

@testset "Mesh refinement" begin
    # The tolerances for precision issues are relative to the size of the domain for domains smaller than 1, and 
    # absolute for larger domains, where they are too small to matter here
    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    nverts = Dict{Float64, Int}()
    for s in SCALES
        tri = triangulate(scale_points(points, s); boundary_nodes = [1, 2, 3, 4, 1], rng = StableRNG(1))
        refine!(tri; max_area = 1.0e-3 * s^2, rng = StableRNG(1))
        @test validate_triangulation(tri; print_result = false)
        nverts[s] = DT.num_solid_vertices(tri)
    end
    for s in SCALES
        @test nverts[s] == nverts[1.0]
    end
end

@testset "Mesh refinement far from the origin" begin
    # A Steiner point used to be rejected if it was within sqrt(eps) of a vertex relative to the size of its coordinates,
    # which for a domain of size 1 at 2^27 ≈ 1.3e8 is 2, so that almost nothing was refined there, and at 2^20 ≈ 1e6
    # the refinement stopped at triangles of size about 0.02. Translating the computed points is not exact, so the
    # results are only compared approximately
    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    tri = triangulate(copy(points); boundary_nodes = [1, 2, 3, 4, 1], rng = StableRNG(1)) # refine! adds to the points it is given
    refine!(tri; max_area = 1.0e-4, rng = StableRNG(1))
    n = DT.num_solid_vertices(tri)
    for t in (2.0^10, 2.0^20, 2.0^27)
        shifted_tri = triangulate([(x + t, y + t) for (x, y) in points]; boundary_nodes = [1, 2, 3, 4, 1], rng = StableRNG(1))
        refine!(shifted_tri; max_area = 1.0e-4, rng = StableRNG(1))
        @test validate_triangulation(shifted_tri; print_result = false)
        @test DT.num_solid_vertices(shifted_tri) ≈ n rtol = 0.05
        @test all(T -> DT.triangle_area(get_point(shifted_tri, triangle_vertices(T)...)...) ≤ 1.0e-4, each_solid_triangle(shifted_tri))
    end
end

@testset "Pole of inaccessibility" begin
    # The pole of inaccessibility of a triangle is its incenter, at a distance of the inradius from the boundary.
    # The default precision is relative to the size of the polygon.
    points = [(0.0, 0.0), (4.0, 0.0), (1.0, 3.0)]
    boundary_nodes = [1, 2, 3, 1]
    inradius = 2DT.triangle_area(points...) / DT.triangle_perimeter(points...)
    for s in SCALES
        pole_s = DT.pole_of_inaccessibility(scale_points(points, s), boundary_nodes)
        dist = DT.distance_to_polygon(pole_s, scale_points(points, s), boundary_nodes)
        @test dist > 0
        @test dist ≥ 0.99s * inradius
    end
end

@testset "Many points at a small scale" begin
    # Every triangle of points at this scale has an area below sqrt(eps), so get_initial_triangle used to shift the
    # insertion order once for each point, recursing until the stack overflowed
    points = 1.0e-5 .* rand(StableRNG(1), 2, 100_000)
    tri = triangulate(points; rng = StableRNG(1))
    @test DT.num_solid_vertices(tri) == 100_000
    @test DT.num_solid_triangles(tri) == 2 * 100_000 - 2 - length(get_convex_hull_vertices(tri)) + 1
end
