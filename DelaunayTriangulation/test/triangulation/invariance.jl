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
    # The precision checks in get_steiner_point are absolute, so refinement stops early for small domains
    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    nverts = Dict{Float64, Int}()
    for s in SCALES
        tri = triangulate(scale_points(points, s); boundary_nodes = [1, 2, 3, 4, 1], rng = StableRNG(1))
        refine!(tri; max_area = 1.0e-3 * s^2, rng = StableRNG(1))
        @test validate_triangulation(tri; print_result = false)
        nverts[s] = DT.num_solid_vertices(tri)
    end
    for s in SCALES
        if s < 1
            @test_broken nverts[s] == nverts[1.0]
        else
            @test nverts[s] == nverts[1.0]
        end
    end
end

@testset "Pole of inaccessibility" begin
    # The pole of inaccessibility of a triangle is its incenter, at a distance of the inradius from the boundary.
    # The default precision is absolute (precision = 1), so for polygons that are not much larger than 1 the
    # search stops before refining any cells and returns the centroid, which is only 90% of the way there.
    points = [(0.0, 0.0), (4.0, 0.0), (1.0, 3.0)]
    boundary_nodes = [1, 2, 3, 1]
    inradius = 2DT.triangle_area(points...) / DT.triangle_perimeter(points...)
    for s in SCALES
        pole_s = DT.pole_of_inaccessibility(scale_points(points, s), boundary_nodes)
        dist = DT.distance_to_polygon(pole_s, scale_points(points, s), boundary_nodes)
        @test dist > 0
        if s ≤ 1
            @test_broken dist ≥ 0.99s * inradius
        else
            @test dist ≥ 0.99s * inradius
        end
    end
end
