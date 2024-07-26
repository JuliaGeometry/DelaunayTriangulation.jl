using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using DataStructures
using CairoMakie
using StaticArrays
using StatsBase
using ElasticArrays
using ..DelaunayTriangulation: Certificate

for PT in subtypes(DT.AbstractPredicateType)
    global x, y = complicated_geometry()
    boundary_nodes, points = convert_boundary_points_to_indices(x, y; existing_points=ElasticMatrix{Float64}(undef, 2, 0))
    _tri = triangulate(points; boundary_nodes, delete_ghosts=false, predicates=PT())
    A = get_area(_tri)
    refine!(_tri; max_area=1e-2A, use_circumcenter=true, predicates=PT())
    _pts = ElasticMatrix(get_points(_tri))
    global tri = DT.Triangulation(_pts, _tri.triangles, _tri.boundary_nodes, _tri.interior_segments,
        _tri.all_segments, _tri.weights, _tri.adjacent, _tri.adjacent2vertex, _tri.graph, _tri.boundary_curves,
        _tri.boundary_edge_map, _tri.ghost_vertex_map, _tri.ghost_vertex_ranges, DT.ConvexHull(_pts, _tri.convex_hull.vertices), _tri.representative_point_list,
        _tri.polygon_hierarchy, _tri.boundary_enricher, _tri.cache)
    DT.compute_representative_points!(tri)
    global rep = DT.get_representative_point_list(tri)
    global ghost_vertex_map = DT.get_ghost_vertex_map(tri)
    global pts = DT.get_points(tri)

    @testset "triangle_orientation" begin
        pts = DT.get_points(tri)
        for T in each_triangle(tri)
            i, j, k = triangle_vertices(T)
            cert3 = DT.triangle_orientation(PT(), tri, T)
            cert4 = DT.triangle_orientation(PT(), tri, i, j, k)
            cert5 = DT.triangle_orientation(PT(), tri, i, pts[:, j], k)
            @test all(DT.is_positively_oriented, (cert3, cert4, cert5))
            @inferred DT.triangle_orientation(PT(), tri, T)
            @inferred DT.triangle_orientation(PT(), tri, i, j, k)
            @inferred DT.triangle_orientation(PT(), tri, i, pts[:, j], k)
        end

        p0 = Float64[5, 5]
        p1 = Float64[4.5, 2.5]
        p2 = Float64[2.5, 1.5]
        p3 = Float64[3, 3.5]
        p4 = Float64[0, 2]
        p5 = Float64[1, 5]
        p6 = Float64[1, 3]
        p7 = Float64[4, -1]
        p8 = Float64[-1, 4]
        points = [p0, p1, p2, p3, p4, p5, p6, p7, p8]
        temptri = triangulate(points; predicates=PT())
        @test DT.is_positively_oriented(DT.triangle_orientation(PT(), temptri, (4, 6, 7)))
        @test DT.is_negatively_oriented(DT.triangle_orientation(PT(), temptri, (4, 7, 6)))
        @test DT.is_negatively_oriented(DT.triangle_orientation(PT(), temptri, (4, 2, 3)))
        @test DT.is_positively_oriented(DT.triangle_orientation(PT(), temptri, (4, 7, 3)))
        @test DT.is_positively_oriented(DT.triangle_orientation(PT(), temptri, (5, 7, 9)))
        @test DT.is_negatively_oriented(DT.triangle_orientation(PT(), temptri, (5, 9, 7)))
        @test DT.is_negatively_oriented(DT.triangle_orientation(PT(), temptri, (3, 8, 5)))
        points = [[1.0, 2.0], [1.0, 5.0], [1.0, 8.0]]
        temptrit = triangulate(points; predicates=PT())
        @test DT.is_degenerate(DT.triangle_orientation(PT(), tri, (1, 2, 3)))
        points = [[0, -3.0], [3.0, 0.0], [0.0, 3.0], [-3.0, 0.0]]
        temptri = triangulate(points; predicates=PT())
        @test DT.is_positively_oriented(DT.triangle_orientation(PT(), temptri, (1, 2, 3))) &&
              DT.is_positively_oriented(DT.triangle_orientation(PT(), temptri, (2, 3, 4))) &&
              DT.is_positively_oriented(DT.triangle_orientation(PT(), temptri, (4, 1, 2)))
    end

    @testset "point_position_relative_to_circumcircle" begin
        for T in each_triangle(tri)
            if !DT.is_ghost_triangle(T)
                i, j, k = triangle_vertices(T)
                for ℓ in (i, j, k)
                    cert3 = DT.point_position_relative_to_circumcircle(PT(), tri, T, ℓ)
                    cert4 = DT.point_position_relative_to_circumcircle(PT(), tri, i, j, k, ℓ)
                    cert5 = DT.point_position_relative_to_circumcircle(PT(), tri, pts[:, i], j, k,
                        pts[:, ℓ])
                    cert6 = DT.point_position_relative_to_circumcircle(PT(), tri, i, pts[:, j], k, ℓ)
                    @test all(DT.is_on, (cert3, cert4, cert5, cert6))
                    @inferred DT.point_position_relative_to_circumcircle(PT(), tri, T, ℓ)
                    @inferred DT.point_position_relative_to_circumcircle(PT(), tri, i, j, k, ℓ)
                    @inferred DT.point_position_relative_to_circumcircle(PT(), tri, pts[:, i], j, k,
                        pts[:, ℓ])
                    @inferred DT.point_position_relative_to_circumcircle(PT(), tri, i, pts[:, j], k, ℓ)
                end
                q = (pts[:, i] .+ pts[:, j] .+ pts[:, k]) ./ 3
                append!(pts, q)
                ℓ = size(pts, 2)
                cert3 = DT.point_position_relative_to_circumcircle(PT(), tri, T, ℓ)
                cert4 = DT.point_position_relative_to_circumcircle(PT(), tri, i, j, k, ℓ)
                @test all(DT.is_inside, (cert3, cert4))
                resize!(pts, 2, ℓ - 1)
            end
        end

        @testset "Testing if a point is in a ghost triangle's circumdisk" begin
            p1 = @SVector[-3.32, 3.53]
            p2 = @SVector[-5.98, 2.17]
            p3 = @SVector[-6.36, -1.55]
            p4 = @SVector[-2.26, -4.31]
            p5 = @SVector[6.34, -3.23]
            p6 = @SVector[-3.24, 1.01]
            p7 = @SVector[0.14, -1.51]
            p8 = @SVector[0.2, 1.25]
            p9 = @SVector[1.0, 4.0]
            p10 = @SVector[4.74, 2.21]
            p11 = @SVector[2.32, -0.27]
            pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
            tri = triangulate(pts; delete_ghosts=false,predicates=PT())
            p12 = @SVector[-1.86, 5.99]
            push!(pts, p12)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 1, 9, DT.𝒢, 12))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 9, DT.𝒢, 1, 12))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 1, 9, 12))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 1, DT.𝒢, 9, 12))
            p13 = @SVector[-1.86, 102.9]
            push!(pts, p13)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 1, 9, DT.𝒢, 13))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 9, DT.𝒢, 1, 12))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 1, 9, 12))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 1, DT.𝒢, 9, 13))
            p14 = @SVector[3.54, 4.684]
            push!(pts, p14)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 9, 10, DT.𝒢, 14))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 10, DT.𝒢, 9, 14))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 9, 10, 14))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, 9, 10, DT.𝒢, 9))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, 10, DT.𝒢, 9, 9))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 9, 10, 9))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 10, 9, 8, 14))
            p15 = @SVector[1.57, 2.514]
            push!(pts, p15)
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 9, 10, DT.𝒢, 15))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 10, DT.𝒢, 9, 15))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 9, 10, 15))
            p16 = @SVector[6.77, 0.269]
            push!(pts, p16)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 10, 5, DT.𝒢, 16))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 5, DT.𝒢, 10, 16))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 10, 5, 16))
            p17 = @SVector[4.21754, 3.00067]
            push!(pts, p17)
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 10, 5, DT.𝒢, 17))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 5, DT.𝒢, 10, 17))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 10, 5, 17))
            p18 = @SVector[4.816, -3.696112]
            push!(pts, p18)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 5, 4, DT.𝒢, 18))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 4, DT.𝒢, 5, 18))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 5, 4, 18))
            p19 = @SVector[6.499685, -2.935]
            push!(pts, p19)
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 5, 4, DT.𝒢, 19))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 4, DT.𝒢, 5, 19))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 5, 4, 19))
            p20 = @SVector[2.79587, -4.020351]
            push!(pts, p20)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 5, 4, DT.𝒢, 20))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 4, DT.𝒢, 5, 20))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 5, 4, 20))
            p21 = @SVector[0.0, -4.0]
            push!(pts, p21)
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 5, 4, DT.𝒢, 21))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 4, DT.𝒢, 5, 21))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 5, 4, 21))
            p22 = @SVector[-2.815953, -4.25729]
            push!(pts, p22)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 4, 3, DT.𝒢, 22))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 3, DT.𝒢, 4, 22))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 4, 3, 22))
            p23 = @SVector[-1.4317, -4.3196]
            push!(pts, p23)
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 4, 3, DT.𝒢, 23))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 3, DT.𝒢, 4, 23))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 4, 3, 23))
            p24 = @SVector[-5.04821, -2.54880]
            push!(pts, p24)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 4, 3, DT.𝒢, 24))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 3, DT.𝒢, 4, 24))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 4, 3, 24))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, 4, 3, DT.𝒢, 4))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, 3, DT.𝒢, 4, 4))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 4, 3, 4))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, 4, 3, DT.𝒢, 3))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, 3, DT.𝒢, 4, 3))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 4, 3, 3))
            p25 = @SVector[-6.3327007, -1.7257]
            push!(pts, p25)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 4, 3, DT.𝒢, 25))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 3, DT.𝒢, 4, 25))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 4, 3, 25))
            p26 = @SVector[-6.444937, -0.54101]
            push!(pts, p26)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 3, 2, DT.𝒢, 26))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 2, DT.𝒢, 3, 26))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 3, 2, 26))
            p27 = @SVector[-5.310, 2.87596]
            push!(pts, p27)
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 2, 1, DT.𝒢, 27))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 1, DT.𝒢, 2, 27))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 2, 1, 27))
            p28 = @SVector[-5.247746, 0.905588]
            push!(pts, p28)
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), PT(), tri, 2, 1, DT.𝒢, 28))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 1, DT.𝒢, 2, 28))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, DT.𝒢, 2, 1, 28))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 8, 7, 11, 28))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 8, 11, 10, 28))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 6, 4, 7, 28))
        end

        @testset "Point on the edge of a ghost triangle's circumdisk" begin
            p8 = (6.0, 6.0)
            p9 = (10.0, -2.0)
            p13 = (8.0, 2.0)
            p14 = (0.0, 0.0)
            pts = [p8, p9, p13, p14]
            tri = triangulate(pts; delete_ghosts=false, predicates=PT())
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, 1, 2, DT.𝒢, 3))
            push!(pts, (2.0, 14.0))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 1, 2, DT.𝒢, 5))
            push!(pts, (12.0, -6.0))
            @test DT.is_outside(DT.point_position_relative_to_circumcircle(PT(), tri, 1, 2, DT.𝒢, 6))
            push!(pts, (34.0, -6.0))
            @test DT.is_inside(DT.point_position_relative_to_circumcircle(PT(), tri, 1, 2, DT.𝒢, 7))
            @test DT.is_on(DT.point_position_relative_to_circumcircle(PT(), tri, 1, 2, DT.𝒢, 1))
        end
    end

    @testset "Operations for picking the correct ghost triangle that a point resides in" begin
        p1 = @SVector[-3.32, 3.53]
        p2 = @SVector[-5.98, 2.17]
        p3 = @SVector[-6.36, -1.55]
        p4 = @SVector[-2.26, -4.31]
        p5 = @SVector[6.34, -3.23]
        p6 = @SVector[-3.24, 1.01]
        p7 = @SVector[0.14, -1.51]
        p8 = @SVector[0.2, 1.25]
        p9 = @SVector[1.0, 4.0]
        p10 = @SVector[4.74, 2.21]
        p11 = @SVector[2.32, -0.27]
        pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
        tri = triangulate(pts; delete_ghosts=false, randomise=false)
        p12 = @SVector[3.538447, 3.99844]
        push!(pts, p12)
        @test DT.is_left(DT.point_position_relative_to_line(PT(), tri, 9, 10, 12))
        @test DT.is_left(DT.point_position_relative_to_line(PT(), tri, 10, DT.𝒢, 12))
        @test DT.is_left(DT.point_position_relative_to_line(PT(), tri, DT.𝒢, 9, 12))
        @test DT.is_inside(DT.point_position_relative_to_triangle(PT(), tri, 9, 10, DT.𝒢, 12))
        @test DT.is_outside(DT.point_position_relative_to_triangle(PT(), tri, 9, 8, 10, 12))
        @test DT.brute_force_search(tri, 12) == (10, DT.𝒢, 9)
        p13 = @SVector[-1.182399, 4.5127]
        push!(pts, p13)
        @test DT.point_position_relative_to_line(PT(), tri, 1, 9, 13) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, 9, 1, 13) |> DT.is_right
        @test DT.point_position_relative_to_line(PT(), tri, 9, DT.𝒢, 13) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, DT.𝒢, 1, 13) |> DT.is_left
        @test DT.point_position_relative_to_triangle(PT(), tri, 1, 9, DT.𝒢, 13) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, 9, DT.𝒢, 1, 13) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, DT.𝒢, 1, 9, 13) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, 1, 6, 8, 13) |> DT.is_outside
        @test DT.brute_force_search(tri, 13; predicates=PT()) == (9, DT.𝒢, 1)
        p14 = @SVector[-4.85877, 3.086]
        push!(pts, p14)
        @test DT.point_position_relative_to_line(PT(), tri, 2, 1, 14) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, 1, DT.𝒢, 14) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, DT.𝒢, 2, 14) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, 10, 5, 14) |> DT.is_right
        @test DT.point_position_relative_to_triangle(PT(), tri, 2, 1, DT.𝒢, 14) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, 1, DT.𝒢, 2, 14) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, DT.𝒢, 2, 1, 14) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, 4, 3, DT.𝒢, 14) |> DT.is_outside
        @test DT.brute_force_search(tri, 14) == (2, 1, DT.𝒢)
        p15 = @SVector[-2.0, -5.0]
        push!(pts, p15)
        @test DT.point_position_relative_to_line(PT(), tri, 5, 4, 15) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, 4, DT.𝒢, 15) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, DT.𝒢, 4, 15) |> DT.is_right
        @test DT.point_position_relative_to_line(PT(), tri, DT.𝒢, 5, 15) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, 5, DT.𝒢, 15) |> DT.is_right
        @test DT.point_position_relative_to_triangle(PT(), tri, 5, 4, DT.𝒢, 15) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, 4, DT.𝒢, 5, 15) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, DT.𝒢, 5, 4, 15) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, 6, 7, 8, 15) |> DT.is_outside
        @test DT.point_position_relative_to_triangle(PT(), tri, 10, 5, DT.𝒢, 15) |> DT.is_outside
        @test DT.point_position_relative_to_triangle(PT(), tri, 1, 9, DT.𝒢, 15) |> DT.is_outside
        @test DT.point_position_relative_to_triangle(PT(), tri, 9, DT.𝒢, 1, 15) |> DT.is_outside
        @test DT.point_position_relative_to_triangle(PT(), tri, DT.𝒢, 1, 9, 15) |> DT.is_outside
        @test brute_force_search(tri, 15) == (5, 4, DT.𝒢)
        p16 = @SVector[16.27, 0.92]
        push!(pts, p16)
        @test DT.point_position_relative_to_line(PT(), tri, 10, 5, 16) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, 5, DT.𝒢, 16) |> DT.is_left
        @test DT.point_position_relative_to_line(PT(), tri, DT.𝒢, 10, 16) |> DT.is_left
        @test DT.point_position_relative_to_triangle(PT(), tri, 10, 5, DT.𝒢, 16) |> DT.is_inside
        @test DT.point_position_relative_to_triangle(PT(), tri, 6, 7, 8, 16) |> DT.is_outside
        @test DT.point_position_relative_to_triangle(PT(), tri, 1, 9, DT.𝒢, 16) |> DT.is_outside
        @test DT.point_position_relative_to_triangle(PT(), tri, 2, 6, 3, 16) |> DT.is_outside
        @test DT.point_position_relative_to_triangle(PT(), tri, 3, 2, DT.𝒢, 16) |> DT.is_outside
        @test DT.point_position_relative_to_line(tri, 10, DT.𝒢, 16) |> DT.is_right
        @test brute_force_search(tri, 16; predicates=PT()) == (10, 5, DT.𝒢)
    end

    @testset "Detailed tests for a simple geometry" begin
        tri, label_map, index_map = simple_geometry()
        pts = get_points(tri)
        orig_length = length(pts)
        DT.compute_representative_points!(tri)
        rep = DT.get_representative_point_list(tri)
        rep[1] = DT.RepresentativeCoordinates(mean(reinterpret(reshape,
                    Float64,
                    pts[unique(reduce(vcat,
                        tri.boundary_nodes[1]))]);
                dims=2)..., 0)
        rep[2] = DT.RepresentativeCoordinates(mean(reinterpret(reshape,
                    Float64,
                    pts[unique(reduce(vcat,
                        tri.boundary_nodes[2]))]);
                dims=2)..., 0)
        rep[3] = DT.RepresentativeCoordinates(mean(reinterpret(reshape,
                    Float64,
                    pts[unique(reduce(vcat,
                        tri.boundary_nodes[3]))]);
                dims=2)..., 0)
        ghost_vertex_map = DT.get_ghost_vertex_map(tri)

        @testset "point_position_relative_to_circumcircle" begin
            i, j, k = index_map["w"], index_map["m"], index_map["v"]
            ℓ = index_map["z"]
            T = (i, j, k)
            cert3 = DT.point_position_relative_to_circumcircle(PT(), tri, T, ℓ)
            cert4 = DT.point_position_relative_to_circumcircle(PT(), tri, i, j, k, ℓ)
            @test all(DT.is_inside, (cert3, cert4))
            ℓ = index_map["q"]
            cert3 = DT.point_position_relative_to_circumcircle(PT(), tri, T, ℓ)
            cert4 = DT.point_position_relative_to_circumcircle(PT(), tri, i, j, k, ℓ)
            @test all(DT.is_outside, (cert3, cert4))
        end

        @testset "point_position_relative_to_line" begin
            f1 = (-5.0, 20.0)
            g1 = (2.2179294, 23.7869)
            push!(pts, f1, g1)
            f1_i = length(pts) - 1
            g1_i = length(pts)
            h1 = (21.517785930112, 22.1693108196262)
            i1 = (22.0683240766855, 21.526199282940)
            push!(pts, h1, i1)
            h1_i = length(pts) - 1
            i1_i = length(pts)
            j1 = (6.8116906951308, 14.9338268341201)
            k1 = (7.3767511907439, 10.1308126214084)
            ℓ1 = (6.0582767009799, 9.3303102526232)
            m1 = (4.3630952141405, 10.3662544945806)
            push!(pts, j1, k1, ℓ1, m1)
            j1_i = length(pts) - 3
            k1_i = length(pts) - 2
            ℓ1_i = length(pts) - 1
            m1_i = length(pts)
            n1 = (15.5995300899034, 9.481188316425)
            o1 = (16.9311921105779, 7.7056389555257)
            p1 = (16.4238922931781, 5.2008461071143)
            q1 = (14.0, 4.0)
            r1 = (14.8385803638037, 3.425296746215)
            s1 = (14.3312805464039, 7.4519890468258)
            push!(pts, n1, o1, p1, q1, r1, s1)
            n1_i = length(pts) - 5
            o1_i = length(pts) - 4
            p1_i = length(pts) - 3
            q1_i = length(pts) - 2
            r1_i = length(pts) - 1
            s1_i = length(pts)
            certs = [("w", "v", "f") => Certificate.Left,
                ("c", "d", "k") => Certificate.Left,
                ("r", "q", "j") => Certificate.Right,
                ("r", "q", "d") => Certificate.Left,
                ("u", "w", "i") => Certificate.Left,
                ("c", "d", "e") => Certificate.Collinear,
                ("g", DT.𝒢, f1_i) => Certificate.Left,
                ("g", DT.𝒢, g1_i) => Certificate.Right,
                (DT.𝒢, "g", f1_i) => Certificate.Right,
                (DT.𝒢, "g", g1_i) => Certificate.Left,
                ("a", DT.𝒢, f1_i) => Certificate.Right,
                ("e", DT.𝒢, h1_i) => Certificate.Left,
                ("e", DT.𝒢, i1_i) => Certificate.Right,
                ("ℓ", DT.𝒢 - 1, j1_i) => Certificate.Right,
                (DT.𝒢 - 1, "ℓ", j1_i) => Certificate.Left,
                ("i", DT.𝒢 - 1, j1_i) => Certificate.Left,
                ("i", DT.𝒢 - 1, m1_i) => Certificate.Right,
                (DT.𝒢 - 1, "i", m1_i) => Certificate.Left,
                ("j", DT.𝒢 - 1, "b1") => Certificate.Left,
                ("j", DT.𝒢 - 1, ℓ1_i) => Certificate.Right,
                ("k", DT.𝒢 - 1, ℓ1_i) => Certificate.Left,
                ("k", DT.𝒢 - 1, k1_i) => Certificate.Right,
                (DT.𝒢 - 1, "k", k1_i) => Certificate.Left,
                ("r", DT.𝒢 - 2, n1_i) => Certificate.Right,
                ("r", DT.𝒢 - 2, o1_i) => Certificate.Left,
                (DT.𝒢 - 3, "r", n1_i) => Certificate.Left,
                (DT.𝒢 - 2, "q", o1_i) => Certificate.Left,
                ("q", DT.𝒢 - 2, o1_i) => Certificate.Right,
                ("p", DT.𝒢 - 3, q1_i) => Certificate.Left,
                ("q", DT.𝒢 - 3, p1_i) => Certificate.Left,
                ("o", DT.𝒢 - 2, q1_i) => Certificate.Right,
                ("m", DT.𝒢 - 2, s1_i) => Certificate.Right,
                (DT.𝒢 - 2, "m", s1_i) => Certificate.Left,
                (DT.𝒢 - 2, "p", r1_i) => Certificate.Left,
                ("p", DT.𝒢 - 2, r1_i) => Certificate.Right]
            for ((i, j, u), cert) in certs
                i = i isa String ? index_map[i] : i
                j = j isa String ? index_map[j] : j
                u = u isa String ? index_map[u] : u
                @test DT.point_position_relative_to_line(PT(), tri, i, j, pts[u]) == cert == DT.point_position_relative_to_line(PT(), tri, i, j, u)
                @inferred DT.point_position_relative_to_line(PT(), tri, i, j, pts[u])
                @inferred DT.point_position_relative_to_line(PT(), tri, i, j, u)
            end
        end

        @testset "point_position_on_line_segment" begin
            t1 = (20.0, 25.0)
            u1 = (20.0, -5.0)
            push!(pts, t1, u1)
            t1_i = length(pts) - 1
            u1_i = length(pts)
            i, j = index_map["c"], index_map["e"]
            k = index_map["d"]
            cert1 = DT.point_position_on_line_segment(PT(), tri, i, j, k)
            cert2 = DT.point_position_on_line_segment(PT(), tri, i, j, t1_i)
            cert3 = DT.point_position_on_line_segment(PT(), tri, i, j, u1_i)
            cert4 = DT.point_position_on_line_segment(PT(), tri, j, i, k)
            cert5 = DT.point_position_on_line_segment(PT(), tri, j, i, t1_i)
            cert6 = DT.point_position_on_line_segment(PT(), tri, j, i, u1_i)
            cert7 = DT.point_position_on_line_segment(PT(), tri, j, i, i)
            cert8 = DT.point_position_on_line_segment(PT(), tri, j, i, j)
            cert9 = DT.point_position_on_line_segment(PT(), tri, i, j, i)
            cert10 = DT.point_position_on_line_segment(PT(), tri, i, j, j)
            cert11 = DT.point_position_on_line_segment(PT(), tri, i, j, pts[j])
            cert12 = DT.point_position_on_line_segment(PT(), tri, pts[i], j, i)
            @test all(DT.is_on, (cert1, cert4))
            @test all(DT.is_left, (cert3, cert5))
            @test all(DT.is_right, (cert2, cert6))
            @test all(DT.is_degenerate, (cert7, cert8, cert9, cert10, cert11, cert12))
            v1 = (-3.0, 0.0)
            w1 = (22.0, 0.0)
            push!(pts, v1, w1)
            v1_i = length(pts) - 1
            w1_i = length(pts)
            i, j = index_map["a"], index_map["c"]
            k = index_map["b"]
            cert1 = DT.point_position_on_line_segment(PT(),tri, i, j, k)
            cert2 = DT.point_position_on_line_segment(PT(),tri, i, j, v1_i)
            @inferred DT.point_position_on_line_segment(PT(),tri, i, j, v1_i)
            cert3 = DT.point_position_on_line_segment(PT(),tri, i, j, w1_i)
            cert4 = DT.point_position_on_line_segment(PT(),tri, j, i, k)
            cert5 = DT.point_position_on_line_segment(PT(),tri, j, i, v1_i)
            cert6 = DT.point_position_on_line_segment(PT(),tri, j, i, w1_i)
            cert7 = DT.point_position_on_line_segment(PT(),tri, j, i, i)
            cert8 = DT.point_position_on_line_segment(PT(),tri, j, i, j)
            cert9 = DT.point_position_on_line_segment(PT(),tri, i, j, i)
            cert10 = DT.point_position_on_line_segment(PT(),tri, i, j, j)
            @inferred DT.point_position_on_line_segment(PT(),tri, i, j, j)
            @inferred DT.point_position_on_line_segment(PT(),tri, i, j, pts[j])
            @inferred DT.point_position_on_line_segment(PT(),tri, i, pts[j], j)
            @test all(DT.is_on, (cert1, cert4))
            @test all(DT.is_left, (cert2, cert6))
            @test all(DT.is_right, (cert3, cert5))
            @test all(DT.is_degenerate, (cert7, cert8, cert9, cert10))
        end

        @testset "line_segment_intersection_type" begin
            z1 = (47.93, 17.33)
            a2 = (53.135, 10.66)
            b2 = (54.8, 18.16)
            d2 = (45.7057, 12.1682)
            push!(pts, z1, a2, b2, d2)
            z1_i = length(pts) - 3
            a2_i = length(pts) - 2
            b2_i = length(pts) - 1
            d2_i = length(pts)
            @test DT.is_single(DT.line_segment_intersection_type(PT(),tri, z1_i, a2_i, b2_i, d2_i))
            @test DT.is_single(DT.line_segment_intersection_type(PT(),tri, b2_i, d2_i, z1_i, a2_i))
            @test DT.is_single(DT.line_segment_intersection_type(PT(),tri, b2_i, d2, z1_i, a2))
            @inferred DT.line_segment_intersection_type(PT(),tri, z1_i, a2_i, b2_i, d2_i)
            @inferred DT.line_segment_intersection_type(PT(),tri, z1, a2, b2_i, d2_i)
            e2 = (48.0, 6.0)
            f2 = (54.0, 6.0)
            g2 = (52.0, 8.0)
            h2 = (52.0, 6.0)
            push!(pts, e2, f2, g2, h2)
            e2_i = length(pts) - 3
            f2_i = length(pts) - 2
            g2_i = length(pts) - 1
            h2_i = length(pts)
            @test DT.is_touching(DT.line_segment_intersection_type(PT(),tri, e2_i, f2_i, g2_i, h2_i))
            @test DT.is_touching(DT.line_segment_intersection_type(PT(),tri, g2_i, h2_i, e2_i, f2_i))
            @test DT.is_touching(DT.line_segment_intersection_type(PT(),tri, e2_i, h2_i, h2_i, g2_i))
            @test DT.is_touching(DT.line_segment_intersection_type(PT(),tri, g2_i, h2_i, h2_i, e2_i))
            @inferred DT.line_segment_intersection_type(PT(),tri, g2_i, h2_i, e2_i, f2_i)
            i2 = (53.58, 35.45)
            j2 = (57.11, 27.3)
            k2 = (49.61, 37.12)
            ℓ2 = (56.39, 40.42)
            push!(pts, i2, j2, k2, ℓ2)
            i2_i = length(pts) - 3
            j2_i = length(pts) - 2
            k2_i = length(pts) - 1
            ℓ2_i = length(pts)
            @test DT.is_none(DT.line_segment_intersection_type(PT(),tri, i2_i, j2_i, k2_i, ℓ2_i))
            @test DT.is_none(DT.line_segment_intersection_type(PT(),tri, k2_i, ℓ2_i, i2_i, j2_i))
            m2 = (50.0, 30.0)
            n2 = (50.0, 24.0)
            o2 = (50.0, 28.0)
            p2 = (50.0, 20.0)
            push!(pts, m2, n2, o2, p2)
            m2_i = length(pts) - 3
            n2_i = length(pts) - 2
            o2_i = length(pts) - 1
            p2_i = length(pts)
            @test DT.is_multiple(DT.line_segment_intersection_type(PT(),tri, m2_i, n2_i, o2_i, p2_i))
            @test DT.is_multiple(DT.line_segment_intersection_type(PT(),tri, n2_i, o2_i, m2_i, p2_i))
            @test DT.is_multiple(DT.line_segment_intersection_type(PT(),tri, o2_i, p2_i, m2_i, n2_i))
            @test DT.is_none(DT.line_segment_intersection_type(PT(),tri, p2_i, n2_i, o2_i, m2_i))
        end

        @testset "point_position_relative_to_triangle" begin
            resize!(pts, orig_length)
            f1 = (0.5742824217282, 13.5106352620416)
            g1 = (1.3454554411888, 9.8289060078422)
            h1 = (3.3355793623777, 12.01804232115)
            i1 = (2.5146532448873, 15.4510060852008)
            j1 = (2.1415050096644, 6.9183497731035)
            k1 = (3.1863200682886, 5.5501395772861)
            ℓ1 = (1.0220603039957, 4.0326700873796)
            m1 = (4, 17.0679817711668)
            n1 = (6.0, 2.0)
            o1 = (4.0, 20.0)
            p1 = (8.0, 20.0)
            q1 = (0.0, 2.0)
            r1 = (2.0, 0.0)
            push!(pts, f1, g1, h1, i1, j1, k1, ℓ1, m1, n1, o1, p1, q1, r1)
            r1_i = length(pts)
            q1_i = length(pts) - 1
            p1_i = length(pts) - 2
            o1_i = length(pts) - 3
            n1_i = length(pts) - 4
            m1_i = length(pts) - 5
            ℓ1_i = length(pts) - 6
            k1_i = length(pts) - 7
            j1_i = length(pts) - 8
            i1_i = length(pts) - 9
            h1_i = length(pts) - 10
            g1_i = length(pts) - 11
            f1_i = length(pts) - 12
            certs = [("g", "h", "b1", f1_i) => Certificate.Inside,
                ("b1", "i", "g", f1_i) => Certificate.Outside,
                ("b1", "h", "j", g1_i) => Certificate.Inside,
                ("i", "b1", "j", h1_i) => Certificate.Inside,
                ("i", "ℓ", "f", h1_i) => Certificate.Outside,
                ("h", "a", "s", ℓ1_i) => Certificate.Inside,
                ("h", "s", "j", k1_i) => Certificate.Inside,
                ("s", "j", "h", j1_i) => Certificate.Inside,
                ("i", "g", "b1", i1_i) => Certificate.Inside,
                ("a", "b", "t", r1_i) => Certificate.On,
                ("b", "t", "a", r1_i) => Certificate.On,
                ("a", "t", "s", n1_i) => Certificate.On,
                ("i", "f", "a1", m1_i) => Certificate.On,
                ("g", "i", "a1", m1_i) => Certificate.On,
                ("g", "a1", "f", o1_i) => Certificate.On,
                ("a1", "f", "g", o1_i) => Certificate.On,
                ("h", "a", "s", q1_i) => Certificate.On,
                ("a1", "f", "g", p1_i) => Certificate.On,
                ("b1", "j", "i", "v") => Certificate.Outside]
            for ((i, j, k, u), cert) in certs
                local cert1, cert2, cert3, cert4, cert5, cert6
                i = i isa String ? index_map[i] : i
                j = j isa String ? index_map[j] : j
                k = k isa String ? index_map[k] : k
                u = u isa String ? index_map[u] : u
                T1 = (i, j, k)
                T2 = [i, j, k]
                if !DT.is_ghost_vertex(i)
                    cert1 = DT.point_position_relative_to_triangle(PT(),tri, i, j, k, u)
                else
                    cert1 = DT.point_position_relative_to_triangle(PT(),tri, pts[i], j, pts[k], u, pts)
                    @inferred DT.point_position_relative_to_triangle(PT(),tri, pts[i], j, pts[k], u)
                end
                cert2 = DT.point_position_relative_to_triangle(PT(),tri, T1, u)
                cert3 = DT.point_position_relative_to_triangle(PT(),tri, T2, u)
                cert4 = DT.point_position_relative_to_triangle(PT(),tri, i, j, k, u)
                cert5 = DT.point_position_relative_to_triangle(PT(),tri, T1, u)
                cert6 = DT.point_position_relative_to_triangle(PT(),tri, T2, pts[u])
                @test all(==(cert), (cert1, cert2, cert3, cert4, cert5, cert6))
                @inferred DT.point_position_relative_to_triangle(PT(),tri, i, j, k, u)
                @inferred DT.point_position_relative_to_triangle(PT(),tri, T1, u)
                @inferred DT.point_position_relative_to_triangle(PT(),tri, T2, u)
                @inferred DT.point_position_relative_to_triangle(PT(),tri, T2, pts[u])
            end
            resize!(pts, orig_length)
            f1 = (-4.0, 20.0)
            g1 = (-2.0, 18.0)
            h1 = (-6.0, 12.0)
            i1 = (-6.0, 22.0)
            j1 = (0.0, 18.0)
            k1 = (0.0, 12.0)
            ℓ1 = (-2.0, 10.0)
            m1 = (-4.0, 24.0)
            n1 = (-4.0, 8.0)
            o1 = (0.0, 4.0)
            p1 = (0.0, 22.0)
            push!(pts, f1, g1, h1, i1, j1, k1, ℓ1, m1, n1, o1, p1)
            p1_i = length(pts)
            o1_i = length(pts) - 1
            n1_i = length(pts) - 2
            m1_i = length(pts) - 3
            ℓ1_i = length(pts) - 4
            k1_i = length(pts) - 5
            j1_i = length(pts) - 6
            i1_i = length(pts) - 7
            h1_i = length(pts) - 8
            g1_i = length(pts) - 9
            f1_i = length(pts) - 10
            certs = [("h", "g", DT.𝒢, f1_i) => Certificate.Inside,
                ("h", "g", DT.𝒢, g1_i) => Certificate.Inside,
                ("h", "g", DT.𝒢, h1_i) => Certificate.Inside,
                (DT.𝒢, "h", "g", i1_i) => Certificate.Inside,
                ("g", DT.𝒢, "h", j1_i) => Certificate.On,
                ("h", "g", DT.𝒢, k1_i) => Certificate.On,
                ("h", "g", DT.𝒢, ℓ1_i) => Certificate.Inside,
                (DT.𝒢, "h", "g", m1_i) => Certificate.Inside,
                ("h", "g", DT.𝒢, n1_i) => Certificate.Outside,
                ("h", "g", DT.𝒢, o1_i) => Certificate.Outside,
                ("h", "g", DT.𝒢, p1_i) => Certificate.Outside,
                ("h", "g", DT.𝒢, "b1") => Certificate.Outside]
            rep[1].x = 10.0
            rep[1].y = 10.0
            for ((i, j, k, u), cert) in certs
                local cert4, cert5, cert6
                i = i isa String ? index_map[i] : i
                j = j isa String ? index_map[j] : j
                k = k isa String ? index_map[k] : k
                u = u isa String ? index_map[u] : u
                T1 = (i, j, k)
                T2 = [i, j, k]
                cert4 = DT.point_position_relative_to_triangle(PT(),tri, i, j, k, u)
                cert5 = DT.point_position_relative_to_triangle(PT(),tri, T1, u)
                cert6 = DT.point_position_relative_to_triangle(PT(),tri, T2, u)
                @test all(==(cert), (cert4, cert5, cert6))
                @inferred DT.point_position_relative_to_triangle(PT(),tri, i, j, k, u)
                @inferred DT.point_position_relative_to_triangle(PT(),tri, T1, u)
                @inferred DT.point_position_relative_to_triangle(PT(),tri, T2, u)
            end
            resize!(pts, orig_length)
            f1 = (2.98004, 22.706)
            g1 = (6.449, 23.343)
            h1 = (10.0, 24.0)
            i1 = (12.4666, 22.777)
            j1 = (18.555, 22.600)
            k1 = (24.0, 24.0)
            push!(pts, f1, g1, h1, i1, j1, k1)
            f1_i = length(pts) - 5
            g1_i = length(pts) - 4
            h1_i = length(pts) - 3
            i1_i = length(pts) - 2
            j1_i = length(pts) - 1
            k1_i = length(pts)
            certs = [("g", "f", DT.𝒢, f1_i) => Certificate.Inside,
                ("f", "e", DT.𝒢, f1_i) => Certificate.Outside,
                ("g", "f", DT.𝒢, g1_i) => Certificate.Inside,
                ("e", DT.𝒢, "f", h1_i) => Certificate.Inside,
                (DT.𝒢, "f", "e", h1_i) => Certificate.Inside,
                ("f", "e", DT.𝒢, i1_i) => Certificate.Inside,
                ("e", DT.𝒢, "f", j1_i) => Certificate.Inside,
                (DT.𝒢, "f", "e", k1_i) => Certificate.Inside,
                ("e", "d", DT.𝒢, k1_i) => Certificate.Inside,
                ("b", "a", DT.𝒢, "s") => Certificate.Outside,
                (DT.𝒢, "b", "a", "v") => Certificate.Outside]
            for ((i, j, k, u), cert) in certs
                local cert1, cert2, cert3, cert4, cert5, cert6
                i = i isa String ? index_map[i] : i
                j = j isa String ? index_map[j] : j
                k = k isa String ? index_map[k] : k
                u = u isa String ? index_map[u] : u
                T1 = (i, j, k)
                T2 = [i, j, k]
                cert4 = DT.point_position_relative_to_triangle(PT(),tri, i, j, k, u)
                cert5 = DT.point_position_relative_to_triangle(PT(),tri, T1, u)
                cert6 = DT.point_position_relative_to_triangle(PT(),tri, T2, u)
                @test all(==(cert), (cert4, cert5, cert6))
            end
            resize!(pts, orig_length)
            f1 = (5.765, 12.7676)
            g1 = (6.0, 15.0)
            h1 = (7.0, 12.0)
            i1 = (6.0, 9.0)
            j1 = (6.0, 8.0)
            k1 = (5.0, 12.0)
            ℓ1 = (4.0, 13.0)
            m1 = (8.0, 13.0)
            n1 = (7.0, 10.0)
            o1 = (6.0, 10.0)
            push!(pts, f1, g1, h1, i1, j1, k1, ℓ1, m1, n1, o1, p1)
            p1_i = length(pts)
            o1_i = length(pts) - 1
            n1_i = length(pts) - 2
            m1_i = length(pts) - 3
            ℓ1_i = length(pts) - 4
            k1_i = length(pts) - 5
            j1_i = length(pts) - 6
            i1_i = length(pts) - 7
            h1_i = length(pts) - 8
            g1_i = length(pts) - 9
            f1_i = length(pts) - 10
            certs = [("ℓ", "i", DT.𝒢 - 1, f1_i) => Certificate.Inside,
                ("i", DT.𝒢 - 1, "ℓ", g1_i) => Certificate.Inside,
                (DT.𝒢 - 1, "ℓ", "i", k1_i) => Certificate.Outside,
                ("i", "j", DT.𝒢 - 1, k1_i) => Certificate.Inside,
                ("i", "j", DT.𝒢 - 1, ℓ1_i) => Certificate.On,
                ("j", "k", DT.𝒢 - 1, o1_i) => Certificate.Inside,
                ("j", "k", DT.𝒢 - 1, h1_i) => Certificate.Outside,
                ("i", DT.𝒢 - 1, "ℓ", i1_i) => Certificate.Outside,
                ("k", DT.𝒢 - 1, "j", j1_i) => Certificate.Inside,
                (DT.𝒢 - 1, "j", "k", i1_i) => Certificate.Inside,
                (DT.𝒢 - 1, "k", "ℓ", k1_i) => Certificate.Outside,
                ("k", "ℓ", DT.𝒢 - 1, n1_i) => Certificate.Inside,
                ("k", "ℓ", DT.𝒢 - 1, h1_i) => Certificate.Inside,
                ("ℓ", DT.𝒢 - 1, "k", m1_i) => Certificate.On,
                ("ℓ", DT.𝒢 - 1, "k", "b1") => Certificate.Outside,
                ("j", "k", DT.𝒢 - 1, "s") => Certificate.Outside]
            for ((i, j, k, u), cert) in certs
                local cert4, cert5, cert6
                i = i isa String ? index_map[i] : i
                j = j isa String ? index_map[j] : j
                k = k isa String ? index_map[k] : k
                u = u isa String ? index_map[u] : u
                T1 = (i, j, k)
                T2 = [i, j, k]
                cert4 = DT.point_position_relative_to_triangle(PT(),tri, i, j, k, u)
                cert5 = DT.point_position_relative_to_triangle(PT(),tri, T1, u)
                cert6 = DT.point_position_relative_to_triangle(PT(),tri, T2, u)
                @test all(==(cert), (cert4, cert5, cert6))
            end
            resize!(pts, orig_length)
            f1 = (15.0, 9.0)
            g1 = (15.5, 8.5)
            h1 = (15.0, 7.5)
            i1 = (16.0, 7.0)
            j1 = (17.5, 8.0)
            k1 = (17.5, 9.5)
            ℓ1 = (17.0, 5.0)
            m1 = (15.5, 3.5)
            n1 = (13.5, 2.5)
            o1 = (12.5, 3.5)
            p1 = (14.21, 8.53)
            q1 = (13.49, 4.22)
            r1 = (15.438, 5.66)
            s1 = (13.26, 5.439)
            t1 = (12.0, 4.5)
            push!(pts, f1, g1, h1, i1, j1, k1, ℓ1, m1, n1, o1, p1, q1, r1, s1, t1)
            t1_i = length(pts)
            s1_i = length(pts) - 1
            r1_i = length(pts) - 2
            q1_i = length(pts) - 3
            p1_i = length(pts) - 4
            o1_i = length(pts) - 5
            n1_i = length(pts) - 6
            m1_i = length(pts) - 7
            ℓ1_i = length(pts) - 8
            k1_i = length(pts) - 9
            j1_i = length(pts) - 10
            i1_i = length(pts) - 11
            h1_i = length(pts) - 12
            g1_i = length(pts) - 13
            f1_i = length(pts) - 14
            certs = [("o", "p", DT.𝒢 - 2, o1_i) => Certificate.On,
                ("p", DT.𝒢 - 2, "o", n1_i) => Certificate.On,
                (DT.𝒢 - 2, "o", "p", q1_i) => Certificate.Inside,
                (DT.𝒢 - 3, "o", "p", t1_i) => Certificate.Outside,
                ("o", DT.𝒢 - 3, "n", s1_i) => Certificate.Outside,
                ("n", "o", DT.𝒢 - 3, t1_i) => Certificate.Outside,
                ("p", "q", DT.𝒢 - 2, p1_i) => Certificate.Outside,
                ("p", "q", DT.𝒢 - 3, r1_i) => Certificate.Inside,
                ("q", DT.𝒢 - 3, "p", m1_i) => Certificate.On,
                (DT.𝒢 - 2, "p", "q", ℓ1_i) => Certificate.On,
                ("r", DT.𝒢 - 3, "q", k1_i) => Certificate.Inside,
                ("r", DT.𝒢 - 3, "q", "d") => Certificate.Outside,
                ("q", "r", DT.𝒢 - 2, j1_i) => Certificate.Inside,
                ("r", DT.𝒢 - 3, "q", i1_i) => Certificate.Inside,
                ("r", DT.𝒢 - 3, "q", "m") => Certificate.Outside,
                ("r", "m", DT.𝒢 - 2, f1_i) => Certificate.Inside,
                ("m", DT.𝒢 - 3, "r", g1_i) => Certificate.Inside,
                ("m", DT.𝒢 - 3, "r", r1_i) => Certificate.Outside,
                ("m", DT.𝒢 - 3, "r", h1_i) => Certificate.Inside,
                ("m", "n", DT.𝒢 - 2, p1_i) => Certificate.Inside,
                (DT.𝒢 - 2, "m", "n", s1_i) => Certificate.Outside]
            for ((i, j, k, u), cert) in certs
                local cert4, cert5, cert6
                i = i isa String ? index_map[i] : i
                j = j isa String ? index_map[j] : j
                k = k isa String ? index_map[k] : k
                u = u isa String ? index_map[u] : u
                T1 = (i, j, k)
                T2 = [i, j, k]
                cert4 = DT.point_position_relative_to_triangle(PT(),tri, i, j, k, u)
                cert5 = DT.point_position_relative_to_triangle(PT(),tri, T1, u)
                cert6 = DT.point_position_relative_to_triangle(PT(),tri, T2, u)
                @test all(==(cert), (cert4, cert5, cert6))
            end
        end
    end
end