using ..DelaunayTriangulation
using Test
using MacroTools
using StableRNGs
const DT = DelaunayTriangulation

@testset "Check that validate_triangulation breaks on triangulations that it should break on" begin
    rng = StableRNG(1)
    for _ in 1:20
        a = [0.0, 0.0]
        b = [3.0, 0.0]
        c = [5.0, 0.0]
        d = [6.0, 0.0]
        e = [2.0, 4.0]
        f = [2.0, 2.0]
        g = [2.0, -2.0]
        h = [6.0, 4.0]
        pts = [a, b, c, d, e, f, g, h]
        tri = triangulate(pts; delete_ghosts=false, rng)
        @test validate_triangulation(tri; print_result=false)
        DT.delete_triangle!(tri, 6, 2, 3)
        @test !validate_triangulation(tri; print_result=false)
        DT.add_triangle!(tri, 6, 3, 2)
        @test !test_state(test_triangle_orientation(tri))
        @test !validate_triangulation(tri; print_result=false)
        tri.points[7] = [2, -0.5]
        @test !validate_triangulation(tri; print_result=false)
        tri.points[7] = [2.0, -2.0]
        tri.adjacent.adjacent[(6, 3)] = DT.∅
        @test !test_state(test_each_edge_has_two_incident_triangles(tri))
        @test !validate_triangulation(tri; print_result=false)
        tri.adjacent.adjacent[(6, 3)] = 2
        tri.adjacent.adjacent[(6, 17)] = 10
        @test !test_state(test_adjacent2vertex_map_matches_adjacent_map(tri))
        @test !test_state(test_adjacent_map_matches_adjacent2vertex_map(tri))
        tri = triangulate(pts; delete_ghosts=false, rng)
        DT.delete_vertex!(tri, 3)
        @test !validate_triangulation(tri; print_result=false)
        tri = triangulate(pts; delete_ghosts=false, rng)
        push!(tri.triangles, (11, 17, 20))
        @test any(!test_state, test_iterators(tri))
        tri = triangulate(pts; delete_ghosts=true, rng)
        push!(tri.triangles, (11, 17, -20))
        @test any(!test_state, test_iterators(tri))
        tri = triangulate(pts; delete_ghosts=false, rng)
        push!(tri.triangles, (11, 17, 20000))
        @test any(!test_state, test_iterators(tri))
    end
end

@testset "Validating constrained triangulations" begin
    rng = StableRNG(2)
    _x, _y = complicated_geometry()
    x = _x
    y = _y
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_1 = triangulate(points; boundary_nodes, rng)
    A = get_area(tri_1)
    refine!(tri_1; max_area=1.0e-3A, use_circumcenter=true, rng)
    @test validate_triangulation(tri_1)

    boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
    tri_2 = triangulate(points; boundary_nodes, delete_ghosts=false, rng)
    A = get_area(tri_2)
    refine!(tri_2; max_area=1.0e-3A, use_circumcenter=true, rng)
    @test validate_triangulation(tri_2)

    x = [0.0, 2.0, 2.0, 0.0, 0.0]
    y = [0.0, 0.0, 2.0, 2.0, 0.0]
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_3 = triangulate(points; boundary_nodes, delete_ghosts=false, rng)
    A = get_area(tri_3)
    refine!(tri_3; max_area=1.0e-3A, use_circumcenter=true, rng)
    @test validate_triangulation(tri_3)

    x = reverse(reverse.(_x[2]))
    y = reverse(reverse.(_y[2]))
    boundary_nodes, points = convert_boundary_points_to_indices(x, y)
    tri_4 = triangulate(points; boundary_nodes, delete_ghosts=false, rng)
    A = get_area(tri_4)
    refine!(tri_4; max_area=1.0e-3A, use_circumcenter=true, rng)
    @test validate_triangulation(tri_4)
end

@testset "is_conformal" begin
    rng = StableRNG(3)
    tri = triangulate(rand(rng, 2, 50); rng)
    @test is_conformal(tri)

    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.48, 0.5)]
    tri = triangulate(points; rng)
    @test is_conformal(tri)
    add_segment!(tri, 1, 3; rng)
    @test !is_conformal(tri)
end

@testset "test_visibility" begin
    rng = StableRNG(4)
    points = [(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (0.0, 6.0), (1.0, 1.0), (3.5, 5.5), (4.5, 4.0)]
    boundary_nodes = [1, 2, 3, 4, 1]
    tri = triangulate(points; boundary_nodes, segments=Set([(2, 4)]), rng)
    @test DT.is_visible(DT.test_visibility(tri, 4, 2, 3))
    @test DT.is_invisible(DT.test_visibility(tri, 5, 1, 7))
end

# A testset that records the results of its tests without failing the enclosing testset, so that tests that fail can be tested
mutable struct CaptureTestSet <: Test.AbstractTestSet
    description::String
    results::Vector{Any}
    CaptureTestSet(description; kwargs...) = new(description, [])
end
Test.record(ts::CaptureTestSet, result) = (push!(ts.results, result); result)
Test.finish(ts::CaptureTestSet) = ts

@testset "test_with_log" begin
    # With an empty message, this is just @test
    rhs = MacroTools.striplines(@macroexpand1 @test_with_log "" f())
    lhs = MacroTools.striplines(:(@test f()))
    @test lhs == rhs

    # Otherwise each test is recorded as usual, and the message is only printed for a test that does not pass
    ts, printed = mktemp() do path, io
        ts = redirect_stdout(io) do
            @testset CaptureTestSet "capture" begin
                @test_with_log "first" 2 + 4 == 6
                @test_with_log "second" 2 + 4 == 7
                @test_with_log "MeC" 2 ≈ 2.0001 rtol = 1e-3 atol = 1e-2
                @test_with_log "third" 2 ≈ 4 rtol = 1e-4 atol = 1e-2
            end
        end
        flush(io)
        ts, read(path, String)
    end
    @test [r isa Test.Pass for r in ts.results] == [true, false, true, false]
    @test !occursin("first", printed) && occursin("Context: second", printed)
    @test !occursin("MeC", printed) && occursin("Context: third", printed)
end

@testset "check" begin
    failures = Dict{String, Any}()
    for x in 1:5
        @check failures x x < 3
        @check failures x x ≈ 2.0 atol = 0.5
        @check failures x x > 0
    end
    # Only the first failure of each check is kept
    @test failures == Dict("x < 3" => 3, "x ≈ 2.0 atol = 0.5" => 1)
end
