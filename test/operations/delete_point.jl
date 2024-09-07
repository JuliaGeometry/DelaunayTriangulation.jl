using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StableRNGs

using StaticArrays

@testset verbose = true "Deleting interior nodes" begin
    @testset "Random point sets" begin
        for PT in subtypes(DT.AbstractPredicateKernel)
            for j in 1:20
                rng1 = StableRNG(j)
                n = rand(rng1, 50:1000)
                points = 20randn(rng1, 2, n)
                tri = triangulate(points; delete_ghosts = false, rng = rng1, predicates = PT())
                deleted_pts = Int[]
                for k in 1:(n ÷ 10)
                    rng2 = StableRNG(j + k * n)
                    i = rand(rng2, each_solid_vertex(tri) |> collect)
                    while DT.is_boundary_node(tri, i)[1]
                        i = rand(rng2, each_solid_vertex(tri) |> collect)
                    end
                    delete_point!(tri, i; rng = rng2, predicates = PT())
                    push!(deleted_pts, i)
                    _tri = triangulate(points; delete_ghosts = false, skip_points = deleted_pts, rng = rng2, predicates = PT())
                    DT.clear_empty_features!(tri)
                    DT.clear_empty_features!(_tri)
                    @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
                    @test get_adjacent(tri) == get_adjacent(_tri)
                    @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
                    @test get_graph(tri) == get_graph(_tri)
                    DT.convex_hull!(tri, predicates = PT())
                    @test validate_triangulation(tri)
                end
            end

            for j in 1:50
                rng = StableRNG(j)
                points = rand(rng, 2, 500)
                tri = triangulate(points; rng, predicates = PT())
                for k in 1:25
                    rng = StableRNG(k)
                    i = rand(rng, each_solid_vertex(tri))
                    while DT.is_boundary_node(tri, i)[1]
                        i = rand(rng, each_solid_vertex(tri))
                    end
                    delete_point!(tri, i, rng = rng, predicates = PT())
                    @test validate_triangulation(tri)
                end
            end
        end
    end

    @testset "Lattice with a single boundary index" begin
        for PT in (DT.ExactKernel, DT.AdaptiveKernel)
            for j in 1:10
                rng1 = StableRNG(j)
                a, b = sort(10randn(rng1, 2))
                c, d = sort(15randn(rng1, 2))
                nx = rand(rng1, 5:25)
                ny = rand(rng1, 5:25)
                tri = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts = false, single_boundary = true, predicates = PT())
                points = get_points(tri)
                n = nx * ny
                for k in 1:(n ÷ 10)
                    rng2 = StableRNG(j + k * n)
                    i = rand(rng2, each_solid_vertex(tri) |> collect)
                    while DT.is_boundary_node(tri, i)[1]
                        i = rand(rng2, each_solid_vertex(tri) |> collect)
                    end
                    delete_point!(tri, i; rng = rng2, predicates = PT())
                    @test validate_triangulation(tri, predicates = PT())
                end
            end
            tri = triangulate_rectangle(0, 1, 0, 1, 25, 25; delete_ghosts = false, single_boundary = true, predicates = PT())
            add_segment!(tri, 7, 7 + 25, predicates = PT())
            @test_throws DT.InvalidVertexDeletionError delete_point!(tri, 2, predicates = PT())
            @test_throws DT.InvalidVertexDeletionError delete_point!(tri, 7, predicates = PT())
            @test_throws DT.InvalidVertexDeletionError delete_point!(tri, 7 + 25, predicates = PT())
            @test_throws DT.InvalidVertexDeletionError delete_point!(tri, -1, predicates = PT())
        end
    end
    @testset "Lattice with multiple boundary indices" begin
        for PT in (DT.ExactKernel, DT.AdaptiveKernel)
            for j in 1:10
                rng1 = StableRNG(j)
                a, b = sort(10randn(rng1, 2))
                c, d = sort(15randn(rng1, 2))
                nx = rand(rng1, 5:25)
                ny = rand(rng1, 5:25)
                tri = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts = false, single_boundary = false, predicates = PT())
                points = get_points(tri)
                n = nx * ny
                for k in 1:(n ÷ 10)
                    rng2 = StableRNG(j + k * n)
                    i = rand(rng2, each_solid_vertex(tri) |> collect)
                    while DT.is_boundary_node(tri, i)[1]
                        i = rand(rng2, each_solid_vertex(tri) |> collect)
                    end
                    delete_point!(tri, i; rng = rng2, predicates = PT())
                    @test validate_triangulation(tri, predicates = PT())
                end
            end
        end
    end
end

@testset "Deleting interior and boundary points for a specific example" begin
    tri = example_with_special_corners()
    rng = StableRNG(292929)
    point = 16
    delete_point!(tri, 16; rng)
    @test validate_triangulation(tri)
    _tri = Triangulation(get_points(tri))
    true_T = [
        10 9 11
        11 9 15
        11 15 12
        10 11 12
        10 12 13
        12 18 13
        13 18 5
        18 4 5
        18 17 4
        14 17 18
        12 14 18
        12 15 14
        14 15 17
        9 8 15
        15 8 7
        15 7 17
        17 6 4
        17 7 6
        6 7 3
        4 6 3
        4 3 1
        1 3 2
        5 4 1
    ]
    for T in eachrow(true_T)
        add_triangle!(_tri, T; update_ghost_edges = true)
    end
    convex_hull!(tri)
    DT.compute_representative_points!(tri)
    DT.clear_empty_features!(_tri)
    DT.clear_empty_features!(tri)
    @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
    @test get_adjacent(tri) == get_adjacent(_tri)
    @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
    @test get_graph(tri) == get_graph(_tri)
    _tri = triangulate(get_points(tri); skip_points = 16, delete_ghosts = false)
    DT.clear_empty_features!(_tri)
    @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
    @test get_adjacent(tri) == get_adjacent(_tri)
    @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
    @test get_graph(tri) == get_graph(_tri)
    @test get_convex_hull(tri) == get_convex_hull(_tri)
    deleted_points = [16]
    points_to_delete = [15, 14, 17, 11, 18, 12]
    local ctr = 5
    for j in 1:100
        _deleted_points = deepcopy(deleted_points)
        _points_to_delete = deepcopy(points_to_delete)
        ctr += 1
        tri = example_with_special_corners()
        tri.points[14] = [1.2, 7.0]
        rng = StableRNG(292929)
        point = 16
        delete_point!(tri, 16; rng)
        while !isempty(_points_to_delete)
            rng = StableRNG(ctr)
            ctr += 1
            i = rand(rng, eachindex(_points_to_delete))
            point = _points_to_delete[i]
            deleteat!(_points_to_delete, i)
            push!(_deleted_points, point)
            rng = StableRNG(ctr)
            ctr += 1
            delete_point!(tri, point; rng)
            convex_hull!(tri)
            DT.compute_representative_points!(tri)
            _tri = triangulate(get_points(tri); skip_points = _deleted_points, delete_ghosts = false, rng)
            @test validate_triangulation(tri)
            rng = StableRNG(ctr)
            ctr += 1
            DT.clear_empty_features!(tri)
            DT.clear_empty_features!(_tri)
            @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
            @test get_adjacent(tri) == get_adjacent(_tri)
            @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
            @test get_graph(tri) == get_graph(_tri)
            @test get_convex_hull(tri) == get_convex_hull(_tri)
        end
    end
end

@testset "Issue #96" begin
    for PT in subtypes(DT.AbstractPredicateKernel)
        j = 96
        rng = StableRNG(j)
        points = rand(rng, 2, 500)
        tri = triangulate(points; rng, predicates = PT())
        orig_tri = deepcopy(tri)
        k = 19
        rng = StableRNG(k)
        i = 183
        delete_point!(tri, i, rng = rng, predicates = PT())
        @test validate_triangulation(tri)
        add_point!(tri, i, rng = rng, predicates = PT())
        @test tri == orig_tri
        @test validate_triangulation(tri)
    end end
