using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using Random
using StableRNGs

@testset "Random constrained Delaunay triangulations" begin
    for PT in subtypes(DT.AbstractPredicateType)
        for i in 1:8
            @info "Testing random constrained Delaunay triangulations. Run: $i; Block: 1."
            rng = StableRNG(i)
            points, edges, mat_edges = get_random_vertices_and_constrained_edges(40, 100, 20, rng)
            # Need to deepcopy edges below, else it gets changed and updated on the first call to tri, which changes the insertion order of the segments and thus comparing tri to _tri might not work
            tri = triangulate(points; segments=deepcopy(edges), rng=StableRNG(i), predicates=PT())
            if i % 5 == 0
                _tri = retriangulate(tri; segments=deepcopy(edges), rng=StableRNG(i), predicates=PT())
                @inferred retriangulate(tri; segments=deepcopy(edges), rng=StableRNG(i), predicates=PT())
                @test tri == _tri
            end
            @test validate_triangulation(tri; predicates=PT())
            empty!(get_all_segments(tri))
            @test !validate_triangulation(tri; predicates=PT(), print_result=false)
        end
        for i in 1:4
            @info "Testing random constrained Delaunay triangulations. Run: $i; Block: 2."
            rng = StableRNG(i^5)
            points, edges, mat_edges = get_random_vertices_and_constrained_edges(200, 500, 100, rng)
            tri = triangulate(points; segments=edges, rng, predicates=PT())
            @test validate_triangulation(tri; predicates=PT())
            empty!(get_all_segments(tri))
            @test !validate_triangulation(tri; predicates=PT(), print_result=false)
        end
    end
end

@testset "Testing Shewchuk's PSLG example" begin
    pts, C = second_shewchuk_example_constrained()
    for PT in subtypes(DT.AbstractPredicateType)
        for i in 1:500
            rng = StableRNG(i^6)
            tri = triangulate(pts; segments=C, rng, predicates=PT())
            @test validate_triangulation(tri; predicates=PT())
        end
    end
end

@testset "Random parabolas" begin
    for PT in (DT.Exact, DT.Adaptive)
        for i in 1:5
            @info "Testing random parabolas. Run: $i. Predicates: $PT"
            rng = StableRNG(i)
            np, nx = 100, 26
            pts = [(2rand(rng) - 1, rand(rng)) for _ in 1:100]
            x = LinRange(-1, 1, 26)
            a = 10.0 .^ (LinRange(0, log10(2), 20)) .- 1
            C = Set{NTuple{2,Int}}()
            for i in eachindex(a)
                y = a[i] * x .^ 2
                append!(pts, zip(x, y))
                push!(C, [(j, j + 1) for j in (np+nx*(i-1)+1):(np+nx*(i-1)+(nx-1))]...)
            end
            tri = triangulate(pts; segments=C, rng, predicates=PT())
            @test validate_triangulation(tri; predicates=PT())
        end
    end
end

@testset "Random collection of straight lines" begin
    for PT in subtypes(DT.AbstractPredicateType)
        for i in 1:4
            @info "Testing random collection of straight lines. Run: $i. Predicates: $PT"
            rng = StableRNG(i)
            pts = NTuple{2,Float64}[]
            C = Set{NTuple{2,Int}}()
            j = 1
            for i in 1:10
                push!(pts, (2i / 11 - 1, 2rand(rng) - 1))
                push!(pts, (2i / 11 - 1, 2rand(rng) - 1))
                push!(C, (j, j + 1))
                j += 2
            end
            x1 = LinRange(-1, 1 - 1e-12, 10)
            y1 = LinRange(-1, -1, 10)
            x2 = LinRange(1, 1, 10)
            y2 = LinRange(-1, 1 - 1e-12, 10)
            x3 = LinRange(1, -1 + 1e-12, 10)
            y3 = LinRange(1, 1, 10)
            x4 = LinRange(-1, -1, 10)
            y4 = LinRange(1, -1 + 1e-12, 10)
            append!(pts, zip(x1, y1), zip(x2, y2), zip(x3, y3), zip(x4, y4))
            push!(C, [(j, j + 1) for j in 21:29]...)
            push!(C, [(j, j + 1) for j in 31:39]...)
            push!(C, [(j, j + 1) for j in 41:49]...)
            push!(C, [(j, j + 1) for j in 51:59]...)
            tri = triangulate(pts; segments=C, rng, predicates=PT())
            @test validate_triangulation(tri; predicates=PT())
        end
    end
end

@testset "Lattice" begin
    for PT in (DT.Exact, DT.Adaptive)
        for m in 1:2
            @info "Testing dense lattice. Run: $m. Predicates: $PT"
            rng = StableRNG(m)
            a = 0.0
            b = 5.0
            c = -3.0
            d = 7.0
            nx = 13
            ny = 20
            tri = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts=false, single_boundary=true, predicates=PT())
            add_segment!(tri, 56, 162; rng, predicates=PT())
            for e in [(1, 249), (1, 250), (1, 251), (1, 26), (1, 39), (1, 52)]
                add_segment!(tri, e; rng, predicates=PT())
            end
            add_segment!(tri, 190, 99; rng, predicates=PT())
            for e in [(99, 113), (113, 101), (101, 115)]
                add_segment!(tri, e; rng, predicates=PT())
            end
            @test validate_triangulation(tri; predicates=PT())

            a = -0.1
            b = 0.1
            c = -0.01
            d = 0.01
            nx = 25
            ny = 25
            tri = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts=false, single_boundary=true, predicates=PT())
            tri = triangulate(get_points(tri); predicates=PT())
            for i in 2:24
                add_segment!(tri, i, 600 + i; rng, predicates=PT())
            end
            @test validate_triangulation(tri; predicates=PT())
            tri = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts=false, single_boundary=true, predicates=PT())
            tri = triangulate(get_points(tri); rng, predicates=PT())
            for e in [(1, 28), (28, 103), (103, 180), (180, 625), (625, 523)]
                add_segment!(tri, e; rng, predicates=PT())
            end
            for e in [(437, 614), (527, 602), (528, 603), (555, 605)]
                add_segment!(tri, e; rng, predicates=PT())
            end
            @test validate_triangulation(tri; predicates=PT())
        end
        for m in 1:10
            rng = StableRNG(m)
            @info "Testing coarse lattice. Run: $m. Predicates: $PT"
            a = 0.0
            b = 1.0
            c = 0.0
            d = 1.0
            nx = 2
            ny = 2
            tri = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts=false, single_boundary=true, predicates=PT())
            add_segment!(tri, 1, 4; rng, predicates=PT())
            @test validate_triangulation(tri; predicates=PT())

            a = 0
            b = 1
            c = 0
            d = 5
            nx = 25
            ny = 3
            tri = triangulate_rectangle(a, b, c, d, nx, ny; predicates=PT(), delete_ghosts=false, single_boundary=true)
            tri = triangulate(get_points(tri); rng, predicates=PT())
            for i in 1:(nx-1)
                u = i
                v = 2nx
                add_segment!(tri, u, v; rng, predicates=PT())
            end
            for i in 51:75
                u = i
                v = 26
                add_segment!(tri, u, v; rng, predicates=PT())
            end
            @test validate_triangulation(tri; predicates=PT())
        end
    end
end

@testset "Triangulating with a deleted exterior" begin
    for PT in subtypes(DT.AbstractPredicateType)
        for i in 1:20
            @info "Testing triangulation of a domain with a hole: Run $i. Predicates: $PT"
            rng = StableRNG(i)
            pts = [(rand(rng), rand(rng)) for _ in 1:50]
            bnd_pts = [(0.3cos(θ), 0.3sin(θ)) .+ 0.5 for θ in LinRange(0, 2π - 1 / 250, 25)]
            bnd_id = [(51:75)..., 51]
            append!(pts, bnd_pts)
            tri = triangulate(pts; boundary_nodes=bnd_id, rng, predicates=PT())
            @test validate_triangulation(tri; predicates=PT())
            _tri = retriangulate(tri; predicates=PT())
            @inferred retriangulate(tri; predicates=PT())
            @test tri == _tri
        end
    end
end

@testset "Triangulation with two curves" begin
    for PT in (DT.Exact, DT.Adaptive)
        for i in 1:5
            @info "Testing triangulation of a domain with two curves: Run $i. Predicates: $PT"
            rng = StableRNG(i)
            pts = [(rand(rng), rand(rng)) for _ in 1:50]
            x1 = LinRange(0, 1, 100)
            y1 = LinRange(0, 0, 100)
            x2 = LinRange(1, 1, 100)
            y2 = LinRange(0, 1, 100)
            x3 = LinRange(1, 0, 100)
            y3 = LinRange(1, 1, 100)
            x4 = LinRange(0, 0, 100)
            y4 = LinRange(1, 0, 100)
            x = [x1..., x2..., x3..., x4...]
            y = [y1..., y2..., y3..., y4...]
            outer_square = [(x, y) for (x, y) in zip(x, y)] |> unique
            push!(outer_square, outer_square[begin])
            outer_square_x = [first.(outer_square)]
            outer_square_y = [last.(outer_square)]
            circ_pts = [(0.3cos(θ), 0.3sin(θ)) .+ 0.5 for θ in LinRange(2π, 0, 50)]
            circ_pts[end] = circ_pts[1]
            inner_circle_x = [first.(circ_pts)]
            inner_circle_y = [last.(circ_pts)]
            x = [outer_square_x, inner_circle_x]
            y = [outer_square_y, inner_circle_y]
            nodes, pts = convert_boundary_points_to_indices(x, y; existing_points=pts)
            tri = triangulate(pts; boundary_nodes=nodes, rng, predicates=PT())
            @test validate_triangulation(tri; predicates=PT())
        end
    end
end

@testset "Adding points into a constrained triangulation; no collinearities" begin
    for PT in subtypes(DT.AbstractPredicateType)
        for L in 1:4
            @info "Testing the addition of points into a constrained triangulation. Run: $L. Predicates: $PT"
            pts, C = example_for_testing_add_point_on_constrained_triangulation()
            tri = triangulate(pts; segments=C, delete_ghosts=false, predicates=PT())
            @test validate_triangulation(tri; predicates=PT())
            DT.push_point!(tri, 2, 1.8)
            add_point!(tri, 15; predicates=PT())
            DT.push_point!(tri, 1.57, 1.778)
            add_point!(tri, 16; predicates=PT())
            T = [
                (1, 10, 11)
                (10, 6, 11)
                (9, 4, 5)
                (4, 6, 10)
                (9, 5, 1)
                (4, 10, 5)
                (13, 1, 14)
                (5, 10, 1)
                (8, 9, 1)
                (13, 8, 1)
                (3, 9, 8)
                (13, 15, 12)
                (12, 2, 7)
                (11, 6, 2)
                (8, 12, 3)
                (8, 13, 12)
                (12, 7, 3)
                (2, 6, 7)
                (15, 2, 12)
                (1, 2, 14)
                (16, 1, 11)
                (14, 15, 13)
                (14, 2, 15)
                (2, 16, 11)
                (2, 1, 16)
                (4, 9, DT.𝒢)
                (9, 3, DT.𝒢)
                (3, 7, DT.𝒢)
                (7, 6, DT.𝒢)
                (6, 4, DT.𝒢)
            ]
            @test DT.compare_triangle_collections(get_triangles(tri), T)
            @test validate_triangulation(tri; predicates=PT())
            (x, y) = (2.3258217552204066, 1.4540267924574883) # a while loop was used to find this point that broke the triangulation (fixed now, obviously)
            DT.push_point!(tri, x, y)
            add_point!(tri, DT.num_points(tri); predicates=PT())
            T = [
                (1, 10, 11)
                (10, 6, 11)
                (9, 4, 5)
                (4, 6, 10)
                (9, 5, 1)
                (4, 10, 5)
                (13, 1, 14)
                (5, 10, 1)
                (8, 9, 1)
                (13, 8, 1)
                (3, 9, 8)
                (13, 15, 17)
                (12, 2, 7)
                (11, 6, 2)
                (8, 12, 3)
                (8, 13, 17)
                (12, 7, 3)
                (2, 6, 7)
                (8, 17, 12)
                (1, 2, 17)
                (16, 1, 11)
                (14, 15, 13)
                (17, 2, 12)
                (2, 16, 11)
                (2, 1, 16)
                (14, 17, 15)
                (14, 1, 17)
                (4, 9, DT.𝒢)
                (9, 3, DT.𝒢)
                (3, 7, DT.𝒢)
                (7, 6, DT.𝒢)
                (6, 4, DT.𝒢)
            ]
            @test DT.compare_triangle_collections(get_triangles(tri), T)
            @test validate_triangulation(tri; predicates=PT())
            for i in 1:25
                x = 4rand()
                y = 5rand()
                DT.push_point!(tri, x, y)
                add_point!(tri, DT.num_points(tri); predicates=PT())
                @test validate_triangulation(tri; predicates=PT())
                @test DT.edge_exists(tri, 1, 2) && DT.edge_exists(tri, 2, 1)
            end
        end
    end
end

@testset "Adding points into a constrained triangulation; interior segment collinearities" begin
    for PT in subtypes(DT.AbstractPredicateType)
        for m in 1:3
            @info "Testing the addition of points into a constrained triangulation with interior segment collinearities. Run: $m. Predicates: $PT"
            pts, C = example_for_testing_add_point_on_constrained_triangulation()
            push!(C, (1, 12))
            tri = triangulate(pts; segments=C, delete_ghosts=false, predicates=PT())
            new_points = [
                (1.0, 3.0),
                (2.0, 3.0),
                (3.0, 3.0),
                (4.0, 3.0)
            ]

            DT.push_point!(tri, new_points[1])
            add_point!(tri, DT.num_points(tri), predicates=PT())
            @test sort_edge_vector(collect(get_interior_segments(tri))) ==
                  sort_edge_vector(collect(Set([(1, 2), (1, 15), (15, 12)]))) ==
                  sort_edge_vector(collect(get_interior_segments(tri)))
            T = [
                (1, 15, 5)
                (10, 6, 11)
                (6, 10, 4)
                (1, 10, 11)
                (1, 5, 10)
                (4, 10, 5)
                (15, 12, 8)
                (15, 9, 5)
                (8, 9, 15)
                (9, 4, 5)
                (14, 13, 15)
                (14, 2, 13)
                (12, 3, 8)
                (6, 2, 11)
                (3, 9, 8)
                (2, 7, 12)
                (1, 14, 15)
                (1, 2, 14)
                (12, 15, 13)
                (13, 2, 12)
                (7, 3, 12)
                (6, 7, 2)
                (2, 1, 11)
                (4, 9, DT.𝒢)
                (9, 3, DT.𝒢)
                (3, 7, DT.𝒢)
                (7, 6, DT.𝒢)
                (6, 4, DT.𝒢)
            ]
            @test DT.compare_triangle_collections(get_triangles(tri), T)
            @test validate_triangulation(tri, predicates=PT())

            DT.push_point!(tri, new_points[2])
            add_point!(tri, DT.num_points(tri), predicates=PT())
            @test sort_edge_vector(collect(get_interior_segments(tri))) ==
                  sort_edge_vector(collect(Set([(1, 2), (1, 15), (15, 16), (16, 12)]))) ==
                  sort_edge_vector(collect(get_interior_segments(tri)))
            T = [
                (16, 9, 15)
                (1, 10, 11)
                (14, 13, 16)
                (5, 1, 15)
                (1, 14, 15)
                (6, 10, 4)
                (6, 11, 10)
                (16, 12, 8)
                (12, 16, 13)
                (15, 9, 5)
                (10, 5, 4)
                (10, 1, 5)
                (14, 16, 15)
                (9, 4, 5)
                (14, 2, 13)
                (16, 8, 9)
                (8, 3, 9)
                (12, 13, 2)
                (12, 2, 7)
                (11, 6, 2)
                (8, 12, 3)
                (1, 2, 14)
                (12, 7, 3)
                (2, 6, 7)
                (2, 1, 11)
                (4, 9, DT.𝒢)
                (9, 3, DT.𝒢)
                (3, 7, DT.𝒢)
                (7, 6, DT.𝒢)
                (6, 4, DT.𝒢)
            ]
            @test DT.compare_triangle_collections(get_triangles(tri), T)
            @test validate_triangulation(tri, predicates=PT())

            DT.push_point!(tri, new_points[3])
            add_point!(tri, DT.num_points(tri), predicates=PT())
            @test sort_edge_vector(collect(get_interior_segments(tri))) ==
                  sort_edge_vector(collect(Set([(1, 2), (1, 15), (15, 16), (16, 17), (17, 12)]))) ==
                  sort_edge_vector(collect(get_interior_segments(tri)))
            T = [
                (13, 16, 14)
                (15, 1, 14)
                (16, 15, 14)
                (5, 1, 15)
                (10, 11, 1)
                (2, 11, 6)
                (16, 13, 17)
                (15, 9, 5)
                (6, 10, 4)
                (6, 11, 10)
                (10, 5, 4)
                (10, 1, 5)
                (5, 9, 4)
                (8, 3, 9)
                (16, 8, 9)
                (17, 2, 12)
                (17, 13, 2)
                (8, 16, 17)
                (9, 15, 16)
                (3, 12, 7)
                (12, 8, 17)
                (7, 12, 2)
                (3, 8, 12)
                (6, 7, 2)
                (14, 2, 13)
                (1, 2, 14)
                (2, 1, 11)
                (4, 9, DT.𝒢)
                (9, 3, DT.𝒢)
                (3, 7, DT.𝒢)
                (7, 6, DT.𝒢)
                (6, 4, DT.𝒢)
            ]
            @test DT.compare_triangle_collections(get_triangles(tri), T)
            @test validate_triangulation(tri, predicates=PT())

            DT.push_point!(tri, new_points[4])
            add_point!(tri, DT.num_points(tri), predicates=PT())
            @test sort_edge_vector(collect(get_interior_segments(tri))) ==
                  sort_edge_vector(collect(Set([(1, 2), (1, 15), (15, 16), (16, 17), (17, 18), (18, 12)]))) ==
                  sort_edge_vector(collect(get_interior_segments(tri)))
            T = [
                (13, 16, 14)
                (15, 1, 14)
                (9, 15, 16)
                (5, 1, 15)
                (10, 11, 1)
                (2, 11, 6)
                (16, 13, 17)
                (15, 9, 5)
                (6, 10, 4)
                (6, 11, 10)
                (10, 5, 4)
                (10, 1, 5)
                (8, 16, 17)
                (15, 14, 16)
                (5, 9, 4)
                (16, 8, 9)
                (17, 13, 2)
                (8, 3, 9)
                (18, 17, 2)
                (3, 12, 7)
                (12, 8, 18)
                (12, 18, 2)
                (8, 17, 18)
                (7, 12, 2)
                (3, 8, 12)
                (6, 7, 2)
                (14, 2, 13)
                (1, 2, 14)
                (2, 1, 11)
                (4, 9, DT.𝒢)
                (9, 3, DT.𝒢)
                (3, 7, DT.𝒢)
                (7, 6, DT.𝒢)
                (6, 4, DT.𝒢)
            ]
            @test DT.compare_triangle_collections(get_triangles(tri), T)
            @test validate_triangulation(tri, predicates=PT())

            add_segment!(tri, 6, 10, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            new_points = LinRange(-1.8, 2.8, 25)
            new_points = collect(new_points)
            shuffle!(new_points)
            foreach(new_points) do p
                DT.push_point!(tri, -3.0, p)
                add_point!(tri, DT.num_points(tri), predicates=PT())
                @test validate_triangulation(tri, predicates=PT())
            end
            @test length(get_interior_segments(tri)) == 32
            @test length(get_all_segments(tri)) == 32

            new_points = LinRange(0.357, 4.8912, 10)
            new_points = collect(new_points)
            shuffle!(new_points)
            foreach(new_points) do p
                DT.push_point!(tri, p, 3.0)
                add_point!(tri, DT.num_points(tri), predicates=PT())
                @test validate_triangulation(tri, predicates=PT())
            end
        end
    end
end

@testset "Adding a point onto a single boundary edge" begin
    for PT in subtypes(DT.AbstractPredicateType)
        for i in 1:8
            tri = triangulate_rectangle(0, 10, 0, 20, 11, 21; delete_ghosts=false, predicates=PT())
            add_point!(tri, 1.5, 0.0; initial_search_point=i, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            @test tri.boundary_nodes[1] == [1, 2, 232, 3, 4, 5, 6, 7, 8, 9, 10, 11]
            @test isempty(tri.interior_segments)
            @test (2, 232) ∈ get_all_segments(tri) || (232, 2) ∈ get_all_segments(tri)
            @test (232, 3) ∈ get_all_segments(tri) || (3, 232) ∈ get_all_segments(tri)
            @test (2, 3) ∉ get_all_segments(tri) && (3, 2) ∉ get_all_segments(tri)
            @test DT.get_boundary_edge_map(tri, 2, 232) == (1, 2)
            @test DT.get_boundary_edge_map(tri, 232, 3) == (1, 3)
            @test_throws KeyError DT.get_boundary_edge_map(tri, 2, 3)
        end
    end
end

@testset "Adding a point onto multiple boundary edges with multiple ghost indices" begin
    for PT in subtypes(DT.AbstractPredicateType)
        for _ in 1:8
            tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; delete_ghosts=false, predicates=PT())
            add_point!(tri, 1.5, 0.0, predicates=PT())
            add_segment!(tri, 1, 45, predicates=PT())
            add_point!(tri, 4.0, 2.5, predicates=PT())
            add_point!(tri, 4.0, 2.6, predicates=PT())
            add_point!(tri, 4.0, 7.3, predicates=PT())
            add_point!(tri, 2.5, 8.0, predicates=PT())
            add_point!(tri, 1.3, 8.0, predicates=PT())
            add_point!(tri, 0.0, 6.7, predicates=PT())
            add_point!(tri, 0.0, 2.5, predicates=PT())
            add_segment!(tri, 11, 43, predicates=PT())
            add_segment!(tri, 4, 34, predicates=PT())
            add_segment!(tri, 4, 23, predicates=PT())
            add_segment!(tri, 45, 27, predicates=PT())
            @test DT.get_boundary_nodes(tri) == [
                [1, 2, 46, 3, 4, 5],
                [5, 10, 15, 47, 48, 20, 25, 30, 35, 40, 49, 45],
                [45, 44, 50, 43, 51, 42, 41],
                [41, 36, 52, 31, 26, 21, 16, 53, 11, 6, 1]
            ]
            @test validate_triangulation(tri, predicates=PT())
        end
    end
end

@testset "Handling only a single boundary index" begin
    for PT in subtypes(DT.AbstractPredicateType)
        for _ in 1:8
            tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; delete_ghosts=false, single_boundary=true, predicates=PT())
            for x in [0.2, 0.3, 1.5, 2.3, 2.8, 3.5]
                add_point!(tri, x, 0.0, predicates=PT())
            end
            for y in [0.2, 0.4, 1.6, 4.5, 6.7, 7.5]
                add_point!(tri, 4.0, y, predicates=PT())
            end
            for x in [0.4, 1.5, 1.8, 2.3, 2.5, 3.9]
                add_point!(tri, x, 8.0, predicates=PT())
            end
            for y in [7.5, 5.9, 3.4, 1.9, 0.1]
                add_point!(tri, 0.0, y, predicates=PT())
            end
            @test DT.get_boundary_nodes(tri) == [1, 46, 47, 2, 48, 3, 49, 50, 4, 51, 5, 52, 53,
                10, 54, 15, 20, 25, 55, 30, 35, 56, 40, 57, 45, 63, 44, 62, 61, 43, 60, 59, 42, 58,
                41, 64, 36, 31, 65, 26, 21, 66, 16, 11, 67, 6, 68, 1]
            @test validate_triangulation(tri, predicates=PT())
        end
    end
end

@testset "Starting with a thin set of boundary nodes, and filling them in with automatic collinearity detection" begin
    @testset "Contiguous boundary" begin
        for PT in subtypes(DT.AbstractPredicateType)
            for _ in 1:8
                tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; delete_ghosts=false, predicates=PT())
                pts = get_points(tri)
                boundary_nodes = [1, 5, 45, 41, 1]
                tri = triangulate(pts; boundary_nodes, randomise=false, delete_ghosts=false, predicates=PT())
                flip_edge!(tri, 1, 7)
                flip_edge!(tri, 2, 8)
                @test sort_edge_vector(collect(get_all_segments(tri))) == sort_edge_vector([
                    (1, 2)
                    (2, 3)
                    (3, 4)
                    (4, 5)
                    (5, 10)
                    (6, 1)
                    (10, 15)
                    (11, 6)
                    (15, 20)
                    (16, 11)
                    (20, 25)
                    (21, 16)
                    (25, 30)
                    (26, 21)
                    (30, 35)
                    (31, 26)
                    (35, 40)
                    (36, 31)
                    (40, 45)
                    (41, 36)
                    (42, 41)
                    (43, 42)
                    (44, 43)
                    (45, 44)
                ])
                T = [
                    (11, 7, 12)
                    (22, 27, 26)
                    (8, 7, 3)
                    (32, 27, 28)
                    (19, 18, 14)
                    (33, 38, 37)
                    (33, 28, 29)
                    (13, 14, 18)
                    (5, 9, 4)
                    (8, 13, 12)
                    (10, 14, 9)
                    (24, 29, 28)
                    (1, 2, 6)
                    (12, 7, 8)
                    (6, 2, 7)
                    (9, 8, 4)
                    (7, 2, 3)
                    (17, 16, 12)
                    (11, 6, 7)
                    (4, 8, 3)
                    (9, 14, 13)
                    (16, 11, 12)
                    (13, 17, 12)
                    (22, 26, 21)
                    (27, 31, 26)
                    (32, 33, 37)
                    (37, 41, 36)
                    (31, 32, 36)
                    (38, 42, 37)
                    (36, 32, 37)
                    (31, 27, 32)
                    (37, 42, 41)
                    (38, 34, 39)
                    (23, 22, 18)
                    (21, 16, 17)
                    (17, 22, 21)
                    (17, 18, 22)
                    (32, 28, 33)
                    (27, 22, 23)
                    (34, 29, 30)
                    (23, 18, 19)
                    (33, 34, 38)
                    (33, 29, 34)
                    (38, 43, 42)
                    (38, 39, 43)
                    (34, 35, 39)
                    (39, 44, 43)
                    (44, 39, 40)
                    (35, 34, 30)
                    (44, 40, 45)
                    (39, 35, 40)
                    (19, 24, 23)
                    (25, 29, 24)
                    (30, 29, 25)
                    (28, 23, 24)
                    (28, 27, 23)
                    (10, 15, 14)
                    (20, 24, 19)
                    (17, 13, 18)
                    (8, 9, 13)
                    (15, 19, 14)
                    (10, 9, 5)
                    (20, 19, 15)
                    (25, 24, 20)
                    (1, 6, DT.𝒢)
                    (6, 11, DT.𝒢)
                    (11, 16, DT.𝒢)
                    (16, 21, DT.𝒢)
                    (21, 26, DT.𝒢)
                    (26, 31, DT.𝒢)
                    (31, 36, DT.𝒢)
                    (36, 41, DT.𝒢)
                    (41, 42, DT.𝒢)
                    (42, 43, DT.𝒢)
                    (43, 44, DT.𝒢)
                    (44, 45, DT.𝒢)
                    (45, 40, DT.𝒢)
                    (40, 35, DT.𝒢)
                    (35, 30, DT.𝒢)
                    (30, 25, DT.𝒢)
                    (25, 20, DT.𝒢)
                    (20, 15, DT.𝒢)
                    (15, 10, DT.𝒢)
                    (10, 5, DT.𝒢)
                    (5, 4, DT.𝒢)
                    (4, 3, DT.𝒢)
                    (3, 2, DT.𝒢)
                    (2, 1, DT.𝒢)
                ]
                @test DT.compare_triangle_collections(get_triangles(tri), T)
                @test get_boundary_nodes(tri) == [1, 2, 3, 4, 5,
                    10, 15, 20, 25, 30, 35, 40, 45,
                    44, 43, 42, 41, 36, 31, 26,
                    21, 16, 11, 6, 1]
                @test validate_triangulation(tri, predicates=PT())
                add_segment!(tri, 1, 45, predicates=PT())
                add_segment!(tri, 6, 44, predicates=PT())
                add_segment!(tri, 34, 4, predicates=PT())
                @test validate_triangulation(tri, predicates=PT())
            end
        end
    end

    @testset "Multiple segments" begin
        for PT in subtypes(DT.AbstractPredicateType)
            for _ in 1:8
                tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; delete_ghosts=false, predicates=PT())
                pts = get_points(tri)
                boundary_nodes = [[1, 5], [5, 45], [45, 41], [41, 1]]
                tri = triangulate(pts; boundary_nodes, randomise=false, delete_ghosts=false, predicates=PT())
                flip_edge!(tri, 1, 7)
                flip_edge!(tri, 2, 8)
                @test sort_edge_vector(collect(get_all_segments(tri))) == sort_edge_vector([
                    (1, 2)
                    (2, 3)
                    (3, 4)
                    (4, 5)
                    (5, 10)
                    (6, 1)
                    (10, 15)
                    (11, 6)
                    (15, 20)
                    (16, 11)
                    (20, 25)
                    (21, 16)
                    (25, 30)
                    (26, 21)
                    (30, 35)
                    (31, 26)
                    (35, 40)
                    (36, 31)
                    (40, 45)
                    (41, 36)
                    (42, 41)
                    (43, 42)
                    (44, 43)
                    (45, 44)
                ])
                T = [
                    (11, 7, 12)
                    (22, 27, 26)
                    (8, 7, 3)
                    (32, 27, 28)
                    (19, 18, 14)
                    (33, 38, 37)
                    (33, 28, 29)
                    (13, 14, 18)
                    (5, 9, 4)
                    (8, 13, 12)
                    (10, 14, 9)
                    (24, 29, 28)
                    (1, 2, 6)
                    (12, 7, 8)
                    (6, 2, 7)
                    (9, 8, 4)
                    (7, 2, 3)
                    (17, 16, 12)
                    (11, 6, 7)
                    (4, 8, 3)
                    (9, 14, 13)
                    (16, 11, 12)
                    (13, 17, 12)
                    (22, 26, 21)
                    (27, 31, 26)
                    (32, 33, 37)
                    (37, 41, 36)
                    (31, 32, 36)
                    (38, 42, 37)
                    (36, 32, 37)
                    (31, 27, 32)
                    (37, 42, 41)
                    (38, 34, 39)
                    (23, 22, 18)
                    (21, 16, 17)
                    (17, 22, 21)
                    (17, 18, 22)
                    (32, 28, 33)
                    (27, 22, 23)
                    (34, 29, 30)
                    (23, 18, 19)
                    (33, 34, 38)
                    (33, 29, 34)
                    (38, 43, 42)
                    (38, 39, 43)
                    (34, 35, 39)
                    (39, 44, 43)
                    (44, 39, 40)
                    (35, 34, 30)
                    (44, 40, 45)
                    (39, 35, 40)
                    (19, 24, 23)
                    (25, 29, 24)
                    (30, 29, 25)
                    (28, 23, 24)
                    (28, 27, 23)
                    (10, 15, 14)
                    (20, 24, 19)
                    (17, 13, 18)
                    (8, 9, 13)
                    (15, 19, 14)
                    (10, 9, 5)
                    (20, 19, 15)
                    (25, 24, 20)
                    (1, 6, DT.𝒢 - 3)
                    (6, 11, DT.𝒢 - 3)
                    (11, 16, DT.𝒢 - 3)
                    (16, 21, DT.𝒢 - 3)
                    (21, 26, DT.𝒢 - 3)
                    (26, 31, DT.𝒢 - 3)
                    (31, 36, DT.𝒢 - 3)
                    (36, 41, DT.𝒢 - 3)
                    (41, 42, DT.𝒢 - 2)
                    (42, 43, DT.𝒢 - 2)
                    (43, 44, DT.𝒢 - 2)
                    (44, 45, DT.𝒢 - 2)
                    (45, 40, DT.𝒢 - 1)
                    (40, 35, DT.𝒢 - 1)
                    (35, 30, DT.𝒢 - 1)
                    (30, 25, DT.𝒢 - 1)
                    (25, 20, DT.𝒢 - 1)
                    (20, 15, DT.𝒢 - 1)
                    (15, 10, DT.𝒢 - 1)
                    (10, 5, DT.𝒢 - 1)
                    (5, 4, DT.𝒢)
                    (4, 3, DT.𝒢)
                    (3, 2, DT.𝒢)
                    (2, 1, DT.𝒢)
                ]
                @test DT.compare_triangle_collections(get_triangles(tri), T)
                @test get_boundary_nodes(tri) == [[1, 2, 3, 4, 5],
                    [5, 10, 15, 20, 25, 30, 35, 40, 45],
                    [45, 44, 43, 42, 41],
                    [41, 36, 31, 26, 21, 16, 11, 6, 1]]
                @test validate_triangulation(tri, predicates=PT())
                add_segment!(tri, 1, 45, predicates=PT())
                add_segment!(tri, 6, 44, predicates=PT())
                add_segment!(tri, 34, 4, predicates=PT())
                @test validate_triangulation(tri, predicates=PT())
            end
        end
    end
end

@testset "Adding points and segments into a multiply-connected domain" begin
    for PT in subtypes(DT.AbstractPredicateType)
        tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; delete_ghosts=false, predicates=PT())
        boundary_nodes = [[
                [1, 5, 45], [45, 41], [41, 1]
            ],
            [[12, 32, 34], [34, 14, 12]]
        ]
        points = get_points(tri)
        rng = StableRNG(19119)
        tri = triangulate(points; boundary_nodes, delete_ghosts=false, randomise=false, rng, predicates=PT())
        flip_edge!(tri, 1, 7)
        flip_edge!(tri, 2, 8)
        T = [
            (11, 7, 12)
            (22, 27, 26)
            (8, 7, 3)
            (33, 38, 37)
            (5, 9, 4)
            (8, 13, 12)
            (10, 14, 9)
            (1, 2, 6)
            (12, 7, 8)
            (6, 2, 7)
            (9, 8, 4)
            (7, 2, 3)
            (17, 16, 12)
            (11, 6, 7)
            (4, 8, 3)
            (9, 14, 13)
            (16, 11, 12)
            (22, 26, 21)
            (27, 31, 26)
            (32, 33, 37)
            (37, 41, 36)
            (31, 32, 36)
            (38, 42, 37)
            (36, 32, 37)
            (31, 27, 32)
            (37, 42, 41)
            (38, 34, 39)
            (21, 16, 17)
            (17, 22, 21)
            (34, 29, 30)
            (33, 34, 38)
            (38, 43, 42)
            (38, 39, 43)
            (34, 35, 39)
            (39, 44, 43)
            (44, 39, 40)
            (35, 34, 30)
            (44, 40, 45)
            (39, 35, 40)
            (25, 29, 24)
            (30, 29, 25)
            (10, 15, 14)
            (20, 24, 19)
            (8, 9, 13)
            (15, 19, 14)
            (10, 9, 5)
            (20, 19, 15)
            (25, 24, 20)
            (2, 1, DT.𝒢)
            (3, 2, DT.𝒢)
            (4, 3, DT.𝒢)
            (5, 4, DT.𝒢)
            (10, 5, DT.𝒢)
            (15, 10, DT.𝒢)
            (20, 15, DT.𝒢)
            (25, 20, DT.𝒢)
            (30, 25, DT.𝒢)
            (35, 30, DT.𝒢)
            (40, 35, DT.𝒢)
            (45, 40, DT.𝒢)
            (44, 45, DT.𝒢 - 1)
            (43, 44, DT.𝒢 - 1)
            (42, 43, DT.𝒢 - 1)
            (41, 42, DT.𝒢 - 1)
            (36, 41, DT.𝒢 - 2)
            (31, 36, DT.𝒢 - 2)
            (26, 31, DT.𝒢 - 2)
            (21, 26, DT.𝒢 - 2)
            (16, 21, DT.𝒢 - 2)
            (11, 16, DT.𝒢 - 2)
            (6, 11, DT.𝒢 - 2)
            (1, 6, DT.𝒢 - 2)
            (17, 12, DT.𝒢 - 3)
            (22, 17, DT.𝒢 - 3)
            (27, 22, DT.𝒢 - 3)
            (32, 27, DT.𝒢 - 3)
            (33, 32, DT.𝒢 - 3)
            (34, 33, DT.𝒢 - 3)
            (29, 34, DT.𝒢 - 4)
            (24, 29, DT.𝒢 - 4)
            (19, 24, DT.𝒢 - 4)
            (14, 19, DT.𝒢 - 4)
            (13, 14, DT.𝒢 - 4)
            (12, 13, DT.𝒢 - 4)
        ]
        C = [
            (1, 2)
            (2, 3)
            (3, 4)
            (4, 5)
            (5, 10)
            (6, 1)
            (10, 15)
            (11, 6)
            (12, 17)
            (13, 12)
            (14, 13)
            (15, 20)
            (16, 11)
            (17, 22)
            (19, 14)
            (20, 25)
            (21, 16)
            (22, 27)
            (24, 19)
            (25, 30)
            (26, 21)
            (27, 32)
            (29, 24)
            (30, 35)
            (31, 26)
            (32, 33)
            (33, 34)
            (34, 29)
            (35, 40)
            (36, 31)
            (40, 45)
            (41, 36)
            (42, 41)
            (43, 42)
            (44, 43)
            (45, 44)
        ]
        @test sort_edge_vector(collect(get_all_segments(tri))) == sort_edge_vector(C)
        @test DT.compare_triangle_collections(get_triangles(tri), T)
        @test validate_triangulation(tri, predicates=PT())
        @test get_boundary_nodes(tri) == [[
                [1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45], [45, 44, 43, 42, 41], [41, 36, 31, 26, 21, 16, 11, 6, 1]
            ],
            [[12, 17, 22, 27, 32, 33, 34], [34, 29, 24, 19, 14, 13, 12]]
        ]
        add_point!(tri, 1.5, 0.0; rng, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        add_point!(tri, 2.5, 0.0; rng, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        add_point!(tri, 4.0, 2.5; rng, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        add_point!(tri, 4.0, 7.5; rng, predicates=PT())
        add_point!(tri, 2.5, 8.0; rng, predicates=PT())
        add_point!(tri, 0.0, 5.5; rng, predicates=PT())
        add_point!(tri, 0.5, 2.2; rng, predicates=PT())
        add_segment!(tri, 21, 27; rng, predicates=PT())
        add_segment!(tri, 14, 45; rng, predicates=PT())
        add_point!(tri, 1.0, 2.5; rng, predicates=PT())
        add_point!(tri, 1.0, 5.5; rng, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        add_point!(tri, 2.5, 2.0; rng, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        add_point!(tri, 3.0, 5.5; rng, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        add_point!(tri, 2.5, 6.0; rng, predicates=PT())
        add_segment!(tri, 12, 5, predicates=PT())
        add_segment!(tri, 2, 52, predicates=PT())
        add_segment!(tri, 53, 21, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())
        @test get_boundary_nodes(tri) == [[
                [1, 2, 46, 3, 47, 4, 5, 10, 15, 48, 20, 25, 30, 35, 40, 49, 45],
                [45, 44, 50, 43, 42, 41],
                [41, 36, 31, 51, 26, 21, 16, 11, 6, 1]
            ],
            [
                [12, 53, 17, 22, 27, 54, 32, 33, 57, 34],
                [34, 56, 29, 24, 19, 14, 55, 13, 12]
            ]
        ]
        @test validate_triangulation(tri, predicates=PT())
        add_point!(tri, 0.5, 4.4, predicates=PT())
        add_point!(tri, 0.5, 4.6, predicates=PT())
        add_point!(tri, 0.5, 4.5, predicates=PT())
        @test validate_triangulation(tri, predicates=PT())

        for L in 1:3
            tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; delete_ghosts=false, predicates=PT())
            boundary_nodes = [[
                    [1, 5, 45], [45, 41], [41, 1]
                ],
                [[12, 32, 34], [34, 14, 12]]
            ]
            points = get_points(tri)
            tri = triangulate(points; boundary_nodes, delete_ghosts=false, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            @test get_boundary_nodes(tri) == [[
                    [1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45], [45, 44, 43, 42, 41], [41, 36, 31, 26, 21, 16, 11, 6, 1]
                ],
                [[12, 17, 22, 27, 32, 33, 34], [34, 29, 24, 19, 14, 13, 12]]
            ]
            add_point!(tri, 1.5, 0.0, predicates=PT())
            add_point!(tri, 2.5, 0.0, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            add_point!(tri, 4.0, 2.5, predicates=PT())
            add_point!(tri, 4.0, 7.5, predicates=PT())
            add_point!(tri, 2.5, 8.0, predicates=PT())
            add_point!(tri, 0.0, 5.5, predicates=PT())
            add_point!(tri, 0.5, 2.2, predicates=PT())
            add_segment!(tri, 21, 27, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            add_segment!(tri, 14, 45, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            add_point!(tri, 1.0, 2.5, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            add_point!(tri, 1.0, 5.5, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            add_point!(tri, 2.5, 2.0, predicates=PT())
            add_point!(tri, 3.0, 5.5, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            add_point!(tri, 2.5, 6.0, predicates=PT())
            add_segment!(tri, 12, 5, predicates=PT())
            add_segment!(tri, 2, 52, predicates=PT())
            add_segment!(tri, 53, 21, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
            @test get_boundary_nodes(tri) == [[
                    [1, 2, 46, 3, 47, 4, 5, 10, 15, 48, 20, 25, 30, 35, 40, 49, 45],
                    [45, 44, 50, 43, 42, 41],
                    [41, 36, 31, 51, 26, 21, 16, 11, 6, 1]
                ],
                [
                    [12, 53, 17, 22, 27, 54, 32, 33, 57, 34],
                    [34, 56, 29, 24, 19, 14, 55, 13, 12]
                ]
            ]
            @test validate_triangulation(tri, predicates=PT())
            add_point!(tri, 0.5, 4.4, predicates=PT())
            add_point!(tri, 0.5, 4.6, predicates=PT())
            add_point!(tri, 0.5, 4.5, predicates=PT())
            @test validate_triangulation(tri, predicates=PT())
        end
    end
end