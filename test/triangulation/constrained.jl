using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using StableRNGs
include("../helper_functions.jl")


@testset "Test random constrained Delaunay triangulations" begin
    rng = StableRNG(191919)
    for i in 1:250
        @show i
        points, edges, mat_edges = get_random_vertices_and_constrained_edges(40, 200, 20, rng)
        tri = triangulate(points; edges, rng)
        @test validate_triangulation(tri)
        empty!(get_all_constrained_edges(tri))
        @test !validate_triangulation(tri)
    end
    for i in 1:25
        @show i
        points, edges, mat_edges = get_random_vertices_and_constrained_edges(200, 1000, 200, rng)
        tri = triangulate(points; edges, rng)
        @test validate_triangulation(tri)
        empty!(get_all_constrained_edges(tri))
        @test !validate_triangulation(tri)
    end
end

@testset "Testing Shewchuk's PSLG example" begin
    pts, C = second_shewchuk_example_constrained()
    for i in 1:500
        @show i
        rng = StableRNG(i^6)
        tri = triangulate(pts; edges=C, rng)
        @test validate_triangulation(tri)
    end
end

@testset "Random parabolas" begin
    for i in 1:50
        @show i
        rng = StableRNG(i)
        pts = [(2rand(rng) - 1, rand(rng)) for _ in 1:500]
        x = LinRange(-1, 1, 250)
        a = LinRange(0.0, 1.0, 8)
        C = Set{NTuple{2,Int64}}()
        for i in eachindex(a)
            y = a[i] * x .^ 2
            append!(pts, zip(x, y))
            push!(C, [(i, i + 1) for i in (500+250*(i-1)+1):(500+250*(i-1)+249)]...)
        end
        tri = triangulate(pts; edges=C, rng)
        @test validate_triangulation(tri)
    end
end

@testset "Random collection of straight lines" begin
    for i in 1:10
        @show i
        rng = StableRNG(i)
        pts = NTuple{2,Float64}[]
        C = Set{NTuple{2,Int64}}()
        j = 1
        for i in 1:100
            push!(pts, (2i / 101 - 1, 2rand(rng) - 1))
            push!(pts, (2i / 101 - 1, 2rand(rng) - 1))
            push!(C, (j, j + 1))
            j += 2
        end
        x1 = LinRange(-1, 1 - 1e-12, 100)
        y1 = LinRange(-1, -1, 100)
        x2 = LinRange(1, 1, 100)
        y2 = LinRange(-1, 1 - 1e-12, 100)
        x3 = LinRange(1, -1 + 1e-12, 100)
        y3 = LinRange(1, 1, 100)
        x4 = LinRange(-1, -1, 100)
        y4 = LinRange(1, -1 + 1e-12, 100)
        append!(pts, zip(x1, y1), zip(x2, y2), zip(x3, y3), zip(x4, y4))
        push!(C, [(j, j + 1) for j in 201:299]...)
        push!(C, [(j, j + 1) for j in 301:399]...)
        push!(C, [(j, j + 1) for j in 401:499]...)
        push!(C, [(j, j + 1) for j in 501:599]...)
        tri = triangulate(pts; edges=C, rng)
        @test validate_triangulation(tri)
    end
end

@testset "Lattice" begin
    for m in 1:20
        rng = StableRNG(m)
        @show m
        a = 0.0
        b = 5.0
        c = -3.0
        d = 7.0
        nx = 13
        ny = 20
        tri = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
        add_edge!(tri, 56, 162; rng)
        for e in [(1, 249), (1, 250), (1, 251), (1, 26), (1, 39), (1, 52)]
            add_edge!(tri, e; rng)
        end
        add_edge!(tri, 190, 99; rng)
        for e in [(99, 113), (113, 101), (101, 115), (115, 127), (127, 141), (141, 177)]
            add_edge!(tri, e; rng)
        end
        @test validate_triangulation(tri)

        a = 0.0
        b = 1.0
        c = 0.0
        d = 1.0
        nx = 2
        ny = 2
        tri = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
        add_edge!(tri, 1, 4; rng)
        @test validate_triangulation(tri)

        a = -0.1
        b = 0.1
        c = -0.01
        d = 0.01
        nx = 25
        ny = 25
        tri = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
        tri = triangulate(get_points(tri))
        for i in 2:24
            add_edge!(tri, i, 600 + i; rng)
        end
        @test validate_triangulation(tri)
        tri = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
        tri = triangulate(get_points(tri); rng)
        for e in [(1, 28), (28, 103), (103, 180), (180, 625), (625, 523), (523, 23), (23, 71), (71, 60), (60, 28)]
            add_edge!(tri, e; rng)
        end
        for e in [(402, 227), (227, 430), (430, 437), (437, 614), (527, 602), (528, 603), (555, 605)]
            add_edge!(tri, e; rng)
        end
        @test validate_triangulation(tri)

        a = 0
        b = 1
        c = 0
        d = 5
        nx = 25
        ny = 3
        tri = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
        tri = triangulate(get_points(tri); rng)
        for i in 1:(nx-1)
            u = i
            v = 2nx
            add_edge!(tri, u, v)
        end
        for i in 51:75
            u = i
            v = 26
            add_edge!(tri, u, v; rng)
        end
        @test validate_triangulation(tri)
    end
end

@testset "Triangulating with a deleted exterior" begin
    for i in 1:500
        @show i
        rng = StableRNG(i)
        pts = [(rand(rng), rand(rng)) for _ in 1:50]
        bnd_pts = [(0.3cos(θ), 0.3sin(θ)) .+ 0.5 for θ in LinRange(0, 2π - 1 / 250, 25)]
        bnd_id = [(51:75)..., 51]
        append!(pts, bnd_pts)
        tri = triangulate(pts; boundary_nodes=bnd_id, rng)
        @test validate_triangulation(tri)
    end
end

@testset "Triangulation with two curves" begin
    for i in 1:50
        @show i
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
        tri = triangulate(pts; boundary_nodes=nodes, rng)
        @test validate_triangulation(tri)
    end
end

if !get(ENV, "CI", false)
    @testset "Triangulating more regions with holes, and non-convex (Gmsh examples)" begin
        _x, _y = complicated_geometry()
        x = _x
        y = _y
        tri_1 = generate_mesh(x, y, 0.1; convert_result=true, add_ghost_triangles=true)
        tri_2 = generate_mesh(x[1], y[1], 0.1; convert_result=true, add_ghost_triangles=true)
        tri_3 = generate_mesh([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0], 0.1;
            convert_result=true, add_ghost_triangles=true)
        tri_4 = generate_mesh(reverse(reverse.(x[2])), reverse(reverse.(y[2])), 0.1; convert_result=true, add_ghost_triangles=true)
        a, b = 0.0, 5.0
        c, d = 3.0, 7.0
        nx = 3
        ny = 3
        tri_5 = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=false)
        tri_6 = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)

        for tri in (tri_1, tri_2, tri_3, tri_4, tri_5, tri_6)
            points = get_points(tri)
            bn_nodes = get_boundary_nodes(tri)
            _tri = triangulate(points; boundary_nodes=bn_nodes, delete_ghosts=false, delete_holes=true)
            @test validate_triangulation(_tri)
            tri ∉ (tri_5, tri_6) && @test DT.compare_triangle_collections(get_triangles(tri), get_triangles(_tri))
        end

        a = 4 / 5
        t = LinRange(0, 2π, 100)
        x = @. a * (2cos(t) + cos(2t))
        y = @. a * (2sin(t) - sin(2t))
        tri = generate_mesh(x, y, 15.0)
        points = get_points(tri)
        bn_nodes = get_boundary_nodes(tri)
        _tri = triangulate(points; boundary_nodes=bn_nodes, delete_ghosts=false)
        @test validate_triangulation(_tri; check_planarity=false, check_ghost_triangle_delaunay=false, check_ghost_triangle_orientation=false)
        validate_statistics(_tri)

        x1 = [collect(LinRange(0, 2, 4)),
            collect(LinRange(2, 2, 4)),
            collect(LinRange(2, 0, 4)),
            collect(LinRange(0, 0, 4))]
        y1 = [collect(LinRange(0, 0, 4)),
            collect(LinRange(0, 6, 4)),
            collect(LinRange(6, 6, 4)),
            collect(LinRange(6, 0, 4))]
        r = 0.5
        h = k = 0.6
        θ = LinRange(2π, 0, 50)
        x2 = [h .+ r .* cos.(θ)]
        y2 = [k .+ r .* sin.(θ)]
        r = 0.2
        h = 1.5
        k = 0.5
        x3 = [h .+ r .* cos.(θ)]
        y3 = [k .+ r .* sin.(θ)]
        x4 = reverse(reverse.([collect(LinRange(1, 1.5, 4)),
            collect(LinRange(1.5, 1.5, 4)),
            collect(LinRange(1.5, 1, 4)),
            collect(LinRange(1, 1, 4))]))
        y4 = reverse(reverse.([collect(LinRange(2, 2, 4)),
            collect(LinRange(2, 5, 4)),
            collect(LinRange(5, 5, 4)),
            collect(LinRange(5, 2, 4))]))
        x5 = [reverse([0.2, 0.5, 0.75, 0.75, 0.2, 0.2])]
        y5 = [reverse([2.0, 2.0, 3.0, 4.0, 5.0, 2.0])]
        x = [x1, x2, x3, x4, x5]
        y = [y1, y2, y3, y4, y5]
        tri = generate_mesh(x, y, 0.2)
        points = get_points(tri)
        bn_nodes = get_boundary_nodes(tri)
        _tri = triangulate(points; boundary_nodes=bn_nodes, delete_ghosts=false)
        @test validate_triangulation(_tri; check_planarity=false, check_ghost_triangle_delaunay=false, check_ghost_triangle_orientation=false)
        validate_statistics(_tri)
    end
end

@testset "Adding points into a constrained triangulation; no collinearities" begin
    for L in 1:10
        pts, C = example_for_testing_add_point_on_constrained_triangulation()
        tri = triangulate(pts; edges=C, delete_ghosts=false)
        @test validate_triangulation(tri)
        DT.push_point!(tri, 2, 1.8)
        add_point!(tri, 15)
        DT.push_point!(tri, 1.57, 1.778)
        add_point!(tri, 16)
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
            (4, 9, DT.BoundaryIndex)
            (9, 3, DT.BoundaryIndex)
            (3, 7, DT.BoundaryIndex)
            (7, 6, DT.BoundaryIndex)
            (6, 4, DT.BoundaryIndex)
        ]
        @test DT.compare_triangle_collections(get_triangles(tri), T)
        @test validate_triangulation(tri)
        (x, y) = (2.3258217552204066, 1.4540267924574883) # a while loop was used to find this point that broke the triangulation (fixed now, obviously)
        DT.push_point!(tri, x, y)
        add_point!(tri, num_points(tri))
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
            (4, 9, DT.BoundaryIndex)
            (9, 3, DT.BoundaryIndex)
            (3, 7, DT.BoundaryIndex)
            (7, 6, DT.BoundaryIndex)
            (6, 4, DT.BoundaryIndex)
        ]
        @test DT.compare_triangle_collections(get_triangles(tri), T)
        @test validate_triangulation(tri)
        for i in 1:250
            @show i, L
            x = 4rand()
            y = 5rand()
            DT.push_point!(tri, x, y)
            add_point!(tri, num_points(tri))
            @test validate_triangulation(tri)
            @test DT.edge_exists(tri, 1, 2) && DT.edge_exists(tri, 2, 1)
        end
    end
end

@testset "Adding points into a constrained triangulation; interior segment collinearities" begin
    for m in 1:10
        @show m
        pts, C = example_for_testing_add_point_on_constrained_triangulation()
        push!(C, (1, 12))
        tri = triangulate(pts; edges=C, delete_ghosts=false)
        new_points = [
            (1.0, 3.0),
            (2.0, 3.0),
            (3.0, 3.0),
            (4.0, 3.0)
        ]

        DT.push_point!(tri, new_points[1])
        add_point!(tri, num_points(tri))
        @test sort_edge_vector(collect(get_constrained_edges(tri))) ==
              sort_edge_vector(collect(Set([(1, 2), (1, 15), (15, 12)]))) ==
              sort_edge_vector(collect(get_constrained_edges(tri)))
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
            (4, 9, DT.BoundaryIndex)
            (9, 3, DT.BoundaryIndex)
            (3, 7, DT.BoundaryIndex)
            (7, 6, DT.BoundaryIndex)
            (6, 4, DT.BoundaryIndex)
        ]
        @test DT.compare_triangle_collections(get_triangles(tri), T)
        @test validate_triangulation(tri)

        DT.push_point!(tri, new_points[2])
        add_point!(tri, num_points(tri))
        @test sort_edge_vector(collect(get_constrained_edges(tri))) ==
              sort_edge_vector(collect(Set([(1, 2), (1, 15), (15, 16), (16, 12)]))) ==
              sort_edge_vector(collect(get_constrained_edges(tri)))
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
            (4, 9, DT.BoundaryIndex)
            (9, 3, DT.BoundaryIndex)
            (3, 7, DT.BoundaryIndex)
            (7, 6, DT.BoundaryIndex)
            (6, 4, DT.BoundaryIndex)
        ]
        @test DT.compare_triangle_collections(get_triangles(tri), T)
        @test validate_triangulation(tri)

        DT.push_point!(tri, new_points[3])
        add_point!(tri, num_points(tri))
        @test sort_edge_vector(collect(get_constrained_edges(tri))) ==
              sort_edge_vector(collect(Set([(1, 2), (1, 15), (15, 16), (16, 17), (17, 12)]))) ==
              sort_edge_vector(collect(get_constrained_edges(tri)))
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
            (4, 9, DT.BoundaryIndex)
            (9, 3, DT.BoundaryIndex)
            (3, 7, DT.BoundaryIndex)
            (7, 6, DT.BoundaryIndex)
            (6, 4, DT.BoundaryIndex)
        ]
        @test DT.compare_triangle_collections(get_triangles(tri), T)
        @test validate_triangulation(tri)

        DT.push_point!(tri, new_points[4])
        add_point!(tri, num_points(tri))
        @test sort_edge_vector(collect(get_constrained_edges(tri))) ==
              sort_edge_vector(collect(Set([(1, 2), (1, 15), (15, 16), (16, 17), (17, 18), (18, 12)]))) ==
              sort_edge_vector(collect(get_constrained_edges(tri)))
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
            (4, 9, DT.BoundaryIndex)
            (9, 3, DT.BoundaryIndex)
            (3, 7, DT.BoundaryIndex)
            (7, 6, DT.BoundaryIndex)
            (6, 4, DT.BoundaryIndex)
        ]
        @test DT.compare_triangle_collections(get_triangles(tri), T)
        @test validate_triangulation(tri)

        add_edge!(tri, 6, 10)
        @test validate_triangulation(tri)
        new_points = LinRange(-1.8, 2.8, 250)
        new_points = collect(new_points)
        shuffle!(new_points)
        foreach(new_points) do p
            DT.push_point!(tri, -3.0, p)
            add_point!(tri, num_points(tri))
            @test validate_triangulation(tri)
        end
        @test length(get_constrained_edges(tri)) == 257
        @test length(get_all_constrained_edges(tri)) == 257

        new_points = LinRange(0.357, 4.8912, 40)
        new_points = collect(new_points)
        shuffle!(new_points)
        foreach(new_points) do p
            DT.push_point!(tri, p, 3.0)
            add_point!(tri, num_points(tri))
            @test validate_triangulation(tri)
        end
    end
end

@testset "Adding a point onto a single boundary edge" begin
    for i in 1:(11*21)
        @show i
        tri = triangulate_rectangle(0, 10, 0, 20, 11, 21; add_ghost_triangles=true)
        add_point!(tri, 1.5, 0.0; initial_search_point=i)
        @test validate_triangulation(tri)
        @test tri.boundary_nodes[1] == [1, 2, 232, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        @test isempty(tri.constrained_edges)
        @test (2, 232) ∈ get_all_constrained_edges(tri) || (232, 2) ∈ get_all_constrained_edges(tri)
        @test (232, 3) ∈ get_all_constrained_edges(tri) || (3, 232) ∈ get_all_constrained_edges(tri)
        @test (2, 3) ∉ get_all_constrained_edges(tri) && (3, 2) ∉ get_all_constrained_edges(tri)
        @test DT.get_boundary_edge_map(tri, 2, 232) == (1, 2)
        @test DT.get_boundary_edge_map(tri, 232, 3) == (1, 3)
        @test_throws KeyError DT.get_boundary_edge_map(tri, 2, 3)
    end
end

@testset "Adding a point onto multiple boundary edges with multiple ghost indices" begin
    for _ in 1:100
        tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; add_ghost_triangles=true)
        add_point!(tri, 1.5, 0.0)
        add_edge!(tri, 1, 45)
        add_point!(tri, 4.0, 2.5)
        add_point!(tri, 4.0, 2.6)
        add_point!(tri, 4.0, 7.3)
        add_point!(tri, 2.5, 8.0)
        add_point!(tri, 1.3, 8.0)
        add_point!(tri, 0.0, 6.7)
        add_point!(tri, 0.0, 2.5)
        add_edge!(tri, 11, 43)
        add_edge!(tri, 4, 34)
        add_edge!(tri, 4, 23)
        add_edge!(tri, 45, 27)
        @test DT.get_boundary_nodes(tri) == [
            [1, 2, 46, 3, 4, 5],
            [5, 10, 15, 47, 48, 20, 25, 30, 35, 40, 49, 45],
            [45, 44, 50, 43, 51, 42, 41],
            [41, 36, 52, 31, 26, 21, 16, 53, 11, 6, 1]
        ]
        @test validate_triangulation(tri)
    end
end

@testset "Handling only a single boundary index" begin
    for _ in 1:100
        tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; add_ghost_triangles=true, single_boundary=true)
        for x in [0.2, 0.3, 1.5, 2.3, 2.8, 3.5]
            add_point!(tri, x, 0.0)
        end
        for y in [0.2, 0.4, 1.6, 4.5, 6.7, 7.5]
            add_point!(tri, 4.0, y)
        end
        for x in [0.4, 1.5, 1.8, 2.3, 2.5, 3.9]
            add_point!(tri, x, 8.0)
        end
        for y in [7.5, 5.9, 3.4, 1.9, 0.1]
            add_point!(tri, 0.0, y)
        end
        @test DT.get_boundary_nodes(tri) == [1, 46, 47, 2, 48, 3, 49, 50, 4, 51, 5, 52, 53,
            10, 54, 15, 20, 25, 55, 30, 35, 56, 40, 57, 45, 63, 44, 62, 61, 43, 60, 59, 42, 58,
            41, 64, 36, 31, 65, 26, 21, 66, 16, 11, 67, 6, 68, 1]
        @test validate_triangulation(tri)
    end
end


@testset "Starting with a thin set of boundary nodes, and filling them in with automatic collinearity detection" begin
    @testset "Contiguous boundary" begin
        for _ in 1:100
            tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; add_ghost_triangles=true)
            pts = get_points(tri)
            boundary_nodes = [1, 5, 45, 41, 1]
            tri = triangulate(pts; boundary_nodes, randomise=false, delete_ghosts=false)
            flip_edge!(tri, 1, 7)
            flip_edge!(tri, 2, 8)
            @test sort_edge_vector(collect(get_all_constrained_edges(tri))) == sort_edge_vector([
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
                (1, 6, DT.BoundaryIndex)
                (6, 11, DT.BoundaryIndex)
                (11, 16, DT.BoundaryIndex)
                (16, 21, DT.BoundaryIndex)
                (21, 26, DT.BoundaryIndex)
                (26, 31, DT.BoundaryIndex)
                (31, 36, DT.BoundaryIndex)
                (36, 41, DT.BoundaryIndex)
                (41, 42, DT.BoundaryIndex)
                (42, 43, DT.BoundaryIndex)
                (43, 44, DT.BoundaryIndex)
                (44, 45, DT.BoundaryIndex)
                (45, 40, DT.BoundaryIndex)
                (40, 35, DT.BoundaryIndex)
                (35, 30, DT.BoundaryIndex)
                (30, 25, DT.BoundaryIndex)
                (25, 20, DT.BoundaryIndex)
                (20, 15, DT.BoundaryIndex)
                (15, 10, DT.BoundaryIndex)
                (10, 5, DT.BoundaryIndex)
                (5, 4, DT.BoundaryIndex)
                (4, 3, DT.BoundaryIndex)
                (3, 2, DT.BoundaryIndex)
                (2, 1, DT.BoundaryIndex)
            ]
            @test DT.compare_triangle_collections(get_triangles(tri), T)
            @test get_boundary_nodes(tri) == [1, 2, 3, 4, 5,
                10, 15, 20, 25, 30, 35, 40, 45,
                44, 43, 42, 41, 36, 31, 26,
                21, 16, 11, 6, 1]
            @test validate_triangulation(tri)
            add_edge!(tri, 1, 45)
            add_edge!(tri, 6, 44)
            add_edge!(tri, 34, 4)
            @test validate_triangulation(tri)
        end
    end

    @testset "Multiple segments" begin
        for _ in 1:100
            tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; add_ghost_triangles=true)
            pts = get_points(tri)
            boundary_nodes = [[1, 5], [5, 45], [45, 41], [41, 1]]
            tri = triangulate(pts; boundary_nodes, randomise=false, delete_ghosts=false)
            flip_edge!(tri, 1, 7)
            flip_edge!(tri, 2, 8)
            @test sort_edge_vector(collect(get_all_constrained_edges(tri))) == sort_edge_vector([
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
                (1, 6, DT.BoundaryIndex - 3)
                (6, 11, DT.BoundaryIndex - 3)
                (11, 16, DT.BoundaryIndex - 3)
                (16, 21, DT.BoundaryIndex - 3)
                (21, 26, DT.BoundaryIndex - 3)
                (26, 31, DT.BoundaryIndex - 3)
                (31, 36, DT.BoundaryIndex - 3)
                (36, 41, DT.BoundaryIndex - 3)
                (41, 42, DT.BoundaryIndex - 2)
                (42, 43, DT.BoundaryIndex - 2)
                (43, 44, DT.BoundaryIndex - 2)
                (44, 45, DT.BoundaryIndex - 2)
                (45, 40, DT.BoundaryIndex - 1)
                (40, 35, DT.BoundaryIndex - 1)
                (35, 30, DT.BoundaryIndex - 1)
                (30, 25, DT.BoundaryIndex - 1)
                (25, 20, DT.BoundaryIndex - 1)
                (20, 15, DT.BoundaryIndex - 1)
                (15, 10, DT.BoundaryIndex - 1)
                (10, 5, DT.BoundaryIndex - 1)
                (5, 4, DT.BoundaryIndex)
                (4, 3, DT.BoundaryIndex)
                (3, 2, DT.BoundaryIndex)
                (2, 1, DT.BoundaryIndex)
            ]
            @test DT.compare_triangle_collections(get_triangles(tri), T)
            @test get_boundary_nodes(tri) == [[1, 2, 3, 4, 5],
                [5, 10, 15, 20, 25, 30, 35, 40, 45],
                [45, 44, 43, 42, 41],
                [41, 36, 31, 26, 21, 16, 11, 6, 1]]
            @test validate_triangulation(tri)
            add_edge!(tri, 1, 45)
            add_edge!(tri, 6, 44)
            add_edge!(tri, 34, 4)
            @test validate_triangulation(tri)
        end
    end
end

@testset "Adding points and segments into a multiply-connected domain" begin
    tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; add_ghost_triangles=true)
    boundary_nodes = [[
            [1, 5, 45], [45, 41], [41, 1]
        ],
        [[12, 32, 34], [34, 14, 12]]
    ]
    points = get_points(tri)
    rng = StableRNG(19119)
    tri = triangulate(points; boundary_nodes, delete_ghosts=false, randomise=false, rng)
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
        (2, 1, DT.BoundaryIndex)
        (3, 2, DT.BoundaryIndex)
        (4, 3, DT.BoundaryIndex)
        (5, 4, DT.BoundaryIndex)
        (10, 5, DT.BoundaryIndex)
        (15, 10, DT.BoundaryIndex)
        (20, 15, DT.BoundaryIndex)
        (25, 20, DT.BoundaryIndex)
        (30, 25, DT.BoundaryIndex)
        (35, 30, DT.BoundaryIndex)
        (40, 35, DT.BoundaryIndex)
        (45, 40, DT.BoundaryIndex)
        (44, 45, DT.BoundaryIndex - 1)
        (43, 44, DT.BoundaryIndex - 1)
        (42, 43, DT.BoundaryIndex - 1)
        (41, 42, DT.BoundaryIndex - 1)
        (36, 41, DT.BoundaryIndex - 2)
        (31, 36, DT.BoundaryIndex - 2)
        (26, 31, DT.BoundaryIndex - 2)
        (21, 26, DT.BoundaryIndex - 2)
        (16, 21, DT.BoundaryIndex - 2)
        (11, 16, DT.BoundaryIndex - 2)
        (6, 11, DT.BoundaryIndex - 2)
        (1, 6, DT.BoundaryIndex - 2)
        (17, 12, DT.BoundaryIndex - 3)
        (22, 17, DT.BoundaryIndex - 3)
        (27, 22, DT.BoundaryIndex - 3)
        (32, 27, DT.BoundaryIndex - 3)
        (33, 32, DT.BoundaryIndex - 3)
        (34, 33, DT.BoundaryIndex - 3)
        (29, 34, DT.BoundaryIndex - 4)
        (24, 29, DT.BoundaryIndex - 4)
        (19, 24, DT.BoundaryIndex - 4)
        (14, 19, DT.BoundaryIndex - 4)
        (13, 14, DT.BoundaryIndex - 4)
        (12, 13, DT.BoundaryIndex - 4)
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
    @test sort_edge_vector(collect(get_all_constrained_edges(tri))) == sort_edge_vector(C)
    @test DT.compare_triangle_collections(get_triangles(tri), T)
    @test validate_triangulation(tri; check_planarity=false)
    @test get_boundary_nodes(tri) == [[
            [1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45], [45, 44, 43, 42, 41], [41, 36, 31, 26, 21, 16, 11, 6, 1]
        ],
        [[12, 17, 22, 27, 32, 33, 34], [34, 29, 24, 19, 14, 13, 12]]
    ]
    add_point!(tri, 1.5, 0.0; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 2.5, 0.0; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 4.0, 2.5; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 4.0, 7.5; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 2.5, 8.0; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 0.0, 5.5; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 0.5, 2.2; rng)
    @test validate_triangulation(tri)
    add_edge!(tri, 21, 27; rng)
    @test validate_triangulation(tri)
    add_edge!(tri, 14, 45; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 1.0, 2.5; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 1.0, 5.5; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 2.5, 2.0; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 3.0, 5.5; rng)
    @test validate_triangulation(tri)
    add_point!(tri, 2.5, 6.0; rng)
    @test validate_triangulation(tri)
    add_edge!(tri, 12, 5)
    @test validate_triangulation(tri)
    add_edge!(tri, 2, 52)
    @test validate_triangulation(tri)
    add_edge!(tri, 53, 21)
    @test validate_triangulation(tri)
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
    @test validate_triangulation(tri)
    add_point!(tri, 0.5, 4.4)
    add_point!(tri, 0.5, 4.6)
    @test validate_triangulation(tri)
    add_point!(tri, 0.5, 4.5)
    @test validate_triangulation(tri)

    for L in 1:100
        @show L
        tri = triangulate_rectangle(0, 4, 0, 8, 5, 9; add_ghost_triangles=true)
        boundary_nodes = [[
                [1, 5, 45], [45, 41], [41, 1]
            ],
            [[12, 32, 34], [34, 14, 12]]
        ]
        points = get_points(tri)
        tri = triangulate(points; boundary_nodes, delete_ghosts=false)
        @test validate_triangulation(tri)
        @test get_boundary_nodes(tri) == [[
                [1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45], [45, 44, 43, 42, 41], [41, 36, 31, 26, 21, 16, 11, 6, 1]
            ],
            [[12, 17, 22, 27, 32, 33, 34], [34, 29, 24, 19, 14, 13, 12]]
        ]
        add_point!(tri, 1.5, 0.0)
        @test validate_triangulation(tri)
        add_point!(tri, 2.5, 0.0)
        @test validate_triangulation(tri)
        add_point!(tri, 4.0, 2.5)
        @test validate_triangulation(tri)
        add_point!(tri, 4.0, 7.5)
        @test validate_triangulation(tri)
        add_point!(tri, 2.5, 8.0)
        @test validate_triangulation(tri)
        add_point!(tri, 0.0, 5.5)
        @test validate_triangulation(tri)
        add_point!(tri, 0.5, 2.2)
        @test validate_triangulation(tri)
        add_edge!(tri, 21, 27)
        @test validate_triangulation(tri)
        add_edge!(tri, 14, 45)
        @test validate_triangulation(tri)
        add_point!(tri, 1.0, 2.5)
        @test validate_triangulation(tri)
        add_point!(tri, 1.0, 5.5)
        @test validate_triangulation(tri)
        add_point!(tri, 2.5, 2.0)
        @test validate_triangulation(tri)
        add_point!(tri, 3.0, 5.5)
        @test validate_triangulation(tri)
        add_point!(tri, 2.5, 6.0)
        @test validate_triangulation(tri)
        add_edge!(tri, 12, 5)
        @test validate_triangulation(tri)
        add_edge!(tri, 2, 52)
        @test validate_triangulation(tri)
        add_edge!(tri, 53, 21)
        @test validate_triangulation(tri)
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
        @test validate_triangulation(tri)
        add_point!(tri, 0.5, 4.4)
        add_point!(tri, 0.5, 4.6)
        @test validate_triangulation(tri)
        add_point!(tri, 0.5, 4.5)
        @test validate_triangulation(tri)
    end
end