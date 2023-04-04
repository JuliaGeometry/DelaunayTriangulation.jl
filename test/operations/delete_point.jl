using ..DelaunayTriangulation
const DT = DelaunayTriangulation
# using CairoMakie
# import SimpleGraphs: relabel, UndirectedGraph
# using DataStructures
using StableRNGs

include("../test_setup.jl")

include("../helper_functions.jl")

#=
_tri = tri
fig, ax, sc = triplot(_tri)
let vert = each_solid_vertex(_tri)
    text!(ax, collect(get_point(_tri, vert...)); text=string.(vert))
end
fig
=#

@testset verbose = true "Deleting interior nodes" begin
    @testset "Random point sets" begin
        for j in 1:20
            rng1 = StableRNG(j)
            n = rand(rng1, 50:1000)
            points = 20randn(rng1, 2, n)
            tri = triangulate(points; delete_ghosts=false, rng=rng1)
            deleted_pts = Int64[]
            for k in 1:(n÷10)
                rng2 = StableRNG(j + k * n)
                @show j, k, n
                i = rand(rng2, each_solid_vertex(tri) |> collect)
                while DT.is_boundary_node(tri, i)[1]
                    i = rand(rng2, each_solid_vertex(tri) |> collect)
                end
                delete_point!(tri, i; rng=rng2)
                push!(deleted_pts, i)
                _tri = triangulate(points; delete_ghosts=false, skip_points=deleted_pts, rng=rng2)
                clear_empty_features!(tri)
                clear_empty_features!(_tri)
                @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
                @test get_adjacent(tri) == get_adjacent(_tri)
                @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
                @test get_graph(tri) == get_graph(_tri)
                convex_hull!(tri)
                @test validate_triangulation(tri)
            end
        end
    end
    @testset "Lattice with a single boundary index" begin
        for j in 1:20
            rng1 = StableRNG(j)
            a, b = sort(10randn(rng1, 2))
            c, d = sort(15randn(rng1, 2))
            nx = rand(rng1, 5:25)
            ny = rand(rng1, 5:25)
            tri = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
            points = get_points(tri)
            n = nx * ny
            for k in 1:(n÷10)
                @show j, k
                rng2 = StableRNG(j + k * n)
                i = rand(rng2, each_solid_vertex(tri) |> collect)
                while DT.is_boundary_node(tri, i)[1]
                    i = rand(rng2, each_solid_vertex(tri) |> collect)
                end
                delete_point!(tri, i; rng=rng2)
                @test validate_triangulation(tri)
            end
        end
    end
    @testset "Lattice with multiple boundary indices" begin
        for j in 1:20
            rng1 = StableRNG(j)
            a, b = sort(10randn(rng1, 2))
            c, d = sort(15randn(rng1, 2))
            nx = rand(rng1, 5:25)
            ny = rand(rng1, 5:25)
            tri = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=false)
            points = get_points(tri)
            n = nx * ny
            for k in 1:(n÷10)
                @show j, k, a, b, c, d, nx, ny
                rng2 = StableRNG(j + k * n)
                i = rand(rng2, each_solid_vertex(tri) |> collect)
                while DT.is_boundary_node(tri, i)[1]
                    i = rand(rng2, each_solid_vertex(tri) |> collect)
                end
                delete_point!(tri, i; rng=rng2)
                @test validate_triangulation(tri)
            end
        end
    end
end

@testset verbose = true "Deleting boundary nodes" begin
    @testset "Random point sets" begin
        for j in 1:20
            rng1 = StableRNG(j)
            n = rand(rng1, 50:1000)
            r = rand(rng1, n)
            θ = rand(rng1, n)
            x = @. sqrt(r) * cos(θ)
            y = @. sqrt(r) * sin(θ)
            points = [(x, y) for (x, y) in zip(x, y)]
            bnd_x = 20cos.(LinRange(0, 2π - eps(Float64), 250))
            bnd_y = 20sin.(LinRange(0, 2π - eps(Float64), 250))
            bnd_points = [(x, y) for (x, y) in zip(bnd_x, bnd_y)]
            append!(points, bnd_points)
            tri = triangulate(points; delete_ghosts=false, rng=rng1)
            deleted_pts = Int64[]
            total_bnd_nodes = DT.num_neighbours(tri, DT.BoundaryIndex)
            for k in 1:(max(1, total_bnd_nodes ÷ 5))
                rng2 = StableRNG(j + k * n)
                @show j, k, n
                i = rand(rng2, each_solid_vertex(tri) |> collect)
                while !DT.is_boundary_node(tri, i)[1]
                    i = rand(rng2, each_solid_vertex(tri) |> collect)
                end
                delete_point!(tri, i; rng=rng2)
                convex_hull!(tri)
                DT.compute_representative_points!(tri)
                push!(deleted_pts, i)
                @test validate_triangulation(tri; ignore_boundary_indices=true)
            end
        end
    end
    @testset "Lattice with a single boundary index" begin
        for j in 1:20
            rng1 = StableRNG(j)
            a, b = sort(10randn(rng1, 2))
            c, d = sort(15randn(rng1, 2))
            nx = rand(rng1, 5:25)
            ny = rand(rng1, 5:25)
            tri = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=true)
            empty!(get_all_constrained_edges(tri))
            points = get_points(tri)
            n = nx * ny
            total_bnd_nodes = DT.num_neighbours(tri, DT.BoundaryIndex)
            for k in 1:(max(1, total_bnd_nodes ÷ 10))
                @show j, k, a, b, c, d, nx, ny
                rng2 = StableRNG(j + k * n)
                i = rand(rng2, each_solid_vertex(tri) |> collect)
                while !DT.is_boundary_node(tri, i)[1]
                    i = rand(rng2, each_solid_vertex(tri) |> collect)
                end
                delete_point!(tri, i; rng=rng2)
                @test validate_triangulation(tri)
            end
        end
    end
    #=
    @testset "Lattice with multiple boundary indices" begin
        for j in 1:50
            rng1 = StableRNG(j)
            a, b = sort(10randn(rng1, 2))
            c, d = sort(15randn(rng1, 2))
            nx = rand(rng1, 5:25)
            ny = rand(rng1, 5:25)
            tri = triangulate_rectangle(a, b, c, d, nx, ny; add_ghost_triangles=true, single_boundary=false)
            points = get_points(tri)
            n = nx * ny
            total_bnd_nodes = DT.num_neighbours(tri, DT.BoundaryIndex) + DT.num_neighbours(tri, DT.BoundaryIndex - 1) + DT.num_neighbours(tri, DT.BoundaryIndex - 2) + DT.num_neighbours(tri, DT.BoundaryIndex - 3)
            for k in 1:(max(1, total_bnd_nodes ÷ 10))
                @show j, k
                rng2 = StableRNG(j + k * n)
                i = rand(rng2, each_solid_vertex(tri) |> collect)
                while !DT.is_boundary_node(tri, i)[1]
                    i = rand(rng2, each_solid_vertex(tri) |> collect)
                end
                @show i
                delete_point!(tri, i; rng=rng2)
                @test validate_triangulation(tri;
                    ignore_boundary_indices=true)
            end
        end
    end
    =#
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
        add_triangle!(_tri, T; update_ghost_edges=true)
    end
    convex_hull!(tri)
    DT.compute_representative_points!(tri)
    clear_empty_features!(_tri)
    clear_empty_features!(tri)
    @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
    @test get_adjacent(tri) == get_adjacent(_tri)
    @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
    @test get_graph(tri) == get_graph(_tri)
    _tri = triangulate(get_points(tri); skip_points=16, delete_ghosts=false)
    clear_empty_features!(_tri)
    @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
    @test get_adjacent(tri) == get_adjacent(_tri)
    @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
    @test get_graph(tri) == get_graph(_tri)
    @test get_convex_hull(tri) == get_convex_hull(_tri)
    deleted_points = [16]
    points_to_delete = [15, 14, 13, 17, 2, 9, 7, 10, 8, 11, 18, 5, 12]
    local ctr = 5
    #TRI = Ref{Any}()
    #_TRI = Ref{Any}()
    for j in 1:100
        _deleted_points = deepcopy(deleted_points)
        _points_to_delete = deepcopy(points_to_delete)
        ctr += 1
        tri = example_with_special_corners()
        tri.points[14] = [1.2, 7.0]
        #TRI[] = tri
        rng = StableRNG(292929)
        point = 16
        delete_point!(tri, 16; rng)
        while !isempty(_points_to_delete)
            rng = StableRNG(ctr)
            ctr += 1
            i = rand(rng, eachindex(_points_to_delete))
            point = _points_to_delete[i]
            #@show ctr - 1, _points_to_delete, _deleted_points, i, point
            deleteat!(_points_to_delete, i)
            push!(_deleted_points, point)
            rng = StableRNG(ctr)
            ctr += 1
            #TRI[] = deepcopy(tri)
            #_TRI[] = _tri
            delete_point!(tri, point; rng)
            convex_hull!(tri)
            DT.compute_representative_points!(tri)
            _tri = triangulate(get_points(tri); skip_points=_deleted_points, delete_ghosts=false, rng)
            @test validate_triangulation(tri)
            rng = StableRNG(ctr)
            ctr += 1
            clear_empty_features!(tri)
            clear_empty_features!(_tri)
            @test DT.compare_triangle_collections(get_triangles(_tri), get_triangles(tri))
            @test get_adjacent(tri) == get_adjacent(_tri)
            @test get_adjacent2vertex(tri) == get_adjacent2vertex(_tri)
            @test get_graph(tri) == get_graph(_tri)
            @test get_convex_hull(tri) == get_convex_hull(_tri)
        end
    end
end