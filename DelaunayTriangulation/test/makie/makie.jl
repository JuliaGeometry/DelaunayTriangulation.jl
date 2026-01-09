using DelaunayTriangulation
using ReferenceTests
using CairoMakie
using Test
using StableRNGs
using GeometryBasics
using Random

const STABLE_RNG = StableRNG(123)
reseed!() = Random.seed!(STABLE_RNG, 123)
_randn(args...) = randn(STABLE_RNG, args...)
_rand(args...) = rand(STABLE_RNG, args...)
const RNG = (; reseed!, randn=_randn, rand=_rand, STABLE_RNG)

@testset "tricontour" begin
    @test_reference "tricontourf.png" begin
        reseed!()
        x = RNG.randn(50)
        y = RNG.randn(50)
        z = -sqrt.(x .^ 2 .+ y .^ 2) .+ 0.1 .* RNG.randn.()

        f, ax, tr = tricontourf(x, y, z)
        scatter!(x, y, color=z, strokewidth=1, strokecolor=:black)
        Colorbar(f[1, 2], tr)
        f
    end

    @test_reference "tricontourf extendhigh extendlow.png" begin
        reseed!()
        x = RNG.randn(50)
        y = RNG.randn(50)
        z = -sqrt.(x .^ 2 .+ y .^ 2) .+ 0.1 .* RNG.randn.()

        f, ax, tr = tricontourf(x, y, z, levels=-1.8:0.2:-0.4, extendhigh=:red, extendlow=:orange)
        scatter!(x, y, color=z, strokewidth=1, strokecolor=:black)
        Colorbar(f[1, 2], tr)
        f
    end

    @test_reference "tricontourf relative mode.png" begin
        reseed!()
        x = RNG.randn(50)
        y = RNG.randn(50)
        z = -sqrt.(x .^ 2 .+ y .^ 2) .+ 0.1 .* RNG.randn.()

        f, ax, tr = tricontourf(x, y, z, mode=:relative, levels=0.2:0.1:1, colormap=:batlow)
        scatter!(x, y, color=z, strokewidth=1, strokecolor=:black, colormap=:batlow)
        Colorbar(f[1, 2], tr)
        f
    end

    @test_reference "tricontourf manual vs delaunay.png" begin
        reseed!()
        n = 20
        angles = range(0, 2pi, length=n + 1)[1:end-1]
        x = [cos.(angles); 2 .* cos.(angles .+ pi / n)]
        y = [sin.(angles); 2 .* sin.(angles .+ pi / n)]
        z = (x .- 0.5) .^ 2 + (y .- 0.5) .^ 2 .+ 0.5 .* RNG.randn.()

        triangulation_inner = reduce(hcat, map(i -> [0, 1, n] .+ i, 1:n))
        triangulation_outer = reduce(hcat, map(i -> [n - 1, n, 0] .+ i, 1:n))
        triangulation = hcat(triangulation_inner, triangulation_outer)

        f, ax, _ = tricontourf(x, y, z, triangulation=triangulation,
            axis=(; aspect=1, title="Manual triangulation"))
        scatter!(x, y, color=z, strokewidth=1, strokecolor=:black)

        tricontourf(f[1, 2], x, y, z, triangulation=Makie.DelaunayTriangulation(),
            axis=(; aspect=1, title="Delaunay triangulation"))
        scatter!(x, y, color=z, strokewidth=1, strokecolor=:black)

        f
    end

    @test_reference "tricontourf with boundary nodes.png" begin
        reseed!()
        n = 20
        angles = range(0, 2pi, length=n + 1)[1:end-1]
        x = [cos.(angles); 2 .* cos.(angles .+ pi / n)]
        y = [sin.(angles); 2 .* sin.(angles .+ pi / n)]
        z = (x .- 0.5) .^ 2 + (y .- 0.5) .^ 2 .+ 0.5 .* RNG.randn.()

        inner = [n:-1:1; n] # clockwise inner
        outer = [(n+1):(2n); n + 1] # counter-clockwise outer
        boundary_nodes = [[outer], [inner]]
        tri = triangulate([x'; y'], boundary_nodes=boundary_nodes)
        f, ax, _ = tricontourf(tri, z)
        scatter!(x, y, color=z, strokewidth=1, strokecolor=:black)
        f
    end

    @test_reference "tricontourf with boundary nodes and edges.png" begin
        reseed!()
        curve_1 = [
            [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0), (15.0, 0.0), (20.0, 0.0), (25.0, 0.0)],
            [(25.0, 0.0), (25.0, 5.0), (25.0, 10.0), (25.0, 15.0), (25.0, 20.0), (25.0, 25.0)],
            [(25.0, 25.0), (20.0, 25.0), (15.0, 25.0), (10.0, 25.0), (5.0, 25.0), (0.0, 25.0)],
            [(0.0, 25.0), (0.0, 20.0), (0.0, 15.0), (0.0, 10.0), (0.0, 5.0), (0.0, 0.0)]
        ]
        curve_2 = [
            [(4.0, 6.0), (4.0, 14.0), (4.0, 20.0), (18.0, 20.0), (20.0, 20.0)],
            [(20.0, 20.0), (20.0, 16.0), (20.0, 12.0), (20.0, 8.0), (20.0, 4.0)],
            [(20.0, 4.0), (16.0, 4.0), (12.0, 4.0), (8.0, 4.0), (4.0, 4.0), (4.0, 6.0)]
        ]
        curve_3 = [
            [(12.906, 10.912), (16.0, 12.0), (16.16, 14.46), (16.29, 17.06),
            (13.13, 16.86), (8.92, 16.4), (8.8, 10.9), (12.906, 10.912)]
        ]
        curves = [curve_1, curve_2, curve_3]
        points = [
            (3.0, 23.0), (9.0, 24.0), (9.2, 22.0), (14.8, 22.8), (16.0, 22.0),
            (23.0, 23.0), (22.6, 19.0), (23.8, 17.8), (22.0, 14.0), (22.0, 11.0),
            (24.0, 6.0), (23.0, 2.0), (19.0, 1.0), (16.0, 3.0), (10.0, 1.0), (11.0, 3.0),
            (6.0, 2.0), (6.2, 3.0), (2.0, 3.0), (2.6, 6.2), (2.0, 8.0), (2.0, 11.0),
            (5.0, 12.0), (2.0, 17.0), (3.0, 19.0), (6.0, 18.0), (6.5, 14.5),
            (13.0, 19.0), (13.0, 12.0), (16.0, 8.0), (9.8, 8.0), (7.5, 6.0),
            (12.0, 13.0), (19.0, 15.0)
        ]
        boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points)
        edges = Set(((1, 19), (19, 12), (46, 4), (45, 12)))

        tri = triangulate(points; boundary_nodes=boundary_nodes, segments=edges, check_arguments=false)
        z = [(x - 1) * (y + 1) for (x, y) in DelaunayTriangulation.each_point(tri)]
        f, ax, _ = tricontourf(tri, z, levels=30)
        f
    end

    @test_reference "tricontourf with provided triangulation.png" begin
        reseed!()
        θ = [LinRange(0, 2π * (1 - 1 / 19), 20); 0]
        xy = Vector{Vector{Vector{NTuple{2,Float64}}}}()
        cx = [0.0, 3.0]
        for i in 1:2
            push!(xy, [[(cx[i] + cos(θ), sin(θ)) for θ in θ]])
            push!(xy, [[(cx[i] + 0.5cos(θ), 0.5sin(θ)) for θ in reverse(θ)]])
        end
        boundary_nodes, points = convert_boundary_points_to_indices(xy)
        tri = triangulate(points; boundary_nodes=boundary_nodes, check_arguments=false)
        z = [(x - 3 / 2)^2 + y^2 for (x, y) in DelaunayTriangulation.each_point(tri)]

        f, ax, tr = tricontourf(tri, z, colormap=:matter)
        f
    end
end

@testset "triplot" begin
    @test_reference "Triplot with points, ghost edges, and convex hull.png" begin
        reseed!()
        pts = RNG.rand(2, 50)
        tri = triangulate(pts; rng=RNG.STABLE_RNG)
        fig, ax, sc = triplot(tri,
            triangle_color=:lightgray, strokewidth=4,
            show_points=true, markersize=20, markercolor=:orange,
            show_ghost_edges=true, ghost_edge_linewidth=4,
            show_convex_hull=true, convex_hull_linewidth=4)
        fig
    end

    @test_reference "Triplot of a constrained triangulation with holes and a custom bounding box.png" begin
        reseed!()
        curve_1 = [[
            (0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (12.0, 0.0), (12.0, 4.0),
            (12.0, 8.0), (14.0, 10.0), (16.0, 12.0), (16.0, 16.0),
            (14.0, 18.0), (12.0, 20.0), (12.0, 24.0), (12.0, 28.0),
            (8.0, 28.0), (4.0, 28.0), (0.0, 28.0), (-2.0, 26.0), (0.0, 22.0),
            (0.0, 18.0), (0.0, 10.0), (0.0, 8.0), (0.0, 4.0), (-4.0, 4.0),
            (-4.0, 0.0), (0.0, 0.0),
        ]]
        curve_2 = [[
            (4.0, 26.0), (8.0, 26.0), (10.0, 26.0), (10.0, 24.0),
            (10.0, 22.0), (10.0, 20.0), (8.0, 20.0), (6.0, 20.0),
            (4.0, 20.0), (4.0, 22.0), (4.0, 24.0), (4.0, 26.0)
        ]]
        curve_3 = [[(4.0, 16.0), (12.0, 16.0), (12.0, 14.0), (4.0, 14.0), (4.0, 16.0)]]
        curve_4 = [[(4.0, 8.0), (10.0, 8.0), (8.0, 6.0), (6.0, 6.0), (4.0, 8.0)]]
        curves = [curve_1, curve_2, curve_3, curve_4]
        points = [
            (2.0, 26.0), (2.0, 24.0), (6.0, 24.0), (6.0, 22.0), (8.0, 24.0), (8.0, 22.0),
            (2.0, 22.0), (0.0, 26.0), (10.0, 18.0), (8.0, 18.0), (4.0, 18.0), (2.0, 16.0),
            (2.0, 12.0), (6.0, 12.0), (2.0, 8.0), (2.0, 4.0), (4.0, 2.0),
            (-2.0, 2.0), (4.0, 6.0), (10.0, 2.0), (10.0, 6.0), (8.0, 10.0), (4.0, 10.0),
            (10.0, 12.0), (12.0, 12.0), (14.0, 26.0), (16.0, 24.0), (18.0, 28.0),
            (16.0, 20.0), (18.0, 12.0), (16.0, 8.0), (14.0, 4.0), (14.0, -2.0),
            (6.0, -2.0), (2.0, -4.0), (-4.0, -2.0), (-2.0, 8.0), (-2.0, 16.0),
            (-4.0, 22.0), (-4.0, 26.0), (-2.0, 28.0), (6.0, 15.0), (7.0, 15.0),
            (8.0, 15.0), (9.0, 15.0), (10.0, 15.0), (6.2, 7.8),
            (5.6, 7.8), (5.6, 7.6), (5.6, 7.4), (6.2, 7.4), (6.0, 7.6),
            (7.0, 7.8), (7.0, 7.4)]
        boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points)
        tri = triangulate(points; randomise=false, boundary_nodes=boundary_nodes, rng=RNG.STABLE_RNG)
        fig, ax, sc = triplot(tri,
            show_points=true,
            show_constrained_edges=true,
            constrained_edge_linewidth=2,
            strokewidth=0.2,
            markersize=15,
            markercolor=:blue,
            show_ghost_edges=true, # not as good because the outer boundary is not convex, but just testing
            marker='x',
            bounding_box=(-5, 20, -5, 35)) # also testing the conversion to Float64 for bbox here
        fig
    end

    @test_reference "Triplot with nonlinear transformation.png" begin
        reseed!()
        f = Figure()
        ax = PolarAxis(f[1, 1])
        points = Point2f[(phi, r) for r in 1:10 for phi in range(0, 2pi, length=36)[1:35]]
        noise = i -> 1.0f-4 * (isodd(i) ? 1 : -1) * i / sqrt(50) # should have small discrepancy
        points = points .+ [Point2f(noise(i), noise(i)) for i in eachindex(points)]
        # The noise forces the triangulation to be unique. Not using RNG to not disrupt the RNG stream later
        tr = triplot!(ax, points)
        f
    end

    @test_reference "Triplot after adding points and make sure the representative_point_list is correctly updated.png" begin
        reseed!()
        points = [(0.0, 0.0), (0.95, 0.0), (1.0, 1.4), (0.0, 1.0)] # not 1 so that we have a unique triangulation
        tri = Observable(triangulate(points; delete_ghosts=false))
        fig, ax, sc = triplot(tri, show_points=true, markersize=14, show_ghost_edges=true, recompute_centers=true, linestyle=:dash)
        for p in [(0.3, 0.5), (-1.5, 2.3), (0.2, 0.2), (0.2, 0.5)]
            add_point!(tri[], p)
        end
        convex_hull!(tri[])
        notify(tri)
        ax = Axis(fig[1, 2])
        triplot!(ax, tri[], show_points=true, markersize=14, show_ghost_edges=true, recompute_centers=true)
        fig
    end

    @test_reference "Triplot Showing ghost edges for a triangulation with disjoint boundaries.png" begin
        reseed!()
        θ = LinRange(0, 2π, 20) |> collect
        θ[end] = 0 # need to make sure that 2π gives the exact same coordinates as 0
        xy = Vector{Vector{Vector{NTuple{2,Float64}}}}()
        cx = 0.0
        for i in 1:2
            ## Make the exterior circle
            push!(xy, [[(cx + cos(θ), sin(θ)) for θ in θ]])
            ## Now the interior circle - clockwise
            push!(xy, [[(cx + 0.5cos(θ), 0.5sin(θ)) for θ in reverse(θ)]])
            cx += 3.0
        end
        boundary_nodes, points = convert_boundary_points_to_indices(xy)
        tri = triangulate(points; boundary_nodes=boundary_nodes, check_arguments=false)
        fig, ax, sc = triplot(tri, show_ghost_edges=true)
        fig
    end
end

@testset "voronoiplot" begin
    @test_reference "Voronoiplot for a tessellation with a custom bounding box.png" begin
        reseed!()
        pts = 25RNG.randn(2, 50)
        tri = triangulate(pts; rng=RNG.STABLE_RNG)
        vorn = voronoi(tri, clip=false)
        fig, ax, sc = voronoiplot(vorn,
            show_generators=true,
            colormap=:RdBu,
            strokecolor=:white,
            strokewidth=4,
            markersize=25,
            marker='x',
            markercolor=:green,
            unbounded_edge_extension_factor=5.0)
        xlims!(ax, -120, 120)
        ylims!(ax, -120, 120)
        fig
    end

    @test_reference "Voronoiplots with clipped tessellation and unbounded polygons.png" begin
        reseed!()
        pts = 25RNG.randn(2, 10)
        tri = triangulate(pts; rng=RNG.STABLE_RNG)
        vorn = voronoi(tri, clip=true)
        fig, ax, sc = voronoiplot(vorn, color=(:blue, 0.2), markersize=20, strokewidth=4)

        # used to be bugged
        points = [(0.0, 1.0), (-1.0, 2.0), (-2.0, -1.0)]
        tri = triangulate(points)
        vorn = voronoi(tri)
        voronoiplot(fig[1, 2], vorn, show_generators=true, strokewidth=4,
            color=[:red, :blue, :green], markercolor=:white, markersize=20)

        fig
    end

    @test_reference "Voronoiplot with a nonlinear transform.png" begin
        reseed!()
        f = Figure()
        ax = PolarAxis(f[1, 1], theta_as_x=false)
        points = Point2d[(r, phi) for r in 1:10 for phi in range(0, 2pi, length=36)[1:35]]
        noise = i -> 1.0f-4 * (isodd(i) ? 1 : -1) * i / sqrt(50) # should have small discrepancy
        points = points .+ [Point2f(noise(i), noise(i)) for i in eachindex(points)] # make triangulation unique
        polygon_color = [r for r in 1:10 for phi in range(0, 2pi, length=36)[1:35]]
        polygon_color_2 = [phi for r in 1:10 for phi in range(0, 2pi, length=36)[1:35]]
        tr = voronoiplot!(ax, points, smooth=false, show_generators=false, color=polygon_color)
        Makie.rlims!(ax, 12) # to make rect clip visible if circular clip doesn't happen
        ax = PolarAxis(f[1, 2], theta_as_x=false)
        tr = voronoiplot!(ax, points, smooth=true, show_generators=false, color=polygon_color_2)
        Makie.rlims!(ax, 12)
        f
    end

    @test_reference "Voronoiplot with some custom bounding boxes may not contain all data sites.png" begin
        reseed!()
        points = [(-3.0, 7.0), (1.0, 6.0), (-1.0, 3.0), (-2.0, 4.0), (3.0, -2.0), (5.0, 5.0), (-4.0, -3.0), (3.0, 8.0)]
        tri = triangulate(points)
        vorn = voronoi(tri)
        color = [:red, :blue, :green, :yellow, :cyan, :magenta, :black, :brown] # the polygon colors should not change even if some are not included (because they're outside of the box)
        fig = Figure()
        ax1 = Axis(fig[1, 1], title="Default")
        voronoiplot!(ax1, vorn, show_generators=true, markersize=14, strokewidth=4, color=color)
        ax2 = Axis(fig[1, 2], title="Some excluded")
        voronoiplot!(ax2, vorn, show_generators=true, markersize=14, strokewidth=4, color=color, clip=BBox(0.0, 5.0, -15.0, 15.0))
        ax3 = Axis(fig[2, 1], title="Bigger range")
        voronoiplot!(ax3, vorn, show_generators=true, markersize=14, strokewidth=4, color=color, clip=(-15.0, 15.0, -15.0, 15.0))
        ax4 = Axis(fig[2, 2], title="Only one polygon")
        voronoiplot!(ax4, vorn, show_generators=true, markersize=14, strokewidth=4, color=color, clip=(10.0, 12.0, 2.0, 5.0))
        for ax in fig.content
            xlims!(ax4, -15, 15)
            ylims!(ax4, -15, 15)
        end
        fig
    end

    @test_reference "Voronoiplot after adding points.png" begin
        reseed!()
        points = Observable([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
        fig, ax, sc = voronoiplot(points, show_generators=true, markersize=36) # make sure any regressions with missing generators are identified, so use 36
        push!(points[], (2.0, 2.0), (0.5, 0.5), (0.25, 0.25), (0.25, 0.75), (0.75, 0.25), (0.75, 0.75))
        notify(points)
        ax2 = Axis(fig[1, 2])
        voronoiplot!(ax2, voronoi(triangulate(points[])), show_generators=true, markersize=36)
        xlims!(ax, -0.5, 2.5)
        ylims!(ax, -0.5, 2.5)
        xlims!(ax2, -0.5, 2.5)
        ylims!(ax2, -0.5, 2.5) # need to make sure all generators are shown, and the bounding box is automatically updated
        fig
    end
end