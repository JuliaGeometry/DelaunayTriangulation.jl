using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
using DataStructures
using Random
using LinearAlgebra
using ReferenceTests
using Test
using StructEquality
using DelimitedFiles

using StableRNGs
@struct_equal DT.InsertionEventHistory
@struct_equal DT.RefinementConstraints
@struct_equal DT.RefinementQueue

const ppoints = rand(2, 50)
const ptri = triangulate(ppoints) # p for placeholder 
const pT = (1, 2, 3)

@testset verbose = true "RefinementConstraints" begin
    @testset "Defaults" begin
        constraints = DT.RefinementConstraints()
        @test constraints.min_angle == 0.0
        @test constraints.max_angle == 180.0
        @test constraints.min_area == 0.0
        @test constraints.max_area == Inf
        @test constraints.max_radius_edge_ratio == Inf
        @test constraints.max_points == typemax(Int)
        @test constraints.seditious_angle == 20.0
        @test !constraints.custom_constraint(ptri, pT)
        @test !DT.violates_custom_constraint(constraints, ptri, pT)
        @test sprint(show, MIME"text/plain"(), constraints) == "RefinementConstraints\n   θₘᵢₙ°: 0.0\n   θₘₐₓ°: 180.0\n   Aₘᵢₙ: 0.0\n   Aₘₐₓ: Inf\n   ρₘₐₓ: Inf\n   max_points: ∞\n   θ_seditious°: 20.0\n   custom_constraint: $(constraints.custom_constraint)"
    end

    @testset "Setting some values" begin
        constraints = DT.RefinementConstraints(min_area=1e-9, max_area=13.7, max_points=505, seditious_angle=60.0)
        @test constraints.min_angle == 0.0
        @test constraints.max_angle == 180.0
        @test constraints.min_area == 1e-9
        @test constraints.max_area == 13.7
        @test constraints.max_radius_edge_ratio == Inf
        @test constraints.max_points == 505
        @test constraints.seditious_angle == 60.0
        @test !constraints.custom_constraint(ptri, pT)
        @test !DT.violates_custom_constraint(constraints, ptri, pT)
        @test !DT.has_max_angle_constraint(constraints)
    end

    @testset "Constraints and setting min_angle" begin
        constraints = DT.RefinementConstraints(min_angle=30.0, custom_constraint=(tri, T) -> T == (1, 2, 3))
        @inferred DT.RefinementConstraints(min_angle=30.0, custom_constraint=(tri, T) -> T == (1, 2, 3))
        @test constraints.min_angle == 30.0
        @test constraints.max_angle == 180.0
        @test constraints.min_area == 0.0
        @test constraints.max_area == Inf
        @test constraints.max_radius_edge_ratio == 1.0
        @test constraints.max_points == typemax(Int)
        @test constraints.seditious_angle == 20.0
        @test constraints.custom_constraint(ptri, pT)
        @test DT.violates_custom_constraint(constraints, ptri, pT)
        @test !constraints.custom_constraint(ptri, (1, 2, 4))
        @test !DT.violates_custom_constraint(constraints, ptri, (1, 2, 4))
        @test !DT.has_max_angle_constraint(constraints)
    end

    @testset "Setting a maximum angle constraint and a non-30 min_angle" begin
        constraints = DT.RefinementConstraints(max_angle=179.9, min_angle=27.59, min_area=1e-9)
        @test constraints.min_angle == 27.59
        @test constraints.max_angle == 179.9
        @test constraints.min_area == 1e-9
        @test constraints.max_area == Inf
        @test constraints.max_radius_edge_ratio == cscd(27.59) / 2
        @test constraints.max_points == typemax(Int)
        @test constraints.seditious_angle == 20.0
        @test !constraints.custom_constraint(ptri, pT)
        @test !DT.violates_custom_constraint(constraints, ptri, pT)
        @test DT.has_max_angle_constraint(constraints)
        @test sprint(show, MIME"text/plain"(), constraints) == "RefinementConstraints\n   θₘᵢₙ°: 27.59\n   θₘₐₓ°: 179.9\n   Aₘᵢₙ: 1.0e-9\n   Aₘₐₓ: Inf\n   ρₘₐₓ: 1.0795840041760634\n   max_points: ∞\n   θ_seditious°: 20.0\n   custom_constraint: $(constraints.custom_constraint)"
    end
end

@testset verbose = true "RefinementQueue" begin
    T = NTuple{3,Int}
    E = NTuple{2,Int}
    F = Float64

    @testset "Initialisation" begin
        queue1 = DT.RefinementQueue{T,E,F}()
        queue2 = DT.RefinementQueue(ptri)
        @inferred DT.RefinementQueue(ptri)
        for queue in (queue1, queue2)
            @test queue.segments == DT.MaxPriorityQueue{E,F}()
            @test typeof(queue.segments) == DT.MaxPriorityQueue{E,F}
            @test queue.triangles == DT.MaxPriorityQueue{T,F}()
            @test typeof(queue.triangles) == DT.MaxPriorityQueue{T,F}
            @test isempty(queue.segments)
            @test isempty(queue.triangles)
            @test !DT.has_segments(queue)
            @test !DT.has_triangles(queue)
            @test isempty(queue)
        end
    end

    @testset "Pushing and popping segments" begin
        queue = DT.RefinementQueue(ptri)
        @test !DT.has_segments(queue)
        @test isempty(queue)
        queue[(1, 2)] = 1.0
        @test DT.has_segments(queue)
        @test !isempty(queue)
        @test !DT.has_triangles(queue)
        queue[(5, 7)] = 0.5
        queue[(10, 3)] = 1.881
        queue[(5, 7)] = 0.6
        queue[(7, 5)] = 5.6
        queue[(11, 12)] = 0.0
        queue[(12, 11)] = 0.0
        queue[(27, 25)] = 150.0
        queue[(17, 13)] = 23.5
        queue[(13, 17)] = 109.591
        @test queue.segments == PriorityQueue(Base.Order.Reverse,
            (1, 2) => 1.0,
            (5, 7) => 5.6,
            (10, 3) => 1.881,
            (11, 12) => 0.0,
            (27, 25) => 150.0,
            (17, 13) => 109.591
        )
        @test sprint(show, MIME"text/plain"(), queue) == "RefinementQueue\n   6 segments: DelaunayTriangulation.MaxPriorityQueue((27, 25) => 150.0, (17, 13) => 109.591, (5, 7) => 5.6, (10, 3) => 1.881, (1, 2) => 1.0, (11, 12) => 0.0)\n   0 triangles: DelaunayTriangulation.MaxPriorityQueue{Tuple{Int64, Int64, Int64}, Float64}()"
        res = Pair{NTuple{2,Int},Float64}[]
        while DT.has_segments(queue)
            push!(res, DT.popfirst_segment!(queue))
        end
        @test res == [
            (27, 25) => 150.0,
            (17, 13) => 109.591,
            (5, 7) => 5.6,
            (10, 3) => 1.881,
            (1, 2) => 1.0,
            (11, 12) => 0.0
        ]
        @test isempty(queue)
        @test !DT.has_segments(queue)
    end

    @testset "Pushing and popping triangles" begin
        queue = DT.RefinementQueue(ptri)
        @test !DT.has_triangles(queue)
        @test isempty(queue)
        queue[(1, 2, 3)] = 1.0
        @test DT.has_triangles(queue)
        @test !isempty(queue)
        @test !DT.has_segments(queue)
        queue[(5, 7, 9)] = 0.5
        queue[(9, 5, 7)] = 0.998
        queue[(10, 3, 4)] = 1.881
        queue[(7, 9, 5)] = 18.375
        queue[(10, 15, 20)] = 0.0
        queue[(10, 17, 21)] = 13818.5
        @test queue.triangles == PriorityQueue(Base.Order.Reverse,
            (1, 2, 3) => 1.0,
            (5, 7, 9) => 18.375,
            (10, 3, 4) => 1.881,
            (10, 17, 21) => 13818.5,
            (10, 15, 20) => 0.0
        )
        @test sprint(show, MIME"text/plain"(), queue) == "RefinementQueue\n   0 segments: DelaunayTriangulation.MaxPriorityQueue{Tuple{Int64, Int64}, Float64}()\n   5 triangles: DelaunayTriangulation.MaxPriorityQueue((10, 17, 21) => 13818.5, (5, 7, 9) => 18.375, (10, 3, 4) => 1.881, (1, 2, 3) => 1.0, (10, 15, 20) => 0.0)"
        res = Pair{NTuple{3,Int},Float64}[]
        while DT.has_triangles(queue)
            push!(res, DT.popfirst_triangle!(queue))
        end
        @test res == [
            (10, 17, 21) => 13818.5,
            (5, 7, 9) => 18.375,
            (10, 3, 4) => 1.881,
            (1, 2, 3) => 1.0,
            (10, 15, 20) => 0.0
        ]
        @test isempty(queue)
        @test !DT.has_triangles(queue)
    end

    @testset "Pushing and popping both segments and triangles / haskey" begin
        queue = DT.RefinementQueue(ptri)
        true_triangle_queue = PriorityQueue{NTuple{3,Int},Float64}(Base.Order.Reverse)
        true_segment_queue = PriorityQueue{NTuple{2,Int},Float64}(Base.Order.Reverse)
        for T in each_solid_triangle(ptri)
            _i, _j, _k = triangle_vertices(T)
            ρ = DT.triangle_radius_edge_ratio(get_point(ptri, T...)...)
            @inferred DT.triangle_radius_edge_ratio(get_point(ptri, T...)...)
            @test !haskey(queue, T) && !haskey(queue, (_j, _k, _i)) && !haskey(queue, (_k, _i, _j))
            queue[T] = ρ
            @test haskey(queue, T) && haskey(queue, (_j, _k, _i)) && haskey(queue, (_k, _i, _j))
            true_triangle_queue[DT.sort_triangle(T)] = ρ
            for e in DT.triangle_edges(T)
                ℓ = norm(get_point(ptri, e[1]) .- get_point(ptri, e[2]))
                if !haskey(queue.segments, e) && !haskey(queue.segments, DT.reverse_edge(e))
                    @test !haskey(queue, e) && !haskey(queue, DT.reverse_edge(e))
                end
                queue[e] = ℓ
                @test haskey(queue, e) && haskey(queue, DT.reverse_edge(e))
                η = e[1] < e[2] ? e : DT.reverse_edge(e)
                true_segment_queue[η] = ℓ
            end
        end
        for e in each_solid_edge(ptri)
            ℓ = norm(get_point(ptri, e[1]) .- get_point(ptri, e[2]))
            η = e[1] < e[2] ? e : DT.reverse_edge(e)
            queue[e] = ℓ
            true_segment_queue[η] = ℓ
            w = get_adjacent(ptri, e)
            T = (e..., w)
            ρ = DT.triangle_radius_edge_ratio(get_point(ptri, T...)...)
            queue[T] = ρ
            true_triangle_queue[DT.sort_triangle(T)] = ρ
            w = get_adjacent(ptri, DT.reverse_edge(e))
            T = (e..., w)
            ρ = DT.triangle_radius_edge_ratio(get_point(ptri, T...)...)
            queue[T] = ρ
            true_triangle_queue[DT.sort_triangle(T)] = ρ
        end
        @test length(true_triangle_queue) == length(queue.triangles)
        @test length(true_segment_queue) == length(queue.segments)
        @test !isempty(queue)
        @test DT.has_segments(queue)
        @test DT.has_triangles(queue)
        triangle_res = Pair{NTuple{3,Int},Float64}[]
        segment_res = Pair{NTuple{2,Int},Float64}[]
        true_triangle_res = Pair{NTuple{3,Int},Float64}[]
        true_segment_res = Pair{NTuple{2,Int},Float64}[]
        for i in 1:length(true_triangle_queue)
            T, ρ = DT.popfirst_triangle!(queue)
            push!(triangle_res, DT.sort_triangle(T) => ρ)
            T, ρ = dequeue_pair!(true_triangle_queue)
            push!(true_triangle_res, T => ρ)
        end
        for i in 1:length(true_segment_queue)
            e, ℓ = DT.popfirst_segment!(queue)
            push!(segment_res, (e[1] < e[2] ? e : DT.reverse_edge(e)) => ℓ)
            e, ℓ = dequeue_pair!(true_segment_queue)
            push!(true_segment_res, e => ℓ)
        end
        @test length(true_triangle_queue) == length(queue.triangles) == 0
        @test length(true_segment_queue) == length(queue.segments) == 0
        @test isempty(queue)
        @test !DT.has_segments(queue)
        @test !DT.has_triangles(queue)
        _compare_pairs(triangle_res, true_triangle_res)
        _compare_pairs(segment_res, true_segment_res)
    end
end

@testset verbose = true "InsertionEventHistory" begin
    @testset "Initialisation" begin
        T = NTuple{3,Int}
        E = NTuple{2,Int}
        F = Float64
        history = DT.InsertionEventHistory(ptri)
        @inferred DT.InsertionEventHistory(ptri)
        @test history.added_boundary_segments ⊢ Set{E}()
        @test history.added_triangles ⊢ Set{T}()
        @test history.deleted_segments ⊢ Set{E}()
        @test history.added_segments ⊢ Set{E}()
        @test history.deleted_boundary_segments ⊢ Set{E}()
        @test history.deleted_triangles ⊢ Set{T}()
        @test DT.triangle_type(history) == T
    end

    @testset "Mutating" begin
        T = NTuple{3,Int}
        E = NTuple{2,Int}
        F = Float64
        history = DT.InsertionEventHistory(ptri)
        @test !DT.has_segment_changes(history)
        add_triangle!(history, (5, 7, 1))
        add_triangle!(history, (10, 3, 5))
        delete_triangle!(history, (1, 2, 3))
        delete_triangle!(history, (10, 75, 1))
        DT.add_edge!(history, (10, 5))
        @test DT.has_segment_changes(history)
        DT.add_edge!(history, (73, 9))
        DT.delete_edge!(history, (17, 3))
        DT.delete_edge!(history, (5, 70))
        DT.split_boundary_edge!(history, 183, 190, 50)
        @test DT.each_added_triangle(history) == history.added_triangles
        @test DT.each_added_segment(history) == history.added_segments
        @test DT.each_added_boundary_segment(history) == history.added_boundary_segments
        @test history.added_boundary_segments == Set(((183, 50), (50, 190)))
        @test history.added_triangles == Set(((5, 7, 1), (10, 3, 5)))
        @test history.deleted_segments == Set(((17, 3), (5, 70)))
        @test history.added_segments == Set(((10, 5), (73, 9)))
        @test history.deleted_boundary_segments == Set(((183, 190),))
        @test history.deleted_triangles == Set(((1, 2, 3), (10, 75, 1)))
        empty!(history)
        @test !DT.has_segment_changes(history)
        @test history.added_boundary_segments ⊢ Set{E}()
        @test history.added_triangles ⊢ Set{T}()
        @test history.deleted_segments ⊢ Set{E}()
        @test history.added_segments ⊢ Set{E}()
        @test history.deleted_boundary_segments ⊢ Set{E}()
        @test history.deleted_triangles ⊢ Set{T}()
    end

    @testset verbose = true "Various operations" begin
        @testset "Adding some points inside" begin
            for _ in 1:10
                points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
                tri = triangulate(points; randomise=false)
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                x, y = 0.9, 0.9
                add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(2, 3, 5), (5, 3, 4), (1, 5, 4), (1, 2, 5)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, 2, 3), (1, 3, 4)]))
                @test !DT.has_segment_changes(history)
                @test isempty(history.added_segments)
                @test isempty(history.deleted_segments)
                @test isempty(history.added_boundary_segments)
                @test isempty(history.deleted_boundary_segments)
                @test tri ≠ orig_tri
                validate_insertion_event_history(tri, orig_tri, history)
                DT.undo_insertion!(tri, history)
                @test tri == orig_tri
                add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                orig_tri = deepcopy(tri)
                new_points = [(0.5, 0.5), (0.2, 0.2), (0.5, 0.75), (0.236, 0.987)]
                for (x, y) in new_points
                    empty!(history)
                    add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                    @test tri ≠ orig_tri # make sure == isn't lying for later 
                    validate_insertion_event_history(tri, orig_tri, history)
                    @test !DT.has_segment_changes(history)
                    DT.undo_insertion!(tri, history)
                    validate_statistics(tri)
                    @test validate_triangulation(tri)
                    @test tri == orig_tri
                    add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                    orig_tri = deepcopy(tri)
                end
                @test tri == triangulate(points; randomise=false)
                new_points = [Tuple(rand(2)) for _ in 1:100]
                for (x, y) in new_points
                    empty!(history)
                    add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                    @test tri ≠ orig_tri # make sure == isn't lying for later 
                    validate_insertion_event_history(tri, orig_tri, history)
                    @test !DT.has_segment_changes(history)
                    DT.undo_insertion!(tri, history)
                    validate_statistics(tri)
                    @test validate_triangulation(tri)
                    @test tri == orig_tri
                    add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                    orig_tri = deepcopy(tri)
                end
                @test tri == triangulate(points; randomise=false)
            end
        end

        @testset "Testing how add_point! handles segments" begin
            for _ in 1:10
                points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.0), (0.5, 1.0)]
                tri = triangulate(points; segments=Set(((5, 6),)), randomise=false)
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                x, y = 0.5, 0.5
                add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(1, 5, 7), (7, 5, 2), (7, 2, 3), (6, 7, 3), (4, 7, 6), (4, 1, 7)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, 5, 4), (5, 2, 3), (5, 3, 6), (4, 5, 6)]))
                @test DT.compare_unoriented_edge_collections(history.added_segments, Set([(5, 7), (7, 6)]))
                @test DT.compare_unoriented_edge_collections(history.deleted_segments, Set([(5, 6)]))
                @test isempty(history.added_boundary_segments)
                @test isempty(history.deleted_boundary_segments)
                @test DT.has_segment_changes(history)
                @test tri ≠ orig_tri
                validate_insertion_event_history(tri, orig_tri, history)
                DT.undo_insertion!(tri, history)
                @test tri == orig_tri
                add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                orig_tri = deepcopy(tri)
                new_points = LinRange(0.0001, 0.9999, 25) |> collect
                setdiff!(new_points, 0.5)
                x = 0.5
                for y in new_points
                    empty!(history)
                    add_point!(tri, x, y, store_event_history=Val(true), event_history=history)
                    @test tri ≠ orig_tri # make sure == isn't lying for later
                    validate_insertion_event_history(tri, orig_tri, history)
                    @test DT.has_segment_changes(history)
                    DT.undo_insertion!(tri, history)
                    @test tri == orig_tri
                    validate_statistics(tri)
                    @test validate_triangulation(tri)
                    add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                    orig_tri = deepcopy(tri)
                end
            end
        end

        @testset "Testing how add_point! handles points on existing boundary segments" begin
            for PT in (DT.Exact, DT.Adaptive)
                for _ in 1:5
                    points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
                    tri = triangulate(points; boundary_nodes=[1, 2, 3, 4, 1], randomise=false, predicates=PT())
                    orig_tri = deepcopy(tri)
                    history = DT.InsertionEventHistory(tri)
                    x, y = 0.5, 0.0
                    add_point!(tri, (x, y), store_event_history=Val(true), event_history=history, predicates=PT())
                    @test validate_triangulation(tri, predicates=PT())
                    @test DT.compare_triangle_collections(history.added_triangles, Set([(1, 5, 4), (5, 2, 3), (4, 5, 3), (5, 1, -1), (2, 5, -1)]))
                    @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, 2, 3), (1, 3, 4), (2, 1, -1)]))
                    @test isempty(history.added_segments)
                    @test isempty(history.deleted_segments)
                    @test DT.has_segment_changes(history)
                    @test DT.compare_unoriented_edge_collections(history.added_boundary_segments, Set([(1, 5), (5, 2)]))
                    @test DT.compare_unoriented_edge_collections(history.deleted_boundary_segments, Set([(1, 2)]))
                    @test tri ≠ orig_tri
                    validate_insertion_event_history(tri, orig_tri, history)
                    DT.undo_insertion!(tri, history)
                    validate_statistics(tri)
                    @test validate_triangulation(tri)
                    @test tri == orig_tri
                    add_point!(tri, (x, y), store_event_history=Val(true), event_history=history, predicates=PT())
                    orig_tri = deepcopy(tri)
                    new_points = LinRange(0.0001, 0.9999, 10) |> collect
                    setdiff!(new_points, 0.5)
                    @test validate_triangulation(tri)
                    for r in new_points
                        for p in ((r, 0.0), (1.0, r), (r, 1.0), (0.0, r))
                            empty!(history)
                            add_point!(tri, p, store_event_history=Val(true), event_history=history, predicates=PT())
                            @test tri ≠ orig_tri # make sure == isn't lying for later
                            validate_insertion_event_history(tri, orig_tri, history)
                            @test DT.has_segment_changes(history)
                            @test tri.boundary_edge_map ≠ orig_tri.boundary_edge_map
                            DT.undo_insertion!(tri, history)
                            @test tri.boundary_edge_map == orig_tri.boundary_edge_map
                            validate_statistics(tri)
                            @test validate_triangulation(tri)
                            @test tri == orig_tri
                            add_point!(tri, p, store_event_history=Val(true), event_history=history, predicates=PT())
                            orig_tri = deepcopy(tri)
                        end
                    end
                end
            end
        end

        @testset "Testing how add_point! handles points on existing boundary segments in a multiply-connected domain" begin
            for _ in 1:10
                points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.3, 0.3), (0.7, 0.3), (0.5, 0.7)]
                boundary_nodes = [[[1, 2], [2, 3], [3, 4], [4, 1]], [[5, 7], [7, 6], [6, 5]]]
                tri = triangulate(points; boundary_nodes, randomise=false)
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                x, y = (0.5, 0.3)
                add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                @test validate_triangulation(tri)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(1, 8, 5), (1, 2, 8), (8, 2, 6), (8, 6, -7), (5, 8, -7)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, 6, 5), (1, 2, 6), (5, 6, -7)]))
                @test isempty(history.added_segments)
                @test isempty(history.deleted_segments)
                @test DT.compare_unoriented_edge_collections(history.added_boundary_segments, Set([(5, 8), (8, 6)]))
                @test DT.compare_unoriented_edge_collections(history.deleted_boundary_segments, Set([(5, 6)]))
                @test DT.has_segment_changes(history)
                @test tri ≠ orig_tri
                validate_insertion_event_history(tri, orig_tri, history)
                DT.undo_insertion!(tri, history)
                validate_statistics(tri)
                @test validate_triangulation(tri)
                @test tri == orig_tri
                add_point!(tri, (x, y), store_event_history=Val(true), event_history=history)
                orig_tri = deepcopy(tri)
                new_points = LinRange(0.31, 0.69, 15) |> collect
                @test validate_triangulation(tri)
                setdiff!(new_points, 0.5)
                for r in new_points
                    for p in ((r, 0.3), (r, 0.0))
                        empty!(history)
                        add_point!(tri, p, store_event_history=Val(true), event_history=history)
                        @test tri ≠ orig_tri # make sure == isn't lying for later
                        validate_insertion_event_history(tri, orig_tri, history)
                        @test DT.has_segment_changes(history)
                        @test tri.boundary_edge_map ≠ orig_tri.boundary_edge_map
                        DT.undo_insertion!(tri, history)
                        @test tri.boundary_edge_map == orig_tri.boundary_edge_map
                        validate_statistics(tri)
                        @test validate_triangulation(tri)
                        @test tri == orig_tri
                        add_point!(tri, p, store_event_history=Val(true), event_history=history)
                        orig_tri = deepcopy(tri)
                    end
                end
            end
        end

        @testset "Testing how split_edge! is handled on non-segments" begin
            for _ in 1:10
                points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.25, 0.25), (0.75, 0.75)]
                tri = triangulate(points; randomise=false)
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                x, y = (0.5, 0.5)
                push!(points, (x, y))
                split_edge!(tri, 5, 6, DT.num_points(tri), Val(true), history)
                validate_insertion_event_history(tri, orig_tri, history)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(5, 7, 4), (7, 6, 4)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(5, 6, 4)]))
                @test isempty(history.added_segments)
                @test isempty(history.deleted_segments)
                @test isempty(history.added_boundary_segments)
                @test isempty(history.deleted_boundary_segments)
                @test !DT.has_segment_changes(history)
                @test tri ≠ orig_tri
                DT.undo_insertion!(tri, history)
                @test tri == orig_tri
                validate_statistics(tri)
                @test validate_triangulation(tri) # don't test this later since split_edge! doesn't restore Delaunay-hood
                push!(points, (x, y))
                split_edge!(tri, 5, 6, DT.num_points(tri), Val(true), history)
                new_points = LinRange(0.26, 0.49, 25) |> collect
                orig_tri = deepcopy(tri)
                for r in new_points
                    p = (r, r)
                    empty!(history)
                    push!(points, p)
                    split_edge!(tri, 5, DT.num_points(tri) - 1, DT.num_points(tri), Val(true), history)
                    validate_insertion_event_history(tri, orig_tri, history)
                    @test DT.compare_triangle_collections(history.added_triangles, Set([(5, DT.num_points(tri), 4), (DT.num_points(tri), DT.num_points(tri) - 1, 4)]))
                    @test DT.compare_triangle_collections(history.deleted_triangles, Set([(5, DT.num_points(tri) - 1, 4)]))
                    @test isempty(history.added_segments)
                    @test isempty(history.deleted_segments)
                    @test isempty(history.added_boundary_segments)
                    @test isempty(history.deleted_boundary_segments)
                    @test !DT.has_segment_changes(history)
                    @test tri ≠ orig_tri
                    DT.undo_insertion!(tri, history)
                    @test tri == orig_tri
                    push!(points, p)
                    split_edge!(tri, 5, DT.num_points(tri) - 1, DT.num_points(tri), Val(true), history)
                    orig_tri = deepcopy(tri)
                end
            end
        end

        @testset "Testing how split_edge! is handled on segments" begin
            for _ in 1:10
                points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
                tri = triangulate(points; randomise=false, segments=Set([(1, 3)]))
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                x, y = (0.5, 0.5)
                push!(points, (x, y))
                split_edge!(tri, 1, 3, DT.num_points(tri), Val(true), history)
                validate_insertion_event_history(tri, orig_tri, history)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(1, 5, 4), (5, 3, 4)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, 3, 4)]))
                @test DT.compare_unoriented_edge_collections(history.added_segments, Set([(1, 5), (5, 3)]))
                @test DT.compare_unoriented_edge_collections(history.deleted_segments, Set([(1, 3)]))
                @test isempty(history.added_boundary_segments)
                @test isempty(history.deleted_boundary_segments)
                @test DT.has_segment_changes(history)
                @test tri ≠ orig_tri
                DT.undo_insertion!(tri, history)
                @test tri == orig_tri
                validate_statistics(tri)
                @test validate_triangulation(tri) # don't test this later since split_edge! doesn't restore Delaunay-hood
                push!(points, (x, y))
                split_edge!(tri, 1, 3, DT.num_points(tri), Val(true), history)
                new_points = LinRange(0.01, 0.49, 25) |> collect
                orig_tri = deepcopy(tri)
                for r in new_points
                    p = (r, r)
                    empty!(history)
                    push!(points, p)
                    split_edge!(tri, 1, DT.num_points(tri) - 1, DT.num_points(tri), Val(true), history)
                    validate_insertion_event_history(tri, orig_tri, history)
                    @test DT.compare_triangle_collections(history.added_triangles, Set([(1, DT.num_points(tri), 4), (DT.num_points(tri), DT.num_points(tri) - 1, 4)]))
                    @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, DT.num_points(tri) - 1, 4)]))
                    @test DT.compare_unoriented_edge_collections(history.added_segments, Set([(1, DT.num_points(tri)), (DT.num_points(tri), DT.num_points(tri) - 1)]))
                    @test DT.compare_unoriented_edge_collections(history.deleted_segments, Set([(1, DT.num_points(tri) - 1)]))
                    @test isempty(history.added_boundary_segments)
                    @test isempty(history.deleted_boundary_segments)
                    @test DT.has_segment_changes(history)
                    @test tri ≠ orig_tri
                    DT.undo_insertion!(tri, history)
                    @test tri == orig_tri
                    push!(points, p)
                    split_edge!(tri, 1, DT.num_points(tri) - 1, DT.num_points(tri), Val(true), history)
                    orig_tri = deepcopy(tri)
                end
            end
        end

        @testset "Testing how split_edge! is handled on boundary segments" begin
            for _ in 1:10
                points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
                tri = triangulate(points; boundary_nodes=[1, 2, 3, 4, 1], randomise=false)
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                x, y = (0.5, 0.0)
                push!(points, (x, y))
                split_edge!(tri, 1, 2, DT.num_points(tri), Val(true), history)
                validate_insertion_event_history(tri, orig_tri, history)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(1, DT.num_points(tri), 3), (DT.num_points(tri), 2, 3)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, 2, 3)]))
                @test isempty(history.added_segments)
                @test isempty(history.deleted_segments)
                @test DT.compare_unoriented_edge_collections(history.added_boundary_segments, Set([(1, 5), (5, 2)]))
                @test DT.compare_unoriented_edge_collections(history.deleted_boundary_segments, Set([(1, 2)]))
                @test DT.has_segment_changes(history)
                @test tri ≠ orig_tri
                DT.undo_insertion!(tri, history)
                @test tri == orig_tri
                validate_statistics(tri)
                @test validate_triangulation(tri)
                push!(points, (x, y))
                split_edge!(tri, 1, 2, DT.num_points(tri), Val(true), history)
                new_points = LinRange(0.01, 0.99, 25) |> collect
                setdiff!(new_points, 0.5)
                orig_tri = deepcopy(tri)
                for r in new_points
                    p = (r, 0.0)
                    empty!(history)
                    push!(points, p)
                    split_edge!(tri, 1, DT.num_points(tri) - 1, DT.num_points(tri), Val(true), history)
                    validate_insertion_event_history(tri, orig_tri, history)
                    @test DT.compare_triangle_collections(history.added_triangles, Set([(1, DT.num_points(tri), 3), (DT.num_points(tri), DT.num_points(tri) - 1, 3)]))
                    @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, DT.num_points(tri) - 1, 3)]))
                    @test isempty(history.added_segments)
                    @test isempty(history.deleted_segments)
                    @test DT.compare_unoriented_edge_collections(history.added_boundary_segments, Set([(1, DT.num_points(tri)), (DT.num_points(tri), DT.num_points(tri) - 1)]))
                    @test DT.compare_unoriented_edge_collections(history.deleted_boundary_segments, Set([(1, DT.num_points(tri) - 1)]))
                    @test DT.has_segment_changes(history)
                    @test tri ≠ orig_tri
                    DT.undo_insertion!(tri, history)
                    @test tri == orig_tri
                    push!(points, p)
                    split_edge!(tri, 1, DT.num_points(tri) - 1, DT.num_points(tri), Val(true), history)
                    orig_tri = deepcopy(tri)
                end
            end
        end

        @testset "Testing how flip_edge! is handled" begin
            for _ in 1:10
                points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
                tri = triangulate(points; randomise=false)
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                flip_edge!(tri, 1, 3, Val(true), history)
                validate_insertion_event_history(tri, orig_tri, history)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(1, 2, 4), (4, 2, 3)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, 2, 3), (4, 1, 3)]))
                @test isempty(history.added_segments)
                @test isempty(history.deleted_segments)
                @test isempty(history.added_boundary_segments)
                @test isempty(history.deleted_boundary_segments)
                @test !DT.has_segment_changes(history)
                @test tri ≠ orig_tri
                DT.undo_insertion!(tri, history, Val(false))
                @test tri == orig_tri
                validate_statistics(tri)
                @test validate_triangulation(tri)
            end
        end

        @testset "Testing how complete_split_edge_and_legalise! is handled for non-segments" begin
            for _ in 1:100
                points = [(0.01, 0.0), (0.99, 0.00027), (1.0, 1.00025), (0.0, 0.994)]
                tri = triangulate(points; randomise=false)
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                push!(points, (0.5, 0.5))
                DT.complete_split_edge_and_legalise!(tri, 2, 4, 5, Val(true), history)
                validate_insertion_event_history(tri, orig_tri, history)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(1, 2, 5), (5, 2, 3), (4, 5, 3), (4, 1, 5)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, 2, 4), (4, 2, 3)]))
                @test isempty(history.added_segments)
                @test isempty(history.deleted_segments)
                @test isempty(history.added_boundary_segments)
                @test isempty(history.deleted_boundary_segments)
                @test !DT.has_segment_changes(history)
                @test tri ≠ orig_tri
                DT.undo_insertion!(tri, history)
                @test tri == orig_tri
                validate_statistics(tri)
                @test validate_triangulation(tri)
                rr = LinRange(0.51, 0.99, 20) |> collect
                push!(points, (0.5, 0.5))
                DT.complete_split_edge_and_legalise!(tri, 2, 4, 5, Val(true), history)
                @test validate_triangulation(tri)
                validate_statistics(tri)
                orig_tri = deepcopy(tri)
                for r in rr
                    empty!(history)
                    push!(points, (r, r))
                    DT.complete_split_edge_and_legalise!(tri, DT.num_points(tri) - 1, 3, DT.num_points(tri), Val(true), history)
                    validate_insertion_event_history(tri, orig_tri, history)
                    @test isempty(history.added_segments)
                    @test isempty(history.deleted_segments)
                    @test isempty(history.added_boundary_segments)
                    @test isempty(history.deleted_boundary_segments)
                    @test !DT.has_segment_changes(history)
                    @test tri ≠ orig_tri
                    @test validate_triangulation(tri)
                    DT.undo_insertion!(tri, history)
                    @test tri == orig_tri
                    validate_statistics(tri)
                    @test validate_triangulation(tri)
                    push!(points, (r, r))
                    DT.complete_split_edge_and_legalise!(tri, DT.num_points(tri) - 1, 3, DT.num_points(tri), Val(true), history)
                    orig_tri = deepcopy(tri)
                end
                @test tri == triangulate(points; randomise=false)
            end
        end

        @testset "Testing how complete_split_edge_and_legalise! is handled for segments" begin
            for _ in 1:10
                points = [(0.01, 0.0), (0.99, 0.00027), (1.0, 1.00025), (0.0, 0.994), (0.5, 0.1), (0.5, 0.9)]
                tri = triangulate(points; randomise=false, segments=Set([(5, 6)]))
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                push!(points, (0.5, 0.5))
                DT.complete_split_edge_and_legalise!(tri, 5, 6, 7, Val(true), history)
                validate_insertion_event_history(tri, orig_tri, history)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(7, 5, 2), (7, 1, 5), (7, 6, 4), (7, 2, 3), (7, 4, 1), (7, 3, 6)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(2, 3, 6), (5, 6, 4), (6, 5, 2), (4, 1, 5)]))
                @test DT.compare_unoriented_edge_collections(history.deleted_segments, Set([(5, 6)]))
                @test DT.compare_unoriented_edge_collections(history.added_segments, Set([(5, 7), (7, 6)]))
                @test isempty(history.added_boundary_segments)
                @test isempty(history.deleted_boundary_segments)
                @test DT.has_segment_changes(history)
                @test tri ≠ orig_tri
                DT.undo_insertion!(tri, history)
                @test tri == orig_tri
                rr = LinRange(0.51, 0.89, 25) |> collect
                push!(points, (0.5, 0.5))
                DT.complete_split_edge_and_legalise!(tri, 5, 6, 7, Val(true), history)
                @test validate_triangulation(tri)
                validate_statistics(tri)
                orig_tri = deepcopy(tri)
                for r in rr
                    empty!(history)
                    push!(points, (0.5, r))
                    DT.complete_split_edge_and_legalise!(tri, DT.num_points(tri) - 1, 6, DT.num_points(tri), Val(true), history)
                    validate_insertion_event_history(tri, orig_tri, history)
                    @test !isempty(history.added_segments)
                    @test !isempty(history.deleted_segments)
                    @test isempty(history.added_boundary_segments)
                    @test isempty(history.deleted_boundary_segments)
                    @test DT.has_segment_changes(history)
                    @test tri ≠ orig_tri
                    @test validate_triangulation(tri)
                    DT.undo_insertion!(tri, history)
                    @test tri == orig_tri
                    validate_statistics(tri)
                    @test validate_triangulation(tri)
                    push!(points, (0.5, r))
                    DT.complete_split_edge_and_legalise!(tri, DT.num_points(tri) - 1, 6, DT.num_points(tri), Val(true), history)
                    orig_tri = deepcopy(tri)
                end
            end
        end

        @testset "Testing how complete_split_edge_and_legalise! is handled for boundary segments" begin
            for _ in 1:10
                points = [(0.01, 0.0), (0.99, 0.00027), (1.0, 1.00025), (0.0, 0.994), (0.5, 0.1), (0.5, 0.9)]
                tri = triangulate(points; randomise=false, boundary_nodes=[1, 2, 3, 4, 1])
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                push!(points, (0.1, 0.0))
                DT.complete_split_edge_and_legalise!(tri, 1, 2, 7, Val(true), history)
                validate_insertion_event_history(tri, orig_tri, history)
                @test DT.compare_triangle_collections(history.added_triangles, Set([(2, 7, -1), (7, 5, 6), (7, 6, 4), (7, 2, 5), (7, 1, -1), (7, 4, 1)]))
                @test DT.compare_triangle_collections(history.deleted_triangles, Set([(1, 2, 5), (2, 1, -1), (5, 6, 4), (5, 4, 1)]))
                @test isempty(history.added_segments)
                @test isempty(history.deleted_segments)
                @test DT.compare_unoriented_edge_collections(history.deleted_boundary_segments, Set([(1, 2)]))
                @test DT.compare_unoriented_edge_collections(history.added_boundary_segments, Set([(1, 7), (7, 2)]))
                @test DT.has_segment_changes(history)
                @test tri ≠ orig_tri
                DT.undo_insertion!(tri, history)
                @test tri == orig_tri
                rr = LinRange(0.2, 0.8, 25) |> collect
                push!(points, (0.1, 0.0))
                DT.complete_split_edge_and_legalise!(tri, 1, 2, 7, Val(true), history)
                @test validate_triangulation(tri)
                validate_statistics(tri)
                orig_tri = deepcopy(tri)
                for r in rr
                    empty!(history)
                    push!(points, (r, 0.0))
                    DT.complete_split_edge_and_legalise!(tri, DT.num_points(tri) - 1, 2, DT.num_points(tri), Val(true), history)
                    validate_insertion_event_history(tri, orig_tri, history)
                    @test isempty(history.added_segments)
                    @test isempty(history.deleted_segments)
                    @test !isempty(history.added_boundary_segments)
                    @test !isempty(history.deleted_boundary_segments)
                    @test DT.has_segment_changes(history)
                    @test tri ≠ orig_tri
                    @test validate_triangulation(tri)
                    DT.undo_insertion!(tri, history)
                    @test tri == orig_tri
                    validate_statistics(tri)
                    @test validate_triangulation(tri)
                    push!(points, (r, 0.0))
                    DT.complete_split_edge_and_legalise!(tri, DT.num_points(tri) - 1, 2, DT.num_points(tri), Val(true), history)
                    orig_tri = deepcopy(tri)
                end
                convex_hull!(tri)
                @test tri == triangulate(points; randomise=false, boundary_nodes=tri.boundary_nodes)
            end
        end

        @testset "Testing delete_point!" begin
            for i in 1:2500
                points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5), (0.2, 0.8), (0.1, 0.785)]
                tri = triangulate(points)
                orig_tri = deepcopy(tri)
                history = DT.InsertionEventHistory(tri)
                delete_point!(tri, 5, store_event_history=Val(true), event_history=history)
                validate_insertion_event_history(tri, orig_tri, history)
                DT.undo_insertion!(tri, history, Val(false))
                validate_statistics(tri)
                @test validate_triangulation(tri)
                @test tri == orig_tri
            end

            for j in 1:250
                rng = StableRNG(j)
                points = rand(rng, 2, 500)
                tri = triangulate(points; rng)
                for k in 1:20
                    rng = StableRNG(k)
                    i = rand(rng, each_solid_vertex(tri))
                    while DT.is_boundary_node(tri, i)[1]
                        i = rand(rng, each_solid_vertex(tri))
                    end
                    orig_tri = deepcopy(tri)
                    history = DT.InsertionEventHistory(tri)
                    delete_point!(tri, i, rng=rng, store_event_history=Val(true), event_history=history)
                    @test validate_triangulation(tri)
                    validate_insertion_event_history(tri, orig_tri, history)
                    DT.undo_insertion!(tri, history, Val(false))
                    @test tri == orig_tri
                end
            end
        end
    end

    @testset "A smaller split_edge! example" begin
        p1 = (0.0, 0.0)
        p2 = (1.0, 0.0)
        p3 = (0.0, 1.0)
        p4 = (1.0, 1.0)
        pts = [p1, p2, p3, p4]
        rng = StableRNG(123)
        tri = triangulate(pts; boundary_nodes=[1, 2, 4, 3, 1], rng, delete_ghosts=false)
        _tri = deepcopy(tri)
        push!(pts, (0.5, 0.5))
        events = DT.InsertionEventHistory(tri)
        DT.complete_split_edge_and_legalise!(tri, 2, 3, 5, Val(true), events)
        @test DT.compare_triangle_collections(events.added_triangles, Set(((2, 5, 1), (3, 5, 4), (5, 3, 1), (5, 2, 4))))
        @test DT.compare_triangle_collections(events.deleted_triangles, Set(((3, 2, 4), (2, 3, 1))))
        for T in events.added_triangles
            @test DT.contains_triangle(tri, T)[2]
            @test !DT.contains_triangle(_tri, T)[2]
        end
        for T in events.deleted_triangles
            @test !DT.contains_triangle(tri, T)[2]
            @test DT.contains_triangle(_tri, T)[2]
        end
        _tri = deepcopy(tri)
        push!(pts, (0.5, 0.0))
        events = DT.InsertionEventHistory(tri)
        DT.complete_split_edge_and_legalise!(tri, 1, 2, 6, Val(true), events)
        @test compare_edge_vectors(collect(events.added_boundary_segments), [(6, 2), (1, 6)])
        @test DT.compare_triangle_collections(events.added_triangles, Set(((2, 6, -1), (6, 2, 5), (1, 6, 5), (6, 1, -1))))
        @test isempty(events.deleted_segments)
        @test isempty(events.added_segments)
        @test compare_edge_vectors(collect(events.deleted_boundary_segments), [(1, 2)])
        @test DT.compare_triangle_collections(events.deleted_triangles, Set(((1, 2, 5), (2, 1, -1))))
        for T in events.added_triangles
            @test DT.contains_triangle(tri, T)[2]
            @test !DT.contains_triangle(_tri, T)[2]
        end
        for T in events.deleted_triangles
            @test !DT.contains_triangle(tri, T)[2]
            @test DT.contains_triangle(_tri, T)[2]
        end
    end
end

@testset verbose = true "RefinementArguments" begin
    @testset "Defaults" begin
        tri = triangulate(rand(2, 50))
        args = DT.RefinementArguments(tri)
        @inferred DT.RefinementArguments(tri)
        @test args.constraints.min_angle == 30.0
        @test args.constraints.max_angle == 180.0
        @test args.constraints.min_area == get_area(tri) / 1e9
        @test args.constraints.max_area == Inf
        @test args.constraints.max_points == 1000^2
        @test args.constraints.seditious_angle == 20.0
        @test !args.constraints.custom_constraint(tri, (1, 2, 3))
        @test args.use_circumcenter
        @test args.use_lens
        @test args.steiner_scale == 0.999
        @test args.rng == Random.default_rng()
        @test args.queue == DT.RefinementQueue(tri)
        @test args.events == DT.InsertionEventHistory(tri)
        @test args.locked_convex_hull
        @test args.had_ghosts
        @test DT.has_ghost_triangles(tri)
        @test DT.has_boundary_nodes(tri) # convex hull 
        @test args.min_steiner_vertex == 51
        segment_list = Set{NTuple{2,Int}}()
        segment_vertices = Set{Int}()
        for i in 1:(length(tri.convex_hull.vertices)-1)
            u = tri.convex_hull.vertices[i]
            v = tri.convex_hull.vertices[i+1]
            push!(segment_list, (u, v))
            push!(segment_vertices, u, v)
            @test !DT.is_subsegment(args, (u, v)) && !DT.is_subsegment(args, (v, u))
            @test DT.is_segment_vertex(args, u) && DT.is_segment_vertex(args, v)
        end
        @test args.segment_list ⊢ segment_list
        @test args.segment_vertices ⊢ segment_vertices
        @test args.midpoint_split_list ⊢ Set{Int}()
        @test args.offcenter_split_list ⊢ Set{Int}()
        @inferred DT.RefinementArguments(tri)
        @test !args.concavity_protection
    end

    @testset "Setting values" begin
        tri = triangulate(rand(2, 50))
        args = DT.RefinementArguments(tri;
            min_angle=30.5,
            max_angle=90.5,
            min_area=0.2,
            max_area=0.9,
            max_points=573,
            custom_constraint=(tri, T) -> T == (1, 2, 3),
            use_circumcenter=true,
            use_lens=false,
            steiner_scale=0.95,
            seditious_angle=25.3,
            rng=StableRNG(123),
            concavity_protection=true
        )
        @test args.constraints.min_angle == 30.5
        @test args.constraints.max_angle == 90.5
        @test args.constraints.min_area == 0.2
        @test args.constraints.max_area == 0.9
        @test args.constraints.max_points == 573
        @test args.constraints.seditious_angle == 25.3
        @test args.constraints.custom_constraint(tri, (1, 2, 3))
        @test !args.constraints.custom_constraint(tri, (1, 2, 4))
        @test args.use_circumcenter
        @test !args.use_lens
        @test args.steiner_scale == 0.95
        @test args.rng == StableRNG(123)
        @test args.queue == DT.RefinementQueue(tri)
        @test args.events == DT.InsertionEventHistory(tri)
        @test args.locked_convex_hull
        @test args.had_ghosts
        @test DT.has_ghost_triangles(tri)
        @test DT.has_boundary_nodes(tri) # convex hull
        @test args.min_steiner_vertex == 51
        segment_list = Set{NTuple{2,Int}}()
        segment_vertices = Set{Int}()
        for i in 1:(length(tri.convex_hull.vertices)-1)
            u = tri.convex_hull.vertices[i]
            v = tri.convex_hull.vertices[i+1]
            push!(segment_list, (u, v))
            push!(segment_vertices, u, v)
            @test !DT.is_subsegment(args, (u, v)) && !DT.is_subsegment(args, (v, u))
            @test DT.is_segment_vertex(args, u) && DT.is_segment_vertex(args, v)
        end
        @test args.segment_list ⊢ segment_list
        @test args.segment_vertices ⊢ segment_vertices
        @test args.midpoint_split_list ⊢ Set{Int}()
        @test args.offcenter_split_list ⊢ Set{Int}()
        @test args.concavity_protection
        @inferred DT.RefinementArguments(tri;
            min_angle=30.5,
            max_angle=90.5,
            min_area=0.2,
            max_area=0.9,
            max_points=573,
            custom_constraint=(tri, T) -> T == (1, 2, 3),
            use_circumcenter=true,
            use_lens=false,
            steiner_scale=0.95,
            rng=StableRNG(123),
            concavity_protection=true
        )
    end

    @testset "Adding ghost triangles" begin
        points = rand(2, 50)
        tri = triangulate(points, delete_ghosts=true)
        args = DT.RefinementArguments(tri)
        @test !args.had_ghosts
        @test DT.has_ghost_triangles(tri)
        @test tri == triangulate(points, delete_ghosts=false, boundary_nodes=tri.convex_hull.vertices)
    end

    @testset "Checking free vertices" begin
        tri = triangulate(rand(2, 50))
        args = DT.RefinementArguments(tri)
        @test DT.is_free(args, 51)
        @test !DT.is_free(args, 50)
        @test all(i -> !DT.is_free(args, i), 1:50)
        @test DT.is_free(args, 52)
        push!(args.midpoint_split_list, 52)
        @test !DT.is_free(args, 52)
        @test DT.is_free(args, 101)
        push!(args.offcenter_split_list, 101)
        @test !DT.is_free(args, 101)
    end

    @testset "Checking midpoint and offcenter splits" begin
        tri = triangulate(rand(2, 50))
        args = DT.RefinementArguments(tri)
        @test !DT.is_midpoint_split(args, 57)
        push!(args.midpoint_split_list, 57)
        @test DT.is_midpoint_split(args, 57)
        @test !DT.is_offcenter_split(args, 101)
        push!(args.offcenter_split_list, 101)
        @test DT.is_offcenter_split(args, 101)
    end

    @testset "Iteration" begin
        tri = triangulate([(rand(), rand()) for _ in 1:50])
        args = DT.RefinementArguments(tri; max_points=55)
        DT.unlock_convex_hull!(tri) # if the random points later are added outside the domain, they don't always get added as a vertex
        @test !DT.keep_iterating(tri, args)
        args.queue[(1, 2, 3)] = 1.0
        @test DT.keep_iterating(tri, args)
        @test !DT.keep_splitting(tri, args)
        args.queue[(1, 2)] = 1.0
        @test DT.keep_iterating(tri, args)
        @test DT.keep_splitting(tri, args)
        popfirst!(args.queue.triangles)
        @test DT.keep_iterating(tri, args)
        @test DT.keep_splitting(tri, args)
        popfirst!(args.queue.segments)
        @test !DT.keep_iterating(tri, args)
        @test !DT.keep_splitting(tri, args)
        args.queue[(1, 2, 3)] = 1.0
        args.queue[(1, 2)] = 1.0
        for _ in 1:4
            p = (rand(), rand())
            add_point!(tri, p)
        end # 54 < 55 still 
        @test DT.keep_iterating(tri, args)
        @test DT.keep_splitting(tri, args)
        add_point!(tri, (rand(), rand())) # now == 55
        @test !DT.keep_iterating(tri, args)
        @test !DT.keep_splitting(tri, args)
        add_point!(tri, (rand(), rand())) # make sure we aren't checking == 
        @test !DT.keep_iterating(tri, args)
        @test !DT.keep_splitting(tri, args)
    end
end

@testset verbose = true "is_encroached/encroaches_upon" begin
    @testset "Testing if point is in diametral circle" begin
        args1 = DT.RefinementArguments(ptri; use_lens=false)
        @inferred DT.RefinementArguments(ptri; use_lens=false)
        args2 = DT.RefinementArguments(ptri; use_lens=true, min_angle=45.0)
        p = Vector{NTuple{2,Float64}}(undef, 7)
        q = Vector{NTuple{2,Float64}}(undef, 7)
        r = Vector{NTuple{2,Float64}}(undef, 7)

        p[1], q[1], r[1] = (2.6, 1.63), (5.68, 1.37), (4.5, 4.63)
        p[2], q[2], r[2] = p[1], q[1], (4.18, 1.96)
        p[3], q[3], r[3] = p[2], q[2], (1.5, 2.0)
        p[4], q[4], r[4] = p[3], q[3], (3.9, 1.324)
        p[5], q[5], r[5] = p[4], q[4], (3.6, -1.25)
        p[6], q[6], r[6] = (2.0, 0.0), (4.0, 0.0), (4.0, 2.0)
        p[7], q[7], r[7] = (2.0, 0.0), (4.0, 2.0), (4.0, 0.0)

        certs = [
            DT.is_outside,
            DT.is_inside,
            DT.is_outside,
            DT.is_inside,
            DT.is_outside,
            DT.is_outside,
            DT.is_on,
        ]

        for (p, q, r, c) in zip(p, q, r, certs)
            @test c(DT.point_position_relative_to_diametral_circle(p, q, r))
            @test c(DT.point_position_relative_to_diametral_circle(q, p, r))
            if c == DT.is_on || c == DT.is_inside
                @test DT.encroaches_upon(p, q, r, args1)
                @inferred DT.encroaches_upon(p, q, r, args1)
                @test DT.encroaches_upon(p, q, r, args2)
            else
                @test !DT.encroaches_upon(p, q, r, args1)
                @test !DT.encroaches_upon(p, q, r, args2)
            end
        end
    end

    @testset "Testing if a point is in a diametral lens" begin
        args = DT.RefinementArguments(ptri; use_lens=true, min_angle=30.0)
        p = (-7.0, 4.0)
        q = (-2.0, 4.0)
        points = Vector{Tuple{NTuple{2,Float64},Function}}(undef, 0)
        push!(points, ((-5.5, 5.0), DT.is_inside))
        push!(points, ((-4.0, 4.5), DT.is_inside))
        push!(points, ((-5.0, 4.5), DT.is_inside))
        push!(points, ((-4.5, 5.0), DT.is_inside))
        push!(points, ((-4.0, 7.0), DT.is_outside))
        push!(points, ((-6.5, 6.0), DT.is_outside))
        push!(points, ((-6.0, 3.5), DT.is_inside))
        push!(points, ((-4.0, 3.0), DT.is_inside))
        push!(points, ((-4.5, 2.0), DT.is_outside))
        push!(points, ((-6.0, 2.0), DT.is_outside))
        push!(points, ((-3.0, 1.0), DT.is_outside))
        push!(points, ((-4.5, 3.0), DT.is_inside))
        push!(points, ((-5.0, 4.0), DT.is_inside))
        push!(points, ((-3.0, 4.0), DT.is_inside))
        push!(points, (p, DT.is_on))
        push!(points, (q, DT.is_on))
        for (pp, c) in points
            @test c(DT.point_position_relative_to_diametral_lens(p, q, pp, 30.0))
            @inferred c(DT.point_position_relative_to_diametral_lens(p, q, pp, 30.0))
            @inferred DT.point_position_relative_to_diametral_lens(p, q, pp, 30.0)
            if c == DT.is_on || c == DT.is_inside
                @test DT.encroaches_upon(p, q, pp, args)
            else
                @test !DT.encroaches_upon(p, q, pp, args)
            end
        end

        p, q = (-7.0, 4.0), (-2.0, 4.0)
        upper_circle, lower_circle, circle = compute_diametral_circle(p, q)
        points_in_diametral_circle, points_outside_diametral_circle = get_points_in_diametral_circle(p, q)
        fig = Figure(fontsize=43)
        ax = Axis(fig[1, 1], width=600, height=600, title="Diametral circle", titlealign=:left)
        scatter!(ax, points_in_diametral_circle, markersize=4, color=:blue)
        scatter!(ax, points_outside_diametral_circle, markersize=4, color=:red)
        #band!(ax, lower_circle, upper_circle, color=(:green, 0.4))
        lines!(ax, circle, color=:black, linewidth=7)
        lines!(ax, [p, q], color=:black, linewidth=7)
        xlims!(ax, -7.5, -1.5)
        ylims!(ax, 1, 7)

        for (i, lens_angle) in enumerate((45.0, 30.0, 20.0, 10.0))
            upper_lens, lower_lens, lens = compute_diametral_lens(p, q, lens_angle)
            points_in_diametral_lens, points_outside_diametral_lens = get_points_in_diametral_lens(p, q, lens_angle)
            ax = Axis(fig[1, 1+i], width=600, height=600, title="Diametral lens (angle $(lens_angle)°)", titlealign=:left)
            scatter!(ax, points_in_diametral_lens, markersize=4, color=:blue)
            scatter!(ax, points_outside_diametral_lens, markersize=4, color=(:red, 0.1))
            #band!(ax, lower_lens, upper_lens, color=(:green, 0.4))
            lines!(ax, lens, color=:black, linewidth=7)
            lines!(ax, [p, q], color=:black, linewidth=7)
            xlims!(ax, -7.5, -1.5)
            ylims!(ax, 1, 7)
        end
        resize_to_layout!(fig)
        fig
        @test_reference "diametral_figures.png" fig
    end

    @testset "Testing if encroached edges are detected" begin
        _x, _y = complicated_geometry()
        x = _x
        y = _y
        boundary_nodes, points = convert_boundary_points_to_indices(x, y)
        tri_1 = triangulate(points; boundary_nodes, delete_ghosts=false)
        boundary_nodes, points = convert_boundary_points_to_indices(x[1], y[1])
        tri_2 = triangulate(points; boundary_nodes, delete_ghosts=false)
        boundary_nodes, points = convert_boundary_points_to_indices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
        tri_3 = triangulate(points; boundary_nodes, delete_ghosts=false)
        boundary_nodes, points = convert_boundary_points_to_indices(reverse(reverse.(x[2])), reverse(reverse.(y[2])))
        tri_4 = triangulate(points; boundary_nodes, delete_ghosts=false)
        a, b = 0.0, 5.0
        c, d = 3.0, 7.0
        nx = 3
        ny = 3
        tri_5 = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts=false, single_boundary=false)
        tri_6 = triangulate_rectangle(a, b, c, d, nx, ny; delete_ghosts=false, single_boundary=true)
        @static if VERSION ≥ v"1.10"
            @inferred triangulate(rand(2, 250))
        end
        for (iii, tri) in enumerate((tri_1, tri_2, tri_3, tri_4, tri_5, tri_6))
            @info "Testing if encroached edges are detected. Run: $iii."
            args = DT.RefinementArguments(tri; use_lens=false)
            in_dt_encroached_edges, not_in_dt_encroached_edges = slow_encroachment_test(tri)
            all_bn = DT.get_all_boundary_nodes(tri)
            for (e, (b, k)) in not_in_dt_encroached_edges
                if DT.initial(e) ∈ all_bn && DT.terminal(e) ∈ all_bn && !DT.contains_segment(tri, e) # e.g. if an edge crosses an interior
                    continue
                end
                @test DT.is_encroached(tri, args, e) == b
            end
            for (e, (b, k)) in in_dt_encroached_edges
                if DT.initial(e) ∈ all_bn && DT.terminal(e) ∈ all_bn && !DT.contains_segment(tri, e)
                    continue
                end
                if tri === tri_1 && e == (11, 12) && k == 112 # k is right on e in this case, so we could adjust the code to check e.g. slightly obtuse, but let's just go past it for now 
                    continue
                end
                @test DT.is_encroached(tri, args, e) == b
            end

            in_dt_encroached_edges_lens_45, not_in_dt_encroached_edges_lens_45 = slow_encroachment_test_diametral_lens(tri, 45.0)
            in_dt_encroached_edges_lens_30, not_in_dt_encroached_edges_lens_30 = slow_encroachment_test_diametral_lens(tri, 30.0)
            in_dt_encroached_edges_lens_20, not_in_dt_encroached_edges_lens_20 = slow_encroachment_test_diametral_lens(tri, 20.0)
            all_bn = DT.get_all_boundary_nodes(tri)
            for (lens_angle, not_in_dt_encroached_edges) in zip((45.0, 30.0, 20.0), (not_in_dt_encroached_edges_lens_45, not_in_dt_encroached_edges_lens_30, not_in_dt_encroached_edges_lens_20))
                @info "Testing encroached edge detection. lens angle: $lens_angle"
                args = DT.RefinementArguments(tri; use_lens=true, min_angle=lens_angle)
                for (e, (b, k)) in not_in_dt_encroached_edges
                    if DT.initial(e) ∈ all_bn && DT.terminal(e) ∈ all_bn && !DT.contains_segment(tri, e) # e.g. if an edge crosses an interior
                        continue
                    end
                    if !b && !DT.unoriented_edge_exists(tri, e) # edge doesn't exist, for collinearities with tri_4 this causes annoying issues since the triangulation isn't unique in that case but the computation of the diametral lens isn't robust enough to distinguish if a point is on the edge or otherwise
                        @test DT.is_encroached(tri, args, e)
                    else
                        @test DT.is_encroached(tri, args, e) == b
                    end
                end
            end
            for (lens_angle, in_dt_encroached_edges) in zip((45.0, 30.0, 20.0), (in_dt_encroached_edges_lens_45, in_dt_encroached_edges_lens_30, in_dt_encroached_edges_lens_20))
                args = DT.RefinementArguments(tri; use_lens=true, min_angle=lens_angle)
                for (e, (b, k)) in in_dt_encroached_edges
                    if DT.initial(e) ∈ all_bn && DT.terminal(e) ∈ all_bn && !DT.contains_segment(tri, e)
                        continue
                    end
                    if tri === tri_1 && e == (11, 12) && (k == 112 || k == 126) # k is not visible from e, so it cannot possibly encroach
                        continue
                    end
                    if (tri ≠ tri_5 && tri ≠ tri_6) || !DT.unoriented_edge_exists(tri, e) # tri_(5/6) has some annoying collinearities, like tri_4 above
                        @test DT.is_encroached(tri, args, e) == b
                    end
                end
            end
        end
    end
end

@testset verbose = true "Subsegment splitting" begin
    @testset "balanced_power_of_two_ternary_split" begin
        x = [0.1992981, 2.9391, 0.0001, 0.891, 7.871, 10.59182, 63.0, 252.0, 0.1]
        nrst = DT.balanced_power_of_two_ternary_split.(x)
        @inferred DT.balanced_power_of_two_ternary_split(x[1])
        @test log2.(nrst) ≈ [-3, 0, -14, -1, 2, 2, 5, 7, -4]

        for _ in 1:1000
            ℓ = 2000rand()^2
            nrst = DT.balanced_power_of_two_ternary_split(ℓ)
            balanced_power_of_two_ternary_split = 1.0
            while ℓ > 3balanced_power_of_two_ternary_split
                balanced_power_of_two_ternary_split = 2balanced_power_of_two_ternary_split
            end
            while ℓ < 1.5balanced_power_of_two_ternary_split
                balanced_power_of_two_ternary_split = 0.5balanced_power_of_two_ternary_split
            end
            x = floor(log2(ℓ / 3))
            y = floor(log2(ℓ / 1.5))
            if abs(x - log2(ℓ)) < abs(y - log2(ℓ))
                @test nrst ≈ exp2(x)
            else
                @test nrst ≈ exp2(y)
            end
            @test nrst ≈ balanced_power_of_two_ternary_split
            _nrst = 2.0^(floor(log2(2ℓ / 3)))
            @test _nrst ≈ nrst
        end
    end

    @testset "balanced_power_of_two_quarternary_split" begin
        for _ in 1:100000
            ℓ = 2000rand()^2
            nrst = DT.balanced_power_of_two_quarternary_split(ℓ)
            @inferred DT.balanced_power_of_two_quarternary_split(ℓ)
            balanced_power_of_two_quarternary_split = 1.0
            while ℓ > 4balanced_power_of_two_quarternary_split
                balanced_power_of_two_quarternary_split = 2balanced_power_of_two_quarternary_split
            end
            while ℓ < 2balanced_power_of_two_quarternary_split
                balanced_power_of_two_quarternary_split = 0.5balanced_power_of_two_quarternary_split
            end
            x = floor(log2(ℓ / 4))
            y = floor(log2(ℓ / 2))
            if abs(x - log2(ℓ)) < abs(y - log2(ℓ))
                @test nrst ≈ exp2(x)
            else
                @test nrst ≈ exp2(y)
            end
            @test nrst ≈ balanced_power_of_two_quarternary_split
        end
    end

    @testset "Testing that ternary concentric shells are split in [1/3, 2/3]ℓ" begin
        for i in 1:500000
            p = 105randn(2)
            q = 237randn(2)
            t = DT.compute_concentric_shell_ternary_split_position(p, q)
            @inferred DT.compute_concentric_shell_ternary_split_position(p, q)
            @test _approx_ispow2(t * norm(p .- q))
            @test 1 / 3 ≤ t ≤ 2 / 3
            @test t ≈ _slow_compute_concentric_shell_ternary_split_position(p, q)
            @test t ≈ _slow_compute_concentric_shell_ternary_split_position(q, p)
        end
    end

    @testset "Testing that quarternary concentric shells are split in [1/4, 1/2]ℓ" begin
        for _ in 1:50000
            p = 105randn(2)
            q = 237randn(2)
            t = DT.compute_concentric_shell_quarternary_split_position(p, q)
            @test _approx_ispow2(t * norm(p .- q))
            @test 1 / 4 ≤ t ≤ 1 / 2
        end
    end

    @testset "segment_vertices_adjoin_other_segments_at_acute_angle/compute_split_position" begin
        for _ in 1:10
            points = [(-10.0, -10.0), (10.0, -10.0), (10.0, 10.0), (-10.0, 10.0), (0.05717272721, 0.00000001), (0.99988881, -0.000000998881)] # don't use (0, 1), else the split positions are all 1/2 since 1/2 is a power of 2 split
            tri = triangulate(points; boundary_nodes=[1, 2, 3, 4, 1], segments=Set([(5, 6)]))
            args = DT.RefinementArguments(tri)
            e = (5, 6) # if e = (6, 5) instead, we would need to be sure that t3 below is replaced by 1 - t3, since we always compute relative to the least vertex. See the comment in compute_split_position.
            adde! = args -> push!(args.segment_list, e)
            dele! = args -> DT.delete_unoriented_edge!(args.segment_list, e)
            segadj = DT.segment_vertices_adjoin_other_segments_at_acute_angle
            p, q = get_point(tri, e...)
            t1 = 1 / 2
            t2 = DT.compute_concentric_shell_ternary_split_position(p, q)
            t3 = DT.compute_concentric_shell_quarternary_split_position(p, q)
            m1_p = p .+ t1 .* (q .- p)
            m1_q = p .+ (1 - t1) .* (q .- p)
            m2_p = p .+ t2 .* (q .- p)
            m2_q = p .+ (1 - t2) .* (q .- p)
            m3_p = p .+ t3 .* (q .- p)
            m3_q = p .+ (1 - t3) .* (q .- p)
            m1 = (p .+ q) ./ 2
            zero_test = (tri, args, p, q, e, mvert) -> begin
                m = DT.compute_split_position(tri, args, e)
                @inferred DT.compute_split_position(tri, args, e)
                @inferred DT.compute_split_position(tri, args, e)
                @test DT.compute_split_position(tri, args, e) == DT.compute_split_position(tri, args, DT.reverse_edge(e))
                @test collect(m) ≈ collect(m1)
                dele!(args)
                m = DT.compute_split_position(tri, args, e)
                DT.compute_split_position(tri, args, e)
                @test DT.compute_split_position(tri, args, e) == DT.compute_split_position(tri, args, DT.reverse_edge(e))
                @test collect(m) ≈ collect(mvert)
                adde!(args)
            end
            one_test = (tri, args, p, q, e, mvert) -> begin
                m = DT.compute_split_position(tri, args, e)
                @inferred DT.compute_split_position(tri, args, e)
                @test DT.compute_split_position(tri, args, e) == DT.compute_split_position(tri, args, DT.reverse_edge(e))
                @test collect(m) ≈ collect(m1)
                dele!(args)
                m = DT.compute_split_position(tri, args, e)
                @test collect(DT.compute_split_position(tri, args, e)) ≈ collect(DT.compute_split_position(tri, args, DT.reverse_edge(e)))
                @test collect(m) ≈ collect(mvert)
                adde!(args)
            end
            two_test = (tri, args, p, q, e, mvert) -> begin
                m = DT.compute_split_position(tri, args, e)
                @inferred DT.compute_split_position(tri, args, e)
                @test collect(DT.compute_split_position(tri, args, e)) ≈ collect(DT.compute_split_position(tri, args, DT.reverse_edge(e)))
                @test collect(m) ≈ collect(m1)
                dele!(args)
                m = DT.compute_split_position(tri, args, e)
                @test collect(DT.compute_split_position(tri, args, e)) ≈ collect(DT.compute_split_position(tri, args, DT.reverse_edge(e)))
                @test collect(m) ≈ collect(mvert)
                adde!(args)
            end

            # case 1: No adjoining segments 
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (0, 0)
            @inferred segadj(tri, e)
            zero_test(tri, args, p, q, e, m1)

            # case 2: A single adjoining segment, but at an obtuse angle
            add_point!(tri, -0.2, 0.6)
            add_segment!(tri, 5, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (0, 0)
            zero_test(tri, args, p, q, e, m1_p)

            add_point!(tri, 1.4, 0.5)
            add_segment!(tri, 6, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (0, 0)
            zero_test(tri, args, p, q, e, m1_p)

            add_point!(tri, -0.6, 0.0)
            add_segment!(tri, 5, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (0, 0)
            zero_test(tri, args, p, q, e, m1_p)

            add_point!(tri, 1.4, 0.0)
            add_segment!(tri, 6, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (0, 0)
            zero_test(tri, args, p, q, e, m1_p)

            add_point!(tri, -0.1, -0.3)
            add_segment!(tri, 5, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (0, 0)
            zero_test(tri, args, p, q, e, m1_p)

            add_point!(tri, 1.2, -0.3)
            add_segment!(tri, 6, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (0, 0)
            zero_test(tri, args, p, q, e, m1_p)

            # case 3: A single adjoining segment at an acute angle taken from below the segment 
            orig_tri = deepcopy(tri)
            add_point!(tri, 0.1, -0.4)
            add_point!(orig_tri, 0.1, -0.4)
            add_segment!(tri, 5, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (1, 5)
            one_test(tri, args, p, q, e, m2_p)

            tri = deepcopy(orig_tri)
            add_point!(tri, 0.5, -0.1)
            add_point!(orig_tri, 0.5, -0.1)
            add_segment!(tri, 6, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (1, 6)
            one_test(tri, args, p, q, e, m2_q)

            # case 4: A single adjoining segment at an acute angle taken from above the segment 
            tri = deepcopy(orig_tri)
            add_point!(tri, 0.2, 0.6)
            add_point!(orig_tri, 0.2, 0.6)
            add_segment!(tri, 5, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (1, 5)
            one_test(tri, args, p, q, e, m2_p)

            tri = deepcopy(orig_tri)
            add_point!(tri, 0.6, 0.2)
            add_point!(orig_tri, 0.6, 0.2)
            add_segment!(tri, 6, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (1, 6)
            one_test(tri, args, p, q, e, m2_q)

            _points = [(0.0, 0.0), (9.0, 0.0), (9.0, 7.0)]
            _tri = triangulate(_points; segments=Set([(1, 2), (1, 3)]))
            @test DT.segment_vertices_adjoin_other_segments_at_acute_angle(_tri, (1, 2)) == (1, 1)
            @test DT.segment_vertices_adjoin_other_segments_at_acute_angle(_tri, (1, 3)) == (1, 1)

            # case 5: Two adjoining segments at an acute angle taken from above the segment from the same vertex
            tri = deepcopy(orig_tri)
            add_point!(tri, 0.4, 0.3)
            add_point!(tri, 0.3, 0.4)
            add_point!(orig_tri, 0.4, 0.3)
            add_point!(orig_tri, 0.3, 0.4)
            add_segment!(tri, 5, DT.num_points(tri) - 1)
            add_segment!(tri, 5, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (1, 5)
            one_test(tri, args, p, q, e, m2_p)

            # case 6: Two adjoining segments at an acute angle
            tri = deepcopy(orig_tri)
            add_point!(tri, 0.4, 0.1)
            add_point!(tri, 0.8, -0.1)
            add_point!(orig_tri, 0.4, 0.1)
            add_point!(orig_tri, 0.8, -0.1)
            add_segment!(tri, 5, DT.num_points(tri) - 1)
            add_segment!(tri, 6, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (2, 0)
            two_test(tri, args, p, q, e, m3_p)

            tri = deepcopy(orig_tri)
            add_point!(tri, 0.1, -0.5)
            add_point!(tri, 0.9, -0.6)
            add_point!(orig_tri, 0.1, -0.5)
            add_point!(orig_tri, 0.9, -0.6)
            add_segment!(tri, 5, DT.num_points(tri) - 1)
            add_segment!(tri, 6, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (2, 0)
            two_test(tri, args, p, q, e, m3_p)

            tri = deepcopy(orig_tri)
            add_point!(tri, 0.1, 0.7)
            add_point!(tri, 0.8, 0.7)
            add_point!(orig_tri, 0.1, 0.7)
            add_point!(orig_tri, 0.8, 0.7)
            add_segment!(tri, 5, DT.num_points(tri) - 1)
            add_segment!(tri, 6, DT.num_points(tri))
            @test segadj(tri, e) == segadj(tri, DT.reverse_edge(e)) == (2, 0)
            two_test(tri, args, p, q, e, m3_p)

            # case 7: Two adjoining segments at an an acute angle with the same shared vertex
            tri = deepcopy(orig_tri)
            add_point!(tri, 0.6, 0.4)
            add_point!(orig_tri, 0.6, 0.4)
            add_segment!(tri, 5, DT.num_points(tri))
            add_segment!(tri, 6, DT.num_points(tri))
            @test segadj(tri, e) == (1, 5) && segadj(tri, DT.reverse_edge(e)) == (1, 6) # should be adjoining unique segments
        end
    end

    @testset "split_subsegment" begin
        for _ in 1:10
            points = [(-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0), (-1.1, 0.0), (1.1, 0.0)]
            orig_points = copy(points)
            tri = triangulate(points; boundary_nodes=[1, 2, 3, 4, 1], segments=Set([(5, 6)]))
            orig_tri = deepcopy(tri)
            args = DT.RefinementArguments(tri)
            e = (5, 6)

            # just a normal split 
            @test DT.is_free(args, 7)
            DT.split_subsegment!(tri, args, e)
            @test tri == triangulate([orig_points; (0.0, 0.0)]; boundary_nodes=[1, 2, 3, 4, 1], segments=Set([(5, 7), (7, 6)]))
            @test DT.is_midpoint_split(args, 7)
            @test !DT.is_midpoint_split(args, 6) && !DT.is_midpoint_split(args, 5)
            @test get_point(tri, 7) == (0.0, 0.0)
            validate_insertion_event_history(tri, orig_tri, args.events)
            @test !DT.is_free(args, 7)

            # abutting segment next to segment on other straight line doesn't constitute a small angle
            @test DT.is_free(args, 8)
            @test DT.is_free(args, 9)
            DT.split_subsegment!(tri, args, (5, 7))
            @test !DT.is_free(args, 8)
            orig_tri = deepcopy(tri)
            orig_points = copy(orig_tri.points)
            DT.split_subsegment!(tri, args, (8, 7))
            @test !DT.is_free(args, 9)
            @test tri == triangulate([orig_points; (-0.275, 0.0)]; boundary_nodes=[1, 2, 3, 4, 1], segments=Set([(8, 7), (7, 6), (5, 8)]))
            @test args.midpoint_split_list ⊢ Set([7, 8, 9])
            @test get_point(tri, 9) == (-0.275, 0.0)
            validate_insertion_event_history(tri, orig_tri, args.events)

            # now adding an adjoining segment not at a small angle
            @test DT.is_free(args, 10)
            @test DT.is_free(args, 11)
            @test DT.is_free(args, 12)
            add_point!(tri, -1.5, 1.5)
            add_segment!(tri, 5, 10)
            DT.split_subsegment!(tri, args, (5, 8))
            @test DT.is_free(args, 10)
            @test !DT.is_free(args, 11)
            orig_tri = deepcopy(tri)
            orig_points = copy(orig_tri.points)
            DT.split_subsegment!(tri, args, (8, 9))
            @test !DT.is_free(args, 12)
            @test collect(get_point(tri, 11)) ≈ [-0.825, 0.0]
            @test collect(get_point(tri, 12)) ≈ [-0.4125, 0.0]
            @test tri == triangulate([orig_points; get_point(tri, 12)]; boundary_nodes=[1, 2, 3, 4, 1], segments=Set([get_interior_segments(tri)..., (8, 9)]))
            @test args.midpoint_split_list ⊢ Set([7, 8, 9, 11, 12])
            validate_insertion_event_history(tri, orig_tri, args.events)

            # a small-input angle
            points = [(0.0, 0.0), (9.0, 0.0), (9.0, 7.0)]
            tri = triangulate(points; segments=Set([(1, 2), (1, 3)]))
            args = DT.RefinementArguments(tri)
            orig_tri = deepcopy(tri)
            orig_points = copy(orig_tri.points)
            @test DT.is_free(args, 4)
            DT.split_subsegment!(tri, args, (1, 2))
            @test !DT.is_free(args, 4)
            @test validate_triangulation(tri)
            @test collect(get_point(tri, 4)) ≈ [4.5, 0.0] # first split is midpoint
            @test DT.is_midpoint_split(args, 4)
            @test DT.is_free(args, 5)
            DT.split_subsegment!(tri, args, (1, 4))
            @test DT.is_offcenter_split(args, 5)
            @test !DT.is_free(args, 5)
            @test validate_triangulation(tri)
            m = DT.compute_split_position(tri, args, (1, 4))
            @test collect(get_point(tri, 5)) ≈ collect(m)
            DT.split_subsegment!(tri, args, (1, 3))
            @test DT.is_offcenter_split(args, 6)
            @test !DT.is_free(args, 6)
            @test validate_triangulation(tri)
            @test collect(get_point(tri, 6)) ≈ [4.5, 3.5] # first split is midpoint 
            DT.split_subsegment!(tri, args, (1, 6))
            @test validate_triangulation(tri)
            m = DT.compute_split_position(tri, args, (1, 6))
            @test collect(get_point(tri, 7)) ≈ collect(m)
            @test !DT.encroaches_upon(get_point(tri, 1), get_point(tri, 5), get_point(tri, 7), args) # edges aligned on the same concentric circular shell should not encroach upon one another 
            @test !DT.encroaches_upon(get_point(tri, 1), get_point(tri, 7), get_point(tri, 5), args) # edges aligned on the same concentric circular shell should not encroach upon one another
            @test norm(get_point(tri, 7)) ≈ norm(get_point(tri, 5))
            DT.split_subsegment!(tri, args, (1, 5))
            DT.split_subsegment!(tri, args, (1, 7))
            @test norm(get_point(tri, 8)) ≈ norm(get_point(tri, 9))
            m1 = DT.compute_split_position(tri, args, (1, 5))
            m2 = DT.compute_split_position(tri, args, (1, 7))
            @test collect(get_point(tri, 8)) ≈ collect(m1)
            @test collect(get_point(tri, 9)) ≈ collect(m2)
            @test norm(get_point(tri, 8)) ≈ norm(get_point(tri, 9))
            @test !DT.encroaches_upon(get_point(tri, 1), get_point(tri, 8), get_point(tri, 9), args) # edges aligned on the same concentric circular shell should not encroach upon one another 
            @test !DT.encroaches_upon(get_point(tri, 1), get_point(tri, 9), get_point(tri, 8), args) # edges aligned on the same concentric circular shell should not encroach upon one another
            @test validate_triangulation(tri)
            DT.unlock_convex_hull!(tri; reconstruct=true)
            @test validate_triangulation(tri)

            # a small-input angle: Same as above, but flipped along vertical axis 
            points = [(0.0, 0.0), (7.0, 0.0), (0.0, 7.0)]
            tri = triangulate(points; segments=Set([(1, 2), (2, 3)]))
            args = DT.RefinementArguments(tri)
            unlock_convex_hull!(tri; reconstruct=true) # so that (1, 3) is not a segment
            DT.split_subsegment!(tri, args, (2, 1))
            @test collect(get_point(tri, 4)) ≈ [3.5, 0.0] # midpoint 
            DT.split_subsegment!(tri, args, (1, 4))
            @test collect(get_point(tri, 5)) ≈ [1.75, 0.0] # should still be midpoint, as doesn't adjoin any segment 
            DT.split_subsegment!(tri, args, (2, 4)) # should be a power of 2 split
            p2, p6, p4 = collect.(get_point(tri, 2, 6, 4))
            @test norm(p2 .- p6) ≈ 2
            @test norm(p4 .- p6) ≈ 1.5
            @test p6 ≈ [5.0, 0.0]
            DT.split_subsegment!(tri, args, (4, 6)) # should be a midpoint 
            @test collect(get_point(tri, 7)) ≈ collect(p4 .+ p6) ./ 2
            DT.split_subsegment!(tri, args, (6, 2)) # should be a midpoint split since we have already had a power of 2 split 
            @test collect(get_point(tri, 8)) ≈ collect(p2 .+ p6) ./ 2
            DT.split_subsegment!(tri, args, (6, 8)) # should be a midpoint split, since there is no adjoining acute segment 
            p8 = get_point(tri, 8)
            @test collect(get_point(tri, 9)) ≈ collect(p6 .+ p8) ./ 2
            DT.split_subsegment!(tri, args, (8, 2)) # should be another midpoint split, since we have already had a power of 2 split
            @test collect(get_point(tri, 10)) ≈ collect(p2 .+ p8) ./ 2
            DT.split_subsegment!(tri, args, (2, 3)) # this is the first time this segment has been split, so it should be a midpoint split 
            @test collect(get_point(tri, 11)) ≈ collect(get_point(tri, 2) .+ get_point(tri, 3)) ./ 2
            DT.split_subsegment!(tri, args, (11, 3)) # this isn't adjoining any segment, so it should be a midpoint split 
            @test collect(get_point(tri, 12)) ≈ collect(get_point(tri, 11) .+ get_point(tri, 3)) ./ 2
            DT.split_subsegment!(tri, args, (11, 2)) # this should be a power of 2 split 
            p2, p11, p13 = collect.(get_point(tri, 2, 11, 13))
            @test norm(p2 .- p13) ≈ 2.0
            @test norm(p11 .- p13) ≈ norm(p2 .- p11) .- 2.0 ≈ 2.949747468305833
            DT.split_subsegment!(tri, args, (13, 2)) # this should be a midpoint split
            @test collect(get_point(tri, 14)) ≈ collect(p2 .+ p13) ./ 2
            DT.split_subsegment!(tri, args, (13, 11)) # this should be a midpoint split since there is no adjoining acute segment 
            @test collect(get_point(tri, 15)) ≈ collect(p11 .+ p13) ./ 2
            @test validate_triangulation(tri)
            @inferred get_point(tri, 5)

            # small-input angle with multiple adjoining segments at acute angles. there will still only ever be one adjoining acute segment for a given vertex, though.
            points = [(0.0, 0.0), (5.0, 7.0), (11.0, 4.0), (3.0, -1.0), (-12.0, -3.0)]
            segments = Set([(1, 2), (1, 3), (1, 4), (1, 5)]) # (1, 5) doesn't adjoin any at an acute angle 
            tri = triangulate(points; segments)
            args = DT.RefinementArguments(tri)
            unlock_convex_hull!(tri; reconstruct=true)
            DT.split_subsegment!(tri, args, (1, 2)) # -> 6
            DT.split_subsegment!(tri, args, (1, 3)) # -> 7
            DT.split_subsegment!(tri, args, (1, 4)) # -> 8
            DT.split_subsegment!(tri, args, (1, 5)) # -> 9. these should all be midpoint splits since they are the initial splits 
            @test collect(get_point(tri, 6)) ≈ collect(get_point(tri, 1) .+ get_point(tri, 2)) / 2
            @test collect(get_point(tri, 7)) ≈ collect(get_point(tri, 1) .+ get_point(tri, 3)) / 2
            @test collect(get_point(tri, 8)) ≈ collect(get_point(tri, 1) .+ get_point(tri, 4)) / 2
            @test collect(get_point(tri, 9)) ≈ collect(get_point(tri, 1) .+ get_point(tri, 5)) / 2
            DT.split_subsegment!(tri, args, (6, 2)) # -> 10. should just be a midpoint split 
            @test collect(get_point(tri, 10)) ≈ collect(get_point(tri, 6) .+ get_point(tri, 2)) / 2
            DT.split_subsegment!(tri, args, (6, 10)) # -> 11. should be another midpoint split 
            @test collect(get_point(tri, 11)) ≈ collect(get_point(tri, 6) .+ get_point(tri, 10)) / 2
            DT.split_subsegment!(tri, args, (1, 6)) # -> 12.
            @test norm(get_point(tri, 12) .- get_point(tri, 6)) ≈ 2.3011626335213133
            @test norm(get_point(tri, 12) .- get_point(tri, 1)) ≈ 2.0
            @test norm(get_point(tri, 1) .- get_point(tri, 6)) ≈ 4.3011626335213133
            DT.split_subsegment!(tri, args, (12, 1)) # -> 13. should be a midpoint split 
            @test collect(get_point(tri, 13)) ≈ collect(get_point(tri, 1) .+ get_point(tri, 12)) / 2
            DT.split_subsegment!(tri, args, (7, 3)) # -> 14. should be a midpoint split since there are no adjoining acute segments 
            @test collect(get_point(tri, 14)) ≈ collect(get_point(tri, 7) .+ get_point(tri, 3)) / 2
            DT.split_subsegment!(tri, args, (14, 7)) # -> 15. another midpoint split 
            @test collect(get_point(tri, 15)) ≈ collect(get_point(tri, 14) .+ get_point(tri, 7)) / 2
            DT.split_subsegment!(tri, args, (7, 1)) # -> 16. there are two adjoining segments, but at the same vertex so this should be a power of 2 split
            @test norm(DT.get_point(tri, 16) .- get_point(tri, 1)) ≈ 2.0
            @test norm(DT.get_point(tri, 16) .- get_point(tri, 7)) ≈ 3.8523499553598124
            @inferred DT.split_subsegment!(tri, args, (1, 16))

            # small-input angle with adjoining segments at each vertex 
            for iii in 1:2
                points = [(0.0, 0.0), (2.2, 0.0), (-0.5, 0.5), (2.5, -0.5)]
                segments = Set([(1, 2), (2, 3), (1, 4)])
                tri = triangulate(points; segments)
                args = DT.RefinementArguments(tri)
                DT.split_subsegment!(tri, args, (1, 2)) # -> 5. midpoint split 
                @test collect(get_point(tri, 5)) ≈ [1.1, 0.0]
                add_segment!(tri, 5, 4)
                if iii == 1
                    e = (5, 2)
                else
                    e = (2, 5)
                end
                DT.split_subsegment!(tri, args, e) # -> 6. Should be a power of 2 split in [1/4, 1/2], relative to vertex 2 since 2 < 5.
                p2, p5, p6 = collect.(get_point(tri, 2, 5, 6))
                @test norm(p2 .- p6) ≤ norm(p5 .- p6)
                @test norm(p2 .- p6) ≈ 1 / 2
                @test norm(p5 .- p6) ≈ 0.6 ≈ norm(p5 .- p2) .- 1 / 2
                DT.split_subsegment!(tri, args, (5, 6)) # -> 7. should be another off-center split 
                p7 = get_point(tri, 7)
                @test norm(p7 .- p5) ≈ 0.25
                @test norm(p7 .- p6) ≈ 0.35 ≈ norm(p5 .- p6) .- 0.25
                orig_ℓ = norm(p2 .- p5) # we want to see if the three new subsegments are ≥ 1/5 original length.
                @test norm(p2 .- p6) ≥ orig_ℓ / 5
                @test norm(p5 .- p6) ≥ orig_ℓ / 5
                @test norm(p7 .- p6) ≥ orig_ℓ / 5
                DT.split_subsegment!(tri, args, (5, 7)) # -> 8. should be a midpoint split
                @test collect(get_point(tri, 8)) ≈ collect(p5 .+ p7) ./ 2
                DT.split_subsegment!(tri, args, (7, 6)) # -> 9. should be a midpoint split
                @test collect(get_point(tri, 9)) ≈ collect(p6 .+ p7) ./ 2
                DT.split_subsegment!(tri, args, (6, 2)) # -> 10. should be a midpoint split
                @test collect(get_point(tri, 10)) ≈ collect(p2 .+ p6) ./ 2
            end
        end
    end
end

@testset "enqueueing and splitting all encroached segments" begin
    for PT in (DT.Exact, DT.Adaptive)
        for iii in 1:5
            for use_lens in (false, true)
                for pass in 1:2
                    @info "Testing enqueuing and splitting all encroached segments. Predicates: $PT; Run: $iii; lens: $use_lens; Block: $pass"
                    if pass == 1
                        p1 = (0.0, 0.0)
                        p2 = (1.0, 0.0)
                        p3 = (0.0, 1.0)
                        p4 = (1.0, 1.0)
                        p5 = (0.5, 0.5)
                        pts = [p1, p2, p3, p4, p5]
                        C = Set{NTuple{2,Int}}()
                        for i in 1:15
                            θ = 2π * rand()
                            r = 0.5sqrt(rand())
                            x = 0.5 + r * cos(θ)
                            y = 0.5 + r * sin(θ)
                            push!(pts, (x, y))
                            push!(C, (5, 5 + i))
                        end
                        tri = triangulate(pts; delete_ghosts=false, boundary_nodes=[1, 2, 4, 3, 1], segments=C, predicates=PT())
                        args = DT.RefinementArguments(tri, use_lens=use_lens, min_angle=25.0, seditious_angle=15.0, max_area=0.1, predicates=PT())
                        DT.enqueue_all_encroached_segments!(args, tri)
                        @inferred DT.enqueue_all_encroached_segments!(args, tri)
                        manual_enqueue = PriorityQueue{NTuple{2,Int},Float64}(Base.Order.Reverse)
                        for e in each_edge(tri)
                            if DT.contains_segment(tri, e...)
                                flag = DT.is_encroached(tri, args, e)
                                if flag
                                    manual_enqueue[e] = DT.edge_length_sqr(tri, e...)
                                end
                            end
                        end
                        compare_encroach_queues(args, manual_enqueue)
                        empty!(args.events)
                        DT.split_all_encroached_segments!(tri, args)
                        @test isempty(args.queue.segments)
                        @test !isempty(args.events.added_triangles)
                        DT.enqueue_all_encroached_segments!(args, tri)
                        @test isempty(args.queue.segments)
                        if use_lens
                            for e in each_segment(tri)
                                @test !DT.is_encroached(tri, args, e)
                            end
                        else
                            @test is_conformal(tri)
                        end
                        @test validate_triangulation(tri)
                    else
                        points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
                        for r in LinRange(0.1, 0.9, 15)
                            push!(points, (r, 0.01))
                            push!(points, (r, 0.99))
                            push!(points, (0.05, r))
                            push!(points, (0.95, r))
                        end
                        for _ in 1:10
                            push!(points, (rand(), rand()))
                        end
                        tri = triangulate(points; boundary_nodes=[1, 2, 3, 4, 1], predicates=PT())
                        args = DT.RefinementArguments(tri, use_lens=use_lens, min_angle=29.0, predicates=PT(), seditious_angle=19.0, max_area=0.05)
                        DT.enqueue_all_encroached_segments!(args, tri)
                        manual_enqueue = PriorityQueue{NTuple{2,Int},Float64}(Base.Order.Reverse)
                        for e in each_edge(tri)
                            if DT.contains_segment(tri, e...)
                                flag = DT.is_encroached(tri, args, e)
                                if flag
                                    manual_enqueue[e] = DT.edge_length_sqr(tri, e...)
                                end
                            end
                        end
                        compare_encroach_queues(args, manual_enqueue)
                        empty!(args.events)
                        DT.split_all_encroached_segments!(tri, args)
                        @test isempty(args.queue.segments)
                        @test !isempty(args.events.added_triangles)
                        DT.enqueue_all_encroached_segments!(args, tri)
                        @test isempty(args.queue.segments)
                        if use_lens
                            for e in each_segment(tri)
                                @test !DT.is_encroached(tri, args, e)
                            end
                        else
                            @test is_conformal(tri)
                        end
                        @test validate_triangulation(tri)
                    end
                end
            end
        end
    end
end

@testset "triangle assessment" begin
    @testset "is_triangle_nestled" begin
        tri = triangulate_rectangle(0, 10, 0, 10, 6, 6, single_boundary=false)
        segments = get_all_segments(tri)
        @inferred get_all_segments(tri)
        for T in each_solid_triangle(tri)
            idx = DT.squared_triangle_lengths_and_smallest_index(get_point(tri, T...)...)[4]
            u, v, w = T
            flag = false
            if idx == 1
                i, j, k = u, v, w
            elseif idx == 2
                i, j, k = v, w, u
            else
                i, j, k = w, u, v
            end
            flag = DT.contains_segment(tri, k, i) && DT.contains_segment(tri, k, j)
            if flag
                @test DT.is_triangle_nestled(tri, (u, v, w), idx) && DT.is_triangle_nestled(tri, (v, w, u), mod1(idx - 1, 3)) && DT.is_triangle_nestled(tri, (w, u, v), mod1(idx - 2, 3))
            else
                @test !DT.is_triangle_nestled(tri, T, idx) && !DT.is_triangle_nestled(tri, (v, w, u), mod1(idx - 1, 3)) && !DT.is_triangle_nestled(tri, (w, u, v), mod1(idx - 2, 3))
            end
        end

        points = [(8.0, 0.0), (2.0, 6.0), (2.0, 2.0)]
        segments = Set([(1, 2)])
        tri = triangulate(points; segments)
        idx = DT.squared_triangle_lengths_and_smallest_index(get_point(tri, (1, 2, 3)...)...)[4]
        @inferred DT.squared_triangle_lengths_and_smallest_index(get_point(tri, (1, 2, 3)...)...)
        @test !DT.is_triangle_nestled(tri, (1, 2, 3), idx) && !DT.is_triangle_nestled(tri, (2, 3, 1), mod1(idx - 1, 3)) && !DT.is_triangle_nestled(tri, (3, 1, 2), mod1(idx - 2, 3))
        add_segment!(tri, 1, 3)
        @test DT.is_triangle_nestled(tri, (1, 2, 3), idx) && DT.is_triangle_nestled(tri, (2, 3, 1), mod1(idx - 1, 3)) && DT.is_triangle_nestled(tri, (3, 1, 2), mod1(idx - 2, 3))

        points = [(0.0, 0.0), (2.0, 0.0), (10.0, 0.0), (1.0, 1.0), (8.0, 2.0), (0.0, 6.0), (10.0, 6.0)]
        segments = Set([(1, 2), (1, 5), (1, 4)])
        tri = triangulate(points; segments)
        T = (1, 4, 6)
        @test !DT.is_triangle_nestled(tri, T, 1) && !DT.is_triangle_nestled(tri, (4, 6, 1), 2) && !DT.is_triangle_nestled(tri, (6, 1, 4), 3)
    end

    @testset "is_triangle_seditious" begin
        for seditious_angle in (20.0, 0.0)
            points = [(0.0, 0.0), (6.0, 0.0), (12.0, 0.0), (11.0, 3.0)]
            midpt = (5.788590116313, 1.5787063953581)
            segments = Set([(1, 2), (2, 3), (1, 4)])
            tri = triangulate(points; segments)
            args = DT.RefinementArguments(tri; seditious_angle)
            @test !DT.is_triangle_seditious(tri, args, 1, 2, 4, 1) && !DT.is_triangle_seditious(tri, args, 2, 4, 1, 2) && !DT.is_triangle_seditious(tri, args, 4, 1, 2, 3)
            @test !DT.is_triangle_seditious(tri, args, 2, 3, 4, 2) && !DT.is_triangle_seditious(tri, args, 3, 4, 2, 3) && !DT.is_triangle_seditious(tri, args, 4, 2, 3, 1)
            push!(points, midpt)
            DT.complete_split_edge_and_legalise!(tri, 1, 4, 5)
            @test DT.is_triangle_seditious(tri, args, 1, 2, 5, 2) && DT.is_triangle_seditious(tri, args, 2, 5, 1, 1) && DT.is_triangle_seditious(tri, args, 5, 1, 2, 3)
            delete!(args.segment_vertices, 1)
            @test !DT.is_triangle_seditious(tri, args, 1, 2, 5, 2) && !DT.is_triangle_seditious(tri, args, 2, 5, 1, 1) && !DT.is_triangle_seditious(tri, args, 5, 1, 2, 3)
            push!(args.segment_vertices, 1)
            points[4] = (10.0, 4.0)
            points[5] = (5.2623546511936, 2.1049418604774)
            @test !DT.is_triangle_seditious(tri, args, 1, 2, 5, 2) && !DT.is_triangle_seditious(tri, args, 2, 5, 1, 1) && !DT.is_triangle_seditious(tri, args, 5, 1, 2, 3)
            points[5] = (5.5708585859039, 2.2283434343616)
            if seditious_angle == 20.0
                @test !DT.is_triangle_seditious(tri, args, 1, 2, 5, 2) && !DT.is_triangle_seditious(tri, args, 2, 5, 1, 1) && !DT.is_triangle_seditious(tri, args, 5, 1, 2, 3)
            else
                @test DT.is_triangle_seditious(tri, args, 1, 2, 5, 2) && DT.is_triangle_seditious(tri, args, 2, 5, 1, 1) && DT.is_triangle_seditious(tri, args, 5, 1, 2, 3)
            end
        end
    end

    @testset "assess_triangle_quality" begin
        points = [(0.0, 0.0), (2.0, 0.0), (10.0, 0.0), (1.0, 1.0), (8.0, 2.0), (0.0, 6.0), (10.0, 6.0)]
        segments = Set([(1, 2), (1, 5), (1, 4)])
        tri = triangulate(points; segments)
        args = DT.RefinementArguments(tri; min_angle=15.0, min_area=3.01, max_area=10.0)
        good_T, bad_T = slow_triangle_assess(tri, args)
        @test DT.contains_triangle((1, 4, 6), good_T)[2] # small angle, but area = 3 < min_area
        @test DT.contains_triangle((6, 4, 5), bad_T)[2] # large area 
        @test DT.contains_triangle((6, 5, 7), bad_T)[2] # large area 
        @test DT.contains_triangle((5, 2, 3), good_T)[2]
        @test DT.contains_triangle((7, 5, 3), good_T)[2]
        @test DT.contains_triangle((1, 5, 4), good_T)[2]
        @test DT.contains_triangle((1, 2, 5), good_T)[2]

        points = [(0.0, 0.0), (6.0, 0.0), (12.0, 0.0), (11.0, 3.0), (5.788590116313, 1.5787063953581)]
        segments = Set([(1, 2), (2, 3), (1, 5), (5, 4)])
        tri = triangulate(points; segments)
        args = DT.RefinementArguments(tri, min_angle=31.0)
        @test !DT.assess_triangle_quality(tri, args, (1, 2, 5))[2] # seditious 
        @test !DT.assess_triangle_quality(tri, args, (5, 2, 4))[2] # seditious 
        @test DT.assess_triangle_quality(tri, args, (4, 2, 3))[2] # small angle
        slow_triangle_assess(tri, args)
    end

    @testset "enqueue_all_bad_triangles!" begin
        for _ in 1:100
            points = [Tuple(rand(2)) for _ in 1:250]
            push!(points, (-1, -1), (2, -1), (2, 2), (-1, 2))
            tri = triangulate(points)
            args = DT.RefinementArguments(tri; min_area=1e-8, max_area=1e-2, min_angle=28.5)
            @test isempty(args.queue.triangles)
            DT.enqueue_all_bad_triangles!(args, tri)
            @test !isempty(args.queue.triangles)
            manual_queue = slow_triangle_assess_queue(tri, args)
            compare_triangle_queues(args, manual_queue)
        end

        points = [(0.0, 0.0), (2.0, 0.0), (10.0, 0.0), (1.0, 1.0), (8.0, 2.0), (0.0, 6.0), (10.0, 6.0)]
        segments = Set([(1, 2), (1, 5), (1, 4)])
        tri = triangulate(points; segments)
        args = DT.RefinementArguments(tri; min_angle=15.0, min_area=3.01, max_area=10.0)
        DT.enqueue_all_bad_triangles!(args, tri)
        queue = PriorityQueue(Base.Order.Reverse,
            (6, 4, 5) => DT.triangle_radius_edge_ratio(get_point(tri, (6, 4, 5)...)...),
            (6, 5, 7) => DT.triangle_radius_edge_ratio(get_point(tri, (6, 5, 7)...)...))
        compare_triangle_queues(args, queue)
    end
end

@testset "steiner_points with circumcenters" begin
    @testset "get_steiner_point" begin
        for _ in 1:10
            tri = triangulate(rand(2, 5000))
            args = DT.RefinementArguments(tri, use_circumcenter=true)
            good_results = Dict{NTuple{3,Int},NTuple{2,Float64}}()
            bad_results = Dict{NTuple{3,Int},NTuple{2,Float64}}()
            for T in each_solid_triangle(tri)
                flag, c = DT.get_steiner_point(tri, args, T)
                DT.is_none(flag) && @test c == DT.triangle_circumcenter(get_point(tri, T...)...)
                @inferred DT.triangle_circumcenter(get_point(tri, T...)...)
                if DT.is_precision_failure(flag)
                    bad_results[T] = c
                else
                    good_results[T] = c
                end
            end
            @test length(bad_results) ≤ 3 # should be quite rare
        end
    end

    @testset "locate_steiner_point" begin
        # random tests
        for _ in 1:10
            tri = triangulate(randn(2, 5000))
            args = DT.RefinementArguments(tri, use_circumcenter=true)
            for T in each_solid_triangle(tri)
                flag, c = DT.get_steiner_point(tri, args, T)
                V, loc_flag = DT.locate_steiner_point(tri, args, T, c)
                _flag = if DT.is_outside(loc_flag)
                    @test DT.is_positively_oriented(DT.triangle_orientation(tri, V))
                    @test DT.is_outside(DT.point_position_relative_to_triangle(tri, V, c))
                else
                    @test DT.is_inside(loc_flag) || DT.is_on(loc_flag)
                    @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, c))
                end
            end
        end

        # circumcenter on hypotenuse
        tri = triangulate([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
        args = DT.RefinementArguments(tri, use_circumcenter=true)
        c = (0.5, 0.5)
        for _ in 1:100
            for T in ((1, 2, 3), (2, 3, 1), (3, 1, 2))
                V, loc_flag = DT.locate_steiner_point(tri, args, T, c)
                @test DT.is_on(loc_flag) && !DT.is_ghost_triangle(V)
                @test !DT.is_outside(DT.point_position_relative_to_triangle(tri, V, c))
            end
        end
    end

    @testset "check_for_invisible_steiner_point" begin
        @testset "many triangles" begin
            for _ in 1:10
                a, b, c, d, e, f, g, h, i, j, k, ℓ = (0.0, 0.0), (8.0, 0.0), (4.0, 0.0),
                (8.0, 4.0), (8.0, 8.0), (4.0, 8.0), (0.0, 8.0),
                (0.0, 4.0), (4.0, 6.0), (2.0, 4.0), (4.0, 2.0), (6.0, 4.0)
                boundary_nodes, points = convert_boundary_points_to_indices([[[a, c, b, d, e, f, g, h, a]], [[i, ℓ, k, j, i]]])
                m, n, o, p, q, r = (1.0, 7.0), (6.0, 7.0), (6.0, 3.0), (2.0, 3.0), (2.5, 2.5), (1.0, 1.0)
                push!(points, m, n, o, p, q, r)
                tri = triangulate(points; boundary_nodes=boundary_nodes)
                args = DT.RefinementArguments(tri, use_circumcenter=true)
                _to_T = ((p, q, r),) -> (findfirst(==(p), points), findfirst(==(q), points), findfirst(==(r), points))
                _ls = T -> DT.locate_steiner_point(tri, args, _to_T(T), DT.triangle_circumcenter(get_point(tri, _to_T(T)...)...))
                _check = (T, V, flag) -> DT.check_for_invisible_steiner_point(tri, V, _to_T(T), flag, DT.triangle_circumcenter(get_point(tri, _to_T(T)...)...))

                T = (g, h, m)
                V, flag = _ls(T)
                @test DT.is_outside(flag)
                c′, V′ = _check(T, V, flag)
                @test V′ == _to_T(T) && c′ == DT.triangle_centroid(T...)

                T = (g, m, f)
                V, flag = _ls(T)
                @test DT.is_outside(flag)
                c′, V′ = _check(T, V, flag)
                @test V′ == _to_T(T) && c′ == DT.triangle_centroid(T...)

                T = (f, m, i)
                V, flag = _ls(T)
                @test DT.is_inside(flag) && DT.compare_triangles(V, _to_T((f, m, i)))
                c′, V′ = _check(T, V, flag)
                @test c′ == DT.triangle_circumcenter(get_point(tri, _to_T(T)...)...) && V′ == _to_T(T)

                T = (f, i, n)
                V, flag = _ls(T)
                @test DT.is_inside(flag) && DT.compare_triangles(V, _to_T((f, i, n)))
                c′, V′ = _check(T, V, flag)
                @test c′ == DT.triangle_circumcenter(get_point(tri, _to_T(T)...)...) && V′ == _to_T(T)

                T = (j, h, p)
                V, flag = _ls(T)
                @test DT.is_on(flag) && (DT.compare_triangles(V, _to_T((h, p, j))) || DT.compare_triangles(V, _to_T((h, r, p))))
                c′, V′ = _check(T, V, flag)
                @test c′ == DT.triangle_circumcenter(get_point(tri, _to_T(T)...)...) && (V′ == _to_T(T) || V′ == _to_T((h, r, p)))

                T = (q, k, j)
                V, flag = _ls(T)
                @test DT.is_outside(flag)
                c′, V′ = _check(T, V, flag)
                @test V′ == _to_T(T) && c′ == DT.triangle_centroid(T...)

                T = (ℓ, k, o)
                V, flag = _ls(T)
                @test DT.is_outside(flag)
                c′, V′ = _check(T, V, flag)
                @test V′ == _to_T(T) && c′ == DT.triangle_centroid(T...)

                T = (k, q, c)
                V, flag = _ls(T)
                @test DT.is_inside(flag) && DT.compare_triangles(V, _to_T((q, r, c)))
                c′, V′ = _check(T, V, flag)
                @test c′ == DT.triangle_circumcenter(get_point(tri, _to_T(T)...)...) && V′ == V

                T = (d, o, b)
                V, flag = _ls(T)
                @test DT.is_inside(flag) && DT.compare_triangles(V, _to_T((d, o, b)))
                c′, V′ = _check(T, V, flag)
                @test c′ == DT.triangle_circumcenter(get_point(tri, _to_T(T)...)...) && V′ == V

                T = (n, i, ℓ)
                V, flag = _ls(T)
                @test DT.is_inside(flag) && DT.compare_triangles(V, _to_T((n, i, ℓ)))
                c′, V′ = _check(T, V, flag)
                @test c′ == DT.triangle_circumcenter(get_point(tri, _to_T(T)...)...) && V′ == V
            end
        end

        @testset "triangle on boundary edge" begin
            for _ in 1:10
                points = [(0.0, 0.0), (2.0, 0.0), (0.0, 2.0)]
                tri = triangulate(points; boundary_nodes=[1, 2, 3, 1])
                args = DT.RefinementArguments(tri, use_circumcenter=true)
                flag, c = DT.get_steiner_point(tri, args, (1, 2, 3))
                @inferred DT.get_steiner_point(tri, args, (1, 2, 3))
                @test DT.is_none(flag) && c == DT.triangle_circumcenter(get_point(tri, 1, 2, 3)...)
                c = (1.0, 1.0) # accuracy
                V, flag = DT.locate_steiner_point(tri, args, (1, 2, 3), c)
                @inferred DT.locate_steiner_point(tri, args, (1, 2, 3), c)
                @test DT.compare_triangles(V, (1, 2, 3)) && DT.is_on(flag)
                c′, V′ = DT.check_for_invisible_steiner_point(tri, V, (1, 2, 3), flag, c)
                @test c′ == (1.0, 1.0) && V′ == V
                @test DT.get_init_for_steiner_point(tri, (1, 2, 3)) == DT.get_init_for_steiner_point(tri, (2, 3, 1)) == DT.get_init_for_steiner_point(tri, (3, 1, 2)) == 1
                @inferred DT.get_init_for_steiner_point(tri, (1, 2, 3))
            end
        end

        @testset "triangle on segment edge" begin
            for _ in 1:10
                points = [(0.0, 0.0), (1.0, 0.0), (0.7, 0.7), (0.0, 1.0)]
                tri = triangulate(points; boundary_nodes=[1, 2, 3, 4, 1], segments=Set([(2, 4)]))
                args = DT.RefinementArguments(tri, use_circumcenter=true)
                T = (1, 2, 4)
                flag, c = DT.get_steiner_point(tri, args, T)
                @test DT.is_none(flag) && c == DT.triangle_circumcenter(get_point(tri, T...)...)
                c = (1 / 2, 1 / 2) # accuracy 
                V, flag = DT.locate_steiner_point(tri, args, T, c)
                @test DT.compare_triangles(V, T) && DT.is_on(flag)
                c′, V′ = DT.check_for_invisible_steiner_point(tri, V, T, flag, c)
                @test c′ == (1 / 2, 1 / 2) && V′ == V
                @test DT.get_init_for_steiner_point(tri, T) == DT.get_init_for_steiner_point(tri, (2, 4, 1)) == DT.get_init_for_steiner_point(tri, (4, 1, 2)) == 1

                T = (2, 3, 4) # this triangle's circumcenter is blocked by a segment
                flag, c = DT.get_steiner_point(tri, args, T)
                @test DT.is_none(flag) && c == DT.triangle_circumcenter(get_point(tri, T...)...)
                V, flag = DT.locate_steiner_point(tri, args, T, c)
                @test DT.compare_triangles(V, T) && DT.is_outside(flag)
                c′, V′ = DT.check_for_invisible_steiner_point(tri, V, T, flag, c)
                @test c′ == DT.triangle_centroid(get_point(tri, T...)...) && V′ == T
                @test DT.get_init_for_steiner_point(tri, T) == DT.get_init_for_steiner_point(tri, (2, 3, 4)) == DT.get_init_for_steiner_point(tri, (3, 4, 2)) == 3
            end
        end
    end

    @testset "split_triangle! / assess_added_triangles!" begin
        for _ in 1:10
            points = [(0.0, 0.0), (0.0, 8.0), (10.0, 8.0), (8.0, 0.0), (4.0, 6.0), (4.5, 5.0), (5.0, 8.0)]
            tri = triangulate(points)
            args = DT.RefinementArguments(tri, use_circumcenter=true)
            orig_tri = deepcopy(tri)
            orig_points = deepcopy(points)

            flag = DT.split_triangle!(tri, args, (2, 1, 5))
            @test tri == orig_tri
            @test DT.is_encroachment_failure(flag)
            segments = args.queue.segments
            @test length(segments) == 1 && segments[(2, 1)] == 8.0^2

            flag = DT.split_triangle!(tri, args, (2, 5, 7))
            @test DT.is_successful_insertion(flag)
            @test !DT.is_free(args, DT.num_points(tri)) # because its circumcenter is on the boundary, it is no longer a free vertex
            validate_insertion_event_history(tri, orig_tri, args.events)
            @test validate_triangulation(tri)
            add_point!(orig_tri, DT.triangle_circumcenter(get_point(orig_tri, (2, 5, 7)...)...))
            @test tri == orig_tri
            orig_tri = deepcopy(tri)
            DT.assess_added_triangles!(args, tri)
            @test length(args.queue.triangles) == 1 && args.queue[(8, 2, 5)] == DT.triangle_radius_edge_ratio(get_point(tri, (8, 2, 5)...)...)

            flag = DT.split_triangle!(tri, args, (6, 1, 4))
            @test tri == orig_tri
            @test DT.is_encroachment_failure(flag)
            segments = args.queue.segments
            @test length(segments) == 2 && segments[(1, 4)] == 8.0^2

            points = [(0.0, 0.0), (0.0, 8.0), (10.0, 8.0), (8.0, 0.0), (4.0, 6.0), (4.0, 1.0), (5.0, 8.0)]
            tri = triangulate(points)
            args = DT.RefinementArguments(tri, use_circumcenter=true)
            orig_tri = deepcopy(tri)
            orig_points = deepcopy(points)
            flag = DT.split_triangle!(tri, args, (6, 1, 4))
            @test DT.is_encroachment_failure(flag) # circumcenter is outside, but the centroid is used
            @test length(args.queue.segments) == 1 && args.queue.segments[(1, 4)] == 8.0^2
            @test tri == orig_tri

            points = [(0.0, 0.0), (0.0, 8.0), (10.0, 8.0), (8.0, 0.0), (4.0, 6.0), (4.0, 4.0), (5.0, 8.0)] # circumcenter is on the boundary edge
            tri = triangulate(points)
            orig_tri = deepcopy(tri)
            orig_points = deepcopy(points)
            args = DT.RefinementArguments(tri, use_circumcenter=true)
            unlock_convex_hull!(tri) # don't want to deal with EncroachmentFailure proper, only if it's outside of the domain
            flag = DT.split_triangle!(tri, args, (6, 1, 4))
            @test DT.is_successful_insertion(flag)
            @test DT.is_free(args, DT.num_points(tri))
            @test length(args.queue.segments) == 0
            @test validate_triangulation(tri)
        end
    end

    @testset "check_for_steiner_point_on_segment" begin
        for PT in subtypes(DT.AbstractPredicateKernel)
            points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            tri = triangulate(points; randomise=false)
            push!(points, (1 / 2, 1 / 2))
            V = (1, 2, 3)
            V′ = (1, 2, 3)
            new_point = 5
            flag = DT.Cert.On
            ison_flag = DT.check_for_steiner_point_on_segment(tri, V, V′, new_point, flag, PT())
            @test !ison_flag
            add_segment!(tri, 1, 3)
            ison_flag = DT.check_for_steiner_point_on_segment(tri, V, V′, new_point, flag, PT())
            @test ison_flag
            ison_flag = DT.check_for_steiner_point_on_segment(tri, V, (2, 3, 1), new_point, flag, PT())
            @test ison_flag
            ison_flag = DT.check_for_steiner_point_on_segment(tri, V, (5, 6, 1), new_point, flag, PT())
            @test !ison_flag
            points[5] = (0.4, 0.555)
            ison_flag = DT.check_for_steiner_point_on_segment(tri, V, V′, new_point, flag, PT())
            @test !ison_flag
        end
    end
end

@testset "finalise!" begin
    tri = triangulate(rand(2, 50))
    args = DT.RefinementArguments(tri)
    @test DT.has_boundary_nodes(tri)
    @test DT.has_ghost_triangles(tri)
    DT.finalise!(tri, args)
    @test !DT.has_boundary_nodes(tri)
    @test DT.has_ghost_triangles(tri)

    tri = triangulate_rectangle(0, 1, 0, 1, 5, 5, delete_ghosts=true)
    args = DT.RefinementArguments(tri)
    @test DT.has_boundary_nodes(tri)
    @test DT.has_ghost_triangles(tri)
    DT.finalise!(tri, args)
    @test DT.has_boundary_nodes(tri)
    @test !DT.has_ghost_triangles(tri)

    tri = triangulate_rectangle(0, 1, 0, 1, 5, 5)
    args = DT.RefinementArguments(tri)
    @test DT.has_boundary_nodes(tri)
    @test DT.has_ghost_triangles(tri)
    DT.finalise!(tri, args)
    @test DT.has_boundary_nodes(tri)
    @test DT.has_ghost_triangles(tri)

    tri = triangulate(rand(2, 50), delete_ghosts=true)
    args = DT.RefinementArguments(tri)
    @test DT.has_boundary_nodes(tri)
    @test DT.has_ghost_triangles(tri)
    DT.finalise!(tri, args)
    @test !DT.has_boundary_nodes(tri)
    @test !DT.has_ghost_triangles(tri)
end

@testset "chew's second algorithm - deleting free vertices from diametral circle" begin
    points = [(0.0, 0.00000006667778), (10.0, 0.0), (10.0, 5.077777777), (0.0000209991, 5.0001991949), (2.0, 3.0), (8.0, 3.0)]
    boundary_nodes = [1, 2, 3, 4, 1]
    segments = Set([(5, 6)])
    tri = triangulate(points; boundary_nodes, segments)
    args = DT.RefinementArguments(tri)
    for p in [(4.0, 2.0), (5.0, 2.5), (3.0, 3.5), (5.0, 4.5), (6.5, 3.5), (7.0, 4.0), (6.5, 1.5), (4.0, 3.5)]
        add_point!(tri, p)
    end
    fig = Figure(fontsize=33)
    ax = Axis(fig[1, 1], width=600, height=400)
    _, _, C = compute_diametral_circle(get_point(tri, 5, 6)...)
    triplot!(ax, tri)
    lines!(ax, C, color=:red)
    DT.split_subsegment!(tri, args, (5, 6))
    ax = Axis(fig[1, 2], width=600, height=400)
    triplot!(ax, tri)
    lines!(ax, C, color=:red)
    resize_to_layout!(fig)
    fig
    @test_reference "delete_free_vertices_around_subsegment.png" fig
    @test validate_triangulation(tri)

    points = [(0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (0.0, 5.0), (2.0, 3.0), (8.0, 3.0), (3.5, 4.0), (9.0, 3.5)]
    boundary_nodes = [1, 2, 3, 4, 1]
    segments = Set([(5, 6), (7, 8)]) # make sure invisible vertices aren't deleted
    tri = triangulate(points; boundary_nodes, segments)
    args = DT.RefinementArguments(tri)
    for p in [(4.0, 2.0), (5.0, 2.5), (3.0, 3.5), (5.0, 4.5), (6.5, 3.5), (7.0, 4.0), (6.5, 1.5), (4.0, 3.5)]
        add_point!(tri, p)
    end
    fig = Figure(fontsize=33)
    ax = Axis(fig[1, 1], width=600, height=400)
    _, _, C = compute_diametral_circle(get_point(tri, 5, 6)...)
    triplot!(ax, tri)
    lines!(ax, C, color=:red)
    DT.split_subsegment!(tri, args, (5, 6))
    ax = Axis(fig[1, 2], width=600, height=400)
    triplot!(ax, tri)
    lines!(ax, C, color=:red)
    resize_to_layout!(fig)
    fig
    @test_reference "delete_free_vertices_around_subsegment2.png" fig
end

@testset verbose = true "refinement testing with circumcenters" begin
    _rng_num(idx1, idx2, idx3, idx4, idx5) = 2^idx1 * 3^idx2 * 5^idx3 * 7^idx4 * 11^idx5

    @testset "A simple convex example" begin
        for (idx1, use_lens) in enumerate((false, true))
            for (idx2, min_angle) in enumerate((20.0, 27.5, 30.0))
                for (idx3, min_area) in enumerate((1e-12,))
                    for (idx4, max_area) in enumerate((1e-1, 1e-2))
                        for (idx5, seditious_angle) in enumerate((10.0, 20.0))
                            @info "Testing refinement of a simple convex example. use_lens: $use_lens; min_angle: $min_angle; min_area: $min_area; max_area: $max_area; seditious_angle: $seditious_angle."
                            rng = StableRNG(_rng_num(idx1, idx2, idx3, idx4, idx5))
                            points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
                            tri = triangulate(points; boundary_nodes=[1, 2, 3, 4, 1], rng)
                            custom_constraint = (_tri, T) -> begin
                                i, j, k = triangle_vertices(T)
                                p, q, r = get_point(_tri, i, j, k)
                                c = (p .+ q .+ r) ./ 3
                                y = gety(c)
                                return y < 1 / 2 && DT.triangle_area(p, q, r) > 1e-3
                            end
                            refine!(tri; min_angle, min_area, max_area, custom_constraint, seditious_angle, use_circumcenter=true, use_lens, rng)
                            args = DT.RefinementArguments(tri; min_angle, min_area, max_area, seditious_angle, custom_constraint, use_circumcenter=true, use_lens)
                            @test validate_refinement(tri, args)
                            if _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(1, 3, 1, 2, 2)
                                fig, ax, sc = triplot(tri)
                                @test_reference "a_simple_convex_example_circle.png" fig by = psnr_equality(15)
                            elseif _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(2, 3, 1, 2, 2)
                                fig, ax, sc = triplot(tri)
                                @test_reference "a_simple_convex_example_lens.png" fig by = psnr_equality(15)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Triangulation with a hole" begin
        for (idx1, use_lens) in enumerate((false, true))
            for (idx2, min_angle) in enumerate((12.9, 27.5, 30.0))
                for (idx3, min_area) in enumerate((1e-12,))
                    for (idx4, max_area) in enumerate((1e-2, 1e-3, 1e-4))
                        for (idx5, seditious_angle) in enumerate((20.0, 40.0))
                            @info "Testing refinement of a triangulation with a hole. use_lens: $use_lens; min_angle: $min_angle; min_area: $min_area; max_area: $max_area; seditious_angle: $seditious_angle."
                            rng = StableRNG(_rng_num(idx1, idx2, idx3, idx4, idx5))
                            points = [(0.0, 0.0), (0.5, 0.1), (1.0, 0.0), (0.9, 0.5), (1.0, 1.0), (0.5, 0.9), (0.0, 1.0),
                                (0.3, 0.3), (0.7, 0.3), (0.7, 0.7), (0.3, 0.7)]
                            boundary_nodes = [[[1, 2, 3, 4, 5, 6, 7, 1]], [[11, 10, 9, 8, 11]]]
                            tri = triangulate(points; boundary_nodes, rng)
                            refine!(tri; min_angle, min_area, max_area, use_circumcenter=true, use_lens, rng, seditious_angle)
                            @test validate_refinement(tri, DT.RefinementArguments(tri; min_angle, min_area, seditious_angle, max_area, use_circumcenter=true, use_lens))
                            if _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(1, 3, 1, 3, 1)
                                fig, ax, sc = triplot(tri)
                                @test_reference "triangulation_with_a_hole_circle.png" fig by = psnr_equality(15)
                            elseif _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(2, 3, 1, 3, 1)
                                fig, ax, sc = triplot(tri)
                                @test_reference "triangulation_with_a_hole_lens.png" fig by = psnr_equality(15)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "A very non-convex example" begin
        for (idx1, use_lens) in enumerate((false, true))
            for (idx2, min_angle) in enumerate((12.9, 15.8, 30.0))
                for (idx3, min_area) in enumerate((1e-12,))
                    for (idx4, max_area) in enumerate((1e-1, 1e-3, 1e-4))
                        for (idx5, seditious_angle) in enumerate((10.0, 20.0, 60.0))
                            @info "Testing refinement of a very non-convex example. use_lens: $use_lens; min_angle: $min_angle; min_area: $min_area; max_area: $max_area; seditious_angle: $seditious_angle."
                            rng = StableRNG(_rng_num(idx1, idx2, idx3, idx4, idx5))
                            points = [(0.0, 0.0), (0.0, 1.0), (0.1, 1.0),
                                (0.1, 0.1), (0.2, 0.1), (0.2, 1.0),
                                (0.3, 1.0), (0.3, 0.0), (0.4, 0.0),
                                (0.4, 0.9), (0.5, 0.9), (0.5, 0.1),
                                (0.6, 0.1), (0.6, 1.0), (0.7, 1.0),
                                (0.7, 0.1), (0.8, 0.1), (0.8, 0.9),
                                (0.9, 0.9), (0.9, 0.1), (0.9, 0.0),
                                (1.0, 0.0), (1.0, -0.1), (0.0, -0.1),
                                (0.0, 0.0)] |> reverse!
                            boundary_nodes, points = convert_boundary_points_to_indices(points)
                            tri = triangulate(points; boundary_nodes, rng)
                            refine!(tri; min_angle, min_area, max_area, use_circumcenter=true, use_lens, seditious_angle, rng)
                            @test validate_refinement(tri; min_angle, min_area, max_area, seditious_angle, use_circumcenter=true, use_lens, check_conformal=false)
                            if _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(1, 3, 1, 3, 2)
                                fig, ax, sc = triplot(tri)
                                @test_reference "a_very_non_convex_example_circle.png" fig
                            elseif _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(2, 3, 1, 3, 2)
                                fig, ax, sc = triplot(tri)
                                @test_reference "a_very_non_convex_example_lens.png" fig
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "A constrained triangulation with multiple holes" begin
        for (idx1, use_lens) in enumerate((false, true))
            for (idx2, min_angle) in enumerate((12.9, 15.8, 30.0))
                for (idx3, min_area) in enumerate((1e-12,))
                    for (idx4, max_area) in enumerate((1e-2, 1e-3, 1e-4))
                        for (idx5, seditious_angle) in enumerate((20.0, 40.0))
                            @info "Testing refinement of a triangulation with multiple holes. use_lens: $use_lens; min_angle: $min_angle; min_area: $min_area; max_area: $max_area; seditious_angle: $seditious_angle."
                            rng = StableRNG(_rng_num(idx1, idx2, idx3, idx4, idx5))
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
                            tri = triangulate(points; boundary_nodes, rng)
                            _max_area = max_area * get_area(tri)
                            refine!(tri; min_angle, min_area, max_area=_max_area, use_circumcenter=true, seditious_angle, use_lens, rng)
                            @test validate_refinement(tri; min_angle, min_area, max_area=_max_area, use_circumcenter=true, seditious_angle, use_lens, check_conformal=false)
                            if _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(1, 3, 1, 3, 1)
                                fig, ax, sc = triplot(tri)
                                @test_reference "a_constrained_triangulation_with_multiple_holes_circle.png" fig by = psnr_equality(15)
                            elseif _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(2, 3, 1, 3, 1)
                                fig, ax, sc = triplot(tri)
                                @test_reference "a_constrained_triangulation_with_multiple_holes_lens.png" fig by = psnr_equality(15)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "another square example, directly testing stats" begin
        for (idx1, use_lens) in enumerate((false, true))
            for (idx2, min_angle) in enumerate((12.9, 15.8, 27.5, 30.0))
                for (idx3, min_area) in enumerate((1e-12,))
                    for (idx4, max_area) in enumerate((1e-1, 1e-2, 1e-4))
                        for (idx5, seditious_angle) in enumerate((10.0, 20.0))
                            @info "Testing refinement of a square. use_lens: $use_lens; min_angle: $min_angle; min_area: $min_area; max_area: $max_area; seditious_angle: $seditious_angle."
                            rng = StableRNG(_rng_num(idx1, idx2, idx3, idx4, idx5))
                            p1 = (0.0, 0.0)
                            p2 = (1.0, 0.0)
                            p3 = (0.0, 1.0)
                            p4 = (1.0, 1.0)
                            pts = [p1, p2, p3, p4]
                            tri = triangulate(pts; rng)
                            refine!(tri; min_angle, min_area, max_area, use_circumcenter=true, seditious_angle, use_lens, rng)
                            stats = statistics(tri)
                            @test DT.get_smallest_angle(stats) ≥ deg2rad(min_angle)
                            @test DT.get_largest_area(stats) ≤ max_area
                            @test DT.get_smallest_area(stats) ≥ min_area
                            @test !DT.is_constrained(tri)
                            @test DT.convex_hull(tri).vertices == DT.convex_hull(tri.points).vertices
                            @test validate_triangulation(tri)
                            validate_statistics(tri)
                            @test validate_refinement(tri; min_angle, min_area, max_area, seditious_angle, use_circumcenter=true, use_lens)
                            if _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(1, 4, 1, 3, 2)
                                fig, ax, sc = triplot(tri)
                                @test_reference "another_square_example_circle.png" fig by = psnr_equality(15)
                            elseif _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(2, 4, 1, 3, 2)
                                fig, ax, sc = triplot(tri)
                                @test_reference "another_square_example_lens.png" fig by = psnr_equality(15)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "algorithm terminates appropriately" begin
        for _ in 1:10
            p1 = (0.0, 0.0)
            p2 = (1.0, 0.0)
            p3 = (0.0, 1.0)
            p4 = (1.0, 1.0)
            pts = [p1, p2, p3, p4]
            tri = triangulate(pts)
            refine!(tri; max_area=1e-6, max_points=5_000, use_circumcenter=true)
            @test DT.num_solid_vertices(tri) == 5000
            @test validate_triangulation(tri)
            @test !validate_refinement(tri; max_area=1e-6, max_points=5_000, use_circumcenter=true, warn=false)
        end
    end

    @testset "avoiding infinite bouncing with concentric circular shells" begin
        for (idx1, use_lens) in enumerate((false, true))
            for (idx2, min_angle) in enumerate((23.0, 27.5, 30.0))
                for (idx3, min_area) in enumerate((1e-12,))
                    for (idx4, max_area) in enumerate((1e-1, 1e-2))
                        for (idx5, seditious_angle) in enumerate((10.0, 20.0))
                            @info "Testing that infinite bouncing is avoided during refinement. use_lens: $use_lens; min_angle: $min_angle; min_area: $min_area; max_area: $max_area; seditious_angle: $seditious_angle."
                            p1 = (0.0, 0.0)
                            p2 = (1.0, 0.0)
                            p3 = (0.0, 0.7)
                            p4 = (1.0, 0.4)
                            pts = [p1, p2, p3, p4]
                            rng = StableRNG(_rng_num(idx1, idx2, idx3, idx4, idx5))
                            tri = triangulate(pts; rng)
                            add_segment!(tri, 1, 4; rng)
                            refine!(tri; min_angle, min_area, max_area, seditious_angle, use_circumcenter=true, use_lens, rng)
                            @test validate_refinement(tri; min_angle, min_area, max_area, seditious_angle, use_circumcenter=true, use_lens)
                            validate_statistics(tri)
                            if _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(1, 3, 1, 2, 2)
                                fig, ax, sc = triplot(tri)
                                @test_reference "avoid_bouncing_concentric_circular_shells_circle.png" fig by = psnr_equality(15)
                            elseif _rng_num(idx1, idx2, idx3, idx4, idx5) == _rng_num(2, 3, 1, 2, 2)
                                fig, ax, sc = triplot(tri)
                                @test_reference "avoid_bouncing_concentric_circular_shells_lens.png" fig by = psnr_equality(15)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "some seditious edge testing and default calls" begin
        for i in 1:10
            @info "Testing seditious edges. Run: $i"
            points = [(0.0, 0.0), (1.0, 0.0), (0.8, 0.2)]
            tri = triangulate(points)
            refine!(tri, use_circumcenter=true)
            @test_reference "seditious_edge_testing_1.png" triplot(tri) by = psnr_equality(15)
            @test validate_refinement(tri, use_circumcenter=true)
            validate_statistics(tri)

            points = [(0.0, 0.0), (1.0, 0.0), (0.8, 0.2)]
            tri = triangulate(points)
            refine!(tri, use_circumcenter=true, use_lens=false)
            @test_reference "seditious_edge_testing_2.png" triplot(tri) by = psnr_equality(15)
            @test validate_refinement(tri, use_circumcenter=true, use_lens=false)
            validate_statistics(tri)

            points = [(0.0, 0.0), (1.0, 0.0), (0.8, 0.2)]
            tri = triangulate(points)
            args = DT.RefinementArguments(tri, use_circumcenter=true, use_lens=false, max_area=1e-3get_area(tri))
            refine!(tri, args)
            @test validate_refinement(tri, args)
            @test_reference "seditious_edge_testing_3.png" triplot(tri) by = psnr_equality(15)
            validate_statistics(tri)
        end
    end

    @testset "Triangulating with an interior hole" begin
        @info "Testing triangulation of an interior hole"
        p1 = (0.0, 0.0)
        p2 = (1.0, 0.0)
        p3 = (1.0, 1.0)
        p4 = (0.0, 1.0)
        p5 = (0.4, 0.4)
        p6 = (0.6, 0.4)
        p7 = (0.6, 0.6)
        p8 = (0.4, 0.6)
        pts = [p1, p2, p3, p4, p5, p6, p7, p8]
        boundary_nodes = [[[1, 2, 3, 4, 1]], [[8, 7, 6, 5, 8]]]
        rng = StableRNG(123)
        tri = triangulate(pts; rng, boundary_nodes=boundary_nodes, delete_ghosts=false)
        add_point!(tri, 0.1, 0.8; rng)
        add_point!(tri, 0.3, 0.2; rng)
        add_point!(tri, 0.7, 0.2; rng)
        add_point!(tri, 0.9, 0.8; rng)
        add_segment!(tri, 9, 10; rng)
        add_segment!(tri, 11, 12; rng)
        add_segment!(tri, 9, 12; rng)
        add_segment!(tri, 10, 11; rng)
        refine!(tri; max_area=0.001, max_points=5000, use_circumcenter=true, use_lens=true, rng)
        stats = statistics(tri)
        @test DT.get_smallest_angle(stats) ≥ deg2rad(30)
        @test DT.get_largest_area(stats) ≤ 0.001
        @test DT.is_constrained(tri)
        @test DT.convex_hull(tri).vertices == DT.convex_hull(tri.points).vertices
        validate_statistics(tri, stats)
        @test validate_triangulation(tri)
        @test validate_refinement(tri; max_area=0.001, max_points=5000, use_circumcenter=true, use_lens=true)
        fig, ax, sc = triplot(tri)
        @test_reference "triangulating_with_an_interior_hole.png" fig by = psnr_equality(15)
    end

    @testset "Refining disjoint sets" begin
        @info "Testing the refinement of disjoint sets"
        θ = LinRange(0, 2π, 30) |> collect
        θ[end] = 0
        xy = Vector{Vector{Vector{NTuple{2,Float64}}}}()
        rng = StableRNG(123)
        cx = 0.0
        for i in 1:2
            push!(xy, [[(cx + cos(θ), sin(θ)) for θ in θ]])
            push!(xy, [[(cx + 0.5cos(θ), 0.5sin(θ)) for θ in reverse(θ)]])
            cx += 3.0
        end
        boundary_nodes, points = convert_boundary_points_to_indices(xy)
        tri = triangulate(points; boundary_nodes=boundary_nodes, rng)
        A = DT.get_area(tri)
        max_area = 0.001A
        min_area = 1e-9A
        refine!(tri; min_area, rng, max_area, use_circumcenter=true)
        stats = statistics(tri)
        validate_statistics(tri)
        @test validate_refinement(tri; min_area, max_area, use_circumcenter=true)
        fig, ax, sc = triplot(tri)
        @test_reference "refining_disjoint_sets.png" fig by = psnr_equality(15)
    end

    @testset "Small angles" begin
        for PT in (DT.Exact, DT.Adaptive)
            ps = 0
            fig = Figure(fontsize=52)
            for i in 1:12
                @info "Testing refinement with small angles. Run: $i; Predicates: $PT"
                if i > 6
                    use_lens = true
                    i -= 6
                    ax = Axis(fig[2, i], title="Lens; $i", width=600, height=600)
                else
                    use_lens = false
                    ax = Axis(fig[1, i], title="Circle; $i", width=600, height=600)
                end
                hidedecorations!(ax)
                hidespines!(ax)
                rng = StableRNG(i)
                p1 = (0.0, 0.0)
                p2 = (1.0, 0.0)
                p3 = (0.0, 1.0)
                p4 = (1.0, 1.0)
                p5 = (0.5, 0.5)
                pts = [p1, p2, p3, p4, p5]
                C = Set{NTuple{2,Int}}()
                for i in 1:20
                    θ = 2π * rand(rng)
                    r = 0.5sqrt(rand(rng))
                    x = 0.5 + r * cos(θ)
                    y = 0.5 + r * sin(θ)
                    push!(pts, (x, y))
                    push!(C, (5, 5 + i))
                end
                tri = triangulate(pts; rng, boundary_nodes=[1, 2, 4, 3, 1], segments=C, predicates=PT())
                refine!(tri; min_angle=27.3, min_area=0.0, use_circumcenter=true, rng, use_lens, predicates=PT())
                stats = statistics(tri)
                ps += DT.get_largest_angle(stats) ≤ max(π - 2 * deg2rad(17.0), 2asin((sqrt(3) - 1) / 2)) # Corollary 8 of "When and Why Ruppert's Algorithm Works. In Twelfth International Meshing Roundtable, pp. 91–102, Santa Fe, NM, Sept 2003."
                validate_statistics(tri)
                @test validate_refinement(tri; min_angle=27.3, min_area=0.0, use_circumcenter=true, warn=false, use_lens, rng)
                triplot!(ax, tri)
            end
            resize_to_layout!(fig)
            @test_reference "refinement_with_small_angles.png" fig
            @test ps == 12
        end
    end

    @testset "Another disjoint domain" begin
        @info "Testing refinement of another disjoint domain"
        rng = StableRNG(123)
        t = (collect ∘ LinRange)(0, 2π, 50)
        t[end] = 0
        r = 2 .* (1 .+ cos.(t))
        x = r .* cos.(t)
        y = r .* sin.(t)
        points = [(x, y) for (x, y) in zip(x, y)]
        x = -(x .+ 3)
        inverse_points = [(x, y) for (x, y) in zip(x, y)]
        reverse!(inverse_points)
        curve_points = [[points], [inverse_points]]
        boundary_nodes, points = convert_boundary_points_to_indices(curve_points)
        C = Set([[(50, i) for i in 98:-1:81]..., [(50, i) for i in 69:-1:51]...,
            [(10, i) for i in 36:49]...])
        C′ = empty(C)
        bem = DT.construct_boundary_edge_map(boundary_nodes)
        for c in C
            c′ = DT.reverse_edge(c)
            c ∉ keys(bem) && c′ ∉ keys(bem) && push!(C′, c)
        end
        tri = triangulate(points; boundary_nodes, rng, segments=C′)
        @test DT.is_disjoint(tri)
        @test validate_triangulation(tri)
        A = DT.get_area(tri)
        refine!(tri, max_area=1e-3A, min_area=1e-12A, rng=rng, use_circumcenter=true, use_lens=true)
        validate_statistics(tri)
        @test validate_refinement(tri, max_area=1e-3A, min_area=0.0, use_circumcenter=true, use_lens=true, warn=false)
        fig, ax, sc = triplot(tri)
        @test_reference "refining_disjoint_sets_2.png" fig by = psnr_equality(15)
    end

    @testset "Tight example with a single triangle boundary interior" begin
        @info "Testing refinement with a single triangle boundary interior"
        A = (0.0, 0.0)
        B = (0.0, 25.0)
        C = (5.0, 25.0)
        D = (5.0, 5.0)
        E = (10.0, 5.0)
        F = (10.0, 10.0)
        G = (25.0, 10.0)
        H = (25.0, 15.0)
        I = (10.0, 15.0)
        J = (10.0, 25.0)
        K = (45.0, 25.0)
        L = (45.0, 20.0)
        M = (40.0, 20.0)
        N = (40.0, 5.0)
        O = (45.0, 5.0)
        P = (45.0, 0.0)
        Q = (10.0, 0.0)
        R = (10.0, -5.0)
        S = (15.0, -5.0)
        T = (15.0, -10.0)
        U = (10.0, -10.0)
        V = (5.0, -10.0)
        W = (5.0, -5.0)
        Z = (5.0, 0.0)
        A1 = (5.0, 2.5)
        B1 = (10.0, 2.5)
        C1 = (38.0, 2.5)
        D1 = (38.0, 20.0)
        E1 = (27.0, 20.0)
        F1 = (27.0, 11.0)
        G1 = (27.0, 4.0)
        H1 = (2.0, 4.0)
        I1 = (2.0, 0.0)
        pts = [A, I1, H1, G1, F1, E1, D1, C1, B1, A1, Z, W, V, U, T, S, R, Q, P, O, N, M, L, K, J, I, H, G, F, E, D, C, B, A]
        J1 = (14.0, 6.0)
        K1 = (14.0, 8.0)
        L1 = (22.0, 8.0)
        inner_pts = [J1, K1, L1, J1]
        boundary_pts = [[pts], [inner_pts]]
        nodes, points = convert_boundary_points_to_indices(boundary_pts)
        rng = StableRNG(123)
        tri = triangulate(points; boundary_nodes=nodes, rng)
        @test validate_triangulation(tri)
        A = get_area(tri)
        refine!(tri; max_area=1e-3A, use_circumcenter=true, rng)
        validate_statistics(tri)
        @test validate_refinement(tri; max_area=1e-3A, rng, use_circumcenter=true, warn=true)
        fig, ax, sc = triplot(tri)
        @test_reference "tight_example_with_a_single_triangle_boundary_interior.png" fig by = psnr_equality(15)
    end

    @testset "A complicated example with tight walls and small angles" begin
        @info "Testing refinement of a complicated example with tight walls and small angles"
        rng = StableRNG(123456)
        A = (0.0, 0.0)
        B = (0.0, 25.0)
        C = (5.0, 25.0)
        D = (5.0, 5.0)
        E = (10.0, 5.0)
        F = (10.0, 10.0)
        G = (25.0, 10.0)
        H = (25.0, 15.0)
        I = (10.0, 15.0)
        J = (10.0, 25.0)
        K = (45.0, 25.0)
        L = (45.0, 20.0)
        M = (40.0, 20.0)
        N = (40.0, 5.0)
        O = (45.0, 5.0)
        P = (45.0, 0.0)
        Q = (10.0, 0.0)
        R = (10.0, -5.0)
        S = (15.0, -5.0)
        T = (15.0, -10.0)
        U = (10.0, -10.0)
        V = (5.0, -10.0)
        W = (5.0, -5.0)
        Z = (5.0, 0.0)
        A1 = (5.0, 2.5)
        B1 = (10.0, 2.5)
        C1 = (38.0, 2.5)
        D1 = (38.0, 20.0)
        E1 = (27.0, 20.0)
        F1 = (27.0, 11.0)
        G1 = (27.0, 4.0)
        H1 = (2.0, 4.0)
        I1 = (2.0, 0.0)
        pts = [A, I1, H1, G1, F1, E1, D1, C1, B1, A1, Z, W, V, U, T, S, R, Q, P, O, N, M, L, K, J, I, H, G, F, E, D, C, B, A]
        J1 = (17.0603265896789, 7.623652007194)
        K1 = (14.8552854162067, 6.5423337394336)
        L1 = (16.6998871670921, 6.9875824379232)
        M1 = (16.0, 6.0)
        N1 = (16.9755173137761, 6.6483453343121)
        O1 = (17.0391242707032, 4.8885528593294)
        P1 = (17.4207660122657, 6.4575244635308)
        Q1 = (17.6327892020226, 4.9945644542079)
        R1 = (22.6789411182379, 6.1818943168468)
        S1 = (21.8096460402344, 6.4787267825065)
        T1 = (26.0, 8.0)
        U1 = (15.0673086059636, 9.086612016517)
        W1 = (15.0, 8.5)
        Z1 = (17.7913089332764, 8.3603005983396)
        inner_pts = [Z1, W1, U1, T1, S1, R1, Q1, P1, O1, N1, M1, L1, K1, J1, Z1]
        boundary_pts = [[pts], [inner_pts]]
        nodes, points = convert_boundary_points_to_indices(boundary_pts)
        push!(points, (20.0, 20.0))
        C = Set{NTuple{2,Int}}()
        for i in 1:50
            θ = 2π * rand(rng)
            r = 4sqrt(rand(rng))
            x = 20 + r * cos(θ)
            y = 20 + r * sin(θ)
            push!(points, (x, y))
            push!(C, (48, 48 + i))
        end
        tri = triangulate(points; boundary_nodes=nodes, segments=C, rng)
        @test validate_triangulation(tri)
        A = get_area(tri)
        refine!(tri; max_area=1.0, rng, use_circumcenter=true)
        validate_statistics(tri)
        fig, ax, sc = triplot(tri)
        @test_reference "complicated_example_with_tight_walls_and_small_angles.png" fig by = psnr_equality(15)
    end

    if !(get(ENV, "CI", "false") == "true")
        @testset "Tasmania" begin
            for PT in (DT.Exact, DT.Adaptive)
                @info "Testing refinement of Tasmnia. Predicates: $PT"
                rng = StableRNG(123)
                tassy_path = joinpath(dirname(dirname(pathof(DelaunayTriangulation))), "test", "tassy.txt")
                tassy = readdlm(tassy_path)
                ymax = @views maximum(tassy[:, 2])
                tassy = [(x, ymax - y) for (x, y) in eachrow(tassy)]
                reverse!(tassy)
                unique!(tassy)
                push!(tassy, tassy[begin])
                boundary_nodes, points = convert_boundary_points_to_indices(tassy)
                tri = triangulate(points; rng, boundary_nodes=boundary_nodes, predicates=PT())
                A = get_area(tri)
                refine!(tri; max_area=1e-2A, use_circumcenter=true, rng, predicates=PT())
                validate_statistics(tri)
                @test validate_refinement(tri; max_area=1e-2A, rng, use_circumcenter=true, warn=false)
                fig, ax, sc = triplot(tri)
                @test_reference "tasmania.png" fig by = psnr_equality(15)
            end
        end
    end

    @testset "Test that nothing is breaking for Float32 inputs" begin
        for PT in (DT.Exact, DT.Adaptive)
            for i in 1:25
                p1 = (0.0f0, 0.0f0)
                p2 = (1.0f0, 0.0f0)
                p3 = (0.0f0, 1.0f0)
                p4 = (1.0f0, 1.0f0)
                pts = [p1, p2, p3, p4]
                tri = triangulate(pts; delete_ghosts=false, predicates=PT())
                validate_statistics(tri)
                refine!(tri; max_area=0.001, max_points=25_000, use_circumcenter=true, predicates=PT())
                stats = statistics(tri)
                @inferred statistics(tri)
                @test DT.get_smallest_angle(stats) ≥ deg2rad(30.0)
                @test DT.get_largest_area(stats) ≤ 0.001f0
                @test !DT.is_constrained(tri)
                @test DT.convex_hull(tri).vertices == DT.convex_hull(tri.points).vertices
                @test validate_triangulation(tri)
                @test validate_refinement(tri; max_area=0.001, max_points=25_000, use_circumcenter=true)
            end

            for _ in 1:5
                rng = StableRNG(123)
                p1 = (0.0f0, 0.0f0)
                p2 = (1.0f0, 0.0f0)
                p3 = (1.0f0, 1.0f0)
                p4 = (0.0f0, 1.0f0)
                p5 = (0.4f0, 0.4f0)
                p6 = (0.6f0, 0.4f0)
                p7 = (0.6f0, 0.6f0)
                p8 = (0.4f0, 0.6f0)
                pts = [p1, p2, p3, p4, p5, p6, p7, p8]
                boundary_nodes = [[[1, 2, 3, 4, 1]], [[8, 7, 6, 5, 8]]]
                tri = triangulate(pts; rng, boundary_nodes=boundary_nodes, delete_ghosts=false, predicates=PT())
                add_point!(tri, 0.1f0, 0.8f0, predicates=PT())
                add_point!(tri, 0.3f0, 0.2f0, predicates=PT())
                add_point!(tri, 0.7f0, 0.2f0, predicates=PT())
                add_point!(tri, 0.9, 0.8, predicates=PT())
                add_segment!(tri, 9, 10, predicates=PT())
                add_segment!(tri, 11, 12, predicates=PT())
                add_segment!(tri, 9, 12, predicates=PT())
                add_segment!(tri, 10, 11, predicates=PT())
                refine!(tri; rng, max_area=0.001f0, max_points=25000, min_angle=24.0f0, min_area=0.0, use_circumcenter=true, predicates=PT())
                stats = statistics(tri)
                @test DT.get_smallest_angle(stats) ≥ deg2rad(24.0f0)
                @test DT.get_largest_area(stats) ≤ 0.001f0
                @test DT.is_constrained(tri)
                @test DT.convex_hull(tri).vertices == DT.convex_hull(tri.points).vertices
                @test validate_refinement(tri; max_area=0.001, max_points=25000, min_angle=24.0f0, min_area=0.0, use_circumcenter=true)
            end
        end
    end

    @testset "Custom refinement example" begin
        @info "Testing refinement with custom constraints"
        rng = StableRNG(123)
        points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        tri = triangulate(points; rng)
        function get_quadrant(p)
            px, py = p
            if px < 0.5 && py < 0.5
                return 1
            elseif px ≥ 0.5 && py < 0.5
                return 2
            elseif px ≥ 0.5 && py ≥ 0.5
                return 3
            else
                return 4
            end
        end
        function custom_refine_ex(_tri, T)
            p, q, r = get_point(_tri, T...)
            c = (p .+ q .+ r) ./ 3
            quad = get_quadrant(c)
            A = DT.triangle_area(p, q, r)
            ρ = DT.triangle_radius_edge_ratio(p, q, r)
            θ = asind(inv(2ρ))
            flag = norm(p .- q) > 0.3 || norm(p .- r) > 0.3 || norm(q .- r) > 0.3
            flag && return true
            if quad == 1
                return A ≥ 0.001 || θ ≤ 15
            elseif quad == 2
                return A ≥ 0.01 || θ ≤ 25
            elseif quad == 3
                return A ≥ 0.005 || θ ≤ 7
            else
                return A ≥ 0.008 || θ ≤ 30.0
            end
        end
        refine!(tri; custom_constraint=custom_refine_ex, rng, use_circumcenter=true)
        stat = Bool[]
        delete_ghost_triangles!(tri)
        for T in each_triangle(tri)
            push!(stat, custom_refine_ex(tri, T))
        end
        @test !any(stat)
        @test validate_refinement(tri; custom_constraint=custom_refine_ex, use_circumcenter=true)
        validate_statistics(tri)
    end

    @testset "Julia logo" begin
        @info "Testing refinement of the Julia logo"
        rng = StableRNG(123)
        C = (15.7109521325776, 33.244486807457)
        D = (14.2705719699703, 32.8530791545746)
        E = (14.3, 27.2)
        F = (14.1, 27.0)
        G = (13.7, 27.2)
        H = (13.4, 27.5)
        I = (13.1, 27.6)
        J = (12.7, 27.4)
        K = (12.5, 27.1)
        L = (12.7, 26.7)
        M = (13.1, 26.5)
        N = (13.6, 26.4)
        O = (14.0, 26.4)
        P = (14.6, 26.5)
        Q = (15.1983491346581, 26.8128534095401)
        R = (15.6, 27.6)
        S = (15.6952958264624, 28.2344688505621)
        T = (17.8088971520274, 33.1192363585346)
        U = (16.3058917649589, 33.0722674401887)
        V = (16.3215480710742, 29.7374742376305)
        W = (16.3841732955354, 29.393035503094)
        Z = (16.6190178872649, 28.9233463196351)
        A1 = (17.0417381523779, 28.5319386667527)
        B1 = (17.5114273358368, 28.3753756055997)
        C1 = (18.1376795804487, 28.3597192994844)
        D1 = (18.7169629067146, 28.5632512789833)
        E1 = (19.2805899268653, 28.8920337074045)
        F1 = (19.26493362075, 28.4536571361762)
        G1 = (20.6426885588962, 28.4223445239456)
        H1 = (20.689657477242, 33.1035800524193)
        I1 = (19.2805899268653, 33.0722674401887)
        J1 = (19.2962462329806, 29.7531305437458)
        K1 = (19.0614016412512, 29.393035503094)
        L1 = (18.7482755189452, 29.236472441941)
        M1 = (18.4508057027546, 29.1425346052493)
        N1 = (18.1689921926793, 29.3147539725175)
        O1 = (17.7932408459121, 29.6278800948235)
        P1 = (22.6466957416542, 35.4207133574833)
        Q1 = (21.2219718851621, 34.9979930923702)
        R1 = (21.2376281912774, 28.4693134422915)
        S1 = (22.6780083538847, 28.4380008300609)
        T1 = (24.5724213938357, 33.1975178891111)
        U1 = (23.3512295168425, 32.8530791545746)
        V1 = (23.3199169046119, 28.4380008300609)
        W1 = (24.6663592305274, 28.3753756055997)
        Z1 = (15.1942940307729, 35.4363696635986)
        A2 = (14.7246048473139, 35.3737444391374)
        B2 = (14.3645098066621, 35.1858687657538)
        C2 = (14.1766341332786, 34.8570863373326)
        D2 = (14.1140089088174, 34.3247719294125)
        E2 = (14.2705719699703, 33.8394264398383)
        F2 = (14.7246048473139, 33.6202381542241)
        G2 = (15.4604512347329, 33.6045818481088)
        H2 = (16.0, 34.0)
        I2 = (15.9771093365377, 34.6848669700643)
        J2 = (15.6170142958859, 35.2328376840997)
        K2 = (24.1653574348379, 35.4520259697138)
        L2 = (23.7739497819555, 35.4363696635986)
        M2 = (23.4608236596496, 35.2641502963303)
        N2 = (23.272947986266, 34.9040552556785)
        O2 = (23.1320412312284, 34.5909291333725)
        P2 = (23.1163849251131, 34.2151777866054)
        Q2 = (23.2886042923813, 33.8081138276077)
        R2 = (23.8209187003014, 33.6045818481088)
        S2 = (24.3062641898756, 33.5576129297629)
        T2 = (24.7602970672192, 33.8550827459536)
        U2 = (25.010797965064, 34.4656786844502)
        V2 = (24.8385785977957, 34.9666804801397)
        W2 = (24.5254524754898, 35.2641502963303)
        Z2 = (25.3708930057158, 37.4716894585871)
        A3 = (24.7916096794498, 37.3464390096648)
        B3 = (24.4471709449133, 36.9550313567823)
        C3 = (24.3062641898756, 36.5636237038999)
        D3 = (24.4941398632592, 35.9999966837492)
        E3 = (25.0264542711793, 35.5929327247515)
        F3 = (25.5587686790994, 35.5929327247515)
        F3 = (25.5587686790994, 35.5929327247515)
        G3 = (26.0, 36.0)
        H3 = (26.1380520053653, 36.5792800100152)
        I3 = (26.0, 37.0)
        J3 = (25.7466443524829, 37.2838137852036)
        K3 = (26.3885529032101, 35.4676822758291)
        L3 = (25.9814889442124, 35.3580881330221)
        M3 = (25.6840191280217, 35.1858687657538)
        N3 = (25.5274560668688, 34.9040552556785)
        O3 = (25.4961434546382, 34.5596165211419)
        P3 = (25.5274560668688, 34.246490398836)
        Q3 = (25.6683628219064, 33.8394264398383)
        R3 = (26.0284578625583, 33.6358944603394)
        S3 = (26.5451159643631, 33.6202381542241)
        T3 = (27.0, 34.0)
        U3 = (27.280962351782, 34.5596165211419)
        V3 = (27.0304614539373, 35.2171813779844)
        W3 = (26.1693646175959, 33.087923746304)
        Z3 = (26.0, 33.0)
        A4 = (25.5274560668688, 32.7278287056522)
        B4 = (25.2612988629087, 32.4147025833463)
        C4 = (25.1830173323322, 32.0702638488098)
        D4 = (25.2299862506781, 31.7727940326191)
        E4 = (25.6527065157911, 31.5222931347744)
        F4 = (26.2946150665183, 31.7258251142732)
        G4 = (26.5607722704784, 32.5086404200381)
        H4 = (27.1557119028596, 32.7434850117675)
        I4 = (27.6097447802033, 32.4929841139228)
        J4 = (27.6410573924338, 32.1015764610403)
        K4 = (27.7193389230103, 31.6005746653509)
        L4 = (27.437525412935, 31.4283552980826)
        M4 = (26.9834925355914, 31.2561359308143)
        N4 = (26.5764285765937, 31.0995728696614)
        O4 = (26.0441141686736, 30.7864467473554)
        P4 = (25.6527065157911, 30.5672584617413)
        Q4 = (25.3239240873699, 30.1915071149741)
        R4 = (25.1673610262169, 29.8783809926682)
        S4 = (25.1047358017558, 29.6122237887082)
        T4 = (25.0890794956405, 29.1895035235952)
        U4 = (25.2926114751393, 28.8294084829433)
        V4 = (25.6840191280217, 28.5632512789833)
        W4 = (26.1537083114806, 28.3753756055997)
        Z4 = (26.8269294744384, 28.391031911715)
        A5 = (27.4844943312809, 28.6102201973292)
        B5 = (27.7342002330051, 28.7239579596219)
        C5 = (27.7264126450755, 28.4202565942047)
        D5 = (29.1825559185446, 28.3922538389457)
        E5 = (29.1545531632856, 32.2146299318021)
        F5 = (29.000538009361, 32.5786657501693)
        G5 = (28.6785063238822, 32.9006974356481)
        H5 = (28.3144705055149, 33.0827153448317)
        I5 = (27.9084305542591, 33.2367304987563)
        J5 = (27.3343740714492, 33.3207387645334)
        K5 = (26.8303244767868, 33.2367304987563)
        L5 = (27.6564057569279, 30.786489413592)
        M5 = (27.6984098898165, 30.3944508399657)
        N5 = (27.6984098898165, 29.7363860913787)
        O5 = (27.5863988687804, 29.4143544059)
        P5 = (27.2643671833016, 29.2043337414573)
        Q5 = (26.9843396307114, 29.1763309861983)
        R5 = (26.6903107004917, 29.3163447624934)
        S5 = (26.5782996794556, 29.7503874690082)
        T5 = (26.7603175886393, 30.3384453294476)
        U5 = (27.3203726938197, 30.7024811478149)
        J_curve = [[C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, C]]
        U_curve = [[T, U, V, W, Z, A1, B1, C1, D1, E1, F1, G1, H1, I1, J1, K1, L1, M1, N1, O1, T]]
        L_curve = [[P1, Q1, R1, S1, P1]]
        I_curve = [[T1, U1, V1, W1, T1]]
        A_curve_outline = [[
            K5, W3, Z3, A4, B4, C4, D4, E4, F4, G4, H4, I4, J4, K4, L4, M4, N4,
            O4, P4, Q4, R4, S4, T4, U4, V4, W4, Z4, A5, B5, C5, D5, E5, F5, G5,
            H5, I5, J5, K5]]
        A_curve_hole = [[L5, M5, N5, O5, P5, Q5, R5, S5, T5, U5, L5]]
        dot_1 = [[Z1, A2, B2, C2, D2, E2, F2, G2, H2, I2, J2, Z1]]
        dot_2 = [[Z2, A3, B3, C3, D3, E3, F3, G3, H3, I3, J3, Z2]]
        dot_3 = [[K2, L2, M2, N2, O2, P2, Q2, R2, S2, T2, U2, V2, W2, K2]]
        dot_4 = [[K3, L3, M3, N3, O3, P3, Q3, R3, S3, T3, U3, V3, K3]]
        curves = [J_curve, U_curve, L_curve, I_curve, A_curve_outline, A_curve_hole, dot_1, dot_2, dot_3, dot_4]
        nodes, points = convert_boundary_points_to_indices(curves)
        tri = triangulate(points; boundary_nodes=nodes, rng)
        A = get_area(tri)
        refine!(tri; min_angle=26.45, max_area=0.005A / 9, rng, use_circumcenter=true)
        validate_statistics(tri)
        @test validate_refinement(tri; min_angle=26.45, max_area=0.005A / 9, rng, use_circumcenter=true)
        fig, ax, sc = triplot(tri)
        @test_reference "julia_logo.png" fig by = psnr_equality(15)
    end
end