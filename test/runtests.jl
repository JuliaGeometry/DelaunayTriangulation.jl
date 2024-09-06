using Dates
using DelaunayTriangulation
using Aqua
using Test
using Random

if isdefined(Docs, :undocumented_names)
    @test isempty(Docs.undocumented_names(DelaunayTriangulation))
end

const ALL_TEST_SCRIPTS = Set{String}()
const NON_TEST_SCRIPTS = Set{String}(["helper_functions.jl", "runtests.jl"])
include("helper_functions.jl")
using .HelperFunctions

ct() = Dates.format(now(), "HH:MM:SS")
function safe_include(filename; name=filename, push=true, verbose=true) # Workaround for not being able to interpolate into SafeTestset test names
    push && push!(ALL_TEST_SCRIPTS, normpath(filename))
    mod = @eval module $(gensym()) end
    @info "[$(ct())] Testing $name"
    @testset verbose = verbose "$(basename(name))" begin
        Random.seed!(0)
        @eval mod using ..HelperFunctions
        @eval mod using ..Test
        Base.include(mod, filename)
    end
end

@testset verbose = true "DelaunayTriangulation.jl" begin
    @testset verbose = true "Aqua" begin
        Aqua.test_all(DelaunayTriangulation; ambiguities=false, project_extras=false) # don't care about julia < 1.2
        Aqua.test_ambiguities(DelaunayTriangulation) # don't pick up Base and Core...
    end

    @testset verbose = true "Triangulation" begin
        safe_include("triangulation/rectangle.jl")
        safe_include("triangulation/bowyer_watson.jl")
        safe_include("triangulation/triangulate.jl")
        safe_include("triangulation/convex_triangulation.jl")
        safe_include("triangulation/constrained.jl")
        safe_include("triangulation/check_args.jl")
        safe_include("triangulation/weighted.jl")
    end

    @testset verbose = true "Interfaces" begin
        safe_include("interfaces/triangles.jl")
        safe_include("interfaces/edges.jl")
        safe_include("interfaces/points.jl")
        safe_include("interfaces/boundary_nodes.jl")
    end

    @testset verbose = true "Data Structures" begin
        safe_include("data_structures/adjacent.jl")
        safe_include("data_structures/adjacent2vertex.jl")
        safe_include("data_structures/graph.jl")
        safe_include("data_structures/convex_hull.jl")
        safe_include("data_structures/triangulation.jl")
        safe_include("data_structures/representative.jl")
        safe_include("data_structures/statistics.jl")
        safe_include("data_structures/shuffled_polygon_linked_list.jl")
        safe_include("data_structures/triangulation_cache.jl")
        safe_include("data_structures/max_priority_queue.jl")
        safe_include("data_structures/queue.jl")
        safe_include("data_structures/curves.jl")
        safe_include("data_structures/rtree.jl")
        safe_include("data_structures/bst.jl")
        safe_include("data_structures/polygon_hierarchy.jl", verbose=false)
    end

    @testset verbose = true "Predicates" begin
        safe_include("predicates/certificate.jl")
        safe_include("predicates/boundaries_and_ghosts.jl")
        safe_include("predicates/general.jl")
        safe_include("predicates/index_and_ghost_handling.jl")
    end

    @testset verbose = true "Utilities" begin
        safe_include("utils.jl")
        safe_include("geo_utils.jl")
        safe_include("helper_function_tests.jl")
    end

    @testset verbose = true "Point Location" begin
        safe_include("point_location/brute_force.jl")
        safe_include("point_location/select_initial_point.jl")
        safe_include("point_location/select_initial_triangle_interior_node.jl")
        safe_include("point_location/interior_edge_intersections.jl")
        safe_include("point_location/ghost_search.jl")
        safe_include("point_location/jump_and_march.jl")
        safe_include("point_location/find_polygon.jl")
    end

    @testset verbose = true "Operations" begin
        safe_include("operations/add_triangle.jl")
        safe_include("operations/delete_triangle.jl")
        safe_include("operations/add_ghost_triangles.jl")
        safe_include("operations/delete_ghost_triangles.jl")
        safe_include("operations/add_point.jl")
        safe_include("operations/flip_edge.jl")
        safe_include("operations/split_triangle.jl")
        safe_include("operations/split_edge.jl")
        safe_include("operations/legalise_edge.jl")
        safe_include("operations/delete_point.jl")
        safe_include("operations/lock_convex_hull.jl")
        safe_include("operations/delete_holes.jl")
    end

    @testset verbose = true "Constrained Triangulation" begin
        safe_include("constrained_triangulation/segment_location.jl")
        safe_include("constrained_triangulation/segment_insertion.jl")
    end

    @testset verbose = true "Refinement" begin
        safe_include("refinement/refine.jl")
        safe_include("refinement/curve_bounded.jl")
    end

    @testset verbose = true "Voronoi" begin
        safe_include("voronoi/voronoi.jl")
    end

    @testset verbose = true "Run the documentation examples" begin
        @testset verbose = true "Check that the applications in the docs run" begin
            app_dir = joinpath(dirname(dirname(pathof(DelaunayTriangulation))), "docs", "src", "literate_applications")
            app_files = readdir(app_dir)
            for file in app_files
                safe_include(joinpath(app_dir, file); push=false)
            end
            mp4_path = joinpath(dirname(dirname(pathof(DelaunayTriangulation))), "cell_simulation.mp4")
            isfile(mp4_path) && rm(mp4_path)
        end

        @testset verbose = true "Test the tutorials" begin
            tut_dir = joinpath(dirname(dirname(pathof(DelaunayTriangulation))), "docs", "src", "literate_tutorials")
            tut_files = readdir(tut_dir)
            for file in tut_files
                safe_include(joinpath(tut_dir, file); push=false)
            end
        end

        @testset verbose = true "Test the readme example" begin
            safe_include("readme_example.jl")
        end
    end

    @testset "All script files are included somewhere" begin
        missing_set = String[]
        test_dir = joinpath(dirname(dirname(pathof(DelaunayTriangulation))), "test", "")
        for (root, dir, files) in walkdir(test_dir)
            for file in files
                filename = normpath(replace(joinpath(root, file), test_dir => ""))
                if endswith(filename, ".jl") && filename ∉ ALL_TEST_SCRIPTS && filename ∉ NON_TEST_SCRIPTS
                    push!(missing_set, filename)
                end
            end
        end
        if !isempty(missing_set)
            @info "There were some test scripts that were not included. These are printed below."
            for script in missing_set
                @info "     $script"
            end
        end
        @test isempty(missing_set)
    end
end

#=
using Random
DT = DelaunayTriangulation
rng = Random.seed!(2)
fi = 4 
tri, submerged, nonsubmerged, weights = get_weighted_example(fi)
points = get_points(tri)
tri = Triangulation(points; weights)
insertion_order = [6, 1, 3, 2, 4]
DT.initialise_bowyer_watson!(tri, insertion_order, predicates)
remaining_points = @view insertion_order[(begin + 3):end]
predicates = AdaptiveKernel()
randomise = true 
try_last_inserted_point = true 
skip_points = ()
num_sample_rule = DT.default_num_samples 
for (num_points, new_point) in enumerate(remaining_points)
    (num_points, new_point) == (3, 5) && continue
    initial_search_point = DT.get_initial_search_point(tri, num_points, new_point, insertion_order, num_sample_rule, rng, try_last_inserted_point)
    DT.add_point_bowyer_watson!(tri, new_point, initial_search_point, rng, predicates)
end
num_points, new_point = 3, 5
initial_search_point = DT.get_initial_search_point(tri, num_points, new_point, insertion_order, num_sample_rule, rng, try_last_inserted_point)
DT.add_point_bowyer_watson!(tri, new_point, initial_search_point, rng, predicates)



tri = triangulate(get_points(tri); weights, insertion_order = [6, 1, 3, 2, 4], rng)

insertion_order = [6, 1, 3, 2, 4, 5]
initial_search_point = DT.get_initial_search_point(tri, num_points, new_point, insertion_order, num_sample_rule, rng, try_last_inserted_point)


new_point = 5
point_indices = each_solid_vertex(tri)
m = DT.default_num_samples(length(point_indices))
try_points = ()
rng = Random.default_rng()
initial_search_point = DT.integer_type(tri)(DT.select_initial_point(tri, new_point; point_indices, m, try_points, rng))
update_representative_point = false
store_event_history = Val(false)
event_history = nothing
concavity_protection = false
V = find_triangle(tri, get_point(tri, new_point); m=nothing, point_indices=nothing, try_points=nothing, k=initial_search_point, concavity_protection, predicates, rng)
peek = Val(false)
int_flag = new_point isa Integer
q = get_point(tri, new_point)
cert = DT.point_position_relative_to_circumcircle(predicates, tri, V, new_point) # redirects to point_position_relative_to_witness_plane
flag = DT.point_position_relative_to_triangle(predicates, tri, V, q)
I = DT.integer_type(tri)
_new_point = DT.is_true(peek) ? new_point : I(new_point)
_new_point = DT.is_true(peek) ? num_points(tri) + 1 : new_point # If we are peeking, then we need to use the number of points in the triangulation as the index for the new point since we don't actually insert the point
i, j, k = triangle_vertices(V)
ℓ₁ = get_adjacent(tri, j, i)
ℓ₂ = get_adjacent(tri, k, j)
ℓ₃ = get_adjacent(tri, i, k)
!DT.is_true(peek) && DT.delete_triangle!(tri, V; protect_boundary=true, update_ghost_edges=false)
i, j, ℓ = i, j, ℓ₁
_r = DT.is_true(peek) ? num_points(tri) + 1 : r


DT.dig_cavity!(tri, new_point, i, j, ℓ₁, flag, V, store_event_history, event_history, peek, predicates)

DT.dig_cavity!(tri, new_point, j, k, ℓ₂, flag, V, store_event_history, event_history, peek, predicates)
DT.dig_cavity!(tri, new_point, k, i, ℓ₃, flag, V, store_event_history, event_history, peek, predicates)


tri = triangulate(get_points(tri); weights=zeros(i * j))
@test validate_triangulation(tri)
=#