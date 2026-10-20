using Dates
using DelaunayTriangulation
using Aqua
using Test
using Random

const ALL_TEST_SCRIPTS = Set{String}()
const NON_TEST_SCRIPTS = Set{String}(["helper_functions.jl", "runtests.jl", (joinpath("helpers", file) for file in readdir(joinpath(@__DIR__, "helpers")))...])
include("helper_functions.jl")
using .HelperFunctions

# To only run some of the test groups below, pass their names as test arguments, e.g.
#   Pkg.test("DelaunayTriangulation"; test_args = ["Voronoi", "Refinement"])
# By default, all the groups are run.
const TEST_GROUPS = ARGS
run_group(name) = isempty(TEST_GROUPS) || name ∈ TEST_GROUPS
let unknown_groups = setdiff(TEST_GROUPS, ["Aqua", "Primitives", "Predicates", "Utilities", "Geometry", "Trees", "Curves", "Triangulation", "Point Location", "Refinement", "Statistics", "Voronoi", "Makie", "Docstring examples", "Run the documentation examples"])
    isempty(unknown_groups) || error("Unknown test groups: $(join(unknown_groups, ", ")).")
end

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
    run_group("Aqua") && @testset verbose = true "Aqua" begin
        Aqua.test_all(DelaunayTriangulation; ambiguities=false, project_extras=false, unbound_args=false) # don't care about julia < 1.2
        Aqua.test_ambiguities(DelaunayTriangulation) # don't pick up Base and Core...
    end

    run_group("Primitives") && @testset verbose = true "Primitives" begin
        safe_include("primitives/triangles.jl")
        safe_include("primitives/edges.jl")
        safe_include("primitives/points.jl")
        safe_include("primitives/boundary_nodes.jl")
    end

    run_group("Predicates") && @testset verbose = true "Predicates" begin
        safe_include("predicates/certificate.jl")
        safe_include("predicates/boundaries_and_ghosts.jl")
        safe_include("predicates/general.jl")
        safe_include("predicates/index_and_ghost_handling.jl")
    end

    run_group("Utilities") && @testset verbose = true "Utilities" begin
        safe_include("utils/utils.jl")
        safe_include("utils/max_priority_queue.jl")
        safe_include("utils/queue.jl")
        safe_include("helper_function_tests.jl")
    end

    run_group("Geometry") && @testset verbose = true "Geometry" begin
        safe_include("geometry/geometry_utils.jl")
        safe_include("geometry/convex_hull.jl")
    end

    run_group("Trees") && @testset verbose = true "Trees" begin
        safe_include("trees/rtree.jl")
        safe_include("trees/polygon_hierarchy.jl", verbose=false)
    end

    run_group("Curves") && @testset verbose = true "Curves" begin
        safe_include("curves/curves.jl")
    end

    run_group("Triangulation") && @testset verbose = true "Triangulation" begin
        safe_include("triangulation/adjacent.jl")
        safe_include("triangulation/adjacent2vertex.jl")
        safe_include("triangulation/graph.jl")
        safe_include("triangulation/triangulation.jl")
        safe_include("triangulation/representative.jl")
        safe_include("triangulation/triangulation_cache.jl")
        safe_include("triangulation/construction/rectangle.jl")
        safe_include("triangulation/construction/bowyer_watson.jl")
        safe_include("triangulation/construction/triangulate.jl")
        safe_include("triangulation/construction/shuffled_polygon_linked_list.jl")
        safe_include("triangulation/construction/convex_triangulation.jl")
        safe_include("triangulation/construction/segment_location.jl")
        safe_include("triangulation/construction/segment_insertion.jl")
        safe_include("triangulation/construction/constrained.jl")
        safe_include("triangulation/construction/check_args.jl")
        safe_include("triangulation/construction/weighted.jl")
        safe_include("triangulation/operations/add_triangle.jl")
        safe_include("triangulation/operations/delete_triangle.jl")
        safe_include("triangulation/operations/add_ghost_triangles.jl")
        safe_include("triangulation/operations/delete_ghost_triangles.jl")
        safe_include("triangulation/operations/add_point.jl")
        safe_include("triangulation/operations/flip_edge.jl")
        safe_include("triangulation/operations/split_triangle.jl")
        safe_include("triangulation/operations/split_edge.jl")
        safe_include("triangulation/operations/legalise_edge.jl")
        safe_include("triangulation/operations/delete_point.jl")
        safe_include("triangulation/operations/lock_convex_hull.jl")
        safe_include("triangulation/operations/delete_holes.jl")
        safe_include("triangulation/invariance.jl")
    end

    run_group("Point Location") && @testset verbose = true "Point Location" begin
        safe_include("point_location/brute_force.jl")
        safe_include("point_location/select_initial_point.jl")
        safe_include("point_location/select_initial_triangle_interior_node.jl")
        safe_include("point_location/interior_edge_intersections.jl")
        safe_include("point_location/ghost_search.jl")
        safe_include("point_location/jump_and_march.jl")
        safe_include("point_location/find_polygon.jl")
    end

    run_group("Refinement") && @testset verbose = true "Refinement" begin
        safe_include("refinement/refine.jl")
        safe_include("refinement/curve_bounded.jl")
    end

    run_group("Statistics") && @testset verbose = true "Statistics" begin
        safe_include("statistics/statistics.jl")
    end

    run_group("Voronoi") && @testset verbose = true "Voronoi" begin
        safe_include("voronoi/voronoi.jl")
        safe_include("voronoi/power.jl")
    end

    run_group("Makie") && @testset verbose = true "Makie" begin
        safe_include("makie/makie.jl")
    end

    run_group("Docstring examples") && @testset verbose = true "Docstring examples" begin
        safe_include("docstring_examples.jl")
    end

    run_group("Run the documentation examples") && @testset verbose = true "Run the documentation examples" begin
        @testset verbose = true "Check that the applications in the docs run" begin
            app_dir = joinpath(dirname(dirname(pathof(DelaunayTriangulation))), "docs", "src", "literate_applications")
            app_files = readdir(app_dir)
            for file in app_files
                safe_include(joinpath(app_dir, file); push=false)
            end
        end

        @testset verbose = true "Test the tutorials" begin
            tut_dir = joinpath(dirname(dirname(pathof(DelaunayTriangulation))), "docs", "src", "literate_tutorials")
            tut_files = readdir(tut_dir)
            for file in tut_files
                safe_include(joinpath(tut_dir, file); push=false)
            end
            # The weighted triangulation and power diagram tutorials record animations into the working directory
            for mp4_file in ("varying_weight.mp4", "varying_weight_power.mp4")
                isfile(mp4_file) && rm(mp4_file)
            end
        end

        @testset verbose = true "Test the readme example" begin
            safe_include("readme_example.jl")
        end
    end

    isempty(TEST_GROUPS) && @testset "All script files are included somewhere" begin
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