# setup LocalPreferences.toml 
using Preferences 
USE_EXACTPREDICATES = ENV["USE_EXACTPREDICATES"]
if USE_EXACTPREDICATES == "true"
    set_preferences!("DelaunayTriangulation", "USE_EXACTPREDICATES" => true)
elseif USE_EXACTPREDICATES == "false" 
    set_preferences!("DelaunayTriangulation", "USE_EXACTPREDICATES" => false)
end # if USE_EXACTPREDICATES == "default", do nothing
@show DelaunayTriangulation.USE_EXACTPREDICATES

throw("...")

# get all the compilation out of the way
using BenchmarkTools
using CairoMakie
using ColorSchemes
using DataStructures
using Dates
using DelaunayTriangulation
using DelimitedFiles
using Documenter
using DocumenterTools
using Downloads
using ElasticArrays
using ExactPredicates
using FastGaussQuadrature
using ForwardDiff
using GeometryBasics
using InteractiveUtils
using LinearAlgebra
using OrderedCollections
using Random
using ReferenceTests
using SafeTestsets
using Setfield
using SimpleGraphs: SimpleGraphs
using SpatialIndexing
using StableRNGs
using StaticArrays
using StaticArraysCore
using StatsBase
using StructEquality
using Aqua
using Test
using Preferences

const ALL_TEST_SCRIPTS = Set{String}()
const NON_TEST_SCRIPTS = Set{String}(["helper_functions.jl", "triangulation_validation.jl", "runtests.jl", "triangulation\\weighted.jl", "triangulation/weighted.jl"])
include("helper_functions.jl")
using .HelperFunctions

ct() = Dates.format(now(), "HH:MM:SS")
function safe_include(filename; name=filename, push=true, verbose = true) # Workaround for not being able to interpolate into SafeTestset test names
    push && push!(ALL_TEST_SCRIPTS, normpath(filename))
    mod = @eval module $(gensym()) end
    @info "[$(ct())] Testing $name"
    @testset verbose = verbose "$name" begin
        @eval mod using ..HelperFunctions
        @eval mod using ..Test
        Base.include(mod, filename)
    end
end

@testset verbose = true "DelaunayTriangulation.jl" begin
    @testset verbose = true "Aqua" begin
        Aqua.test_all(DelaunayTriangulation; ambiguities=false, project_extras=false, stale_deps = load_preference(DelaunayTriangulation, "USE_EXACTPREDICATES", true)) # don't care about julia < 1.2
        Aqua.test_ambiguities(DelaunayTriangulation) # don't pick up Base and Core...
    end    

    @testset verbose = true "Triangulation" begin
        safe_include("triangulation/rectangle.jl")
        safe_include("triangulation/bowyer_watson.jl")
        safe_include("triangulation/triangulate.jl")
        safe_include("triangulation/convex_triangulation.jl")
        safe_include("triangulation/constrained.jl")
        safe_include("triangulation/check_args.jl")
        # Needs to be setup properly before we do more testing of it. Please see the 
        # comments at the end of the weighted.jl file (the one included in the comment below).
        # safe_include("triangulation/weighted.jl")
        # There is also some testing for weighted triangulations at the end of test/predicates/general.jl (search for witness plane).
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
        safe_include("data_structures/polygon_hierarchy.jl", verbose = false)
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
end