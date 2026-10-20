module HelperFunctions

using StaticArrays
using StableRNGs
using LinearAlgebra
using ElasticArrays
using Random
using DataStructures
using DelimitedFiles
using OrderedCollections
using Distributions
using InteractiveUtils
using Test
using DelaunayTriangulation
using DelaunayTriangulation:
    validate_triangulation,
    test_adjacent2vertex_map_matches_adjacent_map,
    test_state,
    test_adjacent_map_matches_adjacent2vertex_map,
    test_each_edge_has_two_incident_triangles,
    test_triangle_orientation,
    test_iterators
import SpatialIndexing as SI
import SimpleGraphs: SimpleGraphs

getxy((0.0, 0.0)) # avoid shadow

const DT = DelaunayTriangulation

# These helpers are split over the files in helpers/ by purpose. To use all of them, e.g. when testing in the REPL, run
#     include(joinpath(pkgdir(DelaunayTriangulation), "test", "helper_functions.jl"))
#     using .HelperFunctions
# with the test environment active.
include("helpers/testing.jl")
include("helpers/comparisons.jl")
include("helpers/examples.jl")
include("helpers/weighted_examples.jl")
include("helpers/validation.jl")
include("helpers/reference_implementations.jl")
include("helpers/curves.jl")
include("helpers/rtree.jl")
include("helpers/sampling.jl")

# Re-exported from DelaunayTriangulation and InteractiveUtils
export validate_triangulation, test_adjacent2vertex_map_matches_adjacent_map, test_state
export test_adjacent_map_matches_adjacent2vertex_map, test_each_edge_has_two_incident_triangles
export test_triangle_orientation, test_iterators, subtypes

end
