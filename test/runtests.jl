using SimpleGraphs
using ExactPredicates
using ExactPredicates.Codegen
using Test
using Random
using DelaunayTriangulation
using CairoMakie
const DT = DelaunayTriangulation

#@testset "DelaunayTriangulation.jl" begin
    @testset "Adjacency data structure" begin
        include("adjacent.jl")
    end
    @testset "Delaunay graph representation" begin
        include("delaunay_graph.jl")
    end
    @testset "History DAG" begin
        include("history_dag.jl")
    end
    @testset "Geometrical and Delaunay predicates" begin
        include("predicates.jl")
    end
    @testset "Representation for a collection of triangles" begin
        include("triangulation_struct.jl")
    end
    @testset "Initialisation for triangulating" begin
        include("initialisation.jl")
    end
    @testset "Flipping an edge in a triangulation" begin
        include("edge_flip.jl")
    end
    @testset "Some examples of Delaunay triangulations" begin
        include("examples.jl")
    end
#end
