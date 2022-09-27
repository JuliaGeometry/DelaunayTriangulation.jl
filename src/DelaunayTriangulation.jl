module DelaunayTriangulation

using SimpleGraphs 
using ExactPredicates
using ExactPredicates.Codegen 
using Random

const TriangleType = NTuple{3, Int64}
const EdgeType = NTuple{2, Int64}
const LargeRightIdx = 0 # p₋₁
const LargeLeftIdx = -1 # p₋₂
const BoundaryIdx = -2  # ∅
const VertexNeighbourVector = Vector{Int64}
const Adjacent2VertexVector = Vector{EdgeType}

include("data_structures.jl")
end
