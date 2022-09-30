using SimpleGraphs
using ExactPredicates
using ExactPredicates.Codegen
using Test
using Random
const TriangleType = NTuple{3,Int64}
const LargeRightIdx = 0 # p₋₁
const LargeLeftIdx = -1 # p₋₂
const EmptyIdx = -2 # ∅
const VertexNeighbourVector = Vector{Int64}
const AdjacentToVertexVector = Vector{NTuple{2,Int64}}


p0 = Float64[5, 5]
p1 = Float64[4.5, 2.5]
p2 = Float64[2.5, 1.5]
p3 = Float64[3, 3.5]
p4 = Float64[0, 2]
p5 = Float64[1, 5]
p6 = Float64[1, 3]
p7 = Float64[4, -1]
p8 = Float64[-1, 4]
pts = [p0, p1, p2, p3, p4, p5, p6, p7, p8]

# @views shuffle!(pts[begin+1:end])


## Add the next point 



## Now let's do a clearer run through


## Do a smaller example 
p1 = Float64[20, 20]
p2 = Float64[0, 6]
p3 = Float64[12, -2]
p4 = Float64[10, 10]
PTS = [p1, p2, p3, p4]
𝒯, 𝒟, 𝒜, 𝒱𝒩 = triangulate(PTS; shuffle_pts=false, trim=false)
@test 𝒯 == TriangleType[
    (1, -1, 2),
    (-1, 0, 2),
    (0, 1, 3),
    (2, 0, 3),
    (1, 2, 4),
    (2, 3, 4),
    (3, 1, 4)
]
@test 𝒜(1, -1) == 2
@test 𝒜(-1, 2) == 1
@test 𝒜(2, 1) == -1
@test 𝒜(-1, 0) == 2
@test 𝒜(0, 2) == -1
@test 𝒜(2, -1) == 0
@test 𝒜(0, 1) == 3
@test 𝒜(1, 3) == 0
@test 𝒜(3, 0) == 1
@test 𝒜(2, 0) == 3
@test 𝒜(0, 3) == 2
@test 𝒜(3, 2) == 0
@test 𝒜(1, 2) == 4
@test 𝒜(2, 4) == 1
@test 𝒜(4, 1) == 2
@test 𝒜(2, 3) == 4
@test 𝒜(3, 4) == 2
@test 𝒜(4, 2) == 3
@test 𝒜(3, 1) == 4
@test 𝒜(1, 4) == 3
@test 𝒜(4, 3) == 1
@test sort(𝒱𝒩(-1)) == [0, 1, 2]
@test sort(𝒱𝒩(0)) == [-1, 1, 2, 3]
@test sort(𝒱𝒩(1)) == [-1, 0, 2, 3, 4]
@test sort(𝒱𝒩(2)) == [-1, 0, 1, 3, 4]
@test sort(𝒱𝒩(3)) == [0, 1, 2, 4]
@test sort(𝒱𝒩(4)) == [1, 2, 3]

w₁ = 𝒜(1, -1)
delete_triangle!(𝒯, TriangleType((1, -1, w₁)))
delete_edge_from_adjacency!(𝒜, 1, -1; protect_boundary=false)
delete_edge_from_adjacency!(𝒜, -1, w; protect_boundary=false)
𝒜[(w₁, 1)] = EmptyIdx

w₂ = 𝒜(-1, 0)
delete_triangle!(𝒯, TriangleType((-1, 0, w₂)))
delete_edge_from_adjacency!(𝒜, -1, 0; protect_boundary=false)
delete_edge_from_adjacency!(𝒜, 0, w₂; protect_boundary=false)

w₃ = 𝒜(0, 1)
delete_triangle!(𝒯, TriangleType((0, 1, w₃)))
delete_edge_from_adjacency!(𝒜, 0, 1; protect_boundary=false)
delete_edge_from_adjacency!(𝒜, 0, w₃; protect_boundary=false)
𝒜[(w₃, w₂)] = EmptyIdx
𝒜[(1, w₃)] = EmptyIdx
𝒯
@test 𝒯 == [(2,)]