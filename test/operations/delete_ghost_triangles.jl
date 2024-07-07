using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StructEquality
using ..DelaunayTriangulation: Triangulation
using StaticArrays
using Preferences

@testset "Deleting ghost triangles" begin
    tri, label_map, index_map = simple_geometry()
    _tri = deepcopy(tri)
    DT.add_ghost_triangles!(tri)
    DT.delete_ghost_triangles!(tri)

    @test tri == _tri
end

@testset "Making sure the ghost triangles are correctly adjusted on disjoint sets" begin
    θ = LinRange(0, 2π, 100) |> collect
    θ[end] = 0
    xy = Vector{Vector{Vector{NTuple{2,Float64}}}}()
    cx = 0.0
    for i in 1:2
        push!(xy, [[(cx + cos(θ) + 1e-3rand(), sin(θ) + 1e-3rand()) for θ in θ]]) # use perturbations to also do this check without ExactPredicates loaded. iszero makes the boundary closed
        push!(xy, [[(cx + 0.5cos(θ) + 1e-3rand(), 0.5sin(θ) + 1e-3rand()) for θ in reverse(θ)]])
        cx += 3.0
    end
    xy[1][1][end] = xy[1][1][1]
    xy[2][1][end] = xy[2][1][1] 
    xy[3][1][end] = xy[3][1][1]
    xy[4][1][end] = xy[4][1][1] 
    boundary_nodes, points = convert_boundary_points_to_indices(xy)
    tri = triangulate(points; boundary_nodes=boundary_nodes)
    @test validate_triangulation(tri)
    validate_statistics(tri)
end