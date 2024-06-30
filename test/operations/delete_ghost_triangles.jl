using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using StructEquality
using ..DelaunayTriangulation: Triangulation

# Base.include(@__MODULE__, "../helper_functions.jl")

@testset "Deleting ghost triangles" begin
    tri, label_map, index_map = simple_geometry()
    _tri = deepcopy(tri)
    DT.add_ghost_triangles!(tri)
    DT.delete_ghost_triangles!(tri)

    @test tri == _tri
end

@testset "Making sure the ghost triangles are correctly adjusted on disjoint sets" begin
    θ = LinRange(0, 2π, 500) |> collect
    θ[end] = 0
    xy = Vector{Vector{Vector{NTuple{2,Float64}}}}()
    cx = 0.0
    for i in 1:2
        push!(xy, [[(cx + cos(θ), sin(θ)) for θ in θ]])
        push!(xy, [[(cx + 0.5cos(θ), 0.5sin(θ)) for θ in reverse(θ)]])
        cx += 3.0
    end
    boundary_nodes, points = convert_boundary_points_to_indices(xy)
    tri = triangulate(points; boundary_nodes=boundary_nodes)
    @test validate_triangulation(tri)
    validate_statistics(tri)
end