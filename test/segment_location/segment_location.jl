using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie 
include("../helper_functions.jl")

tri = shewchuk_example_constrained()

fig, ax, sc = triplot(tri)
let vert = each_solid_vertex(tri)
    text!(ax, collect(get_point(tri, vert...)); text=string.(vert))
end
fig