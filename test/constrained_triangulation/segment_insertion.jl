using ..DelaunayTriangulation
const DT = DelaunayTriangulation
using CairoMakie
include("../helper_functions.jl")

tri = fixed_shewchuk_example_constrained()

fig, ax, sc = triplot(tri)
let vert = each_solid_vertex(tri)
    text!(ax, collect(get_point(tri, vert...)); text=string.(vert))
end
lines!(ax, [get_point(tri, 2, 7)...], color=:blue, linestyle=:dash)
fig

e = (2, 7)
T, C, L, R = DT.locate_intersecting_triangles(tri, e)
DT.delete_intersected_triangles!(tri, T)

P = L 
P̂ = R 

V = L 
m = length(V)
next = zeros(Int64, m-2)
prev = zeros(Int64, m-2)
shuffled_indices = zeros(Int64, m-2)
for i = 2:m-1 
    next[i] = i + 1 
    prev[i] = i - 1
    shuffled_indices[i] = i 
end