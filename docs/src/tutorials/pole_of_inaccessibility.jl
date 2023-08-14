# # Pole of Inaccessibility 
# 

# In this tutorial, we demonstrate how we can compute the pole 
# of inaccessibility of a given polygon, namely the point inside the 
# polygon that is furthest from the boundary, using the algorithm described 
# in [this blogpost](https://blog.mapbox.com/a-new-algorithm-for-finding-a-visual-center-of-a-polygon-7c77e6492fbc) based on recursively subdividing the polygon using quadtree partitioning.
# The polygons should be defined the same way the boundary of a triangulation is defined; see the constrained 
# triangulation tutorials. In particular, the algorithm we use works for polygons that may possibly have 
# interior holes, or even nested holes, as well as for multipolygons. 
# We give a simple example. First, let us define our polygon.
using DelaunayTriangulation
using CairoMakie
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

points = [
    0.0 8.0
    2.0 5.0
    3.0 7.0
    1.81907 8.13422
    3.22963 8.865
    4.24931 7.74335
    4.50423 5.87393
    3.67149 4.3784
    2.73678 2.62795
    5.50691 1.38734
    8.43 2.74691
    9.7046 5.53404
    8.56595 7.79433
    6.71353 9.03494
    4.13034 9.66375
    2.75378 10.3775
    1.0883 10.4965
    -1.138 9.83369
    -2.25965 8.45712
    -2.78649 5.94191
    -1.39292 3.64763
    0.323538 4.97322
    -0.900078 6.6217
    0.98633 9.68074
    0.153591 9.54478
    0.272554 8.66106
    2.90673 8.18521
    2.12497 9.42582
    7.27436 2.7979
    3.0 4.0
    5.33697 1.88019]'
outer_boundary = [ # split into segments for demonstration purposes
    [1, 4, 3, 2],
    [2, 9, 10, 11, 8, 7, 12],
    [12, 6, 13, 5, 14, 15, 16, 17, 16],
    [16, 17, 18, 19, 20, 21, 22, 23, 1]
]
inner_1 = [
    [26, 25, 24], [24, 28, 27, 26]
]
inner_2 = [
    [29, 30, 31, 29]
]
boundary_nodes = [outer_boundary, inner_1, inner_2]
fig = Figure()
ax = Axis(fig[1, 1])
for i in eachindex(boundary_nodes)
    lines!(ax, points[:, reduce(vcat, boundary_nodes[i])], color=:red)
end
fig

# To now find the pole of inaccessibility, use `pole_of_inaccessibility`:
pole = DelaunayTriangulation.pole_of_inaccessibility(points, boundary_nodes)
@test pole[1] ≈ -1.078522500000003 #src 
@test pole[2] ≈ 5.372597499999995 #src

#-
scatter!(ax, pole, color=:blue, markersize=16)
fig
@test_reference joinpath(fig_path, "pole_of_inaccessibility_ex_1.png") fig #src

# This shows that the point inside the red region that is furthest from the boundary 
# is the blue point shown.

# We note that triangulations also store the poles of inaccessibility for each boundary, 
# as these are used to define the orientation of a ghost edge. Here is an example. First, 
# we get the triangulation. 
θ = LinRange(0, 2π, 20) |> collect
θ[end] = 0 # need to make sure that 2π gives the exact same coordinates as 0
xy = Vector{Vector{Vector{NTuple{2,Float64}}}}()
cx = 0.0
for i in 1:2
    global cx
    ## Make the exterior circle
    push!(xy, [[(cx + cos(θ), sin(θ)) for θ in θ]])
    ## Now the interior circle - clockwise
    push!(xy, [[(cx + 0.5cos(θ), 0.5sin(θ)) for θ in reverse(θ)]])
    cx += 3.0
end
boundary_nodes, points = convert_boundary_points_to_indices(xy)
tri = triangulate(points; boundary_nodes=boundary_nodes, check_arguments=false)

# To see the poles, called representative points, we use 
DelaunayTriangulation.get_representative_point_list(tri)

# The keys of the `Dict` refer to the curve index, and the values contain 
# the coordinates. 
fig, ax, sc = triplot(tri, show_ghost_edges=true)
colors = (:red, :blue, :darkgreen, :purple)
for i in eachindex(boundary_nodes)
    lines!(ax, points[reduce(vcat, boundary_nodes[i])], color=colors[i], linewidth=6)
    coords = DelaunayTriangulation.get_representative_point_coordinates(tri, i)
    scatter!(ax, coords, color=colors[i], markersize=16)
end
fig
@test_reference joinpath(fig_path, "pole_of_inaccessibility_ex_2.png") fig #src

# Note that the green and purple boundaries have the same pole of inaccessibility. The 
# first curve, the red curve, is the only one that has the pole of inaccessibility computed 
# with respect to all other boundaries. You can also see that indeed the ghost edges are all 
# oriented relative to the representative points.