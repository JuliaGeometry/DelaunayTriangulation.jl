# # Constrained Triangulations 
# ## Outer Boundary 

# This tutorial now considers the case where, rather than only 
# having constrained segments, we have a constrained outer boundary. This 
# is especially useful as it allows us to, for example, have 
# a non-convex boundary. To start, let us load in the packages we will need.
using DelaunayTriangulation
using CairoMakie
import LinearAlgebra: norm #src
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

# Now, we define some of the points we will be triangulating. 
pts = [
    (-7.36, 12.55), (-9.32, 8.59), (-9.0, 3.0), (-6.32, -0.27),
    (-4.78, -1.53), (2.78, -1.41), (-5.42, 1.45), (7.86, 0.67),
    (10.92, 0.23), (9.9, 7.39), (8.14, 4.77), (13.4, 8.61),
    (7.4, 12.27), (2.2, 13.85), (-3.48, 10.21), (-4.56, 7.35),
    (3.44, 8.99), (3.74, 5.87), (-2.0, 8.0), (-2.52, 4.81),
    (1.34, 6.77), (1.24, 4.15),
]

# To define a boundary, we need to provide a counter-clockwise sequence 
# of indices corresponding to the boundary points, and the first index 
# must match the last index so the boundary is closed. While we could 
# include in `pts` the boundary points that we want to include, and then 
# write down the indices of the points within `pts`, this is cumbersome and often 
# tedious to get correct. So, we instead provide the function 
# [`convert_boundary_points_to_indices`](@ref) which takes in a vector of coordinates, 
# and then returns the correct set of indices. Here is how we use it:
boundary_points = [
    (0.0, 0.0), (2.0, 1.0), (3.98, 2.85), (6.0, 5.0),
    (7.0, 7.0), (7.0, 9.0), (6.0, 11.0), (4.0, 12.0),
    (2.0, 12.0), (1.0, 11.0), (0.0, 9.13), (-1.0, 11.0),
    (-2.0, 12.0), (-4.0, 12.0), (-6.0, 11.0), (-7.0, 9.0),
    (-6.94, 7.13), (-6.0, 5.0), (-4.0, 3.0), (-2.0, 1.0), (0.0, 0.0),
]
boundary_nodes, pts = convert_boundary_points_to_indices(boundary_points; existing_points = pts);

# The keyword argument `existing_points` is so that the points in `boundary_points` get appended (in-place) to 
# `pts`, as we see:
pts

# The `boundary_nodes` is then these indices:
boundary_nodes

# To now triangulate, we use the `boundary_nodes` keyword argument. Like in the last tutorial, 
# we also give a comparison to the unconstrained version.
tri = triangulate(pts)
cons_tri = triangulate(pts; boundary_nodes)

#- 
fig = Figure()
ax1 = Axis(
    fig[1, 1], xlabel = "x", ylabel = L"y",
    title = "(a): Unconstrained", titlealign = :left,
    width = 300, height = 300,
)
ax2 = Axis(
    fig[1, 2], xlabel = "x", ylabel = L"y",
    title = "(b): Constrained", titlealign = :left,
    width = 300, height = 300,
)
triplot!(ax1, tri)
triplot!(ax2, cons_tri, show_constrained_edges = true, show_convex_hull = true)
resize_to_layout!(fig)
fig
@test_reference joinpath(fig_path, "constrained_ex_2.png") fig #src

# Notice now that the boundary in (b) is not convex, as is clear 
# from the convex hull shown in red. You can access the convex hull 
# using [`get_convex_hull(cons_tri)`](@ref get_convex_hull). We also note that the triangulation 
# no longer contains every point in `pts`, as by default all triangles away 
# from the boundary are deleted, so that we do actually have a boundary. If for 
# some reason you do not want this behaviour, use `delete_holes = false`:
full_tri = triangulate(pts; boundary_nodes, delete_holes = false)
fig, ax, sc = triplot(full_tri, show_constrained_edges = true, show_convex_hull = true)
@test_reference joinpath(fig_path, "constrained_ex_3.png") fig #src

# This default behaviour does mean you need to be careful if you use [`DelaunayTriangulation.each_point`](@ref) 
# or [`DelaunayTriangulation.each_point_index`](@ref), as these iterators will contain all points, possibly iterating 
# over points that aren't in the triangulation. For this reason, it is recommended that you 
# use [`each_solid_vertex`](@ref) as a default.

# There are multiple methods available for working directly 
# with the boundary nodes. You can get the boundary nodes 
# using `get_boundary_nodes(tri)`:
get_boundary_nodes(cons_tri)

# Later tutorials also consider other methods for working with the boundary 
# where care needs to be taken with the boundary, or part of the boundary, 
# being considered. For now, here is an example where we use 
# `get_right_boundary_node` to iterate over the boundary in a counter-clockwise 
# order, getting the area of the triangulation using 
# ```math
# A = \dfrac{1}{2}\sum_{i=1}^n \left(y_i + y_{i+1}\right)\left(x_i - x_{i+1}\right).
# ```
# Here is one implementation.
function shoelace_area(tri)
    bn = get_boundary_nodes(tri)
    n = num_boundary_edges(bn) # length(bn) - 1 in this case since bn[1] = bn[end]
    A = 0.0
    for i in 1:n
        vᵢ = get_boundary_nodes(bn, i)
        vᵢ₊₁ = get_boundary_nodes(bn, i + 1)
        pᵢ, pᵢ₊₁ = get_point(tri, vᵢ, vᵢ₊₁)
        xᵢ, yᵢ = getxy(pᵢ)
        xᵢ₊₁, yᵢ₊₁ = getxy(pᵢ₊₁)
        A += (yᵢ + yᵢ₊₁) * (xᵢ - xᵢ₊₁)
    end
    return A / 2
end
shoelace_area(cons_tri)
@test shoelace_area(cons_tri) ≈ get_area(cons_tri) #src

# We also provide a map that contains the edges as the keys (*not* in order), and the values are `Tuple`s `(I, J)` 
# such that `get_boundary_nodes(get_boundary_nodes(cons_tri, I), J)` gives the corresponding 
# edge. The first call, `bn = get_boundary_nodes(cons_tri, I)` is for obtaining the chain of boundary edges 
# containing the boundary edge, and then `get_boundary_nodes(bn, j)` gets the actual edge.
get_boundary_edge_map(cons_tri)

# In our case, the `I` is just `boundary_nodes` since we only have one contiguous boundary. 
# To give an example, take 
bem = get_boundary_edge_map(cons_tri)
e, (I, J) = first(bem)

#-
bn = get_boundary_nodes(cons_tri, I) # same as boundary_nodes for this problem; see the later tutorials 
bn_j = get_boundary_nodes(bn, J)

# This returns `23`, which is the start of the edge `e`. The full edge is given by 
get_boundary_nodes.(Ref(bn), (J, J + 1)) # Ref to not broadcast over bn
@test e == get_boundary_nodes.(Ref(bn), (J, J + 1)) #src

# To give an example, here's how we compute the perimeter of the triangulation. This only 
# needs the edges, so we only consider the `keys` of the map.
function get_perimeter(tri)
    bem = get_boundary_edge_map(tri)
    ℓ = 0.0
    for e in keys(bem)
        u, v = edge_vertices(e)
        p, q = get_point(tri, u, v)
        ℓ += sqrt((getx(p) - getx(q))^2 + (gety(p) - gety(q))^2)
    end
    return ℓ
end
get_perimeter(cons_tri)
@test get_perimeter(cons_tri) ≈ sum(norm(get_point(cons_tri, boundary_nodes[i + 1]) .- get_point(cons_tri, boundary_nodes[i])) for i in 1:(length(boundary_nodes) - 1)) #src
