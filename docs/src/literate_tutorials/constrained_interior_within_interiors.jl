# # Constrained Triangulations 
# ## Domain with Interior Holes within Interior Holes 

# Now we consider triangulating a domain which has not only an 
# interior boundary, but also an interior boundary inside that 
# interior boundary. To start, let us load the packages.
using DelaunayTriangulation
using CairoMakie
using StableRNGs
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

# To represent curves inside of curves, note that we have 
# already had to do this for the outer boundary, where 
# interior curves are clockwise - the opposite orientation 
# of the outer boundary. So, similarly, interiors within 
# other interior curves must be counter-clockwise. Basically,
# the orientation of an interior curve is the opposite 
# orientation to the curve it's inside of. Here are the 
# points we will be triangulating. 
curve_1 = [
    [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0), (15.0, 0.0), (20.0, 0.0), (25.0, 0.0)],
    [(25.0, 0.0), (25.0, 5.0), (25.0, 10.0), (25.0, 15.0), (25.0, 20.0), (25.0, 25.0)],
    [(25.0, 25.0), (20.0, 25.0), (15.0, 25.0), (10.0, 25.0), (5.0, 25.0), (0.0, 25.0)],
    [(0.0, 25.0), (0.0, 20.0), (0.0, 15.0), (0.0, 10.0), (0.0, 5.0), (0.0, 0.0)],
] # outer-most boundary: counter-clockwise  
curve_2 = [
    [(4.0, 6.0), (4.0, 14.0), (4.0, 20.0), (18.0, 20.0), (20.0, 20.0)],
    [(20.0, 20.0), (20.0, 16.0), (20.0, 12.0), (20.0, 8.0), (20.0, 4.0)],
    [(20.0, 4.0), (16.0, 4.0), (12.0, 4.0), (8.0, 4.0), (4.0, 4.0), (4.0, 6.0)],
] # inner boundary: clockwise 
curve_3 = [
    [
        (12.906, 10.912), (16.0, 12.0), (16.16, 14.46), (16.29, 17.06),
        (13.13, 16.86), (8.92, 16.4), (8.8, 10.9), (12.906, 10.912),
    ],
] # this is inside curve_2, so it's counter-clockwise 
curves = [curve_1, curve_2, curve_3]
points = [
    (3.0, 23.0), (9.0, 24.0), (9.2, 22.0), (14.8, 22.8), (16.0, 22.0),
    (23.0, 23.0), (22.6, 19.0), (23.8, 17.8), (22.0, 14.0), (22.0, 11.0),
    (24.0, 6.0), (23.0, 2.0), (19.0, 1.0), (16.0, 3.0), (10.0, 1.0), (11.0, 3.0),
    (6.0, 2.0), (6.2, 3.0), (2.0, 3.0), (2.6, 6.2), (2.0, 8.0), (2.0, 11.0),
    (5.0, 12.0), (2.0, 17.0), (3.0, 19.0), (6.0, 18.0), (6.5, 14.5),
    (13.0, 19.0), (13.0, 12.0), (16.0, 8.0), (9.8, 8.0), (7.5, 6.0),
    (12.0, 13.0), (19.0, 15.0),
]
boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points = points)


# To now triangulate:
rng = StableRNG(123) # triangulation is not unique when there are cocircular points
tri = triangulate(points; boundary_nodes, rng)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "constrained_ex_6.png") fig #src
