# # Constrained Triangulations 
# ## Domain with Interior Holes

# We now consider triangulating a domain which has not only a 
# boundary, but also an interior boundary. To start, let us 
# load the packages. 
using DelaunayTriangulation
using CairoMakie
using StableRNGs 
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

# Let us now define what we will be triangulating. The 
# boundary will be made up of a boundary and three interior 
# holes. To represent the boundaries for this case, 
# we use a vector of vector of vectors. The triple-nested 
# vector structure is necessary to allow for individual 
# boundaries to be segmented (as in the previous tutorial). 
# Moreover, while the outer boundary must be counter-clockwise, 
# the interior boundaries must be clockwise so that the orientations
# of the interiors relative to the boundaries are consistent. Note again 
# that neighbouring segments must connect.
curve_1 = [
    [ # first segment
        (0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (12.0, 0.0), (12.0, 4.0),
        (12.0, 8.0), (14.0, 10.0), (16.0, 12.0), (16.0, 16.0),
        (14.0, 18.0), (12.0, 20.0), (12.0, 24.0), (12.0, 28.0)
    ],
    [ # second segment
        (12.0, 28.0), (8.0, 28.0), (4.0, 28.0), (0.0, 28.0), (-2.0, 26.0), (0.0, 22.0),
        (0.0, 18.0), (0.0, 10.0), (0.0, 8.0), (0.0, 4.0), (-4.0, 4.0),
        (-4.0, 0.0), (0.0, 0.0),
    ]
] # outer: counter-clockwise
curve_2 = [
    [ # first segment
        (4.0, 26.0), (8.0, 26.0), (10.0, 26.0), (10.0, 24.0),
        (10.0, 22.0), (10.0, 20.0)
    ],
    [ # second segment
        (10.0, 20.0), (8.0, 20.0), (6.0, 20.0),
        (4.0, 20.0), (4.0, 22.0), (4.0, 24.0), (4.0, 26.0)
    ]
] # inner: clockwise
curve_3 = [[(4.0, 16.0), (12.0, 16.0), (12.0, 14.0), (4.0, 14.0), (4.0, 16.0)]] # inner: clockwise
curve_4 = [[(4.0, 8.0), (10.0, 8.0), (8.0, 6.0), (6.0, 6.0), (4.0, 8.0)]] # inner: clockwise
curves = [curve_1, curve_2, curve_3, curve_4]
points = [
    (2.0, 26.0), (2.0, 24.0), (6.0, 24.0), (6.0, 22.0), (8.0, 24.0), (8.0, 22.0),
    (2.0, 22.0), (0.0, 26.0), (10.0, 18.0), (8.0, 18.0), (4.0, 18.0), (2.0, 16.0),
    (2.0, 12.0), (6.0, 12.0), (2.0, 8.0), (2.0, 4.0), (4.0, 2.0),
    (-2.0, 2.0), (4.0, 6.0), (10.0, 2.0), (10.0, 6.0), (8.0, 10.0), (4.0, 10.0),
    (10.0, 12.0), (12.0, 12.0), (14.0, 26.0), (16.0, 24.0), (18.0, 28.0),
    (16.0, 20.0), (18.0, 12.0), (16.0, 8.0), (14.0, 4.0), (14.0, -2.0),
    (6.0, -2.0), (2.0, -4.0), (-4.0, -2.0), (-2.0, 8.0), (-2.0, 16.0),
    (-4.0, 22.0), (-4.0, 26.0), (-2.0, 28.0), (6.0, 15.0), (7.0, 15.0),
    (8.0, 15.0), (9.0, 15.0), (10.0, 15.0), (6.2, 7.8),
    (5.6, 7.8), (5.6, 7.6), (5.6, 7.4), (6.2, 7.4), (6.0, 7.6),
    (7.0, 7.8), (7.0, 7.4)]
boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points=points);

# Notice that `curve_1` and `curve_2` are split up into multiple segments. For `curve_3`
# and `curve_4`, note that we have to wrap the entire vector in a vector, essentially 
# treating them as a single segment. 

# Now let us triangulate.
rng = StableRNG(123) # the triangulation is not unique due to cocircular points
tri = triangulate(points; boundary_nodes, rng)
fig, ax, sc = triplot(tri, show_constrained_edges=true, show_convex_hull=true)
fig
@test_reference joinpath(fig_path, "constrained_ex_5.png") fig #src

# Like before, we examine individual segments by referring to them by their ghost vertices, 
# which are still in the order `-1`, `-2`, and so on in the order of the segments provided 
# in `boundary_nodes`. This is a lot more cumbersome to keep track of than the previous
# tutorials since there are many more ghost vertices. This is where the boundary fields 
# become much more useful. For instance, the `boundary_edge_map` in this case is given by:
get_boundary_edge_map(tri)

# The `Tuples` in the `values` are now of the form `((I, J), K)`, with `I` the curve index 
# (with `1` being the outer boundary), `J` the segment index, and `K` the position of the
# edge within the segment. To look at the ghost vertices directly, another useful field is 
# `ghost_vertex_ranges`:
get_ghost_vertex_ranges(tri)

# This field maps a ghost vertex to the complete set of ghost vertices that might be found on the 
# curve corresponding to that ghost vertex. For example, `-3 => -4:-3` means that the ghost vertex 
# `-3` is part of a curve that, in addition to itself, contains the ghost vertex `-4`.
# If you want all the ghost vertex, you can use 
DelaunayTriangulation.all_ghost_vertices(tri)

# which is just `keys(get_ghost_vertex_ranges(tri))`. If you just want to find 
# what curve and what segment a ghost vertex belongs to, you can look at the `ghost_vertex_map`:
get_ghost_vertex_map(tri)

# So that, for example, `-3 => (2, 1)` means that the ghost vertex `-3` corresponds 
# to the first part (from the second `Tuple` element) of the second curve (from the 
# first `Tuple` element). 

# To get all the boundary nodes, you can use
DelaunayTriangulation.get_all_boundary_nodes(tri)

# To give an example of how we can work with this boundary, let us compute the area 
# of the triangulation (a more efficient approach is with [`get_area(tri)`](@ref get_area), but this is just 
# for demonstration). For this, the order of the boundary edges is appropriate, so we must iterate 
# in a way that respects the ordering. 
function get_triangulation_area(tri)
    A = 0.0 
    nc = DelaunayTriangulation.num_curves(tri)
    for curve_index in 1:nc 
        bn = get_boundary_nodes(tri, curve_index)
        ns = DelaunayTriangulation.num_sections(bn)
        for segment_index in 1:ns 
            bnn = get_boundary_nodes(bn, segment_index)
            ne = num_boundary_edges(bnn)
            for i in 1:ne
                vᵢ = get_boundary_nodes(bnn, i)
                vᵢ₊₁ = get_boundary_nodes(bnn, i+1)
                pᵢ, pᵢ₊₁ = get_point(tri, vᵢ, vᵢ₊₁)
                xᵢ, yᵢ = getxy(pᵢ)
                xᵢ₊₁, yᵢ₊₁ = getxy(pᵢ₊₁)
                A += (yᵢ + yᵢ₊₁)*(xᵢ - xᵢ₊₁)
            end
        end
    end
    return A/2
end
A = get_triangulation_area(tri)
@test A ≈ get_area(tri) #src

# This is of course quite a complicated example since we need to take 
# care of the order. If we don't care about order, then the complexity 
# of the code for iterating over a boundary is much simpler. For example, 
# here we compute the perimeter of the boundary, and we also 
# consider the length of each curve and of each segment. 
function get_perimeters(tri)
    total_perimeter = 0.0
    nc = DelaunayTriangulation.num_curves(tri)
    curve_perimeters = zeros(nc) # curve_index => perimeter
    segment_perimeters = Dict{NTuple{2,Int},Float64}() # (curve_index, segment_index) => perimeter
    for (e, ((curve_index, section_index), node_index)) in get_boundary_edge_map(tri)
        u, v = edge_vertices(e)
        p, q = get_point(tri, u, v)
        ℓ = sqrt((getx(p) - getx(q))^2 + (gety(p) - gety(q))^2)
        total_perimeter += ℓ
        curve_perimeters[curve_index] += ℓ
        if haskey(segment_perimeters, (curve_index, section_index))
            segment_perimeters[(curve_index, section_index)] += ℓ
        else
            segment_perimeters[(curve_index, section_index)] = ℓ
        end
    end
    return total_perimeter, curve_perimeters, segment_perimeters
end 
ℓ, cℓ, sℓ = get_perimeters(tri)
@test ℓ ≈ sum(cℓ) #src
@test ℓ ≈ sum(sum.(values(sℓ))) #src