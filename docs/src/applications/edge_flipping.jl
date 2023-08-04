# # Example Applications 
# ## Edge Flipping 

# In this application, we implement the edge flipping algorithm for computing an 
# (unconstrained) Delaunay triangulation. This algorithm, 
# as described by de Berg et al. in the book *Computational Geometry: Algorithms and Applications* (3 ed., 2008),
# is as follows: For a given point set $P$ and some initial triangulation $\mathcal T$ 
# of the point set:
#
# * while $\mathcal T$ contains an illegal edge $e_{ij}$
#   * do
#       * find the triangle $T_{ijk}$ that contains $e_{ij}$* Let $T_{ijk}$ and $T_{ji\ell}$ be the two triangles adjacent to $e_{ij}$.
#       * Remove $e_{ij}$ from $\mathcal T$, and add $e_{k\ell}$ instead. That is, flip the edge $e_{ij}$.
# * return $\mathcal T$
#
# Here, $e_{ij}$ denotes the edge $\overline{p_ip_j}$, where $p_i$ is the $i$th point, and $T_{ijk}$ is the triangle $\triangle p_ip_jp_k$. 
# The edge $e_{ij}$ is *illegal* if the circle through
# $p_i$, $p_j$, and $p_k$ contains $p_\ell$ in its interior, where 
# $p_k$ and $p_\ell$ are the other points on the triangles adjoining $e_{ij}$.
# The idea behind the algorithm is that flipping illegal edges _increases_ the minimum angle 
# of the triangulation, and thus, since the Delaunay triangulation is the triangulation with the
# largest minimum angle, the algorithm converges to the Delaunay triangulation and does 
# so in a finite number of steps. For our implementation of this algorithm,
# we won't be dealing with issues with degeneracies and collinear points just to simplify 
# the exposition.
# 
# One complication with the algorithm is that we need to find a way to 
# construct an initial triangulation $\mathcal T$ for the given set of points. 
# One way to do this is to construct the convex hull, pick one point on the convex hull 
# and connect all other points on the convex hull to that point, and then for all remaining points 
# in the set, find the triangle that contains the point and add the point to the triangulation, connecting 
# that point to the triangle's vertices. Here is one implementation of this. 
using DelaunayTriangulation
using CairoMakie
using StableRNGs
using ReferenceTests #src
using Test #src
fig_path = joinpath(@__DIR__, "../figures") #src

function get_initial_triangulation(points)
    ## Get the convex hull
    ch = convex_hull(points)
    idx = get_indices(ch)

    ## Get an empty triangulation to start
    tri = Triangulation(points)

    ## Now connect each point to the first vertex. Note that the convex hull 
    ## is given in a counter-clockwise order.
    for i in 2:(length(idx)-2)
        add_triangle!(tri, idx[begin], idx[i], idx[i+1])
    end

    ## Now, for all the remaining points, split the triangle that they are inside of 
    interior_vertices = setdiff(each_point_index(points), idx)
    for r in interior_vertices
        T = jump_and_march(tri, get_point(tri, r))
        i, j, k = indices(T)
        split_triangle!(tri, i, j, k, r)
    end
    return tri
end

# As an example:
rng = StableRNG(123)
points = rand(rng, 2, 50)
tri = get_initial_triangulation(points)
fig, ax, sc = triplot(tri)
fig
@test_reference joinpath(fig_path, "edge_flipping_ex_1.png") fig #src

# Now we need to write the function that does the flipping. To keep 
# track of what edges are illegal, we initialise by assuming that all 
# edges are illegal, and we repeatedly check the first edge in the set. 
# When an edge is illegal, we flip it. When an edge $e_{ij}$, with adjacent 
# triangles $T_{ijk}$ and $T_{ji\ell}$, is flipped, the new triangles are 
# $T_{ik\ell}$ and $T_{j\ell k}$, and the edges that need to be checked
# are $e_{ik}$ and $e_{kj}$. With this in mind, here is the implementation.
function edge_flipping_step!(tri, illegal_edges)
    e = first(illegal_edges)
    i, j = edge_indices(e)
    cert = DelaunayTriangulation.is_legal(tri, i, j)
    if DelaunayTriangulation.is_illegal(cert)
        flip_edge!(tri, i, j)
        k = get_adjacent(tri, i, j)
        delete!(illegal_edges, e)
        push!(illegal_edges, (i, k), (k, j))
    else
        delete!(illegal_edges, e)
    end
    return tri
end
function edge_flipping_algorithm!(tri) # tri: an initial triangulation
    illegal_edges = each_edge(tri)
    while !isempty(illegal_edges)
        edge_flipping_step!(tri, illegal_edges)
    end
    return tri
end
function edge_flipping_algorithm(points)
    tri = get_initial_triangulation(points)
    edge_flipping_algorithm!(tri)
    return tri
end

# Here is an example, comparing the triangulation to the triangulation 
# obtained with `triangulate`.
rng = StableRNG(17423)
points = randn(rng, 2, 500)
tri_flip = edge_flipping_algorithm(points)
tri_delaunay = triangulate(points)
fig, ax, sc = triplot(tri_flip; axis=(title="Flipping algorithm",))
triplot!(Axis(fig[1, 2], title="Delaunay algorithm"), tri_delaunay)
fig
@test_reference joinpath(fig_path, "edge_flipping_ex_2.png") fig #src

#-
DelaunayTriangulation.compare_triangle_collections(
    get_triangles(tri_flip),
    get_triangles(tri_delaunay)
)
@test DelaunayTriangulation.compare_triangle_collections(tri_flip.triangles, tri_delaunay.triangles) #src

# Here is an animation of the flipping algorithm in action. 
rng = StableRNG(17423)
points = randn(rng, 2, 100)
tri = get_initial_triangulation(points)
illegal_edges = each_edge(tri)
fig = Figure()
ax = Axis(fig[1, 1])
_tri = Observable(tri)
triplot!(ax, _tri)
record(fig, "flipping_anim.mp4", framerate=24) do io
    ctr = 1
    while !isempty(illegal_edges) 
        _tri[] = edge_flipping_step!(_tri[], illegal_edges)
        ctr % 6 == 0 && recordframe!(io) # don't need to show every single edge flip
        ctr += 1
    end
    for _ in 1:48 # freeze the end for two seconds #src
        recordframe!(io) #src
    end #src
end
# ![Flipping animation](flipping_anim.mp4)

# The next application shows an algorithm that makes this flipping procedure much more efficient.