############################################
##
## UTILITY FUNCTIONS 
##
############################################
"""
    find_root(G::HistoryGraph; method=:brute)

Finds the root of the graph `G`, using one of two methods:
- `method=:brute`: Uses brute force to look for the root, defining the root to be the point with zero out-degree.
- `method=:rng`: Uses randomised searching to look for the root.

See also [`_find_root_brute`](@ref) and [`_find_root_rng`](@ref).

Note that these methods only work when the graph is acyclic.
"""
function find_root(G::HistoryGraph; method=:brute)
    if method == :brute
        return _find_root_brute(G)
    elseif method == :rng
        return _find_root_rng(G)
    end
end
function _find_root_brute(G::HistoryGraph)
    for (k, v) in graph(G).NN
        length(v) == 0 && return k
    end
end
function _find_root_rng(G::HistoryGraph)
    verts = vlist(graph(G))
    num_verts = length(verts)
    starting_node = verts[rand(1:num_verts)]
    for _ in 1:num_verts
        num_choices = in_deg(G, starting_node)
        num_choices == 0 && return starting_node
        starting_node = in_neighbors(G, starting_node)[rand(1:num_choices)]
    end
end

## Accessing data structures directly from a triangulation 
adjacent(DT::Triangulation, uv) = get_edge(adjacent(DT), uv)
adjacent(DT::Triangulation, u, v) = get_edge(adjacent(DT), u, v)
adjacent2vertex(DT::Triangulation, u) = get_edge(adjacent2vertex(DT), u)
neighbours(DT::Triangulation, u) = get_neighbour(graph(DT), u)
get_point(DT::Triangulation, i) = get_point(points(DT), i)
edges(DT::Triangulation) = edges(adjacent(DT))
num_triangles(DT::Triangulation) = length(triangles(DT))
num_points(DT::Triangulation) = length(points(DT))
num_edges(DT::Triangulation) = num_triangles(DT) + num_points(DT) - 1 # Euler's formula; could also use length(edges(DT)) ÷ 2 

"""
    is_point_higher(p, q)
    is_point_lower(p, q)

Tests if `p` is lexicographically higher than `q`, where we say that `p = (xp, yp)`
is lexicographically higher than `q = (xq, yp)` if `yp > yq` or 
`yp = yq` and `xq > xp`.
"""
is_point_higher(p::AbstractPoint, q::AbstractPoint) = (gety(p) > gety(q)) || (gety(p) == gety(q) && getx(q) > getx(p))
is_point_lower(p::AbstractPoint, q::AbstractPoint) = is_point_higher(q, p)

"""
    partial_highest_point_sort!(v, k)

Partially sorts the first `v` so that the first `k` entries are the highest points, with the first 
being the highest (acccording to `is_point_higher`).
"""
partial_highest_point_sort!(v, k) = partialsort!(v, k, lt=is_point_higher)

"""
    num_less(val, v)

Counts the number of values in `v` strictly less than `val`.
"""
num_less(val, v) = count(<(val), v)

"""
    is_delaunay(DT::Triangulation)

Tests if the given triangulation `DT` is Delaunay. This is done by identifying 
if all edges are legal using `islegal`, noting that a legal triangulation is 
necessarily Delaunay.
"""
function is_delaunay(DT::Triangulation)
    adj = adjacent(DT)
    tri_edges = edges(adj)
    for (i, j) in tri_edges
        k = get_edge(adj, i, j) # The convex hull's edges are always included. 
        ℓ = get_edge(adj, j, i)
        k ≠ BoundaryIndex && ℓ ≠ BoundaryIndex && !islegal(i, j, adj, points(DT)) && return false
    end
    return true
end