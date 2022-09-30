"""
    locate_triangle(𝒟::HistoryDAG, pts, p, init=find_root(𝒟; method=:rng))

Given the point location data structure `D` and a set of `pts`, finds the triangle in 
the current triangulation such that `p` is in its interior. The point location starts at `init`.
The function is recursive, and returns a tuple `(tri, flag)`:
    - `tri`: This is the triangle that `p` is in.
    - `flag`: If `flag == 0`, then `p` is on an edge of `tri`. Otherwise, it is in the open interior.
"""
function locate_triangle(𝒟::HistoryDAG, pts, p, init=find_root(𝒟; method=:rng))
    if out_deg(𝒟, init) == 0
        return init, intriangle(init, pts, p)
    end
    out = out_neighbors(𝒟, init)
    for T in out
        intriangle(T, pts, p) ≥ 0 && return locate_triangle(𝒟, pts, p, T)
    end
end

"""
    add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, Tᵢⱼₖ, r)

Given a triangulation `𝒯`, adds the `r`th point of the point set into the triangulation.

# Arguments 
- `𝒯`: The current triangulation.
- `𝒟`: The point location data structure.
- `𝒜`: The adjacency list.
- `𝒜⁻¹`: The adjacent-to-vertex list.
- `𝒟𝒢`: The vertex-neighbour data structure.
-` Tᵢⱼₖ`: The triangle that the `r`th point is inside of. Must be positively oriented.
- `r`: The index of the point in the original point set that is being introduced.

# Outputs 
`𝒯`, `𝒟`, and `𝒜` are all updated in-place.
"""
function add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, Tᵢⱼₖ, r)
    i, j, k = Tᵢⱼₖ # The triangle to be split into three
    delete_triangle!(𝒯, Tᵢⱼₖ) # Now that we've split the triangle, we can remove the triangle
    T₁, T₂, T₃ = TriangleType((i, j, r)), TriangleType((j, k, r)), TriangleType((k, i, r)) # New triangles to add. Note that these triangles are all positively oriented.
    add_triangle!(𝒯, T₁, T₂, T₃) # The three new triangles
    add_triangle!(𝒟, T₁, T₂, T₃) # Add the new triangles into DAG
    add_edge!(𝒟, Tᵢⱼₖ, T₁, T₂, T₃) # Add edges from the old triangle to the new triangles
    update_adjacent!(𝒜, T₁, T₂, T₃) # Add the new edges into the adjacency list
    update_adjacent2vertex_addition!(𝒜⁻¹, i, j, k, r)
    add_neighbour!(𝒟𝒢, r, i, j, k)
    add_neighbour!(𝒟𝒢, i, r)
    add_neighbour!(𝒟𝒢, j, r)
    add_neighbour!(𝒟𝒢, k, r)
    return nothing
end

"""
    flip_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, i, j, k, r)

Performs an edge flip, flipping the edge `(i, j)` into the edge `(k, r)`.

# Arguments
- `𝒯`: The current triangulation.
- `𝒟`: The point location data structure.
- `𝒜`: The adjacency list.
- `𝒜⁻¹`: The adjacent-to-vertex list.
- `𝒟𝒢`: The vertex-neighbour data structure.
- `i, j`: The current edge.
- `k, r`: Indices for the points the edge is flipped onto.

It is assumed that `(i, k, j)` and `(i, j, r)` are positively oriented triangles.

# Outputs 
`𝒯`, `𝒟`, and `𝒜` are all updated in-place.
"""
function flip_edge!(𝒯::Triangles, 𝒟::HistoryDAG,
    𝒜::Adjacent, 𝒜⁻¹::Adjacent2Vertex, 𝒟𝒢::DelaunayGraph, i, j, k, r)
    # The old triangles
    Tᵢₖⱼ = TriangleType((i, k, j))
    Tᵢⱼᵣ = TriangleType((i, j, r))
    delete_triangle!(𝒯, Tᵢₖⱼ, Tᵢⱼᵣ)
    delete_edge!(𝒜, i, j)
    delete_neighbour!(𝒟𝒢, i, j) #delete_neighbour!(𝒟𝒢, j, i)
    # The new triangles 
    Tᵣₖⱼ = TriangleType((r, k, j))
    Tᵣᵢₖ = TriangleType((r, i, k))
    # Add the new triangles to the data structure
    add_triangle!(𝒯, Tᵣₖⱼ, Tᵣᵢₖ)
    add_triangle!(𝒟, Tᵣₖⱼ, Tᵣᵢₖ)
    update_adjacent!(𝒜, Tᵣₖⱼ, Tᵣᵢₖ)
    update_adjacent2vertex_flip!(𝒜⁻¹, i, j, k, r)
    # Connect the new triangles to the replaced triangles in the DAG
    add_edge!(𝒟, Tᵢₖⱼ, Tᵣₖⱼ, Tᵣᵢₖ)
    add_edge!(𝒟, Tᵢⱼᵣ, Tᵣₖⱼ, Tᵣᵢₖ)
    # Add the new neighbours 
    add_neighbour!(𝒟𝒢, r, k) # add_neighbour!(𝒟𝒢, k, r)
    return nothing
end

"""
    legalise_edge!(𝒯, 𝒟, 𝒜, i, j, r, pts)
    
Legalises the edge `(i, j)` if it is illegal.

# Arguments 
- `𝒯`: The current triangulation.
- `𝒟`: The point location data structure.
- `𝒜`: The adjacency list.
- `𝒟𝒢`: The vertex-neighbour data structure.
- `i, j`: The edge to make legal. Nothing happens if `is_legal(i, j, 𝒜, pts)`.
- `r`: The point being added into the triangulation. 
- `pts`: The point set of the triangulation.

# Outputs 
`𝒯`, `𝒟`, `𝒜`, and `𝒟𝒢` are all updated in-place.
"""
function legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, i, j, r, pts)
    @show i, j, r
    if !is_legal(i, j, 𝒜, pts)
        e = 𝒜(j, i)
        @show i, j, e, r
        flip_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, i, j, e, r)
        legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, i, e, r, pts)
        legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, e, j, r, pts)
    end
    return nothing
end

"""
    initialise_triangulation()

This function returns the initial data structures for the Delaunay triangulation:

- `𝒯`: Data structure to contain the list of triangles.
- `𝒟`: The directed acyclic graph storing the history of the triangulation. 
- `𝒜`: The adjacency list.
- `𝒟𝒢`: A dictionary that maps points to their neighbours.
- `root`: The root of `𝒟`, `𝒯[begin]`.
"""
function initialise_triangulation()
    # The data structures
    root = TriangleType((LargeRightIdx + 1, LargeLeftIdx, LargeRightIdx))
    𝒯 = Triangles(Set([root]))
    𝒟 = HistoryDAG()
    𝒜 = Adjacent()
    𝒜⁻¹ = Adjacent2Vertex()
    𝒟𝒢 = DelaunayGraph()
    # Add the root to the DAG
    add_triangle!(𝒟, root)
    # Add the initial adjacencies 
    𝒜[(LargeRightIdx + 1, LargeLeftIdx)] = LargeRightIdx
    𝒜[(LargeLeftIdx, LargeRightIdx)] = LargeRightIdx + 1
    𝒜[(LargeRightIdx, LargeRightIdx + 1)] = LargeLeftIdx
    𝒜[(LargeLeftIdx, LargeRightIdx + 1)] = BoundaryIdx
    𝒜[(LargeRightIdx, LargeLeftIdx)] = BoundaryIdx
    𝒜[(LargeRightIdx + 1, LargeRightIdx)] = BoundaryIdx
    𝒜⁻¹[LargeRightIdx] = Set([(LargeRightIdx + 1, LargeLeftIdx)])
    𝒜⁻¹[LargeLeftIdx] = Set([(LargeRightIdx, LargeRightIdx + 1)])
    𝒜⁻¹[LargeRightIdx+1] = Set([(LargeLeftIdx, LargeRightIdx)])
    # Add the initial neighbours 
    add_neighbour!(𝒟𝒢, LargeRightIdx + 1, LargeLeftIdx, LargeRightIdx)
    add_neighbour!(𝒟𝒢, LargeLeftIdx, LargeRightIdx, LargeRightIdx + 1)
    add_neighbour!(𝒟𝒢, LargeRightIdx, LargeRightIdx + 1, LargeLeftIdx)
    return 𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, root
end

"""
    remove_bounding_triangle!(𝒯::Triangles, 𝒜::Adjacent, 𝒜⁻¹::Adjacent2Vertex, 𝒟𝒢::DelaunayGraph)

Remove the bounding triangle from the triangulation.
"""
function remove_bounding_triangle!(𝒯::Triangles, 𝒜::Adjacent, 𝒜⁻¹::Adjacent2Vertex, 𝒟𝒢::DelaunayGraph)
    for w in (-1, 0)
        neighbours = 𝒜⁻¹[w]
        for (u, v) in neighbours # (u, v, w) is a triangle..
            delete_edge!(𝒜, w, u; protect_boundary=false)
            delete_edge!(𝒜, w, v; protect_boundary=false)
            delete_edge!(𝒜⁻¹, u, v, w)
            delete_edge!(𝒜⁻¹, v, w, u)
            if u ≥ 1 && v ≥ 1 # This can only be a boundary edge
                add_edge!(𝒜⁻¹, BoundaryIdx, u, v)
                𝒜[(u, v)] = BoundaryIdx
            end
            delete_triangle!(𝒯, TriangleType((u, v, w)))
        end
        delete_point!(𝒟𝒢, w)
        delete_point!(𝒜⁻¹, w)
    end
    return nothing
end
function remove_bounding_triangle!(𝒟𝒯::Triangulation)
    remove_bounding_triangle!(triangles(𝒟𝒯), adjacent(𝒟𝒯), adjacent2vertex(𝒟𝒯), graph(𝒟𝒯))
    return nothing
end

function add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, root, pts, r)
    pᵣ = pts[r]
    𝒯ᵢⱼₖ, interior_flag = locate_triangle(𝒟, pts, pᵣ, root)
    i, j, k = 𝒯ᵢⱼₖ
    @show interior_flag
    if interior_flag == 1
        add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, 𝒯ᵢⱼₖ, r)
        legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, i, j, r, pts)
        legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, j, k, r, pts)
        legalise_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, k, i, r, pts)
    else
        eᵢⱼ, pₖ = locate_edge(pᵣ, 𝒯ᵢⱼₖ)
        𝒯ᵢₗⱼ = adjacent_triangles(𝒯, 𝒯ᵢⱼₖ, eᵢⱼ)
        pₗ = select_adjacent_vertex(𝒯, eᵢⱼ, 𝒯ᵢₗⱼ)
        add_edges!(𝒯, 𝒟, pᵣ, pₖ)
        add_edges!(𝒯, 𝒟, pᵣ, pₗ)
        legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pᵢ, pₗ)
        legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pₗ, pⱼ)
        legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pⱼ, pₖ)
        legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pₖ, pᵢ)
    end
end

function triangulate(pts; sort_pts=true, shuffle_pts=true, trim=true)
    Base.require_one_based_indexing(pts)
    𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, root = initialise_triangulation()
    sort_pts && partial_highest_point_sort!(pts, 1) # p0 = pts[begin]
    shuffle_pts && @views shuffle!(pts[begin+1:end])
    for r = (firstindex(pts)+1):lastindex(pts)
        add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢, root, pts, r)
    end
    trim && remove_bounding_triangle!(𝒯, 𝒜, 𝒜⁻¹, 𝒟𝒢)
    return 𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒟𝒢
end