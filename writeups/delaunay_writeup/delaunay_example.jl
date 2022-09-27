using SimpleGraphs
using ExactPredicates
using ExactPredicates.Codegen
using Test
using Random
const TriangleType = NTuple{3,Int64}
const LargeRightIdx = 0 # p₋₁
const LargeLeftIdx = -1 # p₋₂
const EmptyIdx = -2 # ∅
const VertexNeighbourVector = Vector{Int64}
const AdjacentToVertexVector = Vector{NTuple{2, Int64}}
"""
    TriangulationDAG <: AbstractSimpleGraph

Data structure for the Delaunay triangulation. This is a directed acyclic graph 
that is useful for point location.
"""
struct TriangulationDAG <: AbstractSimpleGraph
    graph::DirectedGraph{TriangleType}
    function TriangulationDAG()
        G = DirectedGraph{TriangleType}()
        forbid_loops!(G)
        TDAG = new(G)
        return TDAG
    end
end
const TriDAG = TriangulationDAG
graph(G::TriangulationDAG) = G.graph
for op in (:out_deg, :in_deg, :deg,
    :dual_deg, :in_neighbors,
    :NE, :has, :delete!,
    :elist)
    @eval begin
        SimpleGraphs.$op(G::TriDAG) = SimpleGraphs.$op(graph(G))
        SimpleGraphs.$op(G::TriDAG, u) = SimpleGraphs.$op(graph(G), u)
        SimpleGraphs.$op(G::TriDAG, u, v) = SimpleGraphs.$op(graph(G), u, v)
    end
end
SimpleGraphs.out_neighbors(G::TriDAG, u) = graph(G).N[u] # The version in SimpleGraphs is allocating...
SimpleGraphs.add!(G::TriDAG, u) = SimpleGraphs.add!(graph(G), u)
"""
    SimpleGraphs.add!(G::TriDAG, u)

Add the triangle `u` to the graph `G`. Checks that cyclic permutations of `u`'s indices are not 
in `G` already, returning `false` otherwise.
"""
function SimpleGraphs.add!(G::TriDAG, u)
    for i = 1:3
        i₁, i₂, i₃ = i, (i % 3) + 1, ((i + 1) % 3) + 1
        has(G, TriangleType((u[i₁], u[i₂], u[i₃]))) && return false
    end
    return SimpleGraphs.add!(graph(G), u)
end
"""
    SimpleGraphs.add!(G::TriDAG, u, v)

Adds an edge between triangle `u ∈ G` and `v ∈ G` to the graph `G`. Checks 
are made for cyclic permutations of `u` and `v`'s indices.
"""
function SimpleGraphs.add!(G::TriDAG, u, v)
    for i in 1:3
        for j in 1:3
            i₁, i₂, i₃ = i, (i % 3) + 1, ((i + 1) % 3) + 1
            j₁, j₂, j₃ = j, (j % 3) + 1, ((j + 1) % 3) + 1
            T = TriangleType((u[i₁], u[i₂], u[i₃]))
            V = TriangleType((v[j₁], v[j₂], v[j₃]))
            if has(G, T) && has(G, V)
                SimpleGraphs.add!(graph(G), T, V)
                return true
            end
        end
    end
    return false
end

"""
    find_root(G::TriDAG; method=:brute)

Finds the root of the graph `G`, assuming only one root exists and the only node with in-degree zero is this root. 
There are two methods:

    - `method=:brute`: In this case, all vertices are searched until one is found with in-degree zero.
    - `method=:rng`: In this case, a random integer is selected, defniing the node to start at. 
    We then go upwards until the root is found; note that this only works as the graph is acyclic.
"""
function find_root(G::TriDAG; method=:brute)
    if method == :brute
        return _find_root_brute(G)
    elseif method == :rng
        return _find_root_rng(G)
    end
end
function _find_root_brute(G::TriDAG)
    for (k, v) in graph(G).NN
        length(v) == 0 && return k
    end
end
function _find_root_rng(G::TriDAG)
    verts = vlist(graph(G))
    num_verts = length(verts)
    starting_node = verts[rand(1:num_verts)]
    for _ in 1:num_verts
        num_choices = in_deg(G, starting_node)
        num_choices == 0 && return starting_node
        starting_node = in_neighbors(G, starting_node)[rand(1:num_choices)]
    end
end

"""
    TriangulationAdjacent{A}

This is a data structure for the Delaunay triangulation. It stores a `Dict`, `adjacent`, so that 
`(u, v, adjacent[(u, v)])` is a positively oriented triangle. This struct is also callable, e.g. if 
`𝒜::TriangulationAdjacent`, then `(u, v, 𝒜(u, v))` is a positively oriented triangle.
"""
struct TriangulationAdjacent{A}
    adjacent::A
    function TriangulationAdjacent()
        D = Dict{NTuple{2,Int64},Int64}()
        TADJ = new{typeof(D)}(D) # TADJ = (T)riangulation(ADJ)acent
        return TADJ
    end
end
const TriAdjacent = TriangulationAdjacent
Base.setindex!(adj::TriAdjacent, w, uv) = Base.setindex!(adj.adjacent, w, uv) # (u, v, w) is a positively oriented triangle
Base.getindex(adj::TriAdjacent, uv) = Base.getindex(adj.adjacent, uv) # Find the w so that (u, v, w) is a positively oriented triangle
edges(adj::TriAdjacent) = keys(adj.adjacent) # List of edges - note that there are duplicates, i.e. (u, v) and (v, u) are both stored.
(adj::TriAdjacent)(u, v) = adj[(u, v)] # Find the w so that (u, v, w) is a positively oriented triangle

""" 
    TriangulationAdjacentToVertex{A⁻¹}

This is a data structure for the Delaunay triangulation. It stores a `Dict`, `adjacent2vertex`, so that, 
for each `(u, v) ∈ adjacent2vertex[w]`, `(u, v, w)` is a positively oriented triangle. Note that 
this is the inverse of `TriangulationAdjacent`.
"""
struct TriangulationAdjacentToVertex{A⁻¹}
    adjacent2vertex::A⁻¹
    function TriangulationAdjacentToVertex()
        D = Dict{Int64, AdjacentToVertexVector}()
        TA2V = new{typeof(D)}(D)
        return TA2V 
    end 
end
const TriAdjacent2Vertex = TriangulationAdjacentToVertex
Base.setindex!(adj::TriAdjacent2Vertex, uv, w) = Base.setindex!(adj.adjacent2vertex, uv, w) # (u, v, w) is a positively oriented triangle
Base.getindex(adj::TriAdjacent2Vertex, w) = Base.getindex(adj.adjacent2vertex, w) # Find the (u, v) so that (u, v, w) is a positively oriented triangle

"""
    TriangulationVertexNeighbours{N}

This is a data structure for the Delaunay triangulation. It stores a `Dict`, `neighbours`, so 
that `neighbours[u]` is a list of points that share an edge with the point `u`.
"""
struct TriangulationVertexNeighbours{N}
    neighbours::N
    function TriangulationVertexNeighbours()
        D = Dict{Int64,VertexNeighbourVector}()
        TVN = new{typeof(D)}(D) # TVN = (T)riangulation(V)ertex(N)eighbours 
        return TVN
    end
    TriangulationVertexNeighbours(tvn::N) where {N} = new{N}(tvn)
end
const TriVertexNeighbours = TriangulationVertexNeighbours
Base.setindex!(tvn::TriVertexNeighbours, v::AbstractVector, u) = Base.setindex!(tvn.neighbours, v, u) # This will add v to the list of neighbours in tvn.neighbours[u]
Base.setindex!(tvn::TriVertexNeighbours, v::Integer, u) = Base.setindex!(tvn, [v], u)
Base.getindex(tvn::TriVertexNeighbours, u) = Base.getindex(tvn.neighbours, u) # Find the neighbours of u 
Base.push!(tvn::TriVertexNeighbours, u, v::AbstractVector) = Base.push!(tvn.neighbours[u], v...) # Push v onto the list of neighbours of u 
Base.push!(tvn::TriVertexNeighbours, u, v::Integer) = Base.push!(tvn.neighbours[u], v) # Push v onto the list of neighbours of u 
(tvn::TriVertexNeighbours)(u) = tvn.neighbours[u] # Call tvn to get the neighbours of u
points(tvn::TriVertexNeighbours) = keys(tvn.neighbours)
"""
add_neighbour!(tvn::TriVertexNeighbours, u::Integer, v...; new_point=u ∉ points(tvn))

Adds the point `v` to `tvn.neighbours[u]`, or introduces the key `u` and adds the point. Declare if this 
point is new or not (i.e. an existing key of `tvn.neighbours`) using the keyword `new_point`.
"""
function add_neighbour!(tvn::TriVertexNeighbours, u::Integer, v; new_point=u ∉ points(tvn))
    if !new_point
        push!(tvn, u, v)
        return nothing
    else
        tvn[u] = v
        return nothing
    end
end
function add_neighbour!(tvn::TriVertexNeighbours, u::Integer, v::Integer...; new_point=u ∉ points(tvn))
    if !new_point
        push!(tvn, u, v...)
        return nothing
    else
        tvn[u] = [v...]
    end
end

function ExactPredicates.orient(𝒯::TriangleType, pts)
    u, v, w = 𝒯
    return orient(pts[u], pts[v], pts[w])
end
function ExactPredicates.incircle(pts, i, j, k, ℓ)
    return incircle(pts[i], pts[j], pts[k], pts[ℓ])
end

"""
    leftofline(x₁, y₁, x₂, y₂, x₃, y₃)

Tests if the point (x₁, x₂) is to the left of the half-plane formed 
from the oriented line through (x₂, y₂) to (x₃, y₃).

This is the same as `orient((x₂, y₂), (x₁, y₁), (x₃, y₃))`.
"""
@genpredicate function leftofline(x₁, y₁, x₂, y₂, x₃, y₃)
    Codegen.group!(x₁, y₁, x₂, y₂, x₃, y₃)
    return x₁ * y₂ - x₂ * y₁ - x₁ * y₃ + x₃ * y₁ + x₂ * y₃ - x₃ * y₂
end
"""
    leftofline(pts, p, i, j)

Tests if the point `p` is to the left of the oriented line through 
`pts[i]` to `pts[j]`. Checks are made for non-positive indices.
"""
function leftofline(pts, p, i, j)
    if j == LargeRightIdx && i > LargeRightIdx      # p₋₁ → pᵢ
        return p > pts[i] ? 1 : -1
    elseif j == LargeLeftIdx && i > LargeRightIdx   # p₋₂ → pᵢ
        return p < pts[i] ? 1 : -1
    elseif i == LargeRightIdx && j > LargeRightIdx  # p₋₁ → pᵢ
        return p < pts[j] ? 1 : -1
    elseif i == LargeLeftIdx && j > LargeRightIdx   # p₋₂ → pᵢ
        return p > pts[j] ? 1 : -1
    elseif i == LargeRightIdx && j == LargeLeftIdx  # p₋₁ → p₋₂
        return -1
    elseif i == LargeLeftIdx && j == LargeRightIdx  # p₋₂ → p₋₁
        return 1
    end
    x, y = p
    x₂, y₂ = pts[i]
    x₃, y₃ = pts[j]
    return leftofline(x, y, x₂, y₂, x₃, y₃)
end

function intriangle(e1, e2, e3) # https://stackoverflow.com/a/2049593
    if e1 == 0 || e2 == 0 || e3 == 0
        return 0
    end
    has_neg = e1 < 0 || e2 < 0 || e3 < 0
    has_pos = e1 > 0 || e2 > 0 || e3 > 0
    if has_neg && has_pos
        return -1
    else
        return 1
    end
end
function intriangle(𝒯::TriangleType, pts, p)
    i, j, k = 𝒯
    e1 = leftofline(pts, p, i, j)
    e2 = leftofline(pts, p, j, k)
    e3 = leftofline(pts, p, k, i)
    return intriangle(e1, e2, e3)
end
function intriangle(x₁, y₁, x₂, y₂, x₃, y₃, x, y)
    e1 = leftofline(x, y, x₁, y₁, x₂, y₂)
    e2 = leftofline(x, y, x₂, y₂, x₃, y₃)
    e3 = leftofline(x, y, x₃, y₃, x₁, y₁)
    return intriangle(e1, e2, e3)
end
function intriangle(p₁, p₂, p₃, p)
    x, y = p
    x₁, y₁ = p₁
    x₂, y₂ = p₂
    x₃, y₃ = p₃
    return intriangle(x₁, y₁, x₂, y₂, x₃, y₃, x, y)
end

"""
    locate_triangle(𝒟::TriDAG, pts, p, init=find_root(𝒟; method=:rng))

Given the point location data structure `D` and a set of `pts`, finds the triangle in 
the current triangulation such that `p` is in its interior. The point location starts at `init`.
The function is recursive, and returns a tuple `(tri, flag)`:
    - `tri`: This is the triangle that `p` is in.
    - `flag`: If `flag == 0`, then `p` is on an edge of `tri`. Otherwise, it is in the open interior.
"""
function locate_triangle(𝒟::TriDAG, pts, p, init=find_root(𝒟; method=:rng))
    # Find which triangle in 𝒟, containing points from the point set pts, 
    #   contains the point p, starting from the node root.
    if out_deg(𝒟, init) == 0
        return init, intriangle(init, pts, p)
    end
    out = out_neighbors(𝒟, init)
    for T in out
        intriangle(T, pts, p) ≥ 0 && return locate_triangle(𝒟, pts, p, T)
    end
end

"""
    add_edge!(𝒜::TriAdjacent, u, v, w)

This adds the edge `(u, v)` to the adjacency list `𝒜`, mapping 
`(u, v)` to `w` so that `(u, v, w)` is a positively oriented triangle.
"""
function add_edge!(𝒜::TriAdjacent, u, v, w)
    𝒜[(u, v)] = w
    return nothing
end

"""
    add_edge!(𝒜⁻¹::TriAdjacent2Vertex, w, u, v)

This adds the edge `(u, v)` to `𝒜⁻¹[w]`.
"""
function add_edge!(𝒜⁻¹::TriAdjacent2Vertex, w, u, v)
    uv = get!(AdjacentToVertexVector, 𝒜⁻¹, w) # https://discourse.julialang.org/t/how-do-i-append-add-data-in-dictionary-with-same-key/15891/5?u=legola18
    push!(uv, (u, v))
    return nothing
end

"""
    update_adjacent!(𝒜::TriAdjacent, 𝒜⁻¹::TriAdjacent2Vertex,, 𝒯::TriangleType...)

Updates the adjacency list `𝒜`, and its inverse `𝒜⁻¹`, with the edges from `𝒯`, keying on each edge. 
`𝒯` must be positively oriented.
"""
function update_adjacent!(𝒜::TriAdjacent, 𝒜⁻¹::TriAdjacent2Vertex, 𝒯::TriangleType...)
    for T in 𝒯
        i, j, k = T
        add_edge!(𝒜, i, j, k)
        add_edge!(𝒜, j, k, i)
        add_edge!(𝒜, k, i, j)
        add_edge!(𝒜⁻¹, k, i, j)
        add_edge!(𝒜⁻¹, i, j, k)
        add_edge!(𝒜⁻¹, j, k, i)
    end
    return nothing
end

"""
    delete_triangle!(𝒯, T...)

Deletes the triangle(s) `T` from `𝒯`, and all circular shifts of their indices.
"""
function delete_triangle!(𝒯, T...)
    for T in T
        deleteat!(𝒯, findall(x -> x == T || x == T[[2, 3, 1]] || x == T[[3, 1, 2]], 𝒯))
    end
    return nothing
end

"""
    add_triangle!(𝒯, T...)

Adds the triangle(s) `T` to `𝒯`.
"""
function add_triangle!(𝒯, T...)
    push!(𝒯, T...)
    return nothing
end
function add_triangle!(𝒟::TriDAG, 𝒯...)
    for T in 𝒯
        add!(𝒟, T)
    end
    return nothing
end

function add_edge!(𝒟::TriDAG, 𝒯, 𝒯new...)
    for T in 𝒯new
        add!(𝒟, 𝒯, T)
    end
    return nothing
end

"""
    delete_edge_from_adjacency!(𝒜, i, j; protect_boundary = true)

Deletes the keys `(i, j)` and `(j, i)` from `𝒜`. Use `protect_boundary=true` to avoid 
removing edges on the boundary.
"""
function delete_edge_from_adjacency!(𝒜, i, j; protect_boundary=true)
    if protect_boundary
        𝒜(i, j) ≠ EmptyIdx && delete!(𝒜.adjacent, (i, j))
        𝒜(j, i) ≠ EmptyIdx && delete!(𝒜.adjacent, (j, i))
    else
        delete!(𝒜.adjacent, (i, j))
        delete!(𝒜.adjacent, (j, i))
    end
    wᵢⱼ = 𝒜(i, j)
    wⱼᵢ = 𝒜(j, i)
    return nothing
end

"""
    delete_point_from_neighbour!(𝒱𝒩, u, v)

Removes the point `v` from the neighbourhood of `u`, updating `𝒱𝒩` in-place.
"""
function delete_point_from_neighbour!(𝒱𝒩::TriVertexNeighbours, u, v)
    idx = findfirst(x -> x == v, 𝒱𝒩(u))
    deleteat!(𝒱𝒩[u], idx)
    return nothing
end

"""
    add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, Tᵢⱼₖ, r)

Given a triangulation `𝒯`, adds the `r`th point of the point set into the triangulation.

# Arguments 
- `𝒯`: The current triangulation.
- `𝒟`: The point location data structure.
- `𝒜`: The adjacency list.
- `𝒜⁻¹`: The adjacent-to-vertex list.
- `𝒱𝒩`: The vertex-neighbour data structure.
-` Tᵢⱼₖ`: The triangle that the `r`th point is inside of. Must be positively oriented.
- `r`: The index of the point in the original point set that is being introduced.

# Outputs 
`𝒯`, `𝒟`, and `𝒜` are all updated in-place.
"""
function add_point!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒱𝒩, Tᵢⱼₖ, r)
    i, j, k = Tᵢⱼₖ # The triangle to be split into three
    delete_triangle!(𝒯, Tᵢⱼₖ) # Now that we've split the triangle, we can remove the triangle
    T₁, T₂, T₃ = TriangleType((i, j, r)), TriangleType((j, k, r)), TriangleType((k, i, r)) # New triangles to add. Note that these triangles are all positively oriented.
    add_triangle!(𝒯, T₁, T₂, T₃) # The three new triangles
    add_triangle!(𝒟, T₁, T₂, T₃) # Add the new triangles into DAG
    add_edge!(𝒟, Tᵢⱼₖ, T₁, T₂, T₃) # Add edges from the old triangle to the new triangles
    update_adjacent!(𝒜,𝒜⁻¹, T₁, T₂, T₃) # Add the new edges into the adjacency list
    add_neighbour!(𝒱𝒩, r, i, j, k; new_point=true)
    add_neighbour!(𝒱𝒩, i, r; new_point=false)
    add_neighbour!(𝒱𝒩, j, r; new_point=false)
    add_neighbour!(𝒱𝒩, k, r; new_point=false)
    return nothing
end

"""
    edge_on_large_triangle(i, j)

Returns true if `(i, j)` is an edge of the triangle `(1, -1, 0)`.
"""
function edge_on_large_triangle(i, j)
    if i > (LargeRightIdx + 1) || j > (LargeRightIdx + 1) # = 1 case can be p₀ 
        return false
    elseif (i, j) == (LargeRightIdx + 1, LargeRightIdx) ||
           (i, j) == (LargeRightIdx, LargeLeftIdx) ||
           (i, j) == (LargeLeftIdx, LargeRightIdx + 1) ||
           (i, j) == (LargeRightIdx, LargeRightIdx + 1) ||
           (i, j) == (LargeLeftIdx, LargeRightIdx) ||
           (i, j) == (LargeRightIdx + 1, LargeLeftIdx)
        return true
    else
        return false
    end
end

"""
    is_legal(i, j, 𝒜, pts)

Tests if the edge `(i, j)` is a legal edge. `𝒜` is the adjacency list of the triangulation, and `pts` is the point set.
Returns `true` if the edge is legal.
"""
function is_legal(i, j, k, ℓ, pts)
    if i > LargeRightIdx && j > LargeRightIdx && k > LargeRightIdx && ℓ > LargeRightIdx
        return incircle(pts, i, j, k, ℓ) ≤ 0
    else
        return min(k, ℓ) < min(i, j)
    end
end
function is_legal(i, j, 𝒜, pts)
    edge_on_large_triangle(i, j) && return true
    k, ℓ = 𝒜(i, j), 𝒜(j, i)
    return is_legal(i, j, k, ℓ, pts)
end

"""
    flip_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, i, j, k, r)

Performs an edge flip, flipping the edge `(i, j)` into the edge `(k, r)`.

# Arguments
- `𝒯`: The current triangulation.
- `𝒟`: The point location data structure.
- `𝒜`: The adjacency list.
- `𝒜⁻¹`: The adjacent-to-vertex list.
- `𝒱𝒩`: The vertex-neighbour data structure.
- `i, j`: The current edge.
- `k, r`: Indices for the points the edge is flipped onto.

It is assumed that `(i, k, j)` and `(i, j, r)` are positively oriented triangles.

# Outputs 
`𝒯`, `𝒟`, and `𝒜` are all updated in-place.
"""
function flip_edge!(𝒯, 𝒟, 𝒜, 𝒜⁻¹, 𝒱𝒩, i, j, k, r)
    # The old triangles
    Tᵢₖⱼ = TriangleType((i, k, j))
    Tᵢⱼᵣ = TriangleType((i, j, r))
    delete_triangle!(𝒯, Tᵢₖⱼ, Tᵢⱼᵣ)
    delete_edge_from_adjacency!(𝒜, i, j)
    delete_point_from_neighbour!(𝒱𝒩, i, j)
    delete_point_from_neighbour!(𝒱𝒩, j, i)
    # The new triangles 
    Tᵣₖⱼ = TriangleType((r, k, j))
    Tᵣᵢₖ = TriangleType((r, i, k))
    # Add the new triangles to the data structure
    add_triangle!(𝒯, Tᵣₖⱼ, Tᵣᵢₖ)
    add_triangle!(𝒟, Tᵣₖⱼ, Tᵣᵢₖ)
    update_adjacent!(𝒜, Tᵣₖⱼ)
    update_adjacent!(𝒜, Tᵣᵢₖ)
    # Connect the new triangles to the replaced triangles in the DAG
    add_edge!(𝒟, Tᵢₖⱼ, Tᵣₖⱼ, Tᵣᵢₖ)
    add_edge!(𝒟, Tᵢⱼᵣ, Tᵣₖⱼ, Tᵣᵢₖ)
    # Add the new neighbours 
    add_neighbour!(𝒱𝒩, r, k; new_point=false)
    add_neighbour!(𝒱𝒩, k, r; new_point=false)
    return nothing
end

"""
    legalise_edge!(𝒯, 𝒟, 𝒜, i, j, r, pts)
    
Legalises the edge `(i, j)` if it is illegal.

# Arguments 
- `𝒯`: The current triangulation.
- `𝒟`: The point location data structure.
- `𝒜`: The adjacency list.
- `𝒱𝒩`: The vertex-neighbour data structure.
- `i, j`: The edge to make legal. Nothing happens if `is_legal(i, j, 𝒜, pts)`.
- `r`: The point being added into the triangulation. 
- `pts`: The point set of the triangulation.

# Outputs 
`𝒯`, `𝒟`, `𝒜`, and `𝒱𝒩` are all updated in-place.
"""
function legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, i, j, r, pts)
    if !is_legal(i, j, 𝒜, pts)
        k = 𝒜(j, i)
        flip_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, i, j, k, r)
        legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, i, k, r, pts)
        legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, k, j, r, pts)
    end
    return nothing
end

"""
    initialise_triangulation()

This function returns the initial data structures for the Delaunay triangulation:

- `𝒯`: Data structure to contain the list of triangles.
- `𝒟`: The directed acyclic graph storing the history of the triangulation. 
- `𝒜`: The adjacency list.
- `𝒱𝒩`: A dictionary that maps points to their neighbours.
- `root`: The root of `𝒟`, `𝒯[begin]`.
"""
function initialise_triangulation()
    # The data structures
    𝒯 = TriangleType[(LargeRightIdx + 1, LargeLeftIdx, LargeRightIdx)]
    𝒟 = TriDAG()
    𝒜 = TriAdjacent()
    𝒱𝒩 = TriVertexNeighbours()
    # Add the root to the DAG
    add_triangle!(𝒟, 𝒯[begin])
    root = 𝒯[begin]
    # Add the initial adjacencies 
    𝒜[(LargeRightIdx + 1, LargeLeftIdx)] = LargeRightIdx
    𝒜[(LargeLeftIdx, LargeRightIdx)] = LargeRightIdx + 1
    𝒜[(LargeRightIdx, LargeRightIdx + 1)] = LargeLeftIdx
    𝒜[(LargeLeftIdx, LargeRightIdx + 1)] = EmptyIdx
    𝒜[(LargeRightIdx, LargeLeftIdx)] = EmptyIdx
    𝒜[(LargeRightIdx + 1, LargeRightIdx)] = EmptyIdx
    # Add the initial neighbours 
    add_neighbour!(𝒱𝒩, LargeRightIdx + 1, LargeLeftIdx, LargeRightIdx; new_point=true)
    add_neighbour!(𝒱𝒩, LargeLeftIdx, LargeRightIdx, LargeRightIdx + 1; new_point=true)
    add_neighbour!(𝒱𝒩, LargeRightIdx, LargeRightIdx + 1, LargeLeftIdx; new_point=true)
    return 𝒯, 𝒟, 𝒜, 𝒱𝒩, root
end


function remove_bounding_triangle!()
    # First look at p₋₁
    p₋₁_neighbours = 𝒱𝒩(-1)
    for i in p₋₁_neighbours 
        delete_edge_from_adjacency!(𝒜, -1, i)
        delete_point_from_neighbour!(𝒱𝒩, -1, i)
end

# Test that we can correctly construct and index the vertex neighbour structure 
tvn = Dict(1 => [4, 5, 6, 9], 2 => [11, 9, 8, 15, 16], 3 => [1])
𝒱𝒩 = TriVertexNeighbours(tvn)
𝒱𝒩2 = TriangulationVertexNeighbours(tvn)
@test 𝒱𝒩2.neighbours == 𝒱𝒩.neighbours
@test 𝒱𝒩.neighbours == tvn
@test 𝒱𝒩(1) == [4, 5, 6, 9]
@test 𝒱𝒩(2) == [11, 9, 8, 15, 16]
@test 𝒱𝒩(3) == [1]
@test 𝒱𝒩[1] == [4, 5, 6, 9]
@test 𝒱𝒩[2] == [11, 9, 8, 15, 16]
@test 𝒱𝒩[3] == [1]
𝒱𝒩[4] = 5
@test 𝒱𝒩[4] == [5]
𝒱𝒩[5] = [6]
@test 𝒱𝒩[5] == [6]
𝒱𝒩[7] = [1, 2, 3, 4, 5, 6, 7]
@test 𝒱𝒩[7] == [1, 2, 3, 4, 5, 6, 7]
push!(𝒱𝒩, 3, [1, 2, 3])
@test 𝒱𝒩[3] == [1, 1, 2, 3]
push!(𝒱𝒩, 5, 1)
@test 𝒱𝒩[5] == [6, 1]
push!(𝒱𝒩, 5, 3, 5, 1, 5)
@test 𝒱𝒩[5] == [6, 1, 3, 5, 1, 5]
@test points(𝒱𝒩) == keys(𝒱𝒩.neighbours)
@test collect(points(𝒱𝒩)) == [5, 4, 7, 2, 3, 1]
add_neighbour!(𝒱𝒩, 1, 10)
@test 𝒱𝒩(1) == [4, 5, 6, 9, 10]
add_neighbour!(𝒱𝒩, 1, 11, 12, 13)
@test 𝒱𝒩(1) == [4, 5, 6, 9, 10, 11, 12, 13]
add_neighbour!(𝒱𝒩, 1, [14, 19])
@test 𝒱𝒩(1) == [4, 5, 6, 9, 10, 11, 12, 13, 14, 19]
add_neighbour!(𝒱𝒩, 10, 14)
@test 𝒱𝒩(10) == [14]
add_neighbour!(𝒱𝒩, 11, [15])
@test 𝒱𝒩(11) == [15]
add_neighbour!(𝒱𝒩, 13, 15, 16, 19, 20)
@test 𝒱𝒩(13) == [15, 16, 19, 20]
𝒱𝒩_empty = TriVertexNeighbours()
@test isempty(𝒱𝒩_empty.neighbours)
@test isempty(points(𝒱𝒩_empty))
add_neighbour!(𝒱𝒩_empty, 2, [2, 3, 4, 5, 1])
@test 𝒱𝒩_empty[2] == [2, 3, 4, 5, 1]
@test_throws KeyError push!(𝒱𝒩_empty, 6, 10)
𝒱𝒩_empty[15] = 20
@test 𝒱𝒩_empty(15) == 𝒱𝒩_empty[15] == [20]
add_neighbour!(𝒱𝒩_empty, 15, 27; new_point=false)
@test 𝒱𝒩_empty[15] == [20, 27]
@test_throws KeyError add_neighbour!(𝒱𝒩_empty, 20, 27; new_point=false)
add_neighbour!(𝒱𝒩_empty, 15, 273; new_point=false)
@test 𝒱𝒩_empty[15] == [20, 27, 273]
@test_throws KeyError add_neighbour!(𝒱𝒩_empty, 20, 27, 53, 103; new_point=false)
add_neighbour!(𝒱𝒩_empty, 15, 273, 109, 81; new_point=false)
@test 𝒱𝒩_empty[15] == [20, 27, 273, 273, 109, 81]
add_neighbour!(𝒱𝒩_empty, 29, 31)
@test 𝒱𝒩_empty[29] == [31]
add_neighbour!(𝒱𝒩_empty, 37, 38, 39, 30, 1, 1, 2)
@test 𝒱𝒩_empty[37] == [38, 39, 30, 1, 1, 2]

## Test that we can delete points from a neighbourhood 
tvn = Dict(1 => [4, 5, 6, 9], 2 => [11, 9, 8, 15, 16], 3 => [1])
𝒱𝒩 = TriVertexNeighbours(tvn)
delete_point_from_neighbour!(𝒱𝒩, 1, 6)
@test 𝒱𝒩[1] == [4, 5, 9]
@test 𝒱𝒩[2] == [11, 9, 8, 15, 16]
@test 𝒱𝒩[3] == [1]

# Test that we can correctly add and delete triangles 
𝒯 = TriangleType[(1, 2, 3), (3, 2, 5), (2, 1, 9)]
delete_triangle!(𝒯, (1, 2, 3))
@test 𝒯 == TriangleType[(3, 2, 5), (2, 1, 9)]
delete_triangle!(𝒯, (5, 3, 2))
@test 𝒯 == TriangleType[(2, 1, 9)]
add_triangle!(𝒯, (2, 3, 10))
@test 𝒯 == TriangleType[(2, 1, 9), (2, 3, 10)]
add_triangle!(𝒯, (2, 3, 11), (11, 3, 4), (2, 3, 1))
@test 𝒯 == TriangleType[(2, 1, 9), (2, 3, 10), (2, 3, 11), (11, 3, 4), (2, 3, 1)]
delete_triangle!(𝒯, (2, 1, 9), (10, 2, 3))
@test 𝒯 == TriangleType[(2, 3, 11), (11, 3, 4), (2, 3, 1)]

# Test that we can correctly remove adjacencies 
𝒜n = TriAdjacent()
i, j, k, r = 1, 2, 3, 4
𝒜n[(i, j)] = r
𝒜n[(j, r)] = i
𝒜n[(r, i)] = j
𝒜n[(i, k)] = j
𝒜n[(k, j)] = i
𝒜n[(j, i)] = k
delete_edge_from_adjacency!(𝒜n, i, j)
@test length(𝒜n.adjacent) == 4
@test (i, j) ∉ keys(𝒜n.adjacent)
@test (j, i) ∉ keys(𝒜n.adjacent)
@test 𝒜n(j, r) == i
@test 𝒜n(r, i) == j
@test 𝒜n(i, k) == j
@test 𝒜n(k, j) == i

# Test that we can correctly flip an edge
I = [-2.0, -2.0]
J = [-10.0, 2.0]
K = [-4.0, 4.0]
R = [-8.0, -3.0]
IJKR = [I, J, K, R]
i, j, k, r = 1, 2, 3, 4
𝒟n = TriDAG()
Tᵢₖⱼ = TriangleType((i, k, j))
Tᵢⱼᵣ = TriangleType((i, j, r))
𝒯n = TriangleType[Tᵢₖⱼ, Tᵢⱼᵣ]
add!(𝒟n, Tᵢₖⱼ)
add!(𝒟n, Tᵢⱼᵣ)
𝒜n = TriAdjacent()
𝒜n[(i, j)] = r
𝒜n[(j, r)] = i
𝒜n[(r, i)] = j
𝒜n[(i, k)] = j
𝒜n[(k, j)] = i
𝒜n[(j, i)] = k
𝒱𝒩n = TriVertexNeighbours(Dict(
    i => [k, j, r],
    j => [r, i, k],
    k => [j, i],
    r => [i, j]
))
𝒯𝒯 = deepcopy(𝒯n)
𝒟𝒟 = deepcopy(𝒟n)
𝒜𝒜 = deepcopy(𝒜n)
𝒱𝒱 = deepcopy(𝒱𝒩n)
@test !is_legal(i, j, 𝒜n, IJKR)
flip_edge!(𝒯n, 𝒟n, 𝒜n, 𝒱𝒩n, i, j, k, r)
@test is_legal(r, k, 𝒜n, IJKR)
Tᵣₖⱼ = TriangleType((r, k, j))
Tᵣᵢₖ = TriangleType((r, i, k))
@test length(𝒜n.adjacent) == 6
@test (i, j) ∉ keys(𝒜n.adjacent)
@test (j, i) ∉ keys(𝒜n.adjacent)
@test 𝒯n == TriangleType[Tᵣₖⱼ, Tᵣᵢₖ]
@test 𝒜n(j, r) == k
@test 𝒜n(r, k) == j
@test 𝒜n(k, j) == r
@test 𝒜n(r, i) == k
@test 𝒜n(i, k) == r
@test 𝒜n(k, r) == i
@test 𝒟n.graph.N[(i, k, j)] == Set(TriangleType[(r, k, j), (r, i, k)])
@test 𝒟n.graph.N[(i, j, r)] == Set(TriangleType[(r, k, j), (r, i, k)])
@test all(==(1), orient.(𝒯n, Ref(IJKR)))
@test 𝒱𝒩n(i) == [k, r]
@test 𝒱𝒩n(j) == [r, k]
@test 𝒱𝒩n(r) == [i, j, k]
@test 𝒱𝒩n(k) == [j, i, r]
flip_edge!(𝒯n, 𝒟n, 𝒜n, 𝒱𝒩n, r, k, i, j) # This should go back to the original configuration
@test 𝒯n == TriangleType[(2, 1, 3), (2, 4, 1)]
@test 𝒜𝒜.adjacent == 𝒜n.adjacent
@test all(sort(𝒱𝒱(i)) == sort(𝒱𝒩n(i)) for i in 1:4)
@test !is_legal(i, j, 𝒜n, IJKR)

p0 = Float64[5, 5]
p1 = Float64[4.5, 2.5]
p2 = Float64[2.5, 1.5]
p3 = Float64[3, 3.5]
p4 = Float64[0, 2]
p5 = Float64[1, 5]
p6 = Float64[1, 3]
p7 = Float64[4, -1]
p8 = Float64[-1, 4]
pts = [p0, p1, p2, p3, p4, p5, p6, p7, p8]

# Test triangle orientation
@test orient((4, 6, 7), pts) == 1
@test orient((4, 7, 6), pts) == -1
@test orient((4, 2, 3), pts) == -1
@test orient((4, 7, 3), pts) == 1
@test orient((5, 7, 9), pts) == 1
@test orient((5, 9, 7), pts) == -1
@test orient((3, 8, 5), pts) == -1
@test orient((1, 2, 3), [[1.0, 2.0], [1.0, 5.0], [1.0, 8.0]]) == 0

# Test that we can find points to the left of a line 
G = Float64[-8, -13]
H = Float64[-3, -10]
I = Float64[-4, -13]
J = Float64[-6, -9]
xg, yg = G
xh, yh = H
xi, yi = I
x, y = J
@test leftofline(x, y, xg, yg, xh, yh) == 1 # To the left of line from g to h 
@test leftofline(x, y, xh, yh, xi, yi) == -1 # To the right of line from h to i 
@test leftofline(x, y, xi, yi, xg, yg) == -1 # To the left of line from i to g
J = Float64[-5, -12]
x, y = J
@test leftofline(x, y, xg, yg, xh, yh) == -1 # Point is to the left of all halfplanes, so it is inside
@test leftofline(x, y, xh, yh, xi, yi) == -1 # Point is to the left of all halfplanes, so it is inside
@test leftofline(x, y, xi, yi, xg, yg) == -1 # Point is to the left of all halfplanes, so it is inside

# Test that we can decide whether a point is in a triangle, or on an edge
@test intriangle((2, 3, 4), pts, [3.5, 2.5]) == 1
@test intriangle(p6, p2, p3, [2.0, 3.0]) == 1
x1, y1, x2, y2, x3, y3 = 1 / 10, 1 / 9, 100 / 8, 100 / 3, 100 / 4, 100 / 9
x, y = x1 + (3 / 7) * (x2 - x1), y1 + (3 / 7) * (y2 - y1)
@test intriangle(x1, y1, x2, y2, x3, y3, x, y) == 1
x1, y1, x2, y2, x3, y3 = 2.0, -4.0, 4.0, -4.0, 3.0, -2.0
x, y = 3.0, -4.0
@test intriangle(x1, y1, x2, y2, x3, y3, x, y) == 0
y = y - 1e-12
@test intriangle(x1, y1, x2, y2, x3, y3, x, y) == -1
x, y = 2.5, -3.0
@test intriangle(x1, y1, x2, y2, x3, y3, x, y) == 0
x, y = 3.5, -3.0
@test intriangle(x1, y1, x2, y2, x3, y3, x, y) == 0
@test intriangle((6, 7, 4), pts, [6.0, 6.0]) == -1
@test intriangle((9, 5, 7), pts, [3.0, 3.0]) == -1
@test intriangle((9, 5, 7), pts, [0.0, 3.0]) == 1

## Now let's test what happens for the special points 
@test intriangle((1, 0, -1), pts, [0.0, 0.0]) == 1
p = [0.0, 6.0]
@test leftofline(pts, p, 1, -1) == 1
@test leftofline(pts, p, -1, 1) == -1
@test leftofline(pts, p, -1, 6) == -1
@test intriangle((0, 6, 1), pts, [3.0, 0.0]) == 1

## Test that we can add multiple triangles and edges to the DAG 
𝒟 = TriDAG()
𝒟𝒟 = TriDAG()
T₁ = TriangleType((1, 2, 3))
T₂ = TriangleType((4, 5, 6))
T₃ = TriangleType((7, 8, 9))
add_triangle!(𝒟, T₁, T₂, T₃)
add!(𝒟𝒟, T₁)
add!(𝒟𝒟, T₂)
add!(𝒟𝒟, T₃)
@test graph(𝒟) == graph(𝒟𝒟)
add_triangle!(𝒟, TriangleType((2, 3, 1))) # same triangle as T₁
@test graph(𝒟) == graph(𝒟𝒟)
add_triangle!(𝒟, TriangleType((2, 3, 1)), TriangleType((9, 7, 8))) # same as T₁ and T₃
@test graph(𝒟) == graph(𝒟𝒟)
add_edge!(𝒟, T₁, T₂, T₃)
add!(𝒟𝒟, T₁, T₂)
add!(𝒟𝒟, T₁, T₃)
@test graph(𝒟) == graph(𝒟𝒟)

# partialsort!(pts, 1, rev=true)
𝒯 = TriangleType[(LargeRightIdx + 1, LargeLeftIdx, LargeRightIdx)]
𝒟 = TriDAG()
𝒜 = TriAdjacent()
𝒱𝒩 = TriVertexNeighbours()
add!(𝒟, 𝒯[begin])
root = 𝒯[begin]
@test collect(keys(graph(𝒟).V.dict))[1] == 𝒯[end]
𝒜[(LargeRightIdx + 1, LargeLeftIdx)] = LargeRightIdx
𝒜[(LargeLeftIdx, LargeRightIdx)] = LargeRightIdx + 1
𝒜[(LargeRightIdx, LargeRightIdx + 1)] = LargeLeftIdx
ℬ = TriAdjacent()
update_adjacent!(ℬ, 𝒯[begin])
@test 𝒜.adjacent == ℬ.adjacent
𝒜[(LargeLeftIdx, LargeRightIdx + 1)] = EmptyIdx
𝒜[(LargeRightIdx, LargeLeftIdx)] = EmptyIdx
𝒜[(LargeRightIdx + 1, LargeRightIdx)] = EmptyIdx
add_neighbour!(𝒱𝒩, 1, -1, 0; new_point=true)
add_neighbour!(𝒱𝒩, -1, 0, 1; new_point=true)
add_neighbour!(𝒱𝒩, 0, 1, -1; new_point=true)
T, D, A, VN, _root = initialise_triangulation()
@test T == 𝒯
@test graph(𝒟) == graph(D)
@test 𝒜.adjacent == A.adjacent
@test 𝒱𝒩.neighbours == VN.neighbours
@test root == _root

# @views shuffle!(pts[begin+1:end])

r = 2
pᵣ = pts[r]
Tᵢⱼₖ, flag = locate_triangle(𝒟, pts, pᵣ, root)
_Tᵢⱼₖ, _flag = locate_triangle(𝒟, pts, pᵣ)
@test Tᵢⱼₖ == _Tᵢⱼₖ
@test flag == _flag
# if flag == 1 # on interior 
# Now add the triangles. (i, j, k) becomes (i, j, r), (j, k, r), (k, i, r); these triangles are all positively oriented 
i, j, k = Tᵢⱼₖ
𝒟𝒟 = deepcopy(𝒟)
𝒜𝒜 = deepcopy(𝒜)
𝒯𝒯 = deepcopy(𝒯)
𝒱𝒱 = deepcopy(𝒱𝒩)
push!(𝒯, TriangleType((i, j, r)), TriangleType((j, k, r)), TriangleType((k, i, r)))
deleteat!(𝒯, 1)
add!(𝒟, TriangleType((i, j, r)))
add!(𝒟, TriangleType((j, k, r)))
add!(𝒟, TriangleType((k, i, r)))
add!(𝒟, TriangleType((i, j, k)), TriangleType((i, j, r)))
add!(𝒟, TriangleType((i, j, k)), TriangleType((j, k, r)))
add!(𝒟, TriangleType((i, j, k)), TriangleType((k, i, r)))
update_adjacent!(𝒜, TriangleType((i, j, r)))
update_adjacent!(𝒜, TriangleType((j, k, r)))
update_adjacent!(𝒜, TriangleType((k, i, r)))
add_neighbour!(𝒱𝒩, 1, 2)
add_neighbour!(𝒱𝒩, -1, 2)
add_neighbour!(𝒱𝒩, 0, 2)
add_neighbour!(𝒱𝒩, 2, 1, -1, 0)
add_point!(𝒯𝒯, 𝒟𝒟, 𝒜𝒜, 𝒱𝒱, Tᵢⱼₖ, r)
@test 𝒯 == 𝒯𝒯
@test 𝒟.graph == 𝒟𝒟.graph
@test 𝒜.adjacent == 𝒜𝒜.adjacent
@test 𝒱𝒩.neighbours == 𝒱𝒱.neighbours
@test collect(𝒟.graph.N[(1, -1, 0)]) == TriangleType[(1, -1, 2), (-1, 0, 2), (0, 1, 2)]
@test collect(𝒟.graph.NN[(1, -1, 2)]) == TriangleType[(1, -1, 0)]
@test collect(𝒟.graph.NN[(-1, 0, 2)]) == TriangleType[(1, -1, 0)]
@test collect(𝒟.graph.NN[(0, 1, 2)]) == TriangleType[(1, -1, 0)]
@test !𝒟.graph.looped
@test collect(𝒟.graph.V) == TriangleType[(1, -1, 0), (1, -1, 2), (-1, 0, 2), (0, 1, 2)]
a1, a2, a3 = 𝒯[1]
b1, b2, b3 = 𝒯[2]
c1, c2, c3 = 𝒯[3]
@test TriangleType((1, -1, 0)) ∉ 𝒯
@test (a1, a2, a3) == (1, -1, 2)
@test (b1, b2, b3) == (-1, 0, 2)
@test (c1, c2, c3) == (0, 1, 2)
@test 𝒜(-1, 1) == EmptyIdx
@test 𝒜(0, -1) == EmptyIdx
@test 𝒜(1, 0) == EmptyIdx
@test 𝒜(1, -1) == 2
@test 𝒜(-1, 0) == 2
@test 𝒜(0, 1) == 2
@test 𝒜(-1, 2) == 1
@test 𝒜(2, -1) == 0
@test 𝒜(2, 1) == -1
@test 𝒜(1, 2) == 0
@test 𝒜(2, 0) == 1
@test 𝒜(0, 2) == -1
@test length(𝒜.adjacent) == 12
@test 𝒱𝒩(1) == [-1, 0, 2]
@test 𝒱𝒩(-1) == [0, 1, 2]
@test 𝒱𝒩(0) == [1, -1, 2]
@test 𝒱𝒩(2) == [1, -1, 0]

# Make sure we can detect legal edges
@test edge_on_large_triangle(1, 0)
@test edge_on_large_triangle(0, 1)
@test edge_on_large_triangle(0, -1)
@test edge_on_large_triangle(-1, 0)
@test edge_on_large_triangle(-1, 1)
@test edge_on_large_triangle(1, -1)
@test !edge_on_large_triangle(1, 5)
@test !edge_on_large_triangle(0, 2)
@test !edge_on_large_triangle(-1, 2)
@test !edge_on_large_triangle(0, -2)
@test incircle(pts, 5, 7, 6, 9) == 1
@test incircle(pts, 5, 7, 6, 3) == -1
@test incircle(pts, 5, 7, 6, 3) == -1
@test incircle(pts, 5, 7, 6, 6) == 0
@test incircle(pts, 3, 2, 1, 4) == 1
@test incircle(pts, 3, 2, 1, 6) == 1
@test incircle(pts, 3, 2, 1, 7) == 1
@test incircle(pts, 3, 2, 1, 5) == -1
@test incircle(pts, 3, 2, 1, 8) == -1
@test is_legal(0, -1, 𝒜, pts)
@test is_legal(1, -1, 𝒜, pts)
@test is_legal(1, 0, 𝒜, pts)
@test is_legal(-1, 0, 𝒜, pts)
@test is_legal(-1, 1, 𝒜, pts)
@test is_legal(0, 1, 𝒜, pts)
@test !is_legal(7, 1, 6, 4, pts)
@test is_legal(7, 6, 9, 4, pts)
@test is_legal(7, 3, 4, 5, pts)
@test is_legal(7, 4, 6, 2, pts)
@test is_legal(k, i, 𝒜, pts)
𝒜𝒜 = deepcopy(𝒜)
𝒟𝒟 = deepcopy(𝒟)
𝒯𝒯 = deepcopy(𝒯)
𝒱𝒱 = deepcopy(𝒱𝒩)
legalise_edge!(𝒯𝒯, 𝒟𝒟, 𝒜𝒜, 𝒱𝒱, i, j, r, pts)
legalise_edge!(𝒯𝒯, 𝒟𝒟, 𝒜𝒜, 𝒱𝒱, j, k, r, pts)
legalise_edge!(𝒯𝒯, 𝒟𝒟, 𝒜𝒜, 𝒱𝒱, k, i, r, pts)
@test 𝒯 == 𝒯𝒯
@test 𝒟.graph == 𝒟𝒟.graph
@test 𝒜.adjacent == 𝒜𝒜.adjacent
@test 𝒱𝒩.neighbours == 𝒱𝒱.neighbours

## Add the next point 
r = 3
pᵣ = pts[r]
Tᵢⱼₖ, flag = locate_triangle(𝒟, pts, pᵣ, root)
_Tᵢⱼₖ, _flag = locate_triangle(𝒟, pts, pᵣ)
@test Tᵢⱼₖ == _Tᵢⱼₖ
@test flag == _flag
@test flag == 1
@test Tᵢⱼₖ == TriangleType((LargeLeftIdx, LargeRightIdx, 2))
add_point!(𝒯, 𝒟, 𝒜, 𝒱𝒩, Tᵢⱼₖ, r)
@test length(𝒜.adjacent) == 18
@test 𝒜(-1, 1) == EmptyIdx
@test 𝒜(1, -1) == 2
@test 𝒜(-1, 0) == 3
@test 𝒜(0, -1) == EmptyIdx
@test 𝒜(0, 1) == 2
@test 𝒜(1, 0) == EmptyIdx
@test 𝒜(-1, 2) == 1
@test 𝒜(2, -1) == 3
@test 𝒜(2, 1) == -1
@test 𝒜(1, 2) == 0
@test 𝒜(2, 0) == 1
@test 𝒜(0, 2) == 3
@test 𝒜(2, 3) == 0
@test 𝒜(3, 2) == -1
@test 𝒜(3, -1) == 0
@test 𝒜(-1, 3) == 2
@test 𝒜(3, 0) == 2
@test 𝒜(0, 3) == -1
@test collect(𝒟.graph.N[(1, -1, 0)]) == TriangleType[(1, -1, 2), (-1, 0, 2), (0, 1, 2)]
@test collect(𝒟.graph.N[(-1, 0, 2)]) == TriangleType[(-1, 0, 3), (0, 2, 3), (2, -1, 3)]
@test collect(𝒟.graph.NN[(1, -1, 0)]) == []
@test collect(𝒟.graph.NN[(1, -1, 2)]) == TriangleType[(1, -1, 0)]
@test collect(𝒟.graph.NN[(-1, 0, 2)]) == TriangleType[(1, -1, 0)]
@test collect(𝒟.graph.NN[(-1, 0, 3)]) == TriangleType[(-1, 0, 2)]
@test collect(𝒟.graph.NN[(0, 1, 2)]) == TriangleType[(1, -1, 0)]
@test collect(𝒟.graph.NN[(0, 2, 3)]) == TriangleType[(-1, 0, 2)]
@test collect(𝒟.graph.NN[(2, -1, 3)]) == TriangleType[(-1, 0, 2)]
i, j, k = Tᵢⱼₖ
@test is_legal(i, j, 𝒜, pts)
𝒜𝒜 = deepcopy(𝒜)
𝒟𝒟 = deepcopy(𝒟)
𝒯𝒯 = deepcopy(𝒯)
legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, i, j, r, pts)
@test 𝒯 == 𝒯𝒯
@test 𝒟.graph == 𝒟𝒟.graph
@test 𝒜.adjacent == 𝒜𝒜.adjacent
@test !is_legal(j, k, 𝒜, pts)
@test 𝒜(1, 2) == LargeRightIdx
@test 𝒜(2, 1) == LargeLeftIdx
@test 𝒜(2, 3) == LargeRightIdx
@test 𝒜(3, 2) == LargeLeftIdx
@test 𝒜(LargeLeftIdx, 3) == 2
@test 𝒜(1, LargeLeftIdx) == 2
@test 𝒱𝒩(0) == [1, -1, 2, 3]
@test 𝒱𝒩(1) == [-1, 0, 2]
@test 𝒱𝒩(2) == [1, -1, 0, 3]
@test 𝒱𝒩(3) == [-1, 0, 2]

# Test that we can correctly legalise the edges 
@test is_legal(i, j, 𝒜, pts)
𝒟𝒟 = deepcopy(𝒟)
𝒜𝒜 = deepcopy(𝒜)
𝒯𝒯 = deepcopy(𝒯)
𝒱𝒱 = deepcopy(𝒱𝒩)
legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, i, j, r, pts)
@test 𝒟𝒟.graph == 𝒟.graph
@test 𝒜𝒜.adjacent == 𝒜.adjacent
@test 𝒯𝒯 == 𝒯
@test 𝒱𝒱.neighbours == 𝒱𝒩.neighbours
@test !is_legal(j, k, 𝒜, pts)
_i = 𝒜(k, j)
legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, j, k, r, pts)
@test 𝒟𝒟.graph ≠ 𝒟.graph
@test 𝒜𝒜.adjacent ≠ 𝒜.adjacent
@test 𝒯𝒯 ≠ 𝒯
@test 𝒱𝒩.neighbours ≠ 𝒱𝒱.neighbours
@test is_legal(_i, r, 𝒜, pts)

## Now let's do a clearer run through
p0 = Float64[5, 5]
p1 = Float64[4.5, 2.5]
p2 = Float64[2.5, 1.5]
p3 = Float64[3, 3.5]
p4 = Float64[0, 2]
p5 = Float64[1, 5]
p6 = Float64[1, 3]
p7 = Float64[4, -1]
p8 = Float64[-1, 4]
pts = [p0, p1, p2, p3, p4, p5, p6, p7, p8]

# Initialise 
𝒯, 𝒟, 𝒜, 𝒱𝒩, root = initialise_triangulation()

# Add the second point 
r = 2
pᵣ = pts[r]
𝒯ᵢⱼₖ, interior_flag = locate_triangle(𝒟, pts, pᵣ, root)
i, j, k = 𝒯ᵢⱼₖ
add_point!(𝒯, 𝒟, 𝒜, 𝒱𝒩, 𝒯ᵢⱼₖ, r)
@test is_legal(i, j, 𝒜, pts)
@test is_legal(j, k, 𝒜, pts)
@test is_legal(k, i, 𝒜, pts)
legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, i, j, r, pts)
legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, j, k, r, pts)
legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, k, i, r, pts)
@test 𝒯 == TriangleType[(1, -1, 2), (-1, 0, 2), (0, 1, 2)]

# Add the third point 
r = 3
pᵣ = pts[r]
𝒯ᵢⱼₖ, interior_flag = locate_triangle(𝒟, pts, pᵣ, root)
i, j, k = 𝒯ᵢⱼₖ
add_point!(𝒯, 𝒟, 𝒜, 𝒱𝒩, 𝒯ᵢⱼₖ, r)
@test 𝒯 == TriangleType[(1, -1, 2), (0, 1, 2), (-1, 0, 3), (0, 2, 3), (2, -1, 3)]
@test 𝒜(1, -1) == 2
@test 𝒜(-1, 2) == 1
@test 𝒜(2, 1) == -1
@test 𝒜(0, 1) == 2
@test 𝒜(1, 2) == 0
@test 𝒜(2, 0) == 1
@test 𝒜(-1, 0) == 3
@test 𝒜(0, 3) == -1
@test 𝒜(3, -1) == 0
@test 𝒜(0, 2) == 3
@test 𝒜(2, 3) == 0
@test 𝒜(3, 0) == 2
@test 𝒜(2, -1) == 3
@test 𝒜(-1, 3) == 2
@test 𝒜(3, 2) == -1
@test sort(𝒱𝒩(-1)) == [0, 1, 2, 3]
@test sort(𝒱𝒩(0)) == [-1, 1, 2, 3]
@test sort(𝒱𝒩(1)) == [-1, 0, 2]
@test sort(𝒱𝒩(2)) == [-1, 0, 1, 3]
@test sort(𝒱𝒩(3)) == [-1, 0, 2]
@test is_legal(i, j, 𝒜, pts)
legalise_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, i, j, r, pts)
@test 𝒯 == TriangleType[(1, -1, 2), (0, 1, 2), (-1, 0, 3), (0, 2, 3), (2, -1, 3)]
@test !is_legal(j, k, 𝒜, pts)
@test 𝒜(j, k) == r
@test 𝒜(k, j) == 1
_i, _j, _k = j, k, 1
flip_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, _i, _j, _k, r)
@test 𝒯 == TriangleType[(1, -1, 2), (-1, 0, 3), (2, -1, 3), (3, 1, 2), (3, 0, 1)]
@test is_legal(_i, _k, 𝒜, pts)
@test is_legal(_k, _j, 𝒜, pts)
@test is_legal(j, 1, 𝒜, pts)
@test is_legal(1, k, 𝒜, pts)
@test 𝒜(1, -1) == 2
@test 𝒜(-1, 2) == 1
@test 𝒜(2, 1) == -1
@test 𝒜(-1, 0) == 3
@test 𝒜(0, 3) == -1
@test 𝒜(3, -1) == 0
@test 𝒜(2, -1) == 3
@test 𝒜(-1, 3) == 2
@test 𝒜(3, 2) == -1
@test 𝒜(3, 1) == 2
@test 𝒜(1, 2) == 3
@test 𝒜(2, 3) == 1
@test 𝒜(3, 0) == 1
@test 𝒜(0, 1) == 3
@test 𝒜(1, 3) == 0

# Add the fourth point 
r = 4
pᵣ = pts[r]
𝒯ᵢⱼₖ, interior_flag = locate_triangle(𝒟, pts, pᵣ, root)
@test 𝒯ᵢⱼₖ == TriangleType((2, -1, 3))
@test interior_flag == 1
i, j, k = 𝒯ᵢⱼₖ
add_point!(𝒯, 𝒟, 𝒜, 𝒱𝒩, 𝒯ᵢⱼₖ, r)
@test 𝒯 == TriangleType[(1, -1, 2), (-1, 0, 3), (3, 1, 2), (3, 0, 1), (2, -1, 4), (-1, 3, 4), (3, 2, 4)]
@test 𝒜(1, -1) == 2
@test 𝒜(-1, 2) == 1
@test 𝒜(2, 1) == -1
@test 𝒜(-1, 0) == 3
@test 𝒜(0, 3) == -1
@test 𝒜(3, -1) == 0
@test 𝒜(3, 1) == 2
@test 𝒜(1, 2) == 3
@test 𝒜(2, 3) == 1
@test 𝒜(3, 0) == 1
@test 𝒜(0, 1) == 3
@test 𝒜(1, 3) == 0
@test 𝒜(2, -1) == 4
@test 𝒜(-1, 4) == 2
@test 𝒜(4, 2) == -1
@test 𝒜(-1, 3) == 4
@test 𝒜(3, 4) == -1
@test 𝒜(4, -1) == 3
@test 𝒜(3, 2) == 4
@test 𝒜(2, 4) == 3
@test 𝒜(4, 3) == 2
@test 𝒜(-1, 1) == EmptyIdx
@test 𝒜(0, -1) == EmptyIdx
@test 𝒜(1, 0) == EmptyIdx
@test !is_legal(i, j, 𝒜, pts)
@test 𝒜(j, i) == 1
_i, _j, _k = i, j, 1
flip_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, _i, _j, _k, r)
@test 𝒯 == TriangleType[(-1, 0, 3), (3, 1, 2), (3, 0, 1), (-1, 3, 4), (3, 2, 4), (4, 1, -1), (4, 2, 1)]
@test 𝒜(1, -1) == 4
@test (-1, 2) ∉ keys(𝒜.adjacent)
@test 𝒜(-1, 0) == 3
@test 𝒜(0, 3) == -1
@test 𝒜(3, -1) == 0
@test 𝒜(3, 1) == 2
@test 𝒜(2, 3) == 1
@test 𝒜(3, 0) == 1
@test 𝒜(0, 1) == 3
@test 𝒜(1, 3) == 0
@test (2, -1) ∉ keys(𝒜.adjacent)
@test 𝒜(-1, 4) == 1
@test 𝒜(4, 2) == 1
@test 𝒜(-1, 3) == 4
@test 𝒜(3, 4) == -1
@test 𝒜(4, -1) == 3
@test 𝒜(3, 2) == 4
@test 𝒜(2, 4) == 3
@test 𝒜(4, 3) == 2
@test 𝒜(-1, 1) == EmptyIdx
@test 𝒜(0, -1) == EmptyIdx
@test 𝒜(1, 0) == EmptyIdx
@test 𝒜(4, 1) == -1
@test 𝒜(1, 4) == 2
@test (3, 1) ∈ keys(𝒜.adjacent)
@test is_legal(_i, _k, 𝒜, pts)
@test is_legal(_k, _j, 𝒜, pts)
@test is_legal(4, 1, 𝒜, pts)
@test !is_legal(j, k, 𝒜, pts)
_i, _j, _k = j, k, 0
flip_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, _i, _j, _k, r)
@test 𝒯 == TriangleType[(3, 1, 2), (3, 0, 1), (3, 2, 4), (4, 1, -1), (4, 2, 1), (4, 0, 3), (4, -1, 0)]
@test 𝒜(3, 1) == 2
@test 𝒜(1, 2) == 3
@test 𝒜(2, 3) == 1
@test 𝒜(3, 0) == 1
@test 𝒜(0, 1) == 3
@test 𝒜(1, 3) == 0
@test 𝒜(3, 2) == 4
@test 𝒜(2, 4) == 3
@test 𝒜(4, 3) == 2
@test 𝒜(4, 1) == -1
@test 𝒜(1, -1) == 4
@test 𝒜(-1, 4) == 1
@test 𝒜(4, 2) == 1
@test 𝒜(2, 1) == 4
@test 𝒜(1, 4) == 2
@test 𝒜(4, 0) == 3
@test 𝒜(0, 3) == 4
@test 𝒜(3, 4) == 0
@test 𝒜(4, -1) == 0
@test 𝒜(-1, 0) == 4
@test 𝒜(0, 4) == -1
@test (3, 1) ∈ keys(𝒜.adjacent)
@test is_legal(_i, _k, 𝒜, pts)
@test !is_legal(_k, _j, 𝒜, pts)
_i, _j = _k, _j
_k = 𝒜(_j, _i)
flip_edge!(𝒯, 𝒟, 𝒜, 𝒱𝒩, _i, _j, _k, r)
@test 𝒯 == TriangleType[(3, 1, 2), (3, 2, 4), (4, 1, -1), (4, 2, 1), (4, -1, 0), (4, 1, 3), (4, 0, 1)]
@test 𝒜(-1, 1) == EmptyIdx
@test 𝒜(0, -1) == EmptyIdx
@test 𝒜(1, 0) == EmptyIdx
@test 𝒜(3, 1) == 2
@test 𝒜(1, 2) == 3
@test 𝒜(2, 3) == 1
@test 𝒜(3, 2) == 4
@test 𝒜(2, 4) == 3
@test 𝒜(4, 3) == 2
@test 𝒜(4, 1) == 3
@test 𝒜(1, -1) == 4
@test 𝒜(-1, 4) == 1
@test 𝒜(4, 2) == 1
@test 𝒜(2, 1) == 4
@test 𝒜(1, 4) == 0
@test 𝒜(4, -1) == 0
@test 𝒜(-1, 0) == 4
@test 𝒜(0, 4) == -1
@test 𝒜(4, 1) == 3
@test 𝒜(1, 3) == 4
@test 𝒜(3, 4) == 1
@test 𝒜(4, 0) == 1
@test 𝒜(0, 1) == 4
@test 𝒜(1, 4) == 0
@test is_legal(4, 0, 𝒜, pts)
@test is_legal(k, i, 𝒜, pts)

## Do a smaller example 
p1 = Float64[20, 20]
p2 = Float64[0, 6]
p3 = Float64[12, -2]
p4 = Float64[10, 10]
PTS = [p1, p2, p3, p4]
𝒯, 𝒟, 𝒜, 𝒱𝒩 = triangulate(PTS; shuffle_pts=false, trim=false)
@test 𝒯 == TriangleType[
    (1, -1, 2),
    (-1, 0, 2),
    (0, 1, 3),
    (2, 0, 3),
    (1, 2, 4),
    (2, 3, 4),
    (3, 1, 4)
]
@test 𝒜(1, -1) == 2
@test 𝒜(-1, 2) == 1
@test 𝒜(2, 1) == -1
@test 𝒜(-1, 0) == 2
@test 𝒜(0, 2) == -1
@test 𝒜(2, -1) == 0
@test 𝒜(0, 1) == 3
@test 𝒜(1, 3) == 0
@test 𝒜(3, 0) == 1
@test 𝒜(2, 0) == 3
@test 𝒜(0, 3) == 2
@test 𝒜(3, 2) == 0
@test 𝒜(1, 2) == 4
@test 𝒜(2, 4) == 1
@test 𝒜(4, 1) == 2
@test 𝒜(2, 3) == 4
@test 𝒜(3, 4) == 2
@test 𝒜(4, 2) == 3
@test 𝒜(3, 1) == 4
@test 𝒜(1, 4) == 3
@test 𝒜(4, 3) == 1
@test sort(𝒱𝒩(-1)) == [0, 1, 2]
@test sort(𝒱𝒩(0)) == [-1, 1, 2, 3]
@test sort(𝒱𝒩(1)) == [-1, 0, 2, 3, 4]
@test sort(𝒱𝒩(2)) == [-1, 0, 1, 3, 4]
@test sort(𝒱𝒩(3)) == [0, 1, 2, 4]
@test sort(𝒱𝒩(4)) == [1, 2, 3]

w₁ = 𝒜(1, -1)
delete_triangle!(𝒯, TriangleType((1, -1, w₁)))
delete_edge_from_adjacency!(𝒜, 1, -1; protect_boundary=false)
delete_edge_from_adjacency!(𝒜, -1, w; protect_boundary=false)
𝒜[(w₁, 1)] = EmptyIdx

w₂ = 𝒜(-1, 0)
delete_triangle!(𝒯, TriangleType((-1, 0, w₂)))
delete_edge_from_adjacency!(𝒜, -1, 0; protect_boundary=false)
delete_edge_from_adjacency!(𝒜, 0, w₂; protect_boundary=false)

w₃ = 𝒜(0, 1)
delete_triangle!(𝒯, TriangleType((0, 1, w₃)))
delete_edge_from_adjacency!(𝒜, 0, 1; protect_boundary=false)
delete_edge_from_adjacency!(𝒜, 0, w₃; protect_boundary=false)
𝒜[(w₃, w₂)] = EmptyIdx
𝒜[(1, w₃)] = EmptyIdx
𝒯
@test 𝒯 == [(2,)]