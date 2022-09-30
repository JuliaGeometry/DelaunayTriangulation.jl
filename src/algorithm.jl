############################################
##
## ABSTRACT PRIMITIVES 
##
############################################
"""
    abstract type AbstractEdge{I} end

Abstract type representing an edge with indices of type `I`. 
Typically an `NTuple{2, Int64}`.
"""
abstract type AbstractEdge{I} end

"""
    abstract type AbstractTriangle{I} end

Abstract type representing a triangle with indices of type `I`. 
Typically an `NTuple{3, Int64}`.

`AbstractTriangle`s must have methods for:
- `geti(T)`: The first vertex.
- `getj(T)`: The second vertex. 
- `getk(T)`: The third vertex. 
- `indices(T)`: The indices for all three vertices.
"""
abstract type AbstractTriangle{I} end
shift_triangle_indices_0(T::AbstractTriangle) = (geti(T), getj(T), getk(T))
shift_triangle_indices_1(T::AbstractTriangle) = (getj(T), getk(T), geti(T))
shift_triangle_indices_2(T::AbstractTriangle) = (getk(T), geti(T), getj(T))
shift_triangle_0(T::V) where {I,V<:AbstractTriangle{I}} = V(shift_triangle_indices_0(T))
shift_triangle_1(T::V) where {I,V<:AbstractTriangle{I}} = V(shift_triangle_indices_1(T))
shift_triangle_2(T::V) where {I,V<:AbstractTriangle{I}} = V(shift_triangle_indices_2(T))
shift_triangle_indices(T, m) = m == 0 ? shift_triangle_indices_0(T) : (m == 1 ? shift_triangle_indices_1(T) : shift_triangle_indices_2(T))
shift_triangle(T::V, m) where {I,V<:AbstractTriangle{I}} = V(shift_triangle_indices(T, m))

"""
    abstract type AbstractPoint{T, A} end 

Abstract type representing a point with element type `T` stored in a container of 
type `A`. Typically a `Vector{Float64}` or `SVector{2, Float64}`.

`AbstractPoint`s must have methods for:
- `coords(p)`: The coordinates of the point.
- `getx(p)`: Get the `x`-coordinate of the point.
- `gety(p)`: Get the `y`-coordinate of the point.
"""
abstract type AbstractPoint{T,A} end
Base.Tuple(p::AbstractPoint{T,A}) where {T,A} = (getx(p), gety(p))

############################################
##
## SPECIFIC DATA TYPES 
##
############################################
"""
    struct Edge{I} <: AbstractEdge{I}

A data structure for an edge, with fields

- `initial::I`: The initial point of the edge. Can be accessed via `initial(e::Edge)`.
- `terminal::I`: The terminal point of the edge. Can be accessed via `terminal(e::Edge)`.
"""
struct Edge{I} <: AbstractEdge{I}
    initial::I
    terminal::I
    Edge(i::I, j::I) where {I} = new{I}(i, j)
    Edge(ij::NTuple{2,I}) where {I} = new{I}(ij[1], ij[2])
    Edge{I}(i, j) where {I} = new{I}(convert(I, i), convert(I, j))
end
initial(e::Edge) = e.initial
terminal(e::Edge) = e.terminal

"""
    struct Triangle{I} <: AbstractTriangle{I}

A triangle data structure, with field

- `indices::NTuple{3, I}`: The indices of the triangle's vertices, accessed via `indices`.

The individual indices can be accessed via `geti`, `getj`, and `getk`.
"""
struct Triangle{I} <: AbstractTriangle{I}
    indices::NTuple{3,I}
    Triangle(indices::NTuple{3,I}) where {I} = new{I}(indices)
    Triangle(i::I, j::I, k::I) where {I} = new{I}((i, j, k))
    Triangle{I}(indices) where {I} = new{I}(indices)
end
indices(T::Triangle) = T.indices
geti(T::Triangle) = T.indices[1]
getj(T::Triangle) = T.indices[2]
getk(T::Triangle) = T.indices[3]

"""
    struct Point{T,A} <: AbstractPoint{T,A}

A data structure for a point, with field

- `coords::A`: The coordinates of the point, accessed via `coords`.

The individual coordinates can be expressed via `getx` and `gety`.
"""
struct Point{T,A} <: AbstractPoint{T,A}
    coords::A
    Point(x::T, y::T) where {T} = new{T,Vector{T}}(T[x, y])
    Point{A}(x, y) where {A} = new{eltype(A),A}(convert(A, eltype(A)[x, y]))
    Point(coords::A) where {A} = new{eltype(coords),A}(coords)
end
coords(p::Point) = p.coords
getx(p::Point) = first(coords(p))
gety(p::Point) = last(coords(p))
Base.:(==)(p::Point, q::Point) = isequal(getx(p), getx(q)) && isequal(gety(p), gety(q))

############################################
##
## CONSTANTS 
##
############################################
"""The index used to represent the lower-right vertex of the bounding triangle, p₋₁."""
const LowerRightBoundingIndex = -1
"""The index used to represent the upper vertex of the bounding triangle, p₋₂."""
const UpperBoundingIndex = -2
"""The index used to represent the lower-left of the bounding triangle, p₋₃."""
const LowerLeftBoundingIndex = -3
"""The triangle representing the bounding triangle for any triangulation."""
const BoundingTriangle = Triangle(LowerRightBoundingIndex,
    UpperBoundingIndex,
    LowerLeftBoundingIndex)
"""The index used to represent the unbounded part of the triangulation, i.e. the boundary."""
const BoundaryIndex = 0
"""The first index where physical points begin."""
const FirstPointIndex = 1

############################################
##
## PRIMITIVE COLLECTIONS 
##
############################################
"""
    struct Triangles{I,T<:AbstractTriangle{I}}

Structure for a collection of `Triangle`s, currently implemented using a `Set`. 
The triangles can be accessed via [`triangles`](@ref).
"""
struct Triangles{I,T<:AbstractTriangle{I}}
    triangles::Set{T}
    Triangles(triangles::Set{T}) where {I,T<:AbstractTriangle{I}} = new{I,T}(triangles)
    Triangles(triangles::Base.AbstractVecOrTuple{T}) where {I,T<:AbstractTriangle{I}} = new{I,T}(Set{T}(triangles))
end
triangles(tris::Triangles) = tris.triangles

"""
    add_triangle!(T::Triangles, V::AbstractTriangle...)

Adds the triangles in `V` into the collection of triangles `T`. No 
checks are made for triangles that are equal under circular shifts.
"""
function add_triangle!(T::Triangles, V::AbstractTriangle...)
    push!(triangles(T), V...)
    return nothing
end

"""
    delete_triangle!(T::Triangles, V::AbstractTriangle...)

Deletes the triangle `V` from the collection of triangles `T`. All circular 
shifts of `V`'s indices are also deleted.
"""
function delete_triangle!(T::Triangles, V::AbstractTriangle)
    delete!(triangles(T), V)
    delete!(triangles(T), shift_triangle_1(V))
    delete!(triangles(T), shift_triangle_2(V))
    return nothing
end
@doc (@doc delete_triangle!(::Triangles, ::AbstractTriangle))
function delete_triangle!(T::Triangles, V::AbstractTriangle...)
    for V in V
        delete_triangle!(T, V)
    end
    return nothing
end

"""
    Points{T,A,P<:AbstractPoint{T,A}}

Structure for a collection of `Point`s, currently implemented using a `Vector`.
The points can be accessed using [`points`](@ref).
"""
struct Points{T,A,P<:AbstractPoint{T,A}}
    points::Vector{P}
    Points(points::AbstractVector{P}) where {T,A,P<:AbstractPoint{T,A}} = new{T,A,P}(points)
    Points(points::NTuple{N,P}) where {N,T,A,P<:AbstractPoint{T,A}} = new{T,A,P}(collect(points))
end
points(pts::Points) = pts.points
Base.length(pts::Points) = length(points(pts))
function Points(pts::AbstractVector)
    return Points(Point.(pts))
end

function get_point(pts::Points, i)
    return points(pts)[i]
end

"""
    add_point!(pts::Points, p...)
    
Adds the points in `p` to the collection of points `pts`.
"""
function add_point!(pts::Points, p...)
    push!(points(pts), p...)
    return nothing
end

############################################
##
## CONSTRUCTORS FOR THE TRIANGULATION DATA STRUCTURES 
##
############################################
"""
    struct Adjacent{I,E<:AbstractEdge{I}}

Data structure for the `adjacent` map of a triangulation. The map can be accessed 
via [`adjacent`](@ref). The map is defined using a `Dict{E, I}`, such that the key 
`E(u, v)` is mapped to an `w::I` such that `Triangle(u, v, w)` is positively oriented.
"""
struct Adjacent{I,E<:AbstractEdge{I}}
    adjacent::Dict{E,I}
    function Adjacent{I,E}() where {I,E}
        A = Dict{E,I}()
        adj = new{I,E}(A)
        return adj
    end
    Adjacent() = Adjacent{Int64,Edge{Int64}}()
end
adjacent(adj::Adjacent) = adj.adjacent

"""
    struct Adjacent2Vertex{I,E<:AbstractEdge{I}}

Data structure for the `adjacent2vertex` map of a triangulation. The map can be accessed 
via [`adjacent2vertex`](@ref). The map is defined using a `Dict{I,Set{E}}`, such that the key 
`w` is mapped to a `Set` containing all pairs `(u, v)` such that `Triangle(u, v, w)` is a 
positively oriented triangle.
"""
struct Adjacent2Vertex{I,E<:AbstractEdge{I}}
    adjacent2vertex::Dict{I,Set{E}}
    function Adjacent2Vertex{I,E}() where {I,E}
        D = Dict{I,Set{E}}()
        TA2V = new{I,E}(D)
        return TA2V
    end
    Adjacent2Vertex() = Adjacent2Vertex{Int64,Edge{Int64}}()
end
adjacent2vertex(adj2v::Adjacent2Vertex) = adj2v.adjacent2vertex

"""
    struct DelaunayGraph{I}

Graph representation of a Delaunay triangulation, mapping all points to their neighbours, 
where a point `v` is a neighbour of `u` is `(u, v)` is an edge in the triangulation. The 
implementation is through an `UndirectedGraph`. The graph can be accessed using [`graph`](@ref).
"""
struct DelaunayGraph{I}
    graph::UndirectedGraph{I}
    function DelaunayGraph{I}() where {I}
        DG = UndirectedGraph{I}()
        return new{I}(DG)
    end
    DelaunayGraph() = DelaunayGraph{Int64}()
    DelaunayGraph(DG::UndirectedGraph{I}) where {I} = new{I}(DG)
end
graph(G::DelaunayGraph) = G.graph

"""
    struct HistoryDAG{I,T<:AbstractTriangle{I}} 

The point location graph for the Delaunay triangulation. This is a directed acyclic graph, 
implemented using `DirectedGraph`, that stores the history of the triangulation as points are 
inserted. The graph can be accessed using [`graph`](@ref).
"""
struct HistoryDAG{I,T<:AbstractTriangle{I}}
    graph::DirectedGraph{T}
    function HistoryDAG{I,T}() where {I,T}
        G = DirectedGraph{T}()
        forbid_loops!(G)
        TDAG = new{I,T}(G)
        return TDAG
    end
    HistoryDAG() = HistoryDAG{Int64,Triangle{Int64}}()
    function HistoryDAG(HG::DirectedGraph{T}) where {I,T<:AbstractTriangle{I}}
        forbid_loops!(HG)
        return new{I,T}(HG)
    end
end
graph(G::HistoryDAG) = G.graph

struct Triangulation{A,A2V,DG,H,T,P,R}
    adjacent::A
    adjacent2vertex::A2V
    graph::DG
    history::H
    triangles::T
    points::P
    root::R
end
adjacent(DT::Triangulation) = DT.adjacent
adjacent2vertex(DT::Triangulation) = DT.adjacent2vertex
graph(DT::Triangulation) = DT.graph
history(DT::Triangulation) = DT.history
triangles(DT::Triangulation) = DT.triangles
points(DT::Triangulation) = DT.points
root(DT::Triangulation) = DT.root

############################################
##
## INDEXING AND UPDATING OF THE TRIANGULATION DATA STRUCTURES
##
############################################
## Adjacent  
""" 
    add_edge!(adj::Adjacent, uv::AbstractEdge, w)
    add_edge!(adj::Adjacent{I,E}, u, v, w) where {I,E<:AbstractEdge}

Adds the edge `E(u, v)` to the adjacent map `adj`.
"""
function add_edge!(adj::Adjacent, uv::AbstractEdge, w)
    adjacent(adj)[uv] = w
    return nothing
end
@doc (@doc add_edge!(::Adjacent::AbstractEdge, ::Any))
function add_edge!(adj::Adjacent{I,E}, u, v, w) where {I,E<:AbstractEdge}
    add_edge!(adj, E(u, v), w)
    return nothing
end

"""
    get_edge(adj::Adjacent, uv::AbstractEdge)
    get_edge(adj::Adjacent{I,E}, u, v) where {I,E<:AbstractEdge}

Accessor for the key `E(u, v)` in the adjacent map `adj`, so that `Triangle(u, v, get_edge(adj, u, v))` 
is positively oriented.
"""
function get_edge(adj::Adjacent, uv::AbstractEdge)
    return adjacent(adj)[uv]
end
function get_edge(adj::Adjacent{I,E}, u, v) where {I,E<:AbstractEdge}
    return get_edge(adj, E(u, v))
end

"""
    delete_edge!(adj::Adjacent{I,E}, i, j; protect_boundary=true) where {I,E<:AbstractEdge{I}}

Deletes the edge `E(i, j)` from the map `adj`. This function also deletes `E(j, i)` from `adj`. If 
the boundary needs to be protected, so that `E(i, j)` is not deleted if `E(i, j) = BoundaryIndex`
(and similarly for `E(j, i)`), set `protect_boundary=true`. 
"""
function delete_edge!(adj::Adjacent{I,E}, i, j; protect_boundary=true) where {I,E<:AbstractEdge{I}}
    (!protect_boundary || get_edge(adj, i, j) ≠ BoundaryIndex) && delete!(adjacent(adj), E(i, j))
    (!protect_boundary || get_edge(adj, j, i) ≠ BoundaryIndex) && delete!(adjacent(adj), E(j, i))
    return nothing
end

"""
    add_triangle!(adj::Adjacent, T::AbstractTriangle...)

Adds the edges of the triangle `T` to the adjacent map `adj`.
"""
function add_triangle!(adj::Adjacent, T::AbstractTriangle)
    i, j, k = indices(T)
    add_edge!(adj, i, j, k)
    add_edge!(adj, j, k, i)
    add_edge!(adj, k, i, j)
    return nothing
end
@doc (@doc add_triangle!(::Adjacent, ::AbstractTriangle))
function add_triangle!(adj::Adjacent, T::AbstractTriangle...)
    for T in T
        add_triangle!(adj, T)
    end
    return nothing
end

## Adjacent2Vertex
"""
    add_edge!(adj2v::Adjacent2Vertex{I,E}, w, uv::E) where {I,E<:AbstractEdge{I}}
    add_edge!(adj2v::Adjacent2Vertex{I,E}, w, u, v) where {I,E<:AbstractEdge{I}}

Adds the edge `E(u, v)` to the set of edges in `adjacent2vertex(adj2v)[w]`.
"""
function add_edge!(adj2v::Adjacent2Vertex{I,E}, w, uv::E) where {I,E<:AbstractEdge{I}}
    existing_segments = get!(Set{E}, adjacent2vertex(adj2v), w)
    push!(existing_segments, uv)
    return nothing
end
@doc (@doc add_edge!(::Adjacent2Vertex{I,E}, ::Any, ::E) where {I,E<:AbstractEdge{I}})
function add_edge!(adj2v::Adjacent2Vertex{I,E}, w, u, v) where {I,E<:AbstractEdge{I}}
    add_edge!(adj2v, w, E(u, v))
    return nothing
end

"""
    get_edge(adj2v::Adjacent2Vertex, w)

Accessor for the set of edges `(u, v)` such that `Triangle(u, v, w)` is positively oriented.
"""
function get_edge(adj2v::Adjacent2Vertex, w)
    return adjacent2vertex(adj2v)[w]
end

"""
    delete_edge!(adj2v::Adjacent2Vertex{I,E}, w, uv::E) where {I,E<:AbstractEdge{I}}
    delete_edge!(adj2v::Adjacent2Vertex{I,E}, w, u, v) where {I,E<:AbstractEdge{I}}

Deletes the edge `E(u, v)` from the set of edges around `w` that together define a positively 
oriented triangle.
"""
function delete_edge!(adj2v::Adjacent2Vertex{I,E}, w, uv::E) where {I,E<:AbstractEdge{I}}
    delete!(get_edge(adj2v, w), uv)
    return nothing
end
@doc (@doc delete_edge!(::Adjacent2Vertex{I,E}, ::Any, ::AbstractEdge) where {I,E<:AbstractEdge{I}})
function delete_edge!(adj2v::Adjacent2Vertex{I,E}, w, u, v) where {I,E<:AbstractEdge{I}}
    delete_edge!(adj2v, w, E(u, v))
    return nothing
end

"""
    delete_point!(adj2v::Adjacent2Vertex, w) 

Deletes the point `w` from the map `adj2v`.
"""
function delete_point!(adj2v::Adjacent2Vertex, w)
    delete!(adjacent2vertex(adj2v), w)
    return nothing
end

"""
    update_after_flip!(adj2v::Adjacent2Vertex, i, j, k, r)

Updates the map `adj2v` after the edge `(i, j)` is flipped to `(k, r)`, where 
the triangles `(i, j, r)` and `(k, j, i)` are positively oriented.
"""
function update_after_flip!(adj2v::Adjacent2Vertex, i, j, k, r)
    # Delete the necessary edges
    delete_edge!(adj2v, i, k, j)
    delete_edge!(adj2v, i, j, r)
    delete_edge!(adj2v, j, r, i)
    delete_edge!(adj2v, j, i, k)
    delete_edge!(adj2v, k, j, i)
    delete_edge!(adj2v, r, i, j)
    # Add the new edges 
    add_edge!(adj2v, i, k, r)
    add_edge!(adj2v, j, r, k)
    add_edge!(adj2v, k, j, r)
    add_edge!(adj2v, k, r, i)
    add_edge!(adj2v, r, k, j)
    add_edge!(adj2v, r, i, k)
    return nothing
end

"""
    update_after_insertion!(adj2v::Adjacent2Vertex, i, j, k, r)

Updates the map `adj2v` after the point `r` is put into the positively 
oriented triangle `(i, j, k)`, thus subdividing the triangle into the 
new triangles `(i, j, r)`, `(j, k, r)`, and `(k, i, r)`.`
"""
function update_after_insertion!(adj2v::Adjacent2Vertex, i, j, k, r)
    # Delete old edges 
    delete_edge!(adj2v, i, j, k)
    delete_edge!(adj2v, j, k, i)
    delete_edge!(adj2v, k, i, j)
    # Add the new edges 
    add_edge!(adj2v, i, j, r)
    add_edge!(adj2v, i, r, k)
    add_edge!(adj2v, j, k, r)
    add_edge!(adj2v, j, r, i)
    add_edge!(adj2v, k, i, r)
    add_edge!(adj2v, k, r, j)
    add_edge!(adj2v, r, i, j)
    add_edge!(adj2v, r, k, i)
    add_edge!(adj2v, r, j, k)
    return nothing
end

## DelaunayGraph 
"""
    add_point!(DG::DelaunayGraph, u...)

Adds the point(s) `u` into the vertex set of `DG`.
"""
function add_point!(DG::DelaunayGraph, u)
    add!(graph(DG), u)
    return nothing
end
@doc (@doc add_point!(::DelaunayGraph, ::Any))
function add_point!(DG::DelaunayGraph, u...)
    for v in u
        add_point!(DG, v)
    end
    return nothing
end

"""
    add_neighbour!(DG::DelaunayGraph{I}, u::I, v::I...) where {I}

Adds an edge from `u` to the points in `v` into the graph `DG`. Note that, as 
`DG` is an undirected graph, this also adds an edge from `v` to `u`.
"""
function add_neighbour!(DG::DelaunayGraph{I}, u::I, v::I) where {I}
    add!(graph(DG), u, v)
    return nothing
end
@doc (@doc add_neighbour!(::DelaunayGraph{I}, ::I, ::I) where {I})
function add_neighbour!(DG::DelaunayGraph{I}, u::I, v::I...) where {I}
    for w in v
        add_neighbour!(DG, u, w)
    end
    return nothing
end

"""
    get_neighbour(DG::DelaunayGraph, u)

Extracts the neighbours of `u` from the graph `DG`.
"""
function get_neighbour(DG::DelaunayGraph, u)
    return graph(DG).N[u]
end

"""
    delete_neighbour!(DG::DelaunayGraph, u, v)

Deletes `v` from the neighbourhood of `u`. Note that, as `DG` is an 
undirected graph, this also deletes `u` from the neighbourhood of `v`.
"""
function delete_neighbour!(DG::DelaunayGraph, u, v)
    delete!(graph(DG), u, v)
    return nothing
end

"""
    delete_point!(DG::DelaunayGraph, u...)

Deletes the point(s) `u` from the vertex set of `DG`.
"""
function delete_point!(DG::DelaunayGraph, u)
    delete!(graph(DG), u)
    return nothing
end
@doc (@doc delete_point!(::DelaunayGraph, ::Any))
function delete_point!(DG::DelaunayGraph, u...)
    for v in u
        delete_point!(DG, v)
    end
    return nothing
end

## HistoryDAG 
"""
    SimpleGraphs.out_neighbors(G::HistoryDAG, T)

Returns the set of triangles that `T` was replaced by in the triangulation.
"""
SimpleGraphs.out_neighbors(G::HistoryDAG, T) = graph(G).N[T] # The version in SimpleGraphs is allocating...
SimpleGraphs.in_neighbors(G::HistoryDAG, T) = SimpleGraphs.in_neighbors(graph(G), T) # Need the allocating version for find_root. Only place it's used anyway.

SimpleGraphs.in_deg(G::HistoryDAG, T) = SimpleGraphs.in_deg(graph(G), T)
SimpleGraphs.out_deg(G::HistoryDAG, T) = SimpleGraphs.out_deg(graph(G), T)

"""
    add_triangle!(D::HistoryDAG, T::AbstractTriangle...)

Adds the triangles in `T` into the graph `D`. Checks are made for 
circular shifts of the vertices of `T` to avoid duplicates.
"""
function add_triangle!(D::HistoryDAG, T::AbstractTriangle)
    has(graph(D), T) && return nothing
    has(graph(D), shift_triangle_1(T)) && return nothing
    has(graph(D), shift_triangle_2(T)) && return nothing
    add!(graph(D), T)
    return nothing
end
@doc (@doc add_triangle!(::HistoryDAG, ::AbstractTriangle))
function add_triangle!(D::HistoryDAG, T::AbstractTriangle...)
    for V in T
        add_triangle!(D, V)
    end
end

"""
    add_edge!(G::HistoryDAG, T::AbstractTriangle, V::AbstractTriangle...)

Adds a directed edge from `T` to the triangles in `V`, making checks for circular shifts 
in the vertices of `T` and `V`.
"""
function add_edge!(G::HistoryDAG, T::AbstractTriangle, V::AbstractTriangle)
    for i in 0:2
        Tshift = shift_triangle(T, i)
        if has(graph(G), Tshift)
            for j in 0:2
                Vshift = shift_triangle(V, j)
                if has(graph(G), Vshift)
                    SimpleGraphs.add!(graph(G), Tshift, Vshift)
                    return nothing
                end
            end
        end
    end
    return nothing
end
@doc (@doc add_edge!(::HistoryDAG, ::AbstractTriangle, ::AbstractTriangle))
function add_edge!(G::HistoryDAG, T::AbstractTriangle, V::AbstractTriangle...)
    for V in V
        add_edge!(G, T, V)
    end
    return nothing
end

############################################
##
## PREDICATES 
##
############################################
"""
    ExactPredicates.orient(T::AbstractTriangle, pts::Points)

Tests if the triangle `T` is positively oriented.
"""
function ExactPredicates.orient(T::AbstractTriangle, pts::Points)
    u, v, w = indices(T)
    pu, pv, pw = get_point(pts, u), get_point(pts, v), get_point(pts, w)
    return orient(pu, pv, pw)
end

"""
    ExactPredicates.incircle(pts::Points, i, j, k, ℓ)

Tests if `get_point(pts, ℓ)` is in the circle through `(get_point(pts, i), get_point(pts, j), get_point(pts, k))`.
"""
function ExactPredicates.incircle(pts::Points, i, j, k, ℓ)
    pti, ptj, ptk, ptℓ = get_point(pts, i), get_point(pts, j), get_point(pts, k), get_point(pts, ℓ)
    return incircle(pti, ptj, ptk, ptℓ)
end

"""
    leftofline(p, pᵢ, pⱼ)

Tests if `p` is to the left of the line from `pᵢ` to `pⱼ`.
"""
leftofline(p::AbstractPoint, pᵢ::AbstractPoint, pⱼ::AbstractPoint) = orient(pᵢ, pⱼ, p)
"""
    leftofline(pts, p, i, j)

Tests if the point `p` is to the left of the oriented line through 
`pts[i]` to `pts[j]`. Checks are made for non-positive indices.
"""
function leftofline(pts, p, i, j)
    if i == LowerRightBoundingIndex && j == LowerLeftBoundingIndex
        return -1
    elseif i == LowerRightBoundingIndex && j == UpperBoundingIndex
        return 1
    elseif i == LowerLeftBoundingIndex && j == LowerRightBoundingIndex
        return 1
    elseif i == LowerLeftBoundingIndex && j == UpperBoundingIndex
        return -1
    elseif i == UpperBoundingIndex && j == LowerLeftBoundingIndex
        return 1
    elseif i == UpperBoundingIndex && j == LowerRightBoundingIndex
        return -1
    end
    return leftofline(get_point(pts, i), get_point(pts, j), p) # If a line is from i → j, then k is left of i → j is (i, j, k) is a counter-clockwise triangle
end

"""
    intriangle(T::AbstractTriangle, pts::Points, p::AbstractPoint)

Tests if the point `p` is in the triangle `T`, where the vertices of 
`T = (i, j, k)` are `(pts[i], pts[j], pts[k])`. It is assumed that `𝒯` is 
positively oriented.
"""
function intriangle(e1, e2, e3) # https://stackoverflow.com/a/2049593
    if e1 == 0 || e2 == 0 || e3 == 0
        return 0
    end
    ininterior = e1 > 0 && e2 > 0 && e3 > 0
    if ininterior
        return 1
    else
        return -1
    end
end
function intriangle(T::AbstractTriangle, pts::Points, p::AbstractPoint)
    i, j, k = indices(T)
    e1 = leftofline(pts, p, i, j)
    e2 = leftofline(pts, p, j, k)
    e3 = leftofline(pts, p, k, i)
    return intriangle(e1, e2, e3)
end

"""
    edge_on_bounding_triangle(i, j)

Returns true if `(i, j)` is an edge of bounding triangle $(BoundingTriangle).
"""
function edge_on_bounding_triangle(i, j)
    #=
    if i ≥ FirstPointIndex || j ≥ FirstPointIndex
        return false
    elseif (i, j) == (LowerLeftBoundingIndex, LowerRightBoundingIndex) ||
           (i, j) == (LowerRightBoundingIndex, UpperBoundingIndex) ||
           (i, j) == (UpperBoundingIndex, LowerLeftBoundingIndex) ||
           (i, j) == (LowerRightBoundingIndex, LowerLeftBoundingIndex) ||
           (i, j) == (UpperBoundingIndex, LowerRightBoundingIndex) ||
           (i, j) == (LowerLeftBoundingIndex, UpperBoundingIndex)
        return true
    else
        return false
    end
    =#
    return i < FirstPointIndex && j < FirstPointIndex
end

"""
    islegal(i, j, k, ℓ, pts::Points)
    islegal(i, j, adj::Adjacent, pts::Points)

Returns `true` if the edge `(i, j)` is legal. It is assumed that 
`(i, j, k)` and `(j, i, ℓ)` are positively oriented triangles.
"""
function islegal(i, j, k, ℓ, pts::Points)
    if i ≥ FirstPointIndex && j ≥ FirstPointIndex && k ≥ FirstPointIndex && ℓ ≥ FirstPointIndex
        return incircle(pts, i, j, k, ℓ) ≤ 0
    else
        num_neg = num_less(FirstPointIndex, (i, j, k, ℓ))
        if num_neg == 1
            return i < FirstPointIndex || j < FirstPointIndex
        elseif num_neg == 2
            return min(i, j) < min(k, ℓ)
        else
            throw("Error occured.")
        end
    end
    throw("Error occured.")
end
@doc (@doc islegal(::Any, ::Any, ::Any, ::Any, ::Points))
function islegal(i, j, adj::Adjacent, pts::Points)
    edge_on_bounding_triangle(i, j) && return true
    k = get_edge(adj, i, j)
    ℓ = get_edge(adj, j, i)
    return is_legal(i, j, k, ℓ, pts)
end

############################################
##
## UTILITY FUNCTIONS 
##
############################################
"""
    find_root(G::HistoryDAG; method=:brute)

Finds the root of the graph `G`, using one of two methods:
- `method=:brute`: Uses brute force to look for the root, defining the root to be the point with zero out-degree.
- `method=:rng`: Uses randomised searching to look for the root.

See also [`_find_root_brute`](@ref) and [`_find_root_rng`](@ref).
"""
function find_root(G::HistoryDAG; method=:brute)
    if method == :brute
        return _find_root_brute(G)
    elseif method == :rng
        return _find_root_rng(G)
    end
end
function _find_root_brute(G::HistoryDAG)
    for (k, v) in graph(G).NN
        length(v) == 0 && return k
    end
end
function _find_root_rng(G::HistoryDAG)
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
points(DT::Triangulation, i) = get_point(points(DT), i)

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