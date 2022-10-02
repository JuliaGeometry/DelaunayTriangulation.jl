############################################
##
## ABSTRACT PRIMITIVES 
##
############################################
"""
    abstract type AbstractEdge{I} end

Abstract type representing an edge with indices of type `I`. 
Typically an `NTuple{2, Int64}`. Edges are iterable, writing 

- `i, j = e`,

which returns `(initial(e), terminal(e))`. Note that we define `length(e::AbstractEdge) = 2`.

`AbstractEdge`s must have a method for constucting from a function call, i.e. if `E::AbstractEdge{I}`,
then `E(i, j)` should construct an edge of that type.
"""
abstract type AbstractEdge{I} end
function Base.iterate(e::AbstractEdge, state=1)
    if state == 1
        return (initial(e), state + 1)
    elseif state == 2
        return (terminal(e), state + 1)
    else
        return nothing
    end
end
Base.eltype(::AbstractEdge{I}) where {I} = I
Base.eltype(::Type{AbstractEdge{I}}) where {I} = I
Base.length(e::AbstractEdge) = 2

"""
    abstract type AbstractTriangle{I} end

Abstract type representing a triangle with indices of type `I`. 
Typically an `NTuple{3, Int64}`.

`AbstractTriangle`s must have methods for:
- `geti(T)`: The first vertex.
- `getj(T)`: The second vertex. 
- `getk(T)`: The third vertex. 
- `indices(T)`: The indices for all three vertices.

`AbstractTriangle`s must also have a method for constucting from a function call, i.e. 
if `V::AbstractTriangle{I}`, then `E(i, j, k)` should construct a triangle of that type
with indices `(i, j, k)`.

Circular shifts of a triangle's indices can be obtained via 

- `shift_triangle_0(T) = (geti(T), getj(T), getk(T))`.
- `shift_triangle_1(T) = (getj(T), getk(T), geti(T))`.
- `shift_triangle_2(T) = (getk(T), geti(T), getj(T))`.

Circular shifts of the triangle itself can be obtained via 

- `shift_triangle_0(T::V) where {I,V<:AbstractTriangle{I}} = V(shift_triangle_indices_0(T))`
- `shift_triangle_1(T::V) where {I,V<:AbstractTriangle{I}} = V(shift_triangle_indices_1(T))`
- `shift_triangle_2(T::V) where {I,V<:AbstractTriangle{I}} = V(shift_triangle_indices_2(T))`

See also [`shift_triangle_indices`](@ref) and [`shift_triangle`](@ref).

Triangles are iterable, writing 

- `i, j, k = T`

will return `(geti(T), getj(T), getk(T))`. Note that we define `length(T::AbstractTriangle) = 3`.

Triangles are equal if their nidices are equal (up to a circular shift).
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
function Base.iterate(T::AbstractTriangle, state=1)
    if state == 1
        return (geti(T), state + 1)
    elseif state == 2
        return (getj(T), state + 1)
    elseif state == 3
        return (getk(T), state + 1)
    else
        return nothing
    end
end
Base.eltype(::AbstractTriangle{I}) where {I} = I
Base.eltype(::Type{AbstractTriangle{I}}) where {I} = I
Base.length(T::AbstractTriangle) = 3
function Base.:(==)(T::AbstractTriangle, V::AbstractTriangle)
    return indices(T) == indices(V) ||
           shift_triangle_indices_1(T) == indices(V) ||
           shift_triangle_indices_2(T) == indices(V)
end

"""
    abstract type AbstractPoint{T, A} end 

Abstract type representing a point with element type `T` stored in a container of 
type `A`. Typically a `Vector{Float64}` or `SVector{2, Float64}`.

`AbstractPoint`s must have methods for:
- `coords(p)`: The coordinates of the point.
- `getx(p)`: Get the `x`-coordinate of the point.
- `gety(p)`: Get the `y`-coordinate of the point.

`AbstractPoint`s must also have a method for constucting from a function call, i.e. 
if `P::AbstractPoint{T, A}`, then `P(x, y)` should construct a point of that type
with coordinates `(x, y)`.

Points are iterable, writing 

- `x, y = p`

gives `x = getx(p)` and `y = gety(p)`. Note that we define `length(p) = 2`.
"""
abstract type AbstractPoint{T,A} end
Base.Tuple(p::AbstractPoint{T,A}) where {T,A} = (getx(p), gety(p))
function Base.iterate(p::AbstractPoint, state=1)
    if state == 1
        return (getx(p), state + 1)
    elseif state == 2
        return (gety(p), state + 1)
        return nothing
    end
end
Base.eltype(::AbstractPoint{T,A}) where {T,A} = T
Base.eltype(::Type{AbstractPoint{T,A}}) where {T,A} = T
Base.length(p::AbstractPoint) = 2

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
    Triangle{I}(i, j, k) where {I} = Triangle(convert(I, i), convert(I, j), convert(I, k))
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
    Point{T,A}(x, y) where {T,A} = new{T,A}(A(collect((x, y))))
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
const BoundingTriangleShift = 20
const MinWidthHeight = 0.001

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
    Triangles{I,T}(triangles::Set{F}) where {I,T,F} = new{I,T}(convert(Set{T}, triangles))
    Triangles(triangles::Base.AbstractVecOrTuple{T}) where {I,T<:AbstractTriangle{I}} = new{I,T}(Set{T}(triangles))
    Triangles(triangles...) = Triangles(triangles)
end
triangles(tris::Triangles) = tris.triangles
Base.iterate(T::Triangles) = Base.iterate(triangles(T))
Base.iterate(T::Triangles, state) = Base.iterate(triangles(T), state)
Base.eltype(::Triangles{I,T}) where {I,T<:AbstractTriangle{I}} = T
Base.eltype(::Type{Triangles{I,T}}) where {I,T<:AbstractTriangle{I}} = T
Base.length(T::Triangles) = length(T.triangles)
Base.enumerate(T::Triangles) = Base.enumerate(triangles(T))

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
The points can be accessed using [`points`](@ref). The struct contains the following 
fields: 

- `points::Vector{P}`: The vector of points. 
- `xcentroid::T`: The `x`-coordinate of the centroid, `(xmax + xmin)/2`.
- `ycentroid::T`: The `y`-coordinate of the centroid, `(ymax + ymin)/2`.
- `xmin::T`: The minimum `x`-coordinate in the collection of points.
- `xmax::T`: The maximum `x`-coordinate in the collection of points. 
- `ymin::T`: The minimum `y`-coordinate in the collection of points. 
- `ymax::T`: The maximum `y`-coordinate in the collection of points.
- `width::T`: The width of the collection of points, defined by `xmax - xmin`.
- `height::T`: The height of the collection of points, defined by `ymax - ymin`.
- `max_width_height`: The maximum of `width` and `height`.

(Strictly speaking, the 'centroid' is the centre of the box that tightly bounds the points.)
"""
struct Points{T,A,P<:AbstractPoint{T,A}}
    points::Vector{P}
    xcentroid::T
    ycentroid::T
    xmin::T
    xmax::T
    ymin::T
    ymax::T
    width::T
    height::T
    max_width_height::T
    lower_left_bounding_triangle_coords::P
    lower_right_bounding_triangle_coords::P
    upper_bounding_triangle_coords::P
    function Points(points::AbstractVector{P}) where {T,A,P<:AbstractPoint{T,A}}
        xmin = typemax(T)
        xmax = typemin(T)
        ymin = typemax(T)
        ymax = typemin(T)
        for pt in points
            if getx(pt) < xmin
                xmin = getx(pt)
            end
            if getx(pt) > xmax
                xmax = getx(pt)
            end
            if gety(pt) < ymin
                ymin = gety(pt)
            end
            if gety(pt) > ymax
                ymax = gety(pt)
            end
        end
        width = max(xmax - xmin, MinWidthHeight)
        height = max(ymax - ymin, MinWidthHeight)
        xcentroid = (xmax + xmin) / 2
        ycentroid = (ymax + ymin) / 2
        max_width_height = max(width, height)
        lower_left_bounding_triangle_coords = P(xcentroid - BoundingTriangleShift * max_width_height,
            ycentroid - max_width_height)
        lower_right_bounding_triangle_coords = P(xcentroid + BoundingTriangleShift * max_width_height,
            ycentroid - max_width_height)
        upper_bounding_triangle_coords = P(xcentroid, ycentroid + BoundingTriangleShift * max_width_height)
        new{T,A,P}(points, xcentroid, ycentroid, xmin, xmax, ymin,
            ymax, width, height, max_width_height, lower_left_bounding_triangle_coords,
            lower_right_bounding_triangle_coords, upper_bounding_triangle_coords)
    end
    Points(points::NTuple{N,P}) where {N,T,A,P<:AbstractPoint{T,A}} = Points(collect(points))
end
points(pts::Points) = pts.points
Base.length(pts::Points) = length(points(pts))
Points(pts::AbstractVector) = Points(Point.(pts))
Points(pts::Points) = pts
Points(pts...) = Points(collect(pts))
Base.iterate(pts::Points) = Base.iterate(points(pts))
Base.iterate(pts::Points, state) = Base.iterate(points(pts), state)
Base.eltype(::Points{T,A,P}) where {T,A,P<:AbstractPoint{T,A}} = P
Base.eltype(::Type{Points{T,A,P}}) where {T,A,P<:AbstractPoint{T,A}} = P
width(pts::Points) = pts.width
height(pts::Points) = pts.height
xmin(pts::Points) = pts.xmin
ymin(pts::Points) = pts.ymin
xmax(pts::Points) = pts.xmax
ymax(pts::Points) = pts.ymax
xcentroid(pts::Points) = pts.xcentroid
ycentroid(pts::Points) = pts.ycentroid
max_width_height(pts::Points) = pts.max_width_height
lower_left_bounding_triangle_coords(pts::Points) = pts.lower_left_bounding_triangle_coords
lower_right_bounding_triangle_coords(pts::Points) = pts.lower_right_bounding_triangle_coords
upper_bounding_triangle_coords(pts::Points) = pts.upper_bounding_triangle_coords
Base.size(pts::Points) = Base.size(points(pts))
Base.eachindex(pts::Points) = Base.eachindex(points(pts))
Random.shuffle!(pts::Points) = Random.shuffle!(points(pts))

"""
    get_point(pts::Points, i)

Obtains the `i`th point in `pts`. We also allow for indices 

- `i = LowerLeftBoundingIndex`: This returns the coordinates for the lower-left vertex of the bounding triangle.
- `i = LowerRightBoundingIndex`: This returns the coordinates for the lower-right vertex of the bounding triangle.
- `i = UpperBoundingIndex`: This returns the coordinates for the upper vertex of the bounding triangle.

These latter coordinates are stored in `Points` rather than computed lazily.
"""
function get_point(pts::Points, i)
    if i ≥ FirstPointIndex
        return points(pts)[i]
    elseif i == LowerRightBoundingIndex
        return lower_right_bounding_triangle_coords(pts)
    elseif i == LowerLeftBoundingIndex
        return lower_left_bounding_triangle_coords(pts)
    elseif i == UpperBoundingIndex
        return upper_bounding_triangle_coords(pts)
    end
    throw(BoundsError(pts, i))
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
edges(adj::Adjacent) = keys(adj.adjacent)

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
    i, j, k = T
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
    add_neighbour!(DG::DelaunayGraph{I}, u, v...) 

Adds an edge from `u` to the points in `v` into the graph `DG`. Note that, as 
`DG` is an undirected graph, this also adds an edge from `v` to `u`.
"""
function add_neighbour!(DG::DelaunayGraph, u, v)
    add!(graph(DG), u, v)
    return nothing
end
@doc (@doc add_neighbour!(::DelaunayGraph, ::Any, ::Any))
function add_neighbour!(DG::DelaunayGraph, u, v...)
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
    u, v, w = T
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
function leftofline(pts::Points, p::AbstractPoint, i, j)
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
`T = (i, j, k)` are `(pts[i], pts[j], pts[k])`. It is assumed that `T` is 
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
    i, j, k = T
    if i < FirstPointIndex && j < FirstPointIndex && k < FirstPointIndex
        return 1 # all poiints are in the bounding triangle
    end
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
    return i < FirstPointIndex && j < FirstPointIndex
end

"""
    islegal(i, j, k, ℓ, pts::Points)
    islegal(i, j, adj::Adjacent, pts::Points)

Returns `true` if the edge `(i, j)` is legal. It is assumed that 
`(i, j, k)` and `(j, i, ℓ)` are positively oriented triangles.
"""
function islegal(i, j, k, ℓ, pts::Points)
    #=
    if i ≥ FirstPointIndex && j ≥ FirstPointIndex && k ≥ FirstPointIndex && ℓ ≥ FirstPointIndex
        return incircle(pts, i, j, k, ℓ) ≤ 0
    else
        num_neg = num_less(FirstPointIndex, (i, j, k, ℓ))
        if num_neg == 1
            return !(i < FirstPointIndex || j < FirstPointIndex)
        elseif num_neg == 2
            return min(i, j) < min(k, ℓ)
        else
            throw("Error occured.")
        end
    end
    throw("Error occured.")
    =#
    return incircle(pts, i, j, k, ℓ) ≤ 0
end
@doc (@doc islegal(::Any, ::Any, ::Any, ::Any, ::Points))
function islegal(i, j, adj::Adjacent, pts::Points)
    edge_on_bounding_triangle(i, j) && return true
    k = get_edge(adj, i, j)
    ℓ = get_edge(adj, j, i)
    return islegal(i, j, k, ℓ, pts)
end

############################################
##
## MAIN TRIANGULATION FUNCTIONS
##
############################################
"""
    initialise_triangulation(pts; 
    IntegerType=Int64,
    TriangleType=Triangle{IntegerType},
    EdgeType=Edge{IntegerType})

This function returns the initial form of a [`Triangulation`](@ref) data structure, storing 
the points in `pts` (converted into a `Points` type). You can specify custom integer, 
triangle, and edge types using the keywords `IntegerType`, `TriangleType`, and 
`EdgeType`, respectively.
"""
function initialise_triangulation(pts;
    IntegerType::Type{I}=Int64,
    TriangleType=Triangle{IntegerType},
    EdgeType=Edge{IntegerType}) where {I}
    # The data structures
    root = TriangleType(I(LowerRightBoundingIndex),
        I(UpperBoundingIndex),
        I(LowerLeftBoundingIndex))
    T = Triangles{I,TriangleType}(Set{TriangleType}([root]))
    HG = HistoryDAG{I,TriangleType}()
    adj = Adjacent{I,EdgeType}()
    adj2v = Adjacent2Vertex{I,EdgeType}()
    DG = DelaunayGraph{I}()
    # Add the root to the DAG
    add_triangle!(HG, root)
    # Add the initial adjacencies 
    add_edge!(adj, I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex))
    add_edge!(adj, I(UpperBoundingIndex), I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex))
    add_edge!(adj, I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex), I(UpperBoundingIndex))
    add_edge!(adj, I(LowerRightBoundingIndex), I(LowerLeftBoundingIndex), I(BoundaryIndex))
    add_edge!(adj, I(LowerLeftBoundingIndex), I(UpperBoundingIndex), I(BoundaryIndex))
    add_edge!(adj, I(UpperBoundingIndex), I(LowerRightBoundingIndex), I(BoundaryIndex))
    add_edge!(adj2v, I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex), I(UpperBoundingIndex))
    add_edge!(adj2v, I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex))
    add_edge!(adj2v, I(UpperBoundingIndex), I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex))
    # Add the initial neighbours 
    add_point!(DG, I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex), I(UpperBoundingIndex))
    add_neighbour!(DG, I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex), I(UpperBoundingIndex))
    add_neighbour!(DG, I(LowerRightBoundingIndex), I(LowerLeftBoundingIndex), I(UpperBoundingIndex))
    add_neighbour!(DG, I(UpperBoundingIndex), I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex))
    return Triangulation(adj, adj2v, DG, HG, T, Points(pts), root)
end

"""
    locate_triangle(HG::HistoryDAG, pts, p, init=find_root(HG; method=:rng))

Given the point location data structure `HG` and a set of `pts`, finds the triangle in 
the current triangulation such that `p` is in its interior. The point location starts at `init`.
The function is recursive, and returns a tuple `(tri, flag)`:

    - `tri`: This is the triangle that `p` is in.
    - `flag`: If `flag == 0`, then `p` is on an edge of `tri`. Otherwise, it is in the open interior.
"""
function locate_triangle(HG::HistoryDAG, pts, p, init=find_root(HG; method=:rng))
    if out_deg(HG, init) == 0
        return init, intriangle(init, pts, p)
    end
    out = out_neighbors(HG, init)
    for T in out
        intriangle(T, pts, p) ≥ 0 && return locate_triangle(HG, pts, p, T)
    end
    throw("Failed to find triangle.")
end

"""
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)

Given a triangulation `T`, adds the `r`th point of the point set into the triangulation.

# Arguments 
- `T`: The current triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
-` Tᵢⱼₖ`: The triangle that the `r`th point is inside of. Must be positively oriented.
- `r`: The index of the point in the original point set that is being introduced.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function add_point!(T::Triangles{I,V}, HG::HistoryDAG,
    adj::Adjacent, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, Tᵢⱼₖ::V, r) where {I,V<:AbstractTriangle{I}}
    i, j, k = Tᵢⱼₖ # The triangle to be split into three
    delete_triangle!(T, Tᵢⱼₖ) # Now that we've split the triangle, we can remove the triangle
    T₁, T₂, T₃ = V(i, j, r), V(j, k, r), V(k, i, r) # New triangles to add. Note that these triangles are all positively oriented.
    add_triangle!(T, T₁, T₂, T₃) # The three new triangles
    add_triangle!(HG, T₁, T₂, T₃) # Add the new triangles into DAG
    add_edge!(HG, Tᵢⱼₖ, T₁, T₂, T₃) # Add edges from the old triangle to the new triangles
    add_triangle!(adj, T₁, T₂, T₃) # Add the new edges into the adjacency list
    update_after_insertion!(adj2v, i, j, k, r)
    add_neighbour!(DG, r, i, j, k)
    add_neighbour!(DG, i, r)
    add_neighbour!(DG, j, r)
    add_neighbour!(DG, k, r)
    return nothing
end

"""
    flip_edge!(T, HG, adj, adj2v, DG, i, j, k, r)

Performs an edge flip, flipping the edge `(i, j)` into the edge `(k, r)`.

# Arguments
- `T`: The current triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
- `i, j`: The current edge.
- `k, r`: Indices for the points the edge is flipped onto.

It is assumed that `(i, k, j)` and `(i, j, r)` are positively oriented triangles.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function flip_edge!(T::Triangles{I,V}, HG::HistoryDAG,
    adj::Adjacent, adj2v::Adjacent2Vertex, DG::DelaunayGraph, i, j, k, r) where {I,V<:AbstractTriangle{I}}
    # The old triangles
    Tᵢₖⱼ = V(i, k, j)
    Tᵢⱼᵣ = V(i, j, r)
    delete_triangle!(T, Tᵢₖⱼ, Tᵢⱼᵣ)
    delete_edge!(adj, i, j)
    delete_neighbour!(DG, i, j) #delete_neighbour!(DG, j, i)
    # The new triangles 
    Tᵣₖⱼ = V(r, k, j)
    Tᵣᵢₖ = V(r, i, k)
    # Add the new triangles to the data structure
    add_triangle!(T, Tᵣₖⱼ, Tᵣᵢₖ)
    add_triangle!(HG, Tᵣₖⱼ, Tᵣᵢₖ)
    add_triangle!(adj, Tᵣₖⱼ, Tᵣᵢₖ)
    update_after_flip!(adj2v, i, j, k, r)
    # Connect the new triangles to the replaced triangles in the DAG
    add_edge!(HG, Tᵢₖⱼ, Tᵣₖⱼ, Tᵣᵢₖ)
    add_edge!(HG, Tᵢⱼᵣ, Tᵣₖⱼ, Tᵣᵢₖ)
    # Add the new neighbours 
    add_neighbour!(DG, r, k) # add_neighbour!(𝒟𝒢, k, r)
    return nothing
end

"""
    legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    
Legalises the edge `(i, j)` if it is illegal.

# Arguments 
- `T`: The current triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
- `i, j`: The edge to make legal. Nothing happens if `is_legal(i, j, 𝒜, pts)`.
- `r`: The point being added into the triangulation. 
- `pts`: The point set of the triangulation.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function legalise_edge!(T::Triangles, HG::HistoryDAG,
    adj::Adjacent, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, i, j, r, pts::Points)
    if !islegal(i, j, adj, pts)
        e = get_edge(adj, j, i)
        flip_edge!(T, HG, adj, adj2v, DG, i, j, e, r)
        legalise_edge!(T, HG, adj, adj2v, DG, i, e, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, e, j, r, pts)
    end
    return nothing
end

"""
    remove_bounding_triangle!(DT::Triangulation)
    remove_bounding_triangle!(T, adj, adj2v, DG)

Remove the bounding triangle from the triangulation, where 

- `T`: The current triangulation.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.

These structures are updated in-place.
"""
function remove_bounding_triangle!(T::Triangles{I,V}, adj::Adjacent,
    adj2v::Adjacent2Vertex, DG::DelaunayGraph) where {I,V<:AbstractTriangle{I}}
    for w in BoundingTriangle
        neighbours = get_edge(adj2v, w)
        for (u, v) in neighbours # (u, v, w) is a triangle..
            delete_edge!(adj, w, u; protect_boundary=false)
            delete_edge!(adj, w, v; protect_boundary=false)
            delete_edge!(adj2v, u, v, w)
            delete_edge!(adj2v, v, w, u)
            if u ≥ FirstPointIndex && v ≥ FirstPointIndex # This can only be a boundary edge
                add_edge!(adj2v, BoundaryIndex, u, v)
                add_edge!(adj, u, v, BoundaryIndex)
            end
            delete_triangle!(T, V(u, v, w))
        end
        delete_point!(DG, w)
        delete_point!(adj2v, w)
    end
    return nothing
end
@doc (@doc remove_bounding_triangle!(::Triangles{I,V}, ::Adjacent,
    ::Adjacent2Vertex, ::DelaunayGraph) where {I,V<:AbstractTriangle{I}})
function remove_bounding_triangle!(DT::Triangulation)
    remove_bounding_triangle!(triangles(DT),
        adjacent(DT), adjacent2vertex(DT), graph(DT))
    return nothing
end

"""
    add_point!(DT::Triangulation, r)
    add_point!(T::Triangles, HG::HistoryDAG,
     adj::Adjacent, adj2v::Adjacent2Vertex, 
     DG::DelaunayGraph, root, pts, r)
"""
function add_point!(T::Triangles, HG::HistoryDAG,
    adj::Adjacent, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, root, pts, r)
    pᵣ = get_point(pts, r)
    Tᵢⱼₖ, interior_flag = locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    if interior_flag == 1
        add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
        legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, j, k, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, k, i, r, pts)
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
@doc (@doc add_point!(::Triangles, ::HistoryDAG, ::Adjacent, ::Adjacent2Vertex, ::DelaunayGraph, ::Any, ::Any, ::Any))
function add_point!(DT::Triangulation, r)
    add_point!(triangles(DT), history(DT), adjacent(DT),
        adjacent2vertex(DT), graph(DT), root(DT), points(DT), r)
    return nothing
end

"""
    triangulate(pts; shuffle_pts=true, trim=true)

Computes the Delaunay triangulation of the points in `pts` using randomised incremental 
insertion. The points are shuffled in-place, but this shuffling can be disabled by 
setting `shuffle_pts=false`. The bounding triangle of the triangulation can be retained 
by setting `trim=true`.
"""
function triangulate(pts; shuffle_pts=true, trim=true, method = :berg)
    # Base.require_one_based_indexing(pts)
    DT = initialise_triangulation(pts)
    shuffle_pts && @views shuffle!(points(DT))
    for r in eachindex(points(DT))
        add_point!(DT, r)
    end
    trim && remove_bounding_triangle!(DT)
    return DT
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