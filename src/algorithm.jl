############################################
##
## INTERFACE 
##
############################################
"""
    Interfaces 

This package assumes that a certain set of methods have 
been defined for working with triangles, edges, and points. 
We list the assumed methods below. We also provide predefined 
types that can be readily used, and also provide an example for 
implementing these interfaces. These predefined interfaces can be
seen using `?predefined_types` in the Julia REPL. (This documentation 
here is accessible via `?required_interfaces`.)

# Triangle Interface

Triangles are essentially treated as if they are triples of indices, 
`(i, j, k)`, with each index corresponding to a point. We assume the 
following methods, assuming `TriangleType` is your triangle struct:

- `indices(T::TriangleType)`: Returns the triple of indices for the triangle 
`T`. The returned value should be such that `i, j, k = indices(T)` returns the 
individual indices (in particular, `indices(T)` should be iterable).
- `geti(T::TriangleType)`: If `T` is assumed to act like `(i, j, k)`, then `geti(T)` should return `i`.
- `getj(T::TriangleType)`: Similar to above, except for `j`.
- `getk(T::TriangleType)`: Similar to above, except for `k`.
- `construct_triangle(::Type{Triangle}, i, j, k)`: Constructs a triangle of your type with indices `i`, `j`, and `k`. 

# Edge Interface 

Edges are essentially treated as if they are tuples `(i, j)`, defining 
an oriented edge from the `i`th point to the `j`th point. We assume the 
following method, assuming `EdgeType` is your edge struct:

- `construct_edge(::Type{EdgeType}, i, j)`: Constructs an edge from `i` to `j` of type `EdgeType`.

# Point Interface 

Points are essentially treated as if they are tuples `(x, y)`, defining a point 
with coordinate `(x, y)`. We assume the following methods, assuming `PointType` is 
your point struct:

- `Base.Tuple(p::PointType)` or `ExactPredicates.coord(p::PointType)`: This should return a 
tuple of the coordinates, `(x, y)`. This is needed to work with the predicates defined in 
`ExactPredicates.jl`.
- `getx(p::PointType)`: Returns the `x`-coordinate. 
` `gety(p::PointType)`: Returns the `y`-coordinate.

# Collection of Triangles Interface

We also need an interface for treating a collection of triangles. The collection 
must be iterable, allowing e.g. the use of `push!` and `get!`. We assume the following 
methods, assuming `TrianglesType` is your struct for the collection of triangles:

- `construct_triangles(::Type{TrianglesType})`: Creates an empty collection of type `TrianglesType`. Triangles will be pushed into this using `push!`. We do provide a default implementation of this,
namely `construct_triangles(::Type{T}) where T = T()`, but you could define your own method for this.
- `triangle_type(::Type{TrianglesType}):`: Returns the type of triangle in your collection, allowing for the creation of a triangle that can be pushed into the collection. We do provide a 
default implementation for this, namely `triangle_type(::Type{T}) where T = eltype(T)`, but you could define your own method for this.
- `Base.delete!(T::TrianglesType, V::TriangleType)`: Deletes the triangle `V` from the collection `T`.

# Collection of Edges Interface 

We also need an interface for treating a collection of edges. The collection 
must be iterable, allowing e.g. the use of `push!` and `get!`, and allowing iteration 
of the form `for (i, j) in`.

# Collection of Points Interface

We also need an interface for treating a collection of points. The collection
must be iterable, allowing e.g. the use of `push!` and `get!`. 
We assume the following methods, assuming `PointsType` is the 
struct for your collection of points:

- `get_point(pts::PointsType, i)`: Returns the `i`th point in the collection. 
- `get_point(pts::PointsType, i::Vararg{I, N}) where {I, N}`: If `i::Vararg{I, N}` is, for example, 
`(i, j, k)`, then this should return `get_point(pts, i)`, `get_point(pts, j)`, `get_point(pts, k)`, so that we 
could write e.g. `pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)` instead of writing `get_point` three times.
"""
function required_interfaces end

## Triangle
function indices end            # indices(::Triangle)
function geti end               # geti(::Triangle)
function getj end               # getj(::Triangle)
function getk end               # getk(::Triangle)
function construct_triangle end # construct_triangle(::Type{Triangle}, ::I, ::I, ::I)

shift_triangle_1(T::V) where {V} = construct_triangle(V, getj(T), getk(T), geti(T))
shift_triangle_2(T::V) where {V} = construct_triangle(V, getk(T), geti(T), getj(T))
shift_triangle(T::V, i) where {V} = i == 0 ? T : (i == 1 ? shift_triangle_1(T) : shift_triangle_2(T))

## Edge 
function construct_edge end # construct_edge(::Type{Edge}, ::I, ::I)

## Point 
function getx end
function gety end

## Triangles 
function construct_triangles end # construct_triangles(::Type{Triangles})
function triangle_type end       # triangle_type(::Type{Triangles})

"""
    add_triangle!(T, V::Vararg{Tri, N}) where {Tri, N}

Adds the triangles in `V` into the collection of triangles `T`. No 
checks are made for triangles that are equal under circular shifts.
"""
function add_triangle!(T, V::Vararg{Tri,N}) where {Tri,N}
    push!(T, V...)
    return nothing
end

"""
    delete_triangle!(T, V)
    delete_triangle!(T, V::Vararg{Tri, N}) where {Tri, N}

Deletes the triangle `V` from the collection of triangles `T`. All circular 
shifts of `V`'s indices are also deleted.
"""
function delete_triangle!(T, V)
    delete!(T, V)
    delete!(T, shift_triangle_1(V))
    delete!(T, shift_triangle_2(V))
    return nothing
end
@doc (@doc delete_triangle!(::Any, ::Any))
function delete_triangle!(T, V::Vararg{Tri,N}) where {Tri,N}
    for i in 1:N
        delete_triangle!(T, V[i])
    end
    return nothing
end

## Points 
function get_point end      # get_point(::Points, ::I)
function add_point!(pts, p)
    push!(pts, p)
    return nothing
end
function add_point!(pts, p::Vararg{P,N}) where {P,N}
    for i in 1:N
        add_point!(pts, p[i])
    end
    return nothing
end

############################################
##
## PREDEFINED INTERFACES 
##
############################################
## Triangle
indices(T::NTuple{3,I}) where {I} = T
geti(T::NTuple{3,I}) where {I} = T[1]
getj(T::NTuple{3,I}) where {I} = T[2]
getk(T::NTuple{3,I}) where {I} = T[3]
construct_triangle(::Type{NTuple{3,I}}, i, j, k) where {I} = (i, j, k)

## Edge 
construct_edge(::Type{NTuple{2,I}}, i, j) where {I} = (i, j)

## Point 
getx(p::Base.AbstractVecOrTuple) = p[1]
gety(p::Base.AbstractVecOrTuple) = p[2]

## Triangles 
construct_triangles(::Type{T}) where {T} = T()
triangle_type(::Type{T}) where {T} = eltype(T)

## Points 
function get_point(pts::AbstractVector, i)
    if i ≥ FirstPointIndex
        return pts[i]
    elseif i == LowerRightBoundingIndex
        return lower_right_bounding_triangle_coords(pts)
    elseif i == LowerLeftBoundingIndex
        return lower_left_bounding_triangle_coords(pts)
    elseif i == UpperBoundingIndex
        return upper_bounding_triangle_coords(pts)
    end
    throw(BoundsError(pts, i))
end
function get_point(pts::AbstractVector, i::Vararg{I,N}) where {I,N}
    return ntuple(j -> get_point(pts, i[j]), Val(N))
end
function point_stats(points)
    T = Float64
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
    return xmin, xmax, ymin, ymax, width, height, xcentroid, ycentroid, max_width_height
end
function lower_right_bounding_triangle_coords(points)
    xmin, xmax, ymin, ymax, width, height,
    xcentroid, ycentroid, max_width_height = point_stats(points)
    return (xcentroid + BoundingTriangleShift * max_width_height,
        ycentroid - max_width_height)
end
function lower_left_bounding_triangle_coords(points)
    xmin, xmax, ymin, ymax, width, height,
    xcentroid, ycentroid, max_width_height = point_stats(points)
    return (xcentroid - BoundingTriangleShift * max_width_height,
        ycentroid - max_width_height)
end
function upper_bounding_triangle_coords(points)
    xmin, xmax, ymin, ymax, width, height,
    xcentroid, ycentroid, max_width_height = point_stats(points)
    return (xcentroid,
        ycentroid + BoundingTriangleShift * max_width_height)
end

############################################
##
## PREDICATES 
##
############################################
"""
    isoriented(T, pts)

Tests if the triangle `T` is positively oriented, returning
`1` if the triangle is positively oriented, `0` 
if the triangle is degenerate (i.e. the points are all 
collinear), and `-1` if the triangle is negatively 
oriented.

The argument `pts` should be a collection of points 
such that the `i`th vertex of `T` corresponds to 
`get_point(pts, i)`.

The check is exact, making use of `ExactPredicates.orient` (unless 
you define a new method for it for your point type, obviously).
"""
function isoriented(T, pts)
    i, j, k = indices(T)
    pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
    return orient(pᵢ, pⱼ, pₖ)
end

"""
    isincircle(pts, i, j, k, ℓ)

Tests if the `i`th point is in the circle through the 
`j`th, `k`th, and `ℓ` points (as obtained via `get_point(pts, i)`, etc.).
Returns `1` if the point is inside, `0` if the point is on the circle,
and `-1` if the point is outside. It is assumed that the 
`j`th, `k`th`, and `ℓ`th points are provided counter-clockwise. 

The check is exact, making use of `ExactPredicates.incircle` (unless 
you define a new method for it for your point type, obviously).
"""
function isincircle(pts, i::I, j::I, k::I, ℓ::I) where {I}
    pᵢ, pⱼ, pₖ, pₗ = get_point(pts, i, j, k, ℓ)
    return I(incircle(pᵢ, pⱼ, pₖ, pₗ))
end

"""
    isleftofline(pts, i, j, k)

Tests if the `k`th point if left of the line from the `i`th point 
to the `j`th point (as obtained via `get_point(pts, i)`, etc.). 
Returns `1` if the point is to the left, `0` if the points are 
collinear, and `-1` if the point is to the right.

The check is exact, making use of `ExactPredicates.orient` (unless 
you define a new method for it for your point type, obviously).
"""
function isleftofline(pts, i::I, j::I, k::I) where {I}
    if i == I(LowerRightBoundingIndex) && j == I(LowerLeftBoundingIndex)
        return I(-1)
    elseif i == I(LowerRightBoundingIndex) && j == I(UpperBoundingIndex)
        return I(1)
    elseif i == I(LowerLeftBoundingIndex) && j == I(LowerRightBoundingIndex)
        return I(1)
    elseif i == I(LowerLeftBoundingIndex) && j == I(UpperBoundingIndex)
        return I(-1)
    elseif i == I(UpperBoundingIndex) && j == I(LowerLeftBoundingIndex)
        return I(1)
    elseif i == I(UpperBoundingIndex) && j == I(LowerRightBoundingIndex)
        return I(-1)
    end
    pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
    return I(orient(pᵢ, pⱼ, pₖ))
end

"""
    isintriangle(T, pts, ℓ)

Tests if the `ℓ`th point of `pts` is in the triangle `T`. Returns `1`
if the point is inside, `0` if a point is on an edge of `T`, and `-1` 
if the point is outside. It is assumed that `T` is positively oriented. 

The check is exact, making use of `isleftofline`. 
"""
function isintriangle(T, pts, ℓ::I) where {I}
    i, j, k = indices(T)
    if i < I(FirstPointIndex) && j < I(FirstPointIndex) && k < I(FirstPointIndex)
        return I(1)
    end
    is_left_of_edge_ij = isleftofline(pts, i, j, ℓ)
    is_left_of_edge_jk = isleftofline(pts, j, k, ℓ)
    is_left_of_edge_ki = isleftofline(pts, k, i, ℓ)
    return isintriangle(
        is_left_of_edge_ij,
        is_left_of_edge_jk,
        is_left_of_edge_ki
    )
end
function isintriangle(e1::I, e2::I, e3::I) where {I}
    if e1 == I(0) || e2 == I(0) || e3 == I(0)
        return I(0)
    end
    isininterior = e1 == I(1) && e2 == I(1) && e3 == I(1)
    if isininterior
        return I(1)
    else
        return I(-1)
    end
end

"""
    find_edge(T, pts, ℓ)

Given a triangle `T` and a set of points `pts`, with the `ℓ`th point of `pts` 
on an edge of `T`, find the edge of `T` that the point is on.
"""
function find_edge(T, pts, ℓ)
    i, j, k = indices(T)
    isleftofline(pts, i, j, ℓ) == 0 && return (i, j)
    isleftofline(pts, j, k, ℓ) == 0 && return (j, k)
    isleftofline(pts, k, i, ℓ) == 0 && return (k, i)
    throw("The point $p is not on an edge of $T.")
end

"""
    edge_on_bounding_triangle(i, j)

Returns true if `(i, j)` is an edge of bounding triangle $(BoundingTriangle).
"""
edge_on_bounding_triangle(i, j) = i < FirstPointIndex && j < FirstPointIndex

"""
    islegal(i, j, adj, pts)

Tests if the edge `(i, j)` is legal. `adj` is the adjacent map of the 
triangulation, and `pts` is the point set.
"""
function islegal(i, j, adj, pts)
    edge_on_bounding_triangle(i, j) && return true
    k = get_edge(adj, i, j)
    ℓ = get_edge(adj, j, i)
    return islegal(i, j, k, ℓ, pts)
end
function islegal(i::I, j::I, k::I, ℓ::I, pts) where {I}
    incirc = isincircle(pts, i, j, k, ℓ)
    return incirc == I(0) || incirc == I(-1)
end

############################################
##
## ADJACENT
##
############################################
"""
    Adjacent{I, E}

The adjacent map of the triangulation, mapping an edge `(u, v)` to a vertex `w` 
such that `(u, v, w)` is a positively oriented triangle. The map is stored in the 
field `adjacent` as a `Dict{E, I}`, where `E` is the edge type and `I` is the 
integer type.
"""
struct Adjacent{I,E} # Adjacent{IntegerType, EdgeType}
    adjacent::Dict{E,I}
    function Adjacent{I,E}() where {I,E}
        A = Dict{E,I}()
        adj = new{I,E}(A)
        return adj
    end
end
adjacent(adj::Adjacent) = adj.adjacent
edges(adj::Adjacent) = keys(adjacent(adj))

"""
    get_edge(adj::Adjacent, uv)
    get_edge(adj::Adjacent{I, E}, u, v) where {I, E}

Returns the vertex `w` such that `(u, v, w)` is a positively oriented triangle, indexing 
the adjacent map `adj` at the key `(u, v)`.
"""
get_edge(adj::Adjacent, uv) = adjacent(adj)[uv]
@doc (@doc get_edge(::Adjacent, ::Any))
function get_edge(adj::Adjacent{I,E}, u, v) where {I,E}
    uv = construct_edge(E, u, v)
    w = get_edge(adj, uv)
    return w
end

"""
    add_edge!(adj::Adjacent, uv, w) 
    add_edge!(adj::Adjacent{I, E}, u, v, w) where {I, E}

Adds the edge `(u, v)` into the adjacent map `adj`.
"""
function add_edge!(adj::Adjacent, uv, w)
    adjacent(adj)[uv] = w
    return nothing
end
@doc (@doc add_edge!(::Adjacent, ::Any, ::Any))
function add_edge!(adj::Adjacent{I,E}, u, v, w) where {I,E}
    uv = construct_edge(E, u, v)
    add_edge!(adj, uv, w)
    return nothing
end

"""
    delete_edge!(adj::Adjacent{I, E}, u, v; protect_boundary=true) where {I,E}

Deletes the edge `(u, v)` from the map `adj`. This function also deletes the edge
`(v, u)`. If the boundary needs to be protected, so that `(u, v)` is not deleted if 
the adjacent vertex to `(u, v)` is `$BoundaryIndex` (and similarly for `(v, u)`), 
set `protect_boundary=true`.
"""
function delete_edge!(adj::Adjacent{I,E}, u, v; protect_boundary=true) where {I,E}
    Euv = construct_edge(E, u, v)
    Evu = construct_edge(E, v, u)
    (!protect_boundary || get_edge(adj, Euv) ≠ I(BoundaryIndex)) && delete!(adjacent(adj), Euv)
    (!protect_boundary || get_edge(adj, Evu) ≠ I(BoundaryIndex)) && delete!(adjacent(adj), Evu)
    return nothing
end

"""
    add_triangle!(adj::Adjacent, T...)

Adds the edges of the triangle `T` to the adjacent map `adj`.`
"""
function add_triangle!(adj::Adjacent, T)
    i, j, k = T
    add_edge!(adj, i, j, k)
    add_edge!(adj, j, k, i)
    add_edge!(adj, k, i, j)
end
@doc (@doc add_triangle!(::Adjacent, ::Any))
function add_triangle!(adj::Adjacent, T::Vararg{Tri,N}) where {Tri,N}
    for i in 1:N
        add_triangle!(adj, T[i])
    end
    return nothing
end

############################################
##
## ADJACENT2VERTEX 
##
############################################
"""
    Adjacent2Vertex{I,Es,E}

The adjacent-to-vertex map of the triangulation, mapping a vertex `w` 
to the collection of all edges `(u, v)` such that `(u, v, w)` is a positively 
oriented triangle. The map is stored in `adjacent2vertex` as a `Dict{I, Es}`, 
where `I` is the integer type and `Es` is the type representing the collection 
of edges. The extra parameter `E` is the edge type.
"""
struct Adjacent2Vertex{I,Es,E}
    adjacent2vertex::Dict{I,Es}
    function Adjacent2Vertex{I,Es,E}() where {I,Es,E}
        D = Dict{I,Es}()
        TA2V = new{I,Es,E}(D)
        return TA2V
    end
end
adjacent2vertex(adj2v::Adjacent2Vertex) = adj2v.adjacent2vertex

"""
    get_edge(adj2v::Adjacent2Vertex, w)

Gets the set of edges adjacent to `w` from the map `adj2v`.
"""
get_edge(adj2v::Adjacent2Vertex, w) = adjacent2vertex(adj2v)[w]

"""
    add_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, uv) where {I,Es,E}
    add_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}

Adds the edge `(u, v)` to the set of edges in `adjacent2vertex[w]`.
"""
function add_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, uv) where {I,Es,E}
    existing_segments = get!(Es, adjacent2vertex(adj2v), w)
    push!(existing_segments, uv)
    return nothing
end
@doc (@doc add_edge!(::Adjacent2Vertex{I,Es,E}, w, uv) where {I,Es,E})
function add_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}
    uv = construct_edge(E, u, v)
    add_edge!(adj2v, w, uv)
    return nothing
end

function delete_edge!(adj2v::Adjacent2Vertex, w, uv)
    uvs = get_edge(adj2v, w)
    delete!(uvs, uv)
    return nothing
end
@doc (@doc delete_edge!(::Adjacent2Vertex, ::Any, ::Any))
function delete_edge!(adj2v::Adjacent2Vertex{I,Es,E}, w, u, v) where {I,Es,E}
    uv = construct_edge(E, u, v)
    delete_edge!(adj2v, w, uv)
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

"""
    update_after_split!(adj2v::Adjacent2Vertex, i, j, k, ℓ, r)

Updates the map `adj2v` after the edge `(i, j)`, incident to the two positively 
oriented triangles `(i, j, k)` and `(j, i, ℓ)`, is split into two so that the triangles 
`(i, j, k)` and `(j, i, ℓ`) are split from a point `r` to the points `k` and `ℓ`, 
respectively, assuming `r` is on the edge `(i, j)`.
"""
function update_after_split!(adj2v::Adjacent2Vertex, i, j, k, ℓ, r)
    # Delete old edges 
    delete_edge!(adj2v, i, j, k)
    delete_edge!(adj2v, i, ℓ, j)
    delete_edge!(adj2v, j, k, i)
    delete_edge!(adj2v, j, i, ℓ)
    delete_edge!(adj2v, k, i, j)
    delete_edge!(adj2v, ℓ, j, i)
    # Add new edges 
    add_edge!(adj2v, i, r, k)
    add_edge!(adj2v, i, ℓ, r)
    add_edge!(adj2v, j, k, r)
    add_edge!(adj2v, j, r, ℓ)
    add_edge!(adj2v, k, i, r)
    add_edge!(adj2v, k, r, j)
    add_edge!(adj2v, ℓ, r, i)
    add_edge!(adj2v, ℓ, j, r)
    add_edge!(adj2v, r, k, i)
    add_edge!(adj2v, r, j, k)
    add_edge!(adj2v, r, ℓ, r)
    add_edge!(adj2v, r, ℓ, j)
    return nothing
end

############################################
##
## DELAUNAYGRAPH
##
############################################
"""
    struct DelaunayGraph{I}

Graph representation of a Delaunay triangulation, mapping all points to their neighbours, 
where a point `v` is a neighbour of `u` if `(u, v)` is an edge in the triangulation. The 
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
    add_point!(DG::DelaunayGraph, u)
    add_point!(DG::DelaunayGraph, u::Vararg{I, N})

Addds the point(s) `u` nito the vertex set of `DG`.
"""
function add_point!(DG::DelaunayGraph, u)
    add!(graph(DG), u)
    return nothing
end
@doc (@doc add_point!(::DelaunayGraph, ::Any))
function add_point!(DG::DelaunayGraph, u::Vararg{I,N}) where {I,N}
    for i in 1:N
        add_point!(DG, u[i])
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
function add_neighbour!(DG::DelaunayGraph, u, v::Vararg{I,N}) where {I,N}
    for i in 1:N
        add_neighbour!(DG, u, v[i])
    end
    return nothing
end

"""
    get_neighbour(DG::DelaunayGraph, u)

Gets the neighbours of `u` from the graph `DG`.
"""
function get_neighbour(DG::DelaunayGraph, u)
    return graph(DG).N[u]
end

"""
    delete_neighbour!(DG::DelaunayGraph, u, v)
    delete_neighbour!(DG::DelaunayGraph, u, v::Vararg{I,N}) where {I,N}

Deletes `v` from the neighbourhood of `u`. Note that, as `DG` is an 
undirected graph, this also deletes `u` from the neighbourhood of `v`.
"""
function delete_neighbour!(DG::DelaunayGraph, u, v)
    delete!(graph(DG), u, v)
    return nothing
end
@doc (@doc delete_neighbour!(::DelaunayGraph, ::Any, ::Any))
function delete_neighbour!(DG::DelaunayGraph, u, v::Vararg{I,N}) where {I,N}
    for i in 1:N
        delete_neighbour!(DG, u, v[i])
    end
    return nothing
end

"""
    delete_point!(DG::DelaunayGraph, u)
    delete_point!(DG::DelaunayGraph, u::Vararg{I,N}) where {I,N}

Deletes the point(s) `u` from the vertex set of `DG`.
"""
function delete_point!(DG::DelaunayGraph, u)
    delete!(graph(DG), u)
    return nothing
end
@doc (@doc delete_point!(::DelaunayGraph, ::Any))
function delete_point!(DG::DelaunayGraph, u::Vararg{I,N}) where {I,N}
    for i in 1:N
        delete_point!(DG, u[i])
    end
    return nothing
end

############################################
##
## HISTORYGRAPH
##
############################################
"""
    struct HistoryGraph{I,T<:AbstractTriangle{I}} 

The point location graph for the Delaunay triangulation. This is a directed graph, 
implemented using `DirectedGraph`, that stores the history of the triangulation as points are 
inserted. The graph can be accessed using [`graph`](@ref).
"""
struct HistoryGraph{T}
    graph::DirectedGraph{T}
    function HistoryGraph{T}() where {T}
        G = DirectedGraph{T}()
        forbid_loops!(G)
        TDAG = new{T}(G)
        return TDAG
    end
    HistoryGraph(HG::DirectedGraph{T}) where {T} = new{T}(HG)
end
graph(G::HistoryGraph) = G.graph

"""
    SimpleGraphs.out_neighbors(G::HistoryGraph, T)

Returns the set of triangles that `T` was replaced by in the triangulation.
"""
SimpleGraphs.out_neighbors(G::HistoryGraph, T) = graph(G).N[T] # The version in SimpleGraphs is allocating...

"""
    SimpleGraphs.in_neighbors(G::HistoryGraph, T)

Returns the triangles that were replaced to form the triangle `T`.
"""
SimpleGraphs.in_neighbors(G::HistoryGraph, T) = SimpleGraphs.in_neighbors(graph(G), T) # Need the allocating version for find_root. Only place it's used anyway.

"""
    SimpleGraphs.in_deg(G::HistoryGraph, T)

Returns the number of triangles that were replaced to form the triangle `T`.
"""
SimpleGraphs.in_deg(G::HistoryGraph, T) = SimpleGraphs.in_deg(graph(G), T)

"""
    SimpleGraphs.out_deg(G::HistoryGraph, T)

Returns the number of triangles that `T` was replaced by.
"""
SimpleGraphs.out_deg(G::HistoryGraph, T) = SimpleGraphs.out_deg(graph(G), T)

"""
    add_triangle!(D::HistoryGraph, T)

Adds the triangles in `T` into the graph `D`. Checks are made for 
circular shifts of the vertices of `T` to avoid duplicates.
"""
function add_triangle!(D::HistoryGraph, T)
    has(graph(D), T) && return nothing
    has(graph(D), shift_triangle_1(T)) && return nothing
    has(graph(D), shift_triangle_2(T)) && return nothing
    add!(graph(D), T)
    return nothing
end
@doc (@doc add_triangle!(::HistoryGraph, ::Any))
function add_triangle!(D::HistoryGraph, T::Vararg{Tri,N}) where {Tri,N}
    for i in 1:N
        add_triangle!(D, T[i])
    end
end

"""
    add_edge!(G::HistoryGraph, T, V)
    add_edge!(G::HistoryGraph, T, V::Vararg{Tri, N}) where {Tri, N}

Adds a directed edge from `T` to the triangle(s) in `V`, making checks for circular shifts 
in the vertices of `T` and `V`.
"""
function add_edge!(G::HistoryGraph, T, V)
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
@doc (@doc add_edge!(::HistoryGraph, ::Any, ::Any))
function add_edge!(G::HistoryGraph, T, V::Vararg{Tri,N}) where {Tri,N}
    for V in V
        add_edge!(G, T, V)
    end
    return nothing
end

############################################
##
## TRIANGULATION STRUCTURES
##
############################################
"""
    AbstractTriangulation{T, E, Ts, Es, Ps, I}

Abstract type for a triangulation. The parametric types represent:

- `T`: The triangle type.
- `E`: The edge type.
- `Ts`: The triangle collection type.
- `Es`: The edge collection type.
- `Ps`: The point collection type.
- `I`: The integer type.

See also [`AbstractUnconstrainedTriangulation`](@ref), 
[`AbstractConstrainedTriangulation`](@ref), and 
[`AbstractWeightedTriangulation`](@ref).
"""
abstract type AbstractTriangulation{T,E,Ts,Es,Ps,I} end
adjacent(DT::AbstractTriangulation) = DT.adjacent
adjacent2vertex(DT::AbstractTriangulation) = DT.adjacent2vertex
graph(DT::AbstractTriangulation) = DT.graph
pointlocation(DT::AbstractTriangulation) = DT.pointlocation
triangles(DT::AbstractTriangulation) = DT.triangles
points(DT::AbstractTriangulation) = DT.points
edges(DT::AbstractTriangulation) = edges(adjacent(DT))
adjacent(DT::AbstractTriangulation, uv) = get_edge(adjacent(DT), uv)
adjacent(DT::AbstractTriangulation, u, v) = get_edge(adjacent(DT), u, v)
adjacent2vertex(DT::AbstractTriangulation, u) = get_edge(adjacent2vertex(DT), u)
neighbours(DT::AbstractTriangulation, u) = get_neighbour(graph(DT), u)
get_point(DT::AbstractTriangulation, i) = get_point(points(DT), i)
num_triangles(DT::AbstractTriangulation) = length(triangles(DT))
num_points(DT::AbstractTriangulation) = length(points(DT))
num_edges(DT::AbstractTriangulation) = num_triangles(DT) + num_points(DT) - 1 # Euler's formula

"""
    AbstractUnconstrainedTriangulation{T, E, Ts, Es, Ps, I} <: AbstractTriangulation{T, E, P, Ts, Es, Ps, I}

Abstract type for an unconstrained triangulation. See also [`AbstractTriangulation`](@ref)
and [`AbstractConstrainedTriangulation`](@ref).
"""
abstract type AbstractUnconstrainedTriangulation{T,E,Ts,Es,Ps,I} <: AbstractTriangulation{T,E,Ts,Es,Ps,I} end

"""
    is_delaunay(DT::AbstractUnconstrainedTriangulation)

Tests if the given triangulation `DT` is Delaunay. This is done by identifying 
if all edges are legal using `islegal`, noting that a legal triangulation is 
necessarily Delaunay.
"""
function is_delaunay(DT::AbstractUnconstrainedTriangulation{T,E,Ts,Es,Ps,I}) where {T,E,Ts,Es,Ps,I}
    adj = adjacent(DT)
    pts = points(DT)
    tri_edges = edges(adj)
    for (i, j) in tri_edges
        k = get_edge(adj, i, j) # The convex hull's edges are always included. 
        ℓ = get_edge(adj, j, i)
        k ≠ I(BoundaryIndex) && ℓ ≠ I(BoundaryIndex) && !islegal(i, j, adj, pts) && return false
    end
    return true
end

"""
    AbstractConstrainedTriangulation{T, E, Ts, Es, Ps, I} <: AbstractTriangulation{T, E, Ts, Es, Ps, I}


Abstract type for a constrained triangulation. See also [`AbstractTriangulation`](@ref) and 
[`AbstractUnconstrainedTriangulation`](@ref).
"""
abstract type AbstractConstrainedTriangulation{T,E,Ts,Es,Ps,I} <: AbstractTriangulation{T,E,Ts,Es,Ps,I} end

"""
    AbstractWeightedTriangulation{T, E, P, Ts, Es, Ps, I, W} <: AbstractTriangulation{T, E, P, Ts, Es, Ps, I}

Abstract type for a weighted triangulation. The extra type `W` represents the type used to represent the 
weights. See also [`AbstractTriangulation`](@ref).
"""
abstract type AbstractWeightedTriangulation{T,E,Ts,Es,Ps,I,W} <: AbstractTriangulation{T,E,Ts,Es,Ps,I} end

"""
    UnconstrainedTriangulation{T, E, Ts, Es, Ps, I, PL} <: AbstractUnconstrainedTriangulation{T, E, Ts, Es, Ps, I}

Struct for an unconstrained Delaunay triangulation.

- `adjacent::Adjacent{I, E}`: The adjacent map. See also [`Adjacent`](@ref).
- `adjacent2vertex::Adjacent2Vertex{I, Es, E}`: The adjacent-to-vertex map. See also [`Adjacent2Vertex`](@ref).
- `graph::DelaunayGraph{I}`: The graph representation of the triangulation. See also [`DelaunayGraph`](@ref).
- `pointlocation::PL`: The struct used for point location in the triangulation. For a randomised incremental method, this will likely be a `HistoryGraph{T}`. See also [`HistoryGraph`](@ref).
- `triangles::Ts`: The collection of triangles that together comprise the triangulation. 
- `points::Ps`: The vertex set of the triangulation.

"""
struct UnconstrainedTriangulation{T,E,Ts,Es,Ps,I,PL} <: AbstractUnconstrainedTriangulation{T,E,Ts,Es,Ps,I}
    adjacent::Adjacent{I,E}
    adjacent2vertex::Adjacent2Vertex{I,Es,E}
    graph::DelaunayGraph{I}
    pointlocation::PL
    triangles::Ts
    points::Ps
    function UnconstrainedTriangulation(pts::Ps;
        IntegerType::Type{I}=Int64,
        TriangleType::Type{T}=NTuple{3,IntegerType},
        EdgeType::Type{E}=NTuple{2,IntegerType},
        TrianglesType::Type{Ts}=Set{TriangleType},
        EdgesType::Type{Es}=Set{EdgeType},
        method=:berg) where {I,T,E,Ts,Es,Ps}
        if method == :berg
            # The initial data structures 
            root = construct_triangle(T, I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex))
            tris = construct_triangles(Ts)
            push!(tris, root)
            HG = HistoryGraph{T}()
            adj = Adjacent{I,E}()
            adj2v = Adjacent2Vertex{I,Es,E}()
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
            # Done 
            return new{T,E,Ts,Es,Ps,I,HistoryGraph{T}}(adj, adj2v, DG, HG, tris, pts)
        end
    end
end

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

############################################
##
## TRIANGULATION OPERATIONS
##
############################################
"""
    locate_triangle(HG::HistoryGraph, pts, r, init=find_root(HG; method=:rng))

Given the point location data structure `HG` and a set of `pts`, finds the triangle in 
the current triangulation such that the `r`th point of `pts` is in its interior. The point 
location starts at `init`.

The function is recursive, and returns a tuple `(tri, flag)`:
- `tri`: This is the triangle that `p` is in.
- `flag`: If `flag == 0`, then `p` is on an edge of `tri`. Otherwise, it is in the open interior.
"""
function locate_triangle(HG::HistoryGraph, pts, r::I, init=find_root(HG; method=:rng)) where {I}
    if out_deg(HG, init) == 0
        return init, isintriangle(init, pts, r)
    end
    out = out_neighbors(HG, init)
    for T in out
        isin = isintriangle(T, pts, r)
        (isin == I(0) || isin == I(1)) && return locate_triangle(HG, pts, r, T)
    end
    throw("Failed to find triangle.")
end

"""
    split_triangle!(T, HG, adj, adj2v, DG, i, j, k, ℓ, r)

Given a triangulation `T`, adds the `r`th point of the point set into the triangulation, assumed 
to be on the edge `(i, j)` of the triangulation, by splitting the triangles `(i, j, k)` and `(j, ℓ, i)` 
both into two.

# Arguments 
- `T`: The current set of triangles defining the triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
- `i, j`: The edge that the `r`th point is on.
- `k`: `(i, j, k)` is positively oriented.
- `ℓ`: `(j, i, ℓ`)` is positively oriented. 
- `r`: The point being added, assumed to be on the edge `(i, j)`.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function split_triangle!(T::Ts, HG, adj, adj2v, DG, i, j, k, ℓ, r) where {Ts}
    # The triangles 
    V = triangle_type(Ts)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    Tⱼᵢₗ = construct_triangle(V, j, i, ℓ)
    Tᵢᵣₖ = construct_triangle(V, i, r, k)
    Tᵣⱼₖ = construct_triangle(V, r, j, k)
    Tᵣᵢₗ = construct_triangle(V, r, i, ℓ)
    Tⱼᵣₗ = construct_triangle(V, j, r, ℓ)
    # Delete the old triangles 
    delete_triangle!(T, Tᵢⱼₖ, Tⱼᵢₗ)
    delete_edge!(adj, i, j)
    # Add the new triangles 
    add_triangle!(T, Tᵢᵣₖ, Tᵣⱼₖ, Tᵣᵢₗ, Tⱼᵣₗ)
    add_triangle!(adj, Tᵢᵣₖ, Tᵣⱼₖ, Tᵣᵢₗ, Tⱼᵣₗ)
    update_after_split!(adj2v, i, j, k, ℓ, r)
    add_triangle!(HG, Tᵢᵣₖ, Tᵣⱼₖ, Tᵣᵢₗ, Tⱼᵣₗ)
    add_edge!(HG, Tᵢⱼₖ, Tᵢᵣₖ, Tᵣⱼₖ)
    add_edge!(HG, Tⱼᵢₗ, Tᵣᵢₗ, Tⱼᵣₗ)
    # Update the graph 
    add_point!(DG, r)
    add_neighbour!(DG, r, i, j, k, ℓ)
    delete_neighbour!(DG, i, j)
    return nothing
end

"""
    flip_edge!(T, HG, adj, adj2v, DG, i, j, k, r)

Performs an edge flip, flipping the edge `(i, j)` into the edge `(k, r)`.

# Arguments
- `T`: The current set of triangles defining the triangulation.
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
function flip_edge!(T::Ts, HG, adj, adj2v, DG, i, j, k, r) where {Ts}
    # The old triangles
    V = triangle_type(Ts)
    Tᵢₖⱼ = construct_triangle(V, i, k, j)
    Tᵢⱼᵣ = construct_triangle(V, i, j, r)
    delete_triangle!(T, Tᵢₖⱼ, Tᵢⱼᵣ)
    delete_edge!(adj, i, j)
    delete_neighbour!(DG, i, j) #delete_neighbour!(DG, j, i)
    # The new triangles 
    Tᵣₖⱼ = construct_triangle(V, r, k, j)
    Tᵣᵢₖ = construct_triangle(V, r, i, k)
    # Add the new triangles to the data structure
    add_triangle!(T, Tᵣₖⱼ, Tᵣᵢₖ)
    add_triangle!(HG, Tᵣₖⱼ, Tᵣᵢₖ)
    add_triangle!(adj, Tᵣₖⱼ, Tᵣᵢₖ)
    update_after_flip!(adj2v, i, j, k, r)
    # Connect the new triangles to the replaced triangles in the DAG
    add_edge!(HG, Tᵢₖⱼ, Tᵣₖⱼ, Tᵣᵢₖ)
    add_edge!(HG, Tᵢⱼᵣ, Tᵣₖⱼ, Tᵣᵢₖ)
    # Add the new neighbours 
    add_neighbour!(DG, r, k) # add_neighbour!(DG, k, r)
    return nothing
end

"""
    legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    
Legalises the edge `(i, j)` if it is illegal.

# Arguments 
- `T`: The current set of triangles defining the triangulation.
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
function legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    if !islegal(i, j, adj, pts)
        e = get_edge(adj, j, i)
        flip_edge!(T, HG, adj, adj2v, DG, i, j, e, r)
        legalise_edge!(T, HG, adj, adj2v, DG, i, e, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, e, j, r, pts)
    end
    return nothing
end

"""
    remove_bounding_triangle!(DT::AbstractUnconstrainedTriangulation)
    remove_bounding_triangle!(T, adj, adj2v, DG)

Remove the bounding triangle from the triangulation, where 

- `T`: The current set of triangles defining the triangulation.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.

These structures are updated in-place.
"""
function remove_bounding_triangle!(T::Ts, adj, adj2v, DG::DelaunayGraph{I}) where {I,Ts}
    V = triangle_type(Ts)
    for w in indices(BoundingTriangle)
        w = I(w)
        neighbours = get_edge(adj2v, w)
        for (u, v) in neighbours # (u, v, w) is a triangle..
            delete_edge!(adj, w, u; protect_boundary=false)
            delete_edge!(adj, w, v; protect_boundary=false)
            delete_edge!(adj2v, u, v, w)
            delete_edge!(adj2v, v, w, u)
            if u ≥ I(FirstPointIndex) && v ≥ I(FirstPointIndex) # This can only be a boundary edge
                add_edge!(adj2v, I(BoundaryIndex), u, v)
                add_edge!(adj, u, v, I(BoundaryIndex))
            end
            delete_triangle!(T, construct_triangle(V, u, v, w))
        end
        delete_point!(DG, w)
        delete_point!(adj2v, w)
    end
    return nothing
end
@doc (@doc remove_bounding_triangle!(::Any, ::Any, ::Any, ::DelaunayGraph{I}) where {I})
function remove_bounding_triangle!(DT::AbstractUnconstrainedTriangulation)
    remove_bounding_triangle!(triangles(DT), adjacent(DT), adjacent2vertex(DT), graph(DT))
    return nothing
end

"""
    add_point!(T, HG::HistoryGraph, adj, adj2v, DG, root, pts, r)

Adds `get_point(pts, r)` to the triangulation. 
"""
function add_point!(T, HG::HistoryGraph, adj, adj2v, DG, root, pts, r)
    Tᵢⱼₖ, interior_flag = locate_triangle(HG, pts, r, root)
    if interior_flag == 1
        i, j, k = Tᵢⱼₖ
        split_triangle!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
        legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, j, k, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, k, i, r, pts)
    else
        i, j = find_edge(Tᵢⱼₖ, pts, r)
        k = get_edge(adj, i, j)
        ℓ = get_edge(adj, j, i)
        split_triangle!(T, HG, adj, adj2v, DG, i, j, k, ℓ, r)
        legalise_edge!(T, HG, adj, adj2v, DG, i, ℓ, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, ℓ, j, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, j, k, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, k, i, r, pts)
    end
end

"""
    add_point!(DT::AbstractTriangulation, r::Integer)

Adds the `r`th point of the point set in `DT` into the triangulation.
"""
function add_point!(DT::AbstractTriangulation{T,E,Ts,Es,Ps,I}, r::I) where {T,E,Ts,Es,Ps,I}
    add_point!(triangles(DT), pointlocation(DT), adjacent(DT),
        adjacent2vertex(DT), graph(DT),
        construct_triangle(T, I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex)),
        points(DT), I(r))
    return nothing
end

"""
    add_point!(DT::AbstractTriangulation, p)

Adds the point `p` into the triangulation.
"""
function add_point!(DT::AbstractTriangulation{T,E,Ts,Es,Ps,I}, p) where {T,E,Ts,Es,Ps,I}
    add_point!(points(DT), p)
    add_point!(DT, I(lastindex(points(DT))))
    return nothing
end

############################################
##
## TRIANGULATION METHODS
##
############################################

"""
    triangulate!(DT::AbstractTriangulation; randomise=true, trim=true)

Completes the Delaunay triangulation `DT`. Use the keyword argument 
`randomise=true` to randomise the insertion order, and `trim=true` 
to remove the bounding triangle.
"""
function triangulate!(DT::AbstractTriangulation{T,E,Ts,Es,Ps,I}; randomise=true, trim=true) where {T,E,Ts,Es,Ps,I}
    pt_order = randomise ? shuffle(eachindex(points(DT))) : eachindex(points(DT))
    map(pt_order) do r
        add_point!(DT, I(r))
    end
    trim && remove_bounding_triangle!(DT)
    return nothing
end

"""
    triangulate(pts; randomise=true, trim=true, method=:berg,
        IntegerType::Type{I}=Int64, TriangleType::Type{T}=NTuple{3,I}, EdgeType::Type{E}=NTuple{2,I},
        TrianglesType::Type{Ts}=Set{T}, EdgesType::Type{Es}=Set{E}) where {I,T,E,Ts,Es}

Computes the unconstrained Delaunay triangulation of `pts`.

# Keyword Arguments 
- `randomise=true`: Randomise the insertion order.
- `trim=true`: Remove the bounding triangle.
- `method=:berg`: The method to use. See `?available_methods` for a description of the available methods.
- `IntegerType::Type{I}=Int64`: The integer type.
- `TriangleType::Type{T}=NTuple{3,I}`: The triangle type. 
- `EdgeType::Type{E}=NTuple{2,I}`: The edge type.
- `TrianglesType::Type{Ts}=Set{T}`: The type for a collection of triangles.
- `EdgesType::Type{Es}=Set{E}`: The type for a collection of edges.
"""
function triangulate(pts; randomise=true, trim=true, method=:berg,
    IntegerType::Type{I}=Int64,
    TriangleType::Type{T}=NTuple{3,IntegerType},
    EdgeType::Type{E}=NTuple{2,IntegerType},
    TrianglesType::Type{Ts}=Set{TriangleType},
    EdgesType::Type{Es}=Set{EdgeType}) where {I,T,E,Ts,Es}
    DT = UnconstrainedTriangulation(pts; IntegerType, TriangleType, EdgeType, TrianglesType, EdgesType, method)
    triangulate!(DT; randomise, trim)
    return DT
end