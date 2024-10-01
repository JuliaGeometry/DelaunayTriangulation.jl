"""
    Adjacent2Vertex{IntegerType, EdgesType}

Struct for connectivity information about edges relative to vertices for a triangulation.

# Fields 
- `adjacent2vertex::Dict{IntegerType, EdgesType}`

The map taking `w` to the set of all `(u, v)` such that `(u, v, w)` is a positively oriented triangle in the underlying triangle.

# Constructors 
    Adjacent2Vertex{IntegerType, EdgesType}()
    Adjacent2Vertex(adj2v::Dict{IntegerType, EdgesType})
"""
struct Adjacent2Vertex{I, Es}
    adjacent2vertex::Dict{I, Es}
    Adjacent2Vertex(adj2v::Dict{I, Es}) where {I, Es} = new{I, Es}(adj2v)
end
Adjacent2Vertex{I, Es}() where {I, Es} = Adjacent2Vertex(Dict{I, Es}())
Base.:(==)(adj2v::Adjacent2Vertex, adj2v2::Adjacent2Vertex) = get_adjacent2vertex(adj2v) == get_adjacent2vertex(adj2v2)
function Base.show(io::IO, m::MIME"text/plain", adj2v::Adjacent2Vertex{I, Es}) where {I, Es}
    println(io, "Adjacent2Vertex{", I, ", ", Es, "} with map:")
    show(io, m, get_adjacent2vertex(adj2v))
end
Base.sizehint!(adj2v::Adjacent2Vertex, n) = Base.sizehint!(get_adjacent2vertex(adj2v), n)

Base.copy(adj2v::Adjacent2Vertex) = Adjacent2Vertex(copy(get_adjacent2vertex(adj2v)))

"""
    get_adjacent2vertex(adj2v::Adjacent2Vertex) -> Dict

Returns the `adjacent2vertex` map of `adj2v`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> e1 = Set(((1, 2), (5, 3), (7, 8)));

julia> e2 = Set(((2, 3), (13, 5), (-1, 7)));

julia> d = Dict(9 => e1, 6 => e2);

julia> adj2v = DelaunayTriangulation.Adjacent2Vertex(d)
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 2 entries:
  6 => Set([(13, 5), (-1, 7), (2, 3)])
  9 => Set([(1, 2), (7, 8), (5, 3)])

julia> get_adjacent2vertex(adj2v)
Dict{Int64, Set{Tuple{Int64, Int64}}} with 2 entries:
  6 => Set([(13, 5), (-1, 7), (2, 3)])
  9 => Set([(1, 2), (7, 8), (5, 3)])

julia> get_adjacent2vertex(adj2v) == d
true
```
"""
get_adjacent2vertex(adj2v::Adjacent2Vertex) = adj2v.adjacent2vertex

"""
    get_adjacent2vertex(adj2v::Adjacent2Vertex, w) -> Edges 

Returns the set of edges `E` such that `(u, v, w)` is a positively oriented triangle in the
underlying triangulation for each `(u, v) ∈ E`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj2v = DelaunayTriangulation.Adjacent2Vertex(Dict(1 => Set(((2, 3), (5, 7), (8, 9))), 5 => Set(((1, 2), (7, 9), (8, 3)))))
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 2 entries:
  5 => Set([(1, 2), (8, 3), (7, 9)])
  1 => Set([(8, 9), (5, 7), (2, 3)])

julia> get_adjacent2vertex(adj2v, 1)
Set{Tuple{Int64, Int64}} with 3 elements:
  (8, 9)
  (5, 7)
  (2, 3)

julia> get_adjacent2vertex(adj2v, 5)
Set{Tuple{Int64, Int64}} with 3 elements:
  (1, 2)
  (8, 3)
  (7, 9)
```
"""
function get_adjacent2vertex(adj2v::Adjacent2Vertex, w)
    dict = get_adjacent2vertex(adj2v)
    return dict[w]
end

"""
    add_adjacent2vertex!(adj2v::Adjacent2Vertex, w, uv)
    add_adjacent2vertex!(adj2v::Adjacent2Vertex, w, u, v)

Adds the edge `uv` to the set of edges `E` such that `(u, v, w)` is a positively oriented triangle in the
underlying triangulation for each `(u, v) ∈ E`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj2v = DelaunayTriangulation.Adjacent2Vertex{Int64, Set{NTuple{2, Int64}}}()
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}}()

julia> DelaunayTriangulation.add_adjacent2vertex!(adj2v, 1, (2, 3))
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 1 entry:
  1 => Set([(2, 3)])

julia> DelaunayTriangulation.add_adjacent2vertex!(adj2v, 1, 5, 7)
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 1 entry:
  1 => Set([(5, 7), (2, 3)])

julia> DelaunayTriangulation.add_adjacent2vertex!(adj2v, 17, (5, -1))
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 2 entries:
  17 => Set([(5, -1)])
  1  => Set([(5, 7), (2, 3)])
```
"""
function add_adjacent2vertex!(adj2v::Adjacent2Vertex{I, Es}, w, uv) where {I, Es}
    dict = get_adjacent2vertex(adj2v)
    existing_edges = get!(Es, dict, w)
    add_edge!(existing_edges, uv)
    return adj2v
end
function add_adjacent2vertex!(adj2v::Adjacent2Vertex{I, Es}, w, u, v) where {I, Es}
    E = edge_type(Es)
    uv = construct_edge(E, u, v)
    return add_adjacent2vertex!(adj2v, w, uv)
end

"""
    delete_adjacent2vertex!(adj2v::Adjacent2Vertex, w, uv)
    delete_adjacent2vertex!(adj2v::Adjacent2Vertex, w, u, v)

Deletes the edge `uv` from the set of edges `E` such that `(u, v, w)` is a positively oriented triangle in the
underlying triangulation for each `(u, v) ∈ E`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj2v = DelaunayTriangulation.Adjacent2Vertex(Dict(1 => Set(((2, 3), (5, 7), (8, 9))), 5 => Set(((1, 2), (7, 9), (8, 3)))))
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 2 entries:
  5 => Set([(1, 2), (8, 3), (7, 9)])
  1 => Set([(8, 9), (5, 7), (2, 3)])

julia> DelaunayTriangulation.delete_adjacent2vertex!(adj2v, 5, 8, 3)
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 2 entries:
  5 => Set([(1, 2), (7, 9)])
  1 => Set([(8, 9), (5, 7), (2, 3)])

julia> DelaunayTriangulation.delete_adjacent2vertex!(adj2v, 1, (2, 3))
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 2 entries:
  5 => Set([(1, 2), (7, 9)])
  1 => Set([(8, 9), (5, 7)])
```
"""
function delete_adjacent2vertex!(adj2v::Adjacent2Vertex, w, uv)
    existing_edges = get_adjacent2vertex(adj2v, w)
    delete_edge!(existing_edges, uv)
    return adj2v
end
function delete_adjacent2vertex!(adj2v::Adjacent2Vertex{I, Es}, w, u, v) where {I, Es}
    E = edge_type(Es)
    uv = construct_edge(E, u, v)
    return delete_adjacent2vertex!(adj2v, w, uv)
end

"""
    delete_adjacent2vertex!(adj2v::Adjacent2Vertex, w)

Deletes the vertex `w` from `adj2v`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj2v = DelaunayTriangulation.Adjacent2Vertex(Dict(1 => Set(((2, 3), (5, 7))), 5 => Set(((-1, 2), (2, 3)))))
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 2 entries:
  5 => Set([(-1, 2), (2, 3)])
  1 => Set([(5, 7), (2, 3)])

julia> DelaunayTriangulation.delete_adjacent2vertex!(adj2v, 1)
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 1 entry:
  5 => Set([(-1, 2), (2, 3)])

julia> DelaunayTriangulation.delete_adjacent2vertex!(adj2v, 5)
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}}()
```
"""
function delete_adjacent2vertex!(adj2v::Adjacent2Vertex, w)
    dict = get_adjacent2vertex(adj2v)
    delete!(dict, w)
    return adj2v
end

"""
    add_triangle!(adj2v::Adjacent2Vertex, u, v, w)
    add_triangle!(adj2v::Adjacent2Vertex, T)

Adds the relationships defined by the triangle `T = (u, v, w)` into `adj2v`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj2v = DelaunayTriangulation.Adjacent2Vertex{Int32, Set{NTuple{2, Int32}}}()
Adjacent2Vertex{Int32, Set{Tuple{Int32, Int32}}} with map:
Dict{Int32, Set{Tuple{Int32, Int32}}}()

julia> add_triangle!(adj2v, 17, 5, 8)
Adjacent2Vertex{Int32, Set{Tuple{Int32, Int32}}} with map:
Dict{Int32, Set{Tuple{Int32, Int32}}} with 3 entries:
  5  => Set([(8, 17)])
  8  => Set([(17, 5)])
  17 => Set([(5, 8)])

julia> add_triangle!(adj2v, 1, 5, 13)
Adjacent2Vertex{Int32, Set{Tuple{Int32, Int32}}} with map:
Dict{Int32, Set{Tuple{Int32, Int32}}} with 5 entries:
  5  => Set([(8, 17), (13, 1)])
  13 => Set([(1, 5)])
  8  => Set([(17, 5)])
  17 => Set([(5, 8)])
  1  => Set([(5, 13)])
```
"""
function add_triangle!(adj2v::Adjacent2Vertex, u::Integer, v::Integer, w::Integer)
    add_adjacent2vertex!(adj2v, u, v, w)
    add_adjacent2vertex!(adj2v, v, w, u)
    add_adjacent2vertex!(adj2v, w, u, v)
    return adj2v
end
add_triangle!(adj2v::Adjacent2Vertex, T) = add_triangle!(adj2v, geti(T), getj(T), getk(T))

"""
    delete_triangle!(adj2v::Adjacent2Vertex, u, v, w)
    delete_triangle!(adj2v::Adjacent2Vertex, T)

Deletes the relationships defined by the triangle `T =(u, v, w)` from `adj2v`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj2v = DelaunayTriangulation.Adjacent2Vertex{Int32, Set{NTuple{2, Int32}}}()
Adjacent2Vertex{Int32, Set{Tuple{Int32, Int32}}} with map:
Dict{Int32, Set{Tuple{Int32, Int32}}}()

julia> add_triangle!(adj2v, 1, 2, 3)
Adjacent2Vertex{Int32, Set{Tuple{Int32, Int32}}} with map:
Dict{Int32, Set{Tuple{Int32, Int32}}} with 3 entries:
  2 => Set([(3, 1)])
  3 => Set([(1, 2)])
  1 => Set([(2, 3)])

julia> add_triangle!(adj2v, 17, 5, 2)
Adjacent2Vertex{Int32, Set{Tuple{Int32, Int32}}} with map:
Dict{Int32, Set{Tuple{Int32, Int32}}} with 5 entries:
  5  => Set([(2, 17)])
  2  => Set([(3, 1), (17, 5)])
  17 => Set([(5, 2)])
  3  => Set([(1, 2)])
  1  => Set([(2, 3)])

julia> delete_triangle!(adj2v, 5, 2, 17)
Adjacent2Vertex{Int32, Set{Tuple{Int32, Int32}}} with map:
Dict{Int32, Set{Tuple{Int32, Int32}}} with 5 entries:
  5  => Set()
  2  => Set([(3, 1)])
  17 => Set()
  3  => Set([(1, 2)])
  1  => Set([(2, 3)])

julia> delete_triangle!(adj2v, 2, 3, 1)
Adjacent2Vertex{Int32, Set{Tuple{Int32, Int32}}} with map:
Dict{Int32, Set{Tuple{Int32, Int32}}} with 5 entries:
  5  => Set()
  2  => Set()
  17 => Set()
  3  => Set()
  1  => Set()
```
"""
function delete_triangle!(adj2v::Adjacent2Vertex, u::Integer, v::Integer, w::Integer)
    delete_adjacent2vertex!(adj2v, u, v, w)
    delete_adjacent2vertex!(adj2v, v, w, u)
    delete_adjacent2vertex!(adj2v, w, u, v)
    return adj2v
end
delete_triangle!(adj2v::Adjacent2Vertex, T) = delete_triangle!(adj2v, geti(T), getj(T), getk(T))

"""
    clear_empty_keys!(adj2v::Adjacent2Vertex)

Deletes all vertices `w` from `adj2v` such that `get_adjacent2vertex(adj2v, w)` is empty.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj2v = DelaunayTriangulation.Adjacent2Vertex{Int64, Set{NTuple{2, Int64}}}()
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}}()

julia> add_triangle!(adj2v, 1, 2, 3)
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 3 entries:
  2 => Set([(3, 1)])
  3 => Set([(1, 2)])
  1 => Set([(2, 3)])

julia> delete_triangle!(adj2v, 2, 3, 1)
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}} with 3 entries:
  2 => Set()
  3 => Set()
  1 => Set()

julia> DelaunayTriangulation.clear_empty_keys!(adj2v)
Adjacent2Vertex{Int64, Set{Tuple{Int64, Int64}}} with map:
Dict{Int64, Set{Tuple{Int64, Int64}}}()
```
"""
function clear_empty_keys!(adj2v::Adjacent2Vertex)
    dict = get_adjacent2vertex(adj2v)
    foreach(dict) do (w, S)
        isempty(S) && delete_adjacent2vertex!(adj2v, w)
    end
    return adj2v
end

function Base.empty!(adj2v::Adjacent2Vertex)
    dict = get_adjacent2vertex(adj2v)
    empty!(dict)
    return adj2v
end
