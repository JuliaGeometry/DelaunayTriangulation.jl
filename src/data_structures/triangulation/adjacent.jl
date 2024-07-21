"""
    Adjacent{IntegerType, EdgeType} 

Struct for storing adjacency relationships for a triangulation.

# Fields 
- `adjacent::Dict{EdgeType, IntegerType}`

The map taking edges `(u, v)` to `w` such that `(u, v, w)` is a positively oriented triangle in the underlying triangulation.

# Constructors 
    Adjacent{IntegerType, EdgeType}()
    Adjacent(adjacent::Dict{EdgeType, IntegerType})
"""
struct Adjacent{I, E}
    adjacent::Dict{E, I}
end
Adjacent{I, E}() where {I, E} = Adjacent{I, E}(Dict{E, I}())
Base.:(==)(adj::Adjacent, adj2::Adjacent) = get_adjacent(adj) == get_adjacent(adj2)
function Base.show(io::IO, m::MIME"text/plain", adj::Adjacent{I, E}) where {I, E}
    println(io, "Adjacent{$I, $E}, with map:")
    show(io, m, get_adjacent(adj))
end
Base.sizehint!(adj::Adjacent, n) = sizehint!(get_adjacent(adj), n)

"""
    get_adjacent(adj::Adjacent) -> Dict 

Returns the `adjacent` map of `adj`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> d = Dict((1, 2) => 3, (2, 3) => 1, (3, 1) => 2);

julia> adj = DelaunayTriangulation.Adjacent(d)
Adjacent{Int64, Tuple{Int64, Int64}}, with map:
Dict{Tuple{Int64, Int64}, Int64} with 3 entries:
  (1, 2) => 3
  (3, 1) => 2
  (2, 3) => 1

julia> get_adjacent(adj)
Dict{Tuple{Int64, Int64}, Int64} with 3 entries:
  (1, 2) => 3
  (3, 1) => 2
  (2, 3) => 1

julia> get_adjacent(adj) == d
true
```
"""
get_adjacent(adj::Adjacent) = adj.adjacent


"""
    get_adjacent(adj::Adjacent{I, E}, uv::E) -> Vertex
    get_adjacent(adj::Adjacent{I, E}, u, v) -> Vertex

Returns the vertex `w` such that `(u, v, w)` is a positively oriented triangle in the
underlying triangulation, or `∅` if no such triangle exists.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj = DelaunayTriangulation.Adjacent(Dict((1, 2) => 3, (2, 3) => 1, (3, 1) => 2, (4, 5) => -1))
Adjacent{Int64, Tuple{Int64, Int64}}, with map:
Dict{Tuple{Int64, Int64}, Int64} with 4 entries:
  (4, 5) => -1
  (1, 2) => 3
  (3, 1) => 2
  (2, 3) => 1

julia> get_adjacent(adj, 4, 5)
-1

julia> get_adjacent(adj, (3, 1))
2

julia> get_adjacent(adj, (1, 2))
3

julia> get_adjacent(adj, 17, 5)
0

julia> get_adjacent(adj, (1, 6))
0
```
"""
function get_adjacent(adj::Adjacent{I, E}, uv::E) where {I, E}
    dict = get_adjacent(adj)
    return get(dict, uv, I(∅))
end
function get_adjacent(adj::Adjacent{I, E}, u, v) where {I, E}
    e = construct_edge(E, u, v)
    return get_adjacent(adj, e)
end


"""
    add_adjacent!(adj::Adjacent, uv, w)
    add_adjacent!(adj::Adjacent, u, v, w)

Adds the adjacency relationship `(u, v, w)` to `adj`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj = DelaunayTriangulation.Adjacent{Int64, NTuple{2, Int64}}();

julia> DelaunayTriangulation.add_adjacent!(adj, 1, 2, 3)
Adjacent{Int64, Tuple{Int64, Int64}}, with map:
Dict{Tuple{Int64, Int64}, Int64} with 1 entry:
  (1, 2) => 3

julia> DelaunayTriangulation.add_adjacent!(adj, (2, 3), 1)
Adjacent{Int64, Tuple{Int64, Int64}}, with map:
Dict{Tuple{Int64, Int64}, Int64} with 2 entries:
  (1, 2) => 3
  (2, 3) => 1

julia> DelaunayTriangulation.add_adjacent!(adj, 3, 1, 2)
Adjacent{Int64, Tuple{Int64, Int64}}, with map:
Dict{Tuple{Int64, Int64}, Int64} with 3 entries:
  (1, 2) => 3
  (3, 1) => 2
  (2, 3) => 1
```
"""
function add_adjacent!(adj::Adjacent, uv, w)
    dict = get_adjacent(adj)
    dict[uv] = w
    return adj
end
function add_adjacent!(adj::Adjacent{I, E}, u, v, w) where {I, E}
    e = construct_edge(E, u, v)
    return add_adjacent!(adj, e, w)
end

"""
    delete_adjacent!(adj::Adjacent, uv)
    delete_adjacent!(adj::Adjacent, u, v)

Deletes the edge `uv` from `adj`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj = DelaunayTriangulation.Adjacent(Dict((2, 7) => 6, (7, 6) => 2, (6, 2) => 2, (17, 3) => -1, (-1, 5) => 17, (5, 17) => -1));

julia> DelaunayTriangulation.delete_adjacent!(adj, 2, 7)
Adjacent{Int64, Tuple{Int64, Int64}}, with map:
Dict{Tuple{Int64, Int64}, Int64} with 5 entries:
  (-1, 5) => 17
  (17, 3) => -1
  (6, 2)  => 2
  (5, 17) => -1
  (7, 6)  => 2

julia> DelaunayTriangulation.delete_adjacent!(adj, (6, 2))
Adjacent{Int64, Tuple{Int64, Int64}}, with map:
Dict{Tuple{Int64, Int64}, Int64} with 4 entries:
  (-1, 5) => 17
  (17, 3) => -1
  (5, 17) => -1
  (7, 6)  => 2

julia> DelaunayTriangulation.delete_adjacent!(adj, 5, 17)
Adjacent{Int64, Tuple{Int64, Int64}}, with map:
Dict{Tuple{Int64, Int64}, Int64} with 3 entries:
  (-1, 5) => 17
  (17, 3) => -1
  (7, 6)  => 2
```
"""
function delete_adjacent!(adj::Adjacent, uv)
    dict = get_adjacent(adj)
    delete!(dict, uv)
    return adj
end
function delete_adjacent!(adj::Adjacent{I, E}, u, v) where {I, E}
    e = construct_edge(E, u, v)
    return delete_adjacent!(adj, e)
end

"""
    add_triangle!(adj::Adjacent, u, v, w)
    add_triangle!(adj::Adjacent, T)

Adds the adjacency relationships defined from the triangle `T = (u, v, w)` to `adj`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj = DelaunayTriangulation.Adjacent{Int32, NTuple{2, Int32}}();

julia> add_triangle!(adj, 1, 2, 3)
Adjacent{Int32, Tuple{Int32, Int32}}, with map:
Dict{Tuple{Int32, Int32}, Int32} with 3 entries:
  (1, 2) => 3
  (3, 1) => 2
  (2, 3) => 1

julia> add_triangle!(adj, 6, -1, 7)
Adjacent{Int32, Tuple{Int32, Int32}}, with map:
Dict{Tuple{Int32, Int32}, Int32} with 6 entries:
  (1, 2)  => 3
  (3, 1)  => 2
  (6, -1) => 7
  (-1, 7) => 6
  (2, 3)  => 1
  (7, 6)  => -1
```
"""
function add_triangle!(adj::Adjacent, u::Integer, v::Integer, w::Integer) # method ambiguity
    add_adjacent!(adj, u, v, w)
    add_adjacent!(adj, v, w, u)
    add_adjacent!(adj, w, u, v)
    return adj
end
add_triangle!(adj::Adjacent, T) = add_triangle!(adj, geti(T), getj(T), getk(T))

"""
    delete_triangle!(adj::Adjacent, u, v, w)
    delete_triangle!(adj::Adjacent, T)

Deletes the adjacency relationships defined from the triangle `T = (u, v, w)` from `adj`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> adj = DelaunayTriangulation.Adjacent{Int32, NTuple{2, Int32}}();

julia> add_triangle!(adj, 1, 6, 7);

julia> add_triangle!(adj, 17, 3, 5);

julia> adj
Adjacent{Int32, Tuple{Int32, Int32}}, with map:
Dict{Tuple{Int32, Int32}, Int32} with 6 entries:
  (17, 3) => 5
  (1, 6)  => 7
  (6, 7)  => 1
  (7, 1)  => 6
  (5, 17) => 3
  (3, 5)  => 17

julia> delete_triangle!(adj, 3, 5, 17)
Adjacent{Int32, Tuple{Int32, Int32}}, with map:
Dict{Tuple{Int32, Int32}, Int32} with 3 entries:
  (1, 6) => 7
  (6, 7) => 1
  (7, 1) => 6

julia> delete_triangle!(adj, 7, 1, 6)
Adjacent{Int32, Tuple{Int32, Int32}}, with map:
Dict{Tuple{Int32, Int32}, Int32}()
```
"""
function delete_triangle!(adj::Adjacent, u::Integer, v::Integer, w::Integer) # method ambiguity
    for (i, j) in triangle_edges(u, v, w)
        delete_adjacent!(adj, i, j)
    end
    return adj
end
delete_triangle!(adj::Adjacent, T) = delete_triangle!(adj, geti(T), getj(T), getk(T))

function Base.empty!(adj::Adjacent)
    dict = get_adjacent(adj)
    empty!(dict)
    return adj
end
