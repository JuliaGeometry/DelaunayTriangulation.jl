@doc """
    construct_edge(::Type{E}, i, j) where {E} -> E

Construct an edge of type `E` from vertices `i` and `j`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.construct_edge(NTuple{2,Int}, 2, 5)
(2, 5)

julia> DelaunayTriangulation.construct_edge(Vector{Int32}, 5, 15)
2-element Vector{Int32}:
  5
 15
```
"""
construct_edge
construct_edge(::Type{NTuple{2,I}}, i, j) where {I} = (I(i), I(j))
construct_edge(::Type{Vector{I}}, i, j) where {I} = I[i, j]
construct_edge(::Type{A}, i, j) where {I,A<:AbstractVector{I}} = A(I[i, j])

"""
    initial(e) -> Vertex

Get the initial vertex of `e`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> e = (1, 3);

julia> DelaunayTriangulation.initial(e)
1

julia> e = [2, 5];

julia> DelaunayTriangulation.initial(e)
2
```
"""
initial(e) = e[1]

"""
    terminal(e) -> Vertex

Get the terminal vertex of `e`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> e = (1, 7);

julia> DelaunayTriangulation.terminal(e)
7

julia> e = [2, 13];

julia> DelaunayTriangulation.terminal(e)
13
```
"""
terminal(e) = e[2]

"""
    edge_vertices(e) -> NTuple{2, Vertex}

Get the vertices of `e`

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> e = (1, 5);

julia> edge_vertices(e)
(1, 5)

julia> e = [23, 50];

julia> edge_vertices(e)
(23, 50)
```
"""
edge_vertices(e) = (initial(e), terminal(e))

"""
    reverse_edge(e) -> Edge

Get the edge with the vertices of `e` in reverse order.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> e = (17, 3);

julia> DelaunayTriangulation.reverse_edge(e)
(3, 17)

julia> e = [1, 2];

julia> DelaunayTriangulation.reverse_edge(e)
2-element Vector{Int64}:
 2
 1
```
"""
reverse_edge(e) = construct_edge(typeof(e), terminal(e), initial(e))

"""
    compare_unoriented_edges(u, v) -> Bool

Compare the unoriented edges `u` and `v`, i.e. compare the vertices of `u` and `v`
in any order.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> u = (1, 3);

julia> v = (5, 3);

julia> DelaunayTriangulation.compare_unoriented_edges(u, v)
false

julia> v = (1, 3);

julia> DelaunayTriangulation.compare_unoriented_edges(u, v)
true

julia> v = (3, 1);

julia> DelaunayTriangulation.compare_unoriented_edges(u, v)
true
```
"""
function compare_unoriented_edges(u, v)
  i, j = edge_vertices(u)
  k, ℓ = edge_vertices(v)
  return (i, j) == (k, ℓ) || (i, j) == (ℓ, k)
end

"""
    num_edges(E) -> Integer 

Get the number of edges in `E`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> e = [(1, 2), (3, 4), (1, 5)];

julia> num_edges(e)
3
```
"""
num_edges(E) = length(E)

"""
    edge_type(E) -> DataType

Get the type of edges in `E`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> e = Set(((1,2),(2,3),(17,5)))
Set{Tuple{Int64, Int64}} with 3 elements:
  (1, 2)
  (17, 5)
  (2, 3)

julia> DelaunayTriangulation.edge_type(e)
Tuple{Int64, Int64}

julia> e = [[1,2],[3,4],[17,3]]
3-element Vector{Vector{Int64}}:
 [1, 2]
 [3, 4]
 [17, 3]

julia> DelaunayTriangulation.edge_type(e)
Vector{Int64} (alias for Array{Int64, 1})
```
"""
edge_type(E) = eltype(E)

@doc """
    contains_edge(i, j, E) -> Bool
    contains_edge(e, E) -> Bool 

Check if `E` contains the edge `e = (i, j)`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> E = Set(((1,3),(17,3),(1,-1)))
Set{Tuple{Int64, Int64}} with 3 elements:
  (1, -1)
  (17, 3)
  (1, 3)

julia> DelaunayTriangulation.contains_edge((1, 2), E)
false

julia> DelaunayTriangulation.contains_edge((17, 3), E)
true

julia> DelaunayTriangulation.contains_edge(3, 17, E) # order
false

julia> E = [[1,2],[5,13],[-1,1]]
3-element Vector{Vector{Int64}}:
 [1, 2]
 [5, 13]
 [-1, 1]

julia> DelaunayTriangulation.contains_edge(1, 2, E)
true
```
"""
contains_edge
contains_edge(e, E) = e ∈ E
contains_edge(i, j, E) = contains_edge(construct_edge(edge_type(E), i, j), E)

"""
    contains_unoriented_edge(e, E) -> Bool
  
Check if `E` contains the unoriented edge `e`, i.e. check if `E` contains the edge `e` or `reverse_edge(e)`.
"""
contains_unoriented_edge(e, E) = contains_edge(e, E) || contains_edge(reverse_edge(e), E)

"""
    add_to_edges!(E, e)

Add the edge `e` to `E`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> E = Set(((1, 2),(3,5)))
Set{Tuple{Int64, Int64}} with 2 elements:
  (1, 2)
  (3, 5)

julia> DelaunayTriangulation.add_to_edges!(E, (1, 5))
Set{Tuple{Int64, Int64}} with 3 elements:
  (1, 2)
  (3, 5)
  (1, 5)
```
"""
add_to_edges!(E, e) = push!(E, e)

"""
    add_edge!(E, e...)

Add the edges `e...` to `E`.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> E = Set(((1,5),(17,10),(5,3)))
Set{Tuple{Int64, Int64}} with 3 elements:
  (5, 3)
  (17, 10)
  (1, 5)

julia> DelaunayTriangulation.add_edge!(E, (3, 2))

julia> E
Set{Tuple{Int64, Int64}} with 4 elements:
  (3, 2)
  (5, 3)
  (17, 10)
  (1, 5)

julia> DelaunayTriangulation.add_edge!(E, (1, -3), (5, 10), (1, -1))

julia> E
Set{Tuple{Int64, Int64}} with 7 elements:
  (3, 2)
  (5, 10)
  (1, -3)
  (1, -1)
  (5, 3)
  (17, 10)
  (1, 5)
```
"""
function add_edge!(E, e...)
  for v in e
    add_to_edges!(E, v)
  end
end

"""
    delete_from_edges!(E, e)

Delete the edge `e` from `E`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> E = Set(([1,2],[5,15],[17,10],[5,-1]))
Set{Vector{Int64}} with 4 elements:
  [17, 10]
  [5, 15]
  [5, -1]
  [1, 2]

julia> DelaunayTriangulation.delete_from_edges!(E, [5, 15])
Set{Vector{Int64}} with 3 elements:
  [17, 10]
  [5, -1]
  [1, 2]
```
"""
delete_from_edges!(E, e) = delete!(E, e)

"""
    delete_edge!(E, e...)

Delete the edges `e...` from `E`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> E = Set(([1,2],[10,15],[1,-1],[13,23],[1,5]))
Set{Vector{Int64}} with 5 elements:
  [10, 15]
  [1, 5]
  [1, 2]
  [1, -1]
  [13, 23]

julia> DelaunayTriangulation.delete_edge!(E, [10,15])

julia> E
Set{Vector{Int64}} with 4 elements:
  [1, 5]
  [1, 2]
  [1, -1]
  [13, 23]

julia> DelaunayTriangulation.delete_edge!(E, [1,5], [1, -1])

julia> E
Set{Vector{Int64}} with 2 elements:
  [1, 2]
  [13, 23]
```
"""
function delete_edge!(E, e...)
  for v in e
    delete_from_edges!(E, v)
  end
end

"""
    delete_unoriented_edge!(E, e)

Delete the unoriented edge `e` from `E`, i.e. delete both the edge `e` and `reverse_edge(e)`.
"""
delete_unoriented_edge!(E, e) = delete_edge!(E, e, reverse_edge(e))

@doc """
    each_edge(E) -> Iterator

Get an iterator over the edges in `E`.

# Examples 
```jldoctest 
julia> using DelaunayTriangulation

julia> E = Set(((1,2),(1,3),(2,-1)))
Set{Tuple{Int64, Int64}} with 3 elements:
  (1, 2)
  (1, 3)
  (2, -1)

julia> each_edge(E)
Set{Tuple{Int64, Int64}} with 3 elements:
  (1, 2)
  (1, 3)
  (2, -1)
```
"""
each_edge
each_edge(E) = E
each_edge(E::AbstractMatrix) = eachcol(E)

@doc """
    random_edge([rng], E) -> E

Get a random edge from `E`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation, StableRNGs

julia> E = Set(((1,2),(10,15),(23,20)))
Set{Tuple{Int64, Int64}} with 3 elements:
  (1, 2)
  (23, 20)
  (10, 15)

julia> rng = StableRNG(123);

julia> DelaunayTriangulation.random_edge(rng, E)
(10, 15)

julia> DelaunayTriangulation.random_edge(rng, E)
(10, 15)

julia> DelaunayTriangulation.random_edge(rng, E)
(23, 20)
```
```julia-repl 
julia> DelaunayTriangulation.random_edge(E)
(1, 2)
```
"""
random_edge
random_edge(E) = rand(E)
random_edge(rng, E) = rand(rng, E)

"""
  compare_unoriented_edge_collections(E, F) -> Bool

Tests if the edge collections `E` and `F` are equal, ignoring edge orientation.
"""
function compare_unoriented_edge_collections(E, F)
  num_edges(E) ≠ num_edges(F) && return false
  for e in each_edge(E)
    contains_unoriented_edge(e, F) || return false
  end
  return true
end

"""
    edges_are_disjoint(e, e′) -> Bool 

Returns `true` if `e` and `e′` have any shared vertex, and `false` otherwise.
"""
function edges_are_disjoint(e, e′)
    w = get_shared_vertex(e, e′)
    return w == ∅
end