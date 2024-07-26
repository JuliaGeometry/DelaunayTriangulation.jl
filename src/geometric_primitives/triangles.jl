@doc """
    construct_triangle(::Type{T}, i, j, k) where {T} -> Triangle

Construct a triangle of type `T` from vertices `i`, `j`, and `k`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.construct_triangle(NTuple{3,Int}, 1, 2, 3)
(1, 2, 3)

julia> DelaunayTriangulation.construct_triangle(Vector{Int32}, 1, 2, 3)
3-element Vector{Int32}:
 1
 2
 3
```
"""
construct_triangle
construct_triangle(::Type{NTuple{3,I}}, i, j, k) where {I} = (I(i), I(j), I(k))
construct_triangle(::Type{Vector{I}}, i, j, k) where {I} = I[i, j, k]
construct_triangle(::Type{A}, i, j, k) where {I,A<:AbstractVector{I}} = A(I[i, j, k])

"""
    geti(T) -> Vertex

Get the first vertex of `T`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.geti((1, 2, 3))
1

julia> DelaunayTriangulation.geti([2, 5, 1])
2
```
"""
geti(T) = T[1]

"""
    getj(T) -> Vertex

Get the second vertex of `T`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.getj((5, 6, 13))
6

julia> DelaunayTriangulation.getj([10, 19, 21])
19
```
"""
getj(T) = T[2]

"""
    getk(T) -> Vertex

Get the third vertex of `T`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.getk((1,2,3))
3

julia> DelaunayTriangulation.getk([1,2,3])
3
```
"""
getk(T) = T[3]

"""
    triangle_vertices(T) -> NTuple{3, Vertex}

Returns the vertices of `T` as a `Tuple`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> triangle_vertices((1, 5, 17))
(1, 5, 17)

julia> triangle_vertices([5, 18, 23]) # -> tuple
(5, 18, 23)
```
"""
triangle_vertices(T) = (geti(T), getj(T), getk(T))

"""
    triangle_type(::Type{T}) -> DataType

Get the triangle type of `T`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.triangle_type(Set{NTuple{3,Int64}})
Tuple{Int64, Int64, Int64}

julia> DelaunayTriangulation.triangle_type(Vector{NTuple{3,Int32}})
Tuple{Int32, Int32, Int32}

julia> DelaunayTriangulation.triangle_type(Vector{Vector{Int64}})
Vector{Int64} (alias for Array{Int64, 1})
```
"""
triangle_type(::Type{T}) where {T} = eltype(T)

"""
    num_triangles(T) -> Integer

Get the number of triangles in `T`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> T1, T2, T3 = (1, 5, 10), (17, 23, 10), (-1, 10, 5);

julia> T = Set((T1, T2, T3));

julia> num_triangles(T)
3
```
"""
num_triangles(T) = length(T)

@doc """
    triangle_edges(T) -> NTuple{3, Edge}
    triangle_edges(i, j, k) -> NTuple{3, Edge}

Get the edges of `T = (i, j, k)` as a `Tuple`, in particular 

    ((i, j), (j, k), (k, i)).

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> T = (1, 2, 3);

julia> DelaunayTriangulation.triangle_edges(T)
((1, 2), (2, 3), (3, 1))

julia> DelaunayTriangulation.triangle_edges(1, 2, 3)
((1, 2), (2, 3), (3, 1))
```
"""
triangle_edges
triangle_edges(i, j, k) = ((i, j), (j, k), (k, i))
triangle_edges(T) = triangle_edges(geti(T), getj(T), getk(T))

"""
    rotate_triangle(T, rotation) -> Triangle

Rotate the vertices of `T` by `rotation`. In particular, if 
`T = (i, j, k)`:

- `rotation = 0`: `(i, j, k)`
- `rotation = 1`: `(j, k, i)`
- `rotation = 2`: `(k, i, j)`
- Otherwise, return `rotate_triangle(T, rotation % 3)`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> T = (1, 2, 3)
(1, 2, 3)

julia> DelaunayTriangulation.rotate_triangle(T, 0)
(1, 2, 3)

julia> DelaunayTriangulation.rotate_triangle(T, 1)
(2, 3, 1)

julia> DelaunayTriangulation.rotate_triangle(T, 2)
(3, 1, 2)

julia> DelaunayTriangulation.rotate_triangle(T, 3)
(1, 2, 3)
```
"""
function rotate_triangle(T::V, ::Val{N}) where {N,V}
    i, j, k = triangle_vertices(T)
    N < 0 && throw(ArgumentError("Cannot rotate triangle $T by a negative amount."))
    if N == 0
        return T
    elseif N == 1
        return construct_triangle(V, j, k, i)
    elseif N == 2
        return construct_triangle(V, k, i, j)
    else
        return rotate_triangle(T, Val(N % 3))
    end
end
rotate_triangle(T::V, N::Integer) where {V} = rotate_triangle(T, Val(N))

"""
    construct_positively_oriented_triangle(::Type{V}, i, j, k, points) where {V} -> V

Construct a triangle of type `V` from vertices `i`, `j`, and `k` such that the
triangle is positively oriented, using `points` for the coordinates.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> points = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)];

julia> DelaunayTriangulation.construct_positively_oriented_triangle(NTuple{3, Int}, 1, 2, 3, points)
(2, 1, 3)

julia> DelaunayTriangulation.construct_positively_oriented_triangle(NTuple{3, Int}, 2, 3, 1, points)
(3, 2, 1)

julia> DelaunayTriangulation.construct_positively_oriented_triangle(NTuple{3, Int}, 2, 1, 3, points)
(2, 1, 3)

julia> DelaunayTriangulation.construct_positively_oriented_triangle(NTuple{3, Int}, 3, 2, 1, points)
(3, 2, 1)

julia> points = [(1.0, 1.0), (2.5, 2.3), (17.5, 23.0), (50.3, 0.0), (-1.0, 2.0), (0.0, 0.0), (5.0, 13.33)];

julia> DelaunayTriangulation.construct_positively_oriented_triangle(Vector{Int}, 5, 3, 2, points)
3-element Vector{Int64}:
 3
 5
 2

julia> DelaunayTriangulation.construct_positively_oriented_triangle(Vector{Int}, 7, 1, 2, points)
3-element Vector{Int64}:
 7
 1
 2

julia> DelaunayTriangulation.construct_positively_oriented_triangle(Vector{Int}, 7, 2, 1, points)
3-element Vector{Int64}:
 2
 7
 1

julia> DelaunayTriangulation.construct_positively_oriented_triangle(Vector{Int}, 5, 4, 3, points)
3-element Vector{Int64}:
 5
 4
 3
```
"""
function construct_positively_oriented_triangle(::Type{V}, i, j, k, points, predicates::AbstractPredicateType=def_alg222()) where {V}
    p, q, r = get_point(points, i, j, k)
    orientation = triangle_orientation(predicates, p, q, r)
    if is_negatively_oriented(orientation)
        return construct_triangle(V, j, i, k)
    else
        return construct_triangle(V, i, j, k)
    end
end

"""
    compare_triangles(T, V) -> Bool

Compare the triangles `T` and `V` by comparing their vertices up to rotation.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> T1 = (1, 5, 10);

julia> T2 = (17, 23, 20);

julia> DelaunayTriangulation.compare_triangles(T1, T2)
false

julia> T2 = (5, 10, 1);

julia> DelaunayTriangulation.compare_triangles(T1, T2)
true

julia> T2 = (10, 1, 5);

julia> DelaunayTriangulation.compare_triangles(T1, T2)
true

julia> T2 = (10, 5, 1);

julia> DelaunayTriangulation.compare_triangles(T1, T2)
false
```
"""
function compare_triangles(T, V)
    i, j, k = triangle_vertices(T)
    u, v, w = triangle_vertices(V)
    return (i, j, k) == (u, v, w) ||
           (i, j, k) == (v, w, u) ||
           (i, j, k) == (w, u, v)
end

@doc """
    contains_triangle(T, V) -> (Triangle, Bool)

Check if the collection of triangles `V` contains the triangle `T` up to rotation.
The `Triangle` returned is the triangle in `V` that is equal to `T` up to rotation, 
or `T` if no such triangle exists. The `Bool` is `true` if `V` contains `T`,
and `false` otherwise.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> V = Set(((1, 2, 3), (4, 5, 6), (7, 8, 9)))
Set{Tuple{Int64, Int64, Int64}} with 3 elements:
  (7, 8, 9)
  (4, 5, 6)
  (1, 2, 3)

julia> DelaunayTriangulation.contains_triangle((1, 2, 3), V)
((1, 2, 3), true)

julia> DelaunayTriangulation.contains_triangle((2, 3, 1), V)
((1, 2, 3), true)

julia> DelaunayTriangulation.contains_triangle((10, 18, 9), V)
((10, 18, 9), false)

julia> DelaunayTriangulation.contains_triangle(9, 7, 8, V)
((7, 8, 9), true)
```
"""
contains_triangle
function contains_triangle(T::F, V::S) where {F,S}
    if F ≠ triangle_type(S)
        i, j, k = triangle_vertices(T)
        Tfix = construct_triangle(triangle_type(S), i, j, k)
    else
        Tfix = T
    end
    if Tfix ∈ V
        return Tfix, true
    end
    T1 = rotate_triangle(Tfix, Val(1))
    if T1 ∈ V
        return T1, true
    end
    T2 = rotate_triangle(Tfix, Val(2))
    if T2 ∈ V
        return T2, true
    end
    return Tfix, false
end
function contains_triangle(i, j, k, V::Vs) where {Vs}
    T = construct_triangle(triangle_type(Vs), i, j, k)
    return contains_triangle(T, V)
end

@doc """
    sort_triangle(T) -> Triangle
    sort_triangle(i, j, k) -> Triangle

Sort the triangle `T = (i, j, k)` so that its last vertex is the smallest, 
respecting the orientation of `T`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.sort_triangle((1, 5, 3))
(5, 3, 1)

julia> DelaunayTriangulation.sort_triangle((1, -1, 2))
(2, 1, -1)

julia> DelaunayTriangulation.sort_triangle((3, 2, 1))
(3, 2, 1)
```
"""
sort_triangle
function sort_triangle(T::V) where {V}
    i, j, k = triangle_vertices(T)
    minijk = min(i, j, k)
    if minijk == i
        return construct_triangle(V, j, k, i)
    elseif minijk == j
        return construct_triangle(V, k, i, j)
    else
        return construct_triangle(V, i, j, k)
    end
end
sort_triangle(i::Integer, j::Integer, k::Integer) = sort_triangle((i, j, k))
@inline function sort_triangle(i, j, k) # in case one of the vertices is a point rather than an integer 
    # Note that the return type is a 3-Union 
    if is_ghost_vertex(i)
        return j, k, i
    elseif is_ghost_vertex(j)
        return k, i, j
    else
        return i, j, k
    end
end

"""
    add_to_triangles!(T, V)

Add the triangle `V` to the collection of triangles `T`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> T = Set(((1, 2, 3), (17, 8, 9)));

julia> DelaunayTriangulation.add_to_triangles!(T, (1, 5, 12))
Set{Tuple{Int64, Int64, Int64}} with 3 elements:
  (1, 5, 12)
  (1, 2, 3)
  (17, 8, 9)

julia> DelaunayTriangulation.add_to_triangles!(T, (-1, 3, 6))
Set{Tuple{Int64, Int64, Int64}} with 4 elements:
  (1, 5, 12)
  (1, 2, 3)
  (17, 8, 9)
  (-1, 3, 6)
```
"""
function add_to_triangles!(T::Ts, V::VT) where {Ts,VT}
    if VT ≠ triangle_type(Ts)
        i, j, k = triangle_vertices(V)
        U = construct_triangle(triangle_type(Ts), i, j, k)
        return push!(T, U)
    else
        return push!(T, V)
    end
end

@doc """
    add_triangle!(T, V...)
    add_triangle!(T, i, j, k)

Add the triangles `V...` or `V = (i, j, k)` to the collection of triangles `T`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> T = Set(((1, 2, 3), (4, 5, 6)))
Set{Tuple{Int64, Int64, Int64}} with 2 elements:
  (4, 5, 6)
  (1, 2, 3)

julia> add_triangle!(T, (7, 8, 9));

julia> add_triangle!(T, (10, 11, 12), (13, 14, 15));

julia> add_triangle!(T, 16, 17, 18);

julia> T
Set{Tuple{Int64, Int64, Int64}} with 6 elements:
  (7, 8, 9)
  (10, 11, 12)
  (4, 5, 6)
  (13, 14, 15)
  (16, 17, 18)
  (1, 2, 3)
```
"""
add_triangle!
function add_triangle!(T, V...)
    for W in V
        add_to_triangles!(T, W)
    end
    return T
end
function add_triangle!(T::Ts, i::Integer, j::Integer, k::Integer) where {Ts} # method ambiguity 
    V = triangle_type(Ts)
    I = number_type(V)
    F = construct_triangle(V, I(i), I(j), I(k))
    return add_to_triangles!(T, F)
end

"""
    delete_from_triangles!(T, V)

Delete the triangle `V` from the collection of triangles `T`.
Only deletes `V` if `V` is in `T` up to rotation.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> V = Set(((1, 2, 3), (4, 5, 6), (7, 8, 9)))
Set{Tuple{Int64, Int64, Int64}} with 3 elements:
  (7, 8, 9)
  (4, 5, 6)
  (1, 2, 3)

julia> DelaunayTriangulation.delete_from_triangles!(V, (4, 5, 6))
Set{Tuple{Int64, Int64, Int64}} with 2 elements:
  (7, 8, 9)
  (1, 2, 3)

julia> DelaunayTriangulation.delete_from_triangles!(V, (9, 7, 8))
Set{Tuple{Int64, Int64, Int64}} with 1 element:
  (1, 2, 3)
```
"""
function delete_from_triangles!(V, T)
    rotated_T, pred = contains_triangle(T, V)
    pred && delete!(V, rotated_T)
    return V
end

@doc """
    delete_triangle!(V, T...)
    delete_triangle!(V, i, j, k)

Delete the triangles `T...` from the collection of triangles `V`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> V = Set(((1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12), (13, 14, 15)))
Set{Tuple{Int64, Int64, Int64}} with 5 elements:
  (7, 8, 9)
  (10, 11, 12)
  (4, 5, 6)
  (13, 14, 15)
  (1, 2, 3)

julia> delete_triangle!(V, (6, 4, 5))
Set{Tuple{Int64, Int64, Int64}} with 4 elements:
  (7, 8, 9)
  (10, 11, 12)
  (13, 14, 15)
  (1, 2, 3)

julia> delete_triangle!(V, (10, 11, 12), (1, 2, 3))
Set{Tuple{Int64, Int64, Int64}} with 2 elements:
  (7, 8, 9)
  (13, 14, 15)

julia> delete_triangle!(V, 8, 9, 7)
Set{Tuple{Int64, Int64, Int64}} with 1 element:
  (13, 14, 15)
```
"""
delete_triangle!
function delete_triangle!(V, T...)
    for W in T
        delete_from_triangles!(V, W)
    end
    return V
end
function delete_triangle!(V::Vs, i::Integer, j::Integer, k::Integer) where {Vs} # method ambiguity 
    T = triangle_type(Vs)
    I = number_type(T)
    F = construct_triangle(T, I(i), I(j), I(k))
    return delete_from_triangles!(V, F)
end

@doc """
    each_triangle(T) -> Iterator

Return an iterator over the triangles in `T`.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> T = Set(((1, 2, 3), (-1, 5, 10), (17, 13, 18)));

julia> each_triangle(T)
Set{Tuple{Int64, Int64, Int64}} with 3 elements:
  (-1, 5, 10)
  (1, 2, 3)
  (17, 13, 18)

julia> T = [[1, 2, 3], [10, 15, 18], [1, 5, 6]];

julia> each_triangle(T)
3-element Vector{Vector{Int64}}:
 [1, 2, 3]
 [10, 15, 18]
 [1, 5, 6]
```
"""
each_triangle
each_triangle(T) = T
each_triangle(T::AbstractMatrix) = eachcol(T)

"""
    compare_triangle_collections(T, V) -> Bool

Compare the collections of triangles `T` and `V` by comparing their triangles
according to [`compare_triangles`](@ref).

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> T = Set(((1, 2, 3), (4, 5, 6), (7, 8, 9)));

julia> V = [[2, 3, 1], [4, 5, 6], [9, 7, 8]];

julia> DelaunayTriangulation.compare_triangle_collections(T, V)
true

julia> V[1] = [17, 19, 20];

julia> DelaunayTriangulation.compare_triangle_collections(T, V)
false

julia> V = [[1, 2, 3], [8, 9, 7]];

julia> DelaunayTriangulation.compare_triangle_collections(T, V)
false
```
"""
function compare_triangle_collections(T, V)
    num_triangles(T) == num_triangles(V) || return false
    for τ in each_triangle(T)
        _, flag = contains_triangle(τ, V)
        !flag && return false
    end
    return true
end

"""
    sort_triangles(T) -> Triangle

Sort the triangles in `T` so that the first vertex of each triangle is the largest,
respecting the orientation of the triangles. See [`sort_triangle`](@ref).

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> T = Set(((1, 3, 2), (5, 2, 3), (10, 1, 13), (-1, 10, 12), (10, 1, 17), (5, 8, 2)))
Set{Tuple{Int64, Int64, Int64}} with 6 elements:
  (5, 8, 2)
  (10, 1, 13)
  (10, 1, 17)
  (5, 2, 3)
  (1, 3, 2)
  (-1, 10, 12)

julia> DelaunayTriangulation.sort_triangles(T)
Set{Tuple{Int64, Int64, Int64}} with 6 elements:
  (13, 10, 1)
  (3, 5, 2)
  (10, 12, -1)
  (5, 8, 2)
  (17, 10, 1)
  (3, 2, 1)
```
"""
function sort_triangles(T::Ts) where {Ts}
    V = Ts()
    for τ in each_triangle(T)
        σ = sort_triangle(τ)
        add_to_triangles!(V, σ)
    end
    return V
end