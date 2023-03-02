"""
    construct_triangle(::Type{T}, i, j, k) where {T}

Constructs a triangle with indices `(i, j, k)` with the 
type `T`. The following methods are currently defined:

    construct_triangle(::Type{NTuple{3, I}}, i, j, k) where {I} 
    construct_triangle(::Type{A}, i, j, k) where {I, A <: AbstractVector{I}}

You can extend this function as you need, ensuring you define it 
on the type rather than on an instance of the type.
"""
function construct_triangle end
function construct_triangle(::Type{F}, i, j, k) where {F}
    return error("The construct_triangle function has not been defined for the type $F.")
end
function construct_triangle(::Type{NTuple{3,I}}, i, j, k) where {I}
    return (I(i), I(j), I(k))
end
function construct_triangle(::Type{A}, i, j, k) where {I,A<:AbstractVector{I}}
    return A([i, j, k])
end

"""
    geti(T::F) where {F}

From a triangle `T`, extract the `i`th index, i.e. the first. 
The following methods are currently defined: 

    geti(T::NTuple{3, I})
    geti(T::A) where {I,A<:AbstractVector{I}} 

You can extend this function as you need.
"""
function geti end
function geti(::F) where {F}
    return error("The geti function has not been defined for the type $F.")
end
geti(T::NTuple{3,I}) where {I} = T[1]
geti(T::A) where {I,A<:AbstractVector{I}} = T[1]

"""
    getj(T::F) where {F}

From a triangle `T`, extract the `j`th index, i.e. the second. 
The following methods are currently defined: 

    getj(T::NTuple{3, I})
    getj(T::A) where {I,A<:AbstractVector{I}} 

You can extend this function as you need.
"""
function getj end
function getj(::F) where {F}
    return error("The getj function has not been defined for the type $F.")
end
getj(T::NTuple{3,I}) where {I} = T[2]
getj(T::A) where {I,A<:AbstractVector{I}} = T[2]

"""
    getk(T::F) where {F}

From a triangle `T`, extract the `k`th index, i.e. the third. 
The following methods are currently defined: 

    getk(T::NTuple{3, I})
    getk(T::A) where {I,A<:AbstractVector{I}} 

You can extend this function as you need.
"""
function getk end
function getk(::F) where {F}
    return error("The getk function has not been defined for the type $F.")
end
getk(T::NTuple{3,I}) where {I} = T[3]
getk(T::A) where {I,A<:AbstractVector{I}} = T[3]

"""
    indices(T)

Returns a tuple `(i, j, k)` containing the indices of the triangle `T`. The 
indices are obtained using [`geti`](@ref), [`getj`](@ref), and [`getk`](@ref).
"""
function indices end
indices(T) = (geti(T), getj(T), getk(T))

"""
    integer_type(::Type{T}) where {T} 

Returns the integer type used for representnig a triangle with indices 
`(i, j, k)` with the type `T`. The following methods are currently defined:

    integer_type(::Type{NTuple{N, I}}) where {N, I} 
    integer_type(::Type{A}) where {I, A <: AbstractVector{I}}

You can extend this function as you need, ensuring you define it 
on the type rather than on an instance of the type.
"""
function integer_type end
function integer_type(::Type{F}) where {F}
    return error("The integer_type function has not been defined for the type $F.")
end
integer_type(::Type{NTuple{N,I}}) where {N,I} = I
integer_type(::Type{A}) where {I,A<:AbstractVector{I}} = I

"""
    triangle_edges(T)
    triangle_edges(i, j, k)

Returns an iterator over each edge of the triangle `T`. In particular, 
returns `((geti(T), getj(T)), (getj(T), getk(T)), (getk(T), geti(T)))`.
The latter method uses the indices directly.
"""
function triangle_edges end
triangle_edges(i, j, k) = ((i, j), (j, k), (k, i))
triangle_edges(T) = triangle_edges(geti(T), getj(T), getk(T))

"""
    rotate_triangle(T::V, i) where {V}

Given a triangle `T`, rotates the indices an amount `i`. In particular,
if `T = (i, j, k)`:

- `i = 0`: Returns `(i, j, k)`.
- `i = 1`: Returns `(j, k, i)`.
- `i = 2`: Returns `(k, i, j)`.
- Otherwise, return `rotate_triangle(T, i % 3)`.
"""
function rotate_triangle end
function rotate_triangle(T::V, i) where {V}
    if i == 0
        return T
    elseif i == 1
        return construct_triangle(V, getj(T), getk(T), geti(T))
    elseif i == 2
        return construct_triangle(V, getk(T), geti(T), getj(T))
    else
        return rotate_triangle(T, i % 3)
    end
end

"""
    construct_positively_oriented_triangle(::Type{V}, i, j, k, points) where {V}

Given a triangle type `V`, indices `(i, j, k)` corresponding to points in `points`, 
returns either `construct_triangle(V, i, j, k)` or `construct_triangle(V, j, i, k)`,
whichever is not negatively oriented.
"""
function construct_positively_oriented_triangle(::Type{V}, i, j, k, points) where {V}
    p, q, r = get_point(points, i, j, k)
    orientation = triangle_orientation(p, q, r)
    if is_negatively_oriented(orientation)
        return construct_triangle(V, j, i, k)
    else
        return construct_triangle(V, i, j, k)
    end
end

"""
    compare_triangles(T, V) 

Tests if the triangle `T` is equal to the triangle `V`, with equality 
defined on the indices. In particular, `compare_triangles((i, j, k), (i, j, k))`,
`compare_triangles((i, j, k), (j, k, i))`, and `compare_triangles((i, j, k), (k, i, j))`
are all true.
"""
function compare_triangles end
function compare_triangles(T, V)
    i, j, k = indices(T)
    u, v, w = indices(V)
    return (i, j, k) == (u, v, w) ||
           (i, j, k) == (v, w, u) ||
           (i, j, k) == (w, u, v)
end

"""
    initialise_triangles(::Type{S})

For a given type `S` for some collection (e.g. a `Set`), returns an
empty instance of that collection. The only method defined is

    initialise_triangles(::Type{S}) where {T, S <: Set{T}}

which returns a `Set{T}()`. You can extend this function as you need, making sure 
you extend it for the type rather than for instances of that type.
"""
function initialise_triangles end
function initialise_triangles(::Type{F}) where {F}
    return error("The initialise_triangles function has not been defined for the type $F.")
end
initialise_triangles(::Type{S}) where {T,S<:Set{T}} = S()

"""
    triangle_type(::Type{S}) where {S}

For a given type `S` representing a collection of triangles, 
returns the type of triangle used inside `S`, e.g. `NTuple{3, Int64}`
if `S = Set{NTuple{3, Int64}}`. The only methods defined are

    triangle_type(::Type{S}) where {T, S <: Set{T}}
    triangle_type(::Type{A}) where {T, A <: AbstractVector{E}}

which returns `T`. You can extend this function as you need, making sure 
you extend it for the type rather than for instances of that type.
"""
function triangle_type end
function triangle_type(::Type{F}) where {F}
    return error("The triangle_type function has not been defined for the type $F.")
end
triangle_type(::Type{S}) where {T,S<:Set{T}} = T
triangle_type(::Type{A}) where {T,A<:AbstractVector{T}} = T

"""
    num_triangles(T::S) where {S}

Given a collection of triangles `T`, returns the number of triangles 
in `T`. The only method currently defined is 

    num_triangles(T::Set)

which returns `length(T)`. You can extend this function as you need.
"""
function num_triangles end
function num_triangles(::F) where {F}
    return error("The num_triangles function has not been defined for the type $F.")
end
num_triangles(T::Set) = length(T)

"""
    contains_triangle(T::F, V::S) where {F, S}

Given a collection of triangles `V` of type `S`, containing triangles 
of type `F`, checks if `V` includes the triangle `T`. In particular, 
the returned value is 

    (rotate_triangle(T, m), true)

if `V` includes the triangle `T`, and `m` is the integer such that 
`rotate_triangle(T, m)` is the form of `T` that `V` contains. If there 
is no such `m`, the returned value is simply 

    (T, false).

To use this function, your type needs to only have a definition for `T ∈ V` for 
testing specific rotations of a triangle. This function then performs the checks 
by checking `T`, then `rotate_triangle(T, 1)`,
then `rotate_triangle(T, 2)`, and if none of those succeed just returns 
`(T, false)`. You can extend this function as you need. We also define

    contains_triangle(i, j, k, V::Ts) where {Ts},

and this method just calls the two-argument method after constructing 
the triangle with [`construct_triangle`](@ref).
"""
function contains_triangle end
function contains_triangle(T::F, V::S) where {F,S}
    if T ∈ V
        return (T, true)
    end
    T1 = rotate_triangle(T, 1)
    if T1 ∈ V
        return (T1, true)
    end
    T2 = rotate_triangle(T, 2)
    if T2 ∈ V
        return (T2, true)
    end
    return (T, false)
end
function contains_triangle(i, j, k, V::Ts) where {Ts}
    F = triangle_type(Ts)
    T = construct_triangle(F, i, j, k)
    return contains_triangle(T, V)
end

"""
    sort_triangle(T::V) where {V}
    
Given a triangle `T = (i, j, k)`, sorts it so that the first index is the smallest, maintaining 
the orientation of `T`.
"""
function sort_triangle(T::V) where {V}
    i, j, k = indices(T)
    minijk = min(i, j, k)
    if minijk == i
        return construct_triangle(V, i, j, k)
    elseif minijk == j
        return construct_triangle(V, j, k, i)
    else
        return construct_triangle(V, k, i, j)
    end
end

"""
    add_to_triangles!(T::S, V) where {S}

Given a collection of triangles `T`, pushes `V` into it. The only 
methods currently defined are

    add_to_triangles!(T::Set{F}, V::F) where {F}
    add_to_triangles!(T::Set{F}, V) where {F}

which simply call `push!(T, V)`. The latter method reconstructs `V` 
using [`indices`] and [`construct_triangle`](@ref). You can extend this function  
as you need. 
"""
function add_to_triangles! end
function add_to_triangles!(::S, V) where {S}
    return error("The add_to_triangles! function has not been defined for the type $S.")
end
function add_to_triangles!(T::Set{F}, V::F) where {F}
    push!(T, V)
    return nothing
end
function add_to_triangles!(T::Set{F}, V::G) where {F,G}
    i, j, k = indices(V)
    V = construct_triangle(F, i, j, k)
    add_to_triangles!(T, V)
    return nothing
end

"""
    add_triangle!(T, V...)

Given a collection of triangles `T`, adds all the triangles `V...` into it. 
To extend this method to other collections, see [`add_to_triangles!`](@ref).
"""
function add_triangle!(T, V::Vararg{F,N}) where {F,N}
    for i in 1:N
        add_to_triangles!(T, V[i])
    end
    return nothing
end

"""
    add_triangle!(T::Ts, i::Integer, j::Integer, k::Integer) where {Ts}

Given a collection of triangles `T`, adds the triangle `(i, j, k)` into it.
"""
function add_triangle!(T::Ts, i::Integer, j::Integer, k::Integer) where {Ts}
    V = triangle_type(Ts)
    I = integer_type(V)
    F = construct_triangle(V, I(i), I(j), I(k))
    add_triangle!(T, F)
    return nothing
end

"""
    delete_from_triangles!(V::S, T::F) where {S}

Given a collection of triangles `V` of type `S`, containing 
triangles of type `F`, deletes the triangle `T` from `V`. 
The only method currently defined is 

    delete_from_triangles!(V::Set{F}, T::F) where {F}.
    
which just calls `delete!` on `V`. The function already assumes that `T` 
is already in `V`, and that `T` doesn't need to be rotated at all.
"""
function delete_from_triangles! end
function delete_from_triangles!(::S, T) where {S}
    return error("The delete_from_triangles! function has not been defined for the type $S.")
end
function delete_from_triangles!(V::Set{F}, T::F) where {F}
    delete!(V, T)
    return nothing
end

"""
    delete_triangle!(V, T...)

Given a collection of triangles `V`, deletes all the triangles `T...` from it. 
Checks are made for rotated forms of `T` that `V` includes. For example, if 
`T = (i, j, k)` but `V` contains `(j, k, i)`, then `(j, k, i)` will be deleted.
The function also checks to see if the triangle is in `V` at all prior to deleting it.
To extend this method to other collections, see [`delete_from_triangles!`](@ref).
"""
function delete_triangle!(V, T::Vararg{F,N}) where {F,N}
    for i in 1:N
        rotated_T, pred = contains_triangle(T[i], V)
        pred && delete_from_triangles!(V, rotated_T)
    end
    return nothing
end

"""
    delete_triangle!(T::Ts, i::Integer, j::Integer, k::Integer) where {Ts}

Given a collection of triangles `T`, deletes the triangle `(i, j, k)` from it. 
Checks are made to see if `T` needs to be rotated, or if `(i, j, k)` is in `T` 
at all.
"""
function delete_triangle!(T::Ts, i::Integer, j::Integer, k::Integer) where {Ts}
    V = triangle_type(Ts)
    I = integer_type(V)
    F = construct_triangle(V, I(i), I(j), I(k))
    rotated_F, pred = contains_triangle(F, T)
    pred && delete_triangle!(T, rotated_F)
    return nothing
end

"""
    each_triangle(V::F) where {F}

For a given collection of triangles `V`, returns an iterator that 
goes over each triangle in the collection. The methods currently 
defined are 

    each_triangle(V::Set)
    each_triangle(V::AbstractMatrix)

with the first method simply returning `V`, and the second returning 
`eachcol(V)`. You can extend this function as you need.
"""
function each_triangle end
function each_triangle(::F) where {F}
    return error("The each_triangle function has not been defined for the type $F.")
end
each_triangle(V::Set) = V
each_triangle(V::AbstractMatrix) = eachcol(V)

"""
    compare_triangle_collections(T, V)

Given two collections of triangles `T` and `V`, tests if they are equal, with 
equality defined according to [`compare_triangles`](@ref).
"""
function compare_triangle_collections(T, V)
    num_triangles(T) ≠ num_triangles(V) && return false
    for τ in each_triangle(T)
        _, flag = contains_triangle(τ, V)
        !flag && return false
    end
    return true
end

"""
    sort_triangles(T::Ts) where {Ts}

Sorts the triangles in the collection `T` so that each triangle's first vertex 
has the smallest value. The orientation of each triangle is preserved.
"""
function sort_triangles(T::Ts) where {Ts}
    tris = initialise_triangles(Ts)
    for τ in each_triangle(T)
        σ = sort_triangle(τ)
        add_triangle!(tris, σ)
    end
    return tris
end

"""
    remove_duplicate_triangles(T::Ts) where {Ts} 

Removes duplicate triangles from `T`. This procedure also sorts the triangles 
so that the first index of each triangle is the smallest. Orientations are 
preserved.
"""
function remove_duplicate_triangles(T::Ts) where {Ts}
    V = sort_triangles(T)
    Ts <: Set || unique!(V)
    return V
end