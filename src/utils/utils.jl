@doc """
    number_type(x) -> DataType

Given a container `x`, returns the number type used for storing coordinates.

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.number_type([1, 2, 3])
Int64

julia> DelaunayTriangulation.number_type((1, 2, 3))
Int64

julia> DelaunayTriangulation.number_type([1.0 2.0 3.0; 4.0 5.0 6.0])
Float64

julia> DelaunayTriangulation.number_type([[[1, 2, 3, 4, 5, 1]], [[6, 8, 9], [9, 10, 11], [11, 12, 6]]])
Int64

julia> DelaunayTriangulation.number_type((1.0f0, 2.0f0))
Float32

julia> DelaunayTriangulation.number_type(Vector{Float64})
Float64

julia> DelaunayTriangulation.number_type(Vector{Vector{Float64}})
Float64

julia> DelaunayTriangulation.number_type(NTuple{2, Float64})
Float64
```
"""
number_type
number_type(x) = number_type(typeof(x))
number_type(::Type{T}) where {T <: AbstractArray} = number_type(eltype(T))
number_type(::Type{<:NTuple{N, T}}) where {N, T} = number_type(T)
number_type(::Type{<:NTuple{0}}) = Any
number_type(::Type{Tuple{}}) = Any
number_type(::Type{T}) where {T} = T

@doc """
    get_ghost_vertex(i, j, k) -> Vertex
    get_ghost_vertex(i, j) -> Vertex

Given three vertices `i`, `j`, and `k`, returns the ghost vertex among them. If none of them are ghost vertices, returns `k`.
The two-argument version is equivalent to `get_ghost_vertex(i, j, j)`.

# Examples 
```jldoctest 
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.get_ghost_vertex(1, 7, -2)
-2

julia> DelaunayTriangulation.get_ghost_vertex(-1, 2, 3)
-1

julia> DelaunayTriangulation.get_ghost_vertex(1, 5, 10)
10

julia> DelaunayTriangulation.get_ghost_vertex(1, -1)
-1

julia> DelaunayTriangulation.get_ghost_vertex(-5, 2)
-5
```
"""
get_ghost_vertex
function get_ghost_vertex(i, j, k)
    is_ghost_vertex(i) && return i
    is_ghost_vertex(j) && return j
    return k
end
function get_ghost_vertex(i, j)
    is_ghost_vertex(i) && return i
    return j
end

@doc """
    is_true(b) -> Bool

Returns `b` represents a `true` value, and `false` otherwise.

# Examples
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.is_true(true)
true

julia> DelaunayTriangulation.is_true(false)
false

julia> DelaunayTriangulation.is_true(Val(true))
true

julia> DelaunayTriangulation.is_true(Val(false))
false
```
"""
is_true
is_true(b::Bool) = b
is_true(b::Val{true}) = true
is_true(b::Val{false}) = false

"""
    choose_uvw(e1, e2, e3, u, v, w) -> (Vertex, Vertex, Vertex)

Choose values for `(u, v, w)` based on the Booleans `(e1, e2, e3)`, 
assuming only one is true. The three cases are: 

- If `e1`, returns `(u, v, w)`.
- If `e2`, returns `(v, w, u)`.
- If `e3`, returns `(w, u, v)`.
"""
function choose_uvw(e1, e2, e3, u, v, w)
    e1 && return (u, v, w)
    e2 && return (v, w, u)
    return (w, u, v)
end

"""
    is_circular(A) -> Bool

Tests if `A` is circular, meaning that `A[begin] == A[end]`.
"""
is_circular(A) = isempty(A) || (A[begin] == A[end])

"""
    circular_equality(A, B, by=isequal) -> Bool

Compares the two circular vectors `A` and `B` for equality up to circular rotation, 
using `by` to compare individual elements.
"""
function circular_equality(A, B, by::F = isequal) where {F}
    @assert is_circular(A) && is_circular(B) "A and B must be circular"
    length(A) ≠ length(B) && return false
    isempty(A) && return true # isempty(B) is true as well because of the previous assertion 
    _A = @views A[begin:(end - 1)]
    _B = @views B[begin:(end - 1)]
    same_idx = findfirst(Base.Fix1(by, _A[begin]), _B)
    same_idx === nothing && return false
    n = length(_A)
    for (i, a) in pairs(_A)
        j = mod1(i + same_idx - 1, n)
        b = _B[j]
        !by(a, b) && return false
    end
    return true
end

"""
    get_ordinal_suffix(i) -> String

Returns the ordinal suffix for the integer `i`. 

# Examples 
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.get_ordinal_suffix(1)
"st"

julia> DelaunayTriangulation.get_ordinal_suffix(2)
"nd"

julia> DelaunayTriangulation.get_ordinal_suffix(3)
"rd"

julia> DelaunayTriangulation.get_ordinal_suffix(4)
"th"

julia> DelaunayTriangulation.get_ordinal_suffix(5)
"th"

julia> DelaunayTriangulation.get_ordinal_suffix(6)
"th"

julia> DelaunayTriangulation.get_ordinal_suffix(11)
"th"

julia> DelaunayTriangulation.get_ordinal_suffix(15)
"th"

julia> DelaunayTriangulation.get_ordinal_suffix(100)
"th"
```
"""
function get_ordinal_suffix(i) # https://stackoverflow.com/a/13627586
    let j = i % 10, k = i % 100
        if j == 1 && k ≠ 11
            return "st"
        elseif j == 2 && k ≠ 12
            return "nd"
        elseif j == 3 && k ≠ 13
            return "rd"
        else
            return "th"
        end
    end
end

"""
    nextindex_circular(C, i) -> Integer 

Returns the next index after `i` in the circular vector `C`. 
"""
nextindex_circular(C, i) = i == lastindex(C) ? firstindex(C) : i + 1

"""
    previndex_circular(C, i) -> Integer

Returns the previous index before `i` in the circular vector `C`.
"""
previndex_circular(C, i) = i == firstindex(C) ? lastindex(C) - 1 : i - 1

"""
    replace_boundary_triangle_with_ghost_triangle(tri::Triangulation, V) -> Triangle

Given a boundary triangle `V` of `tri`, returns the adjacent ghost triangle. Note that 
for triangles in a corner of a domain, like a lattice triangulation, there are two choices 
of ghost triangle.
"""
function replace_boundary_triangle_with_ghost_triangle(tri::Triangulation, V)
    u, v, w = triangle_vertices(V)
    T = triangle_type(tri)
    is_boundary_edge(tri, u, v) && return construct_triangle(T, v, u, get_adjacent(tri, v, u))
    is_boundary_edge(tri, v, w) && return construct_triangle(T, w, v, get_adjacent(tri, w, v))
    return construct_triangle(T, u, w, get_adjacent(tri, u, w))
end

"""
    replace_ghost_triangle_with_boundary_triangle(tri::Triangulation, V) -> Triangle

Given a ghost triangle `V` of `tri`, returns the adjacent boundary triangle.
"""
function replace_ghost_triangle_with_boundary_triangle(tri::Triangulation, V)
    T = sort_triangle(V)
    u, v, w = triangle_vertices(T) # w is the ghost triangle 
    return construct_triangle(triangle_type(tri), v, u, get_adjacent(tri, v, u))
end

@doc raw"""
    iterated_neighbourhood(tri::Triangulation, i, d) -> Set{Vertex}

Returns the set of vertices of `tri` in the iterated neighbourhood of the vertex `i` of depth `d`. 

# Extended help 
The ``d``-times iterated neighbourhood is defined by 

```math 
N_i^d = \bigcup_{j \in N_i^{d-1}} N_j \setminus \{i\},
```
where ``N_i^1 = N_i`` is the set of neighbours of ``i``.
"""
function iterated_neighbourhood(tri::Triangulation, i, d)
    I = integer_type(tri)
    neighbours = Set{I}()
    sizehint!(neighbours, ceil(I, 6^(d / 2)))
    return iterated_neighbourhood!(neighbours, tri, i, d)
end
function iterated_neighbourhood!(neighbours, tri, i, d)
    empty!(neighbours)
    i_neighbours = get_neighbours(tri, i)
    I = integer_type(tri)
    for j in i_neighbours
        if !is_ghost_vertex(j)
            push!(neighbours, j)
        end
    end
    for _ in 2:d
        new_neighbours = Set{I}() # don't want to mutate the iterator while iterating
        for j in neighbours
            for k in get_neighbours(tri, j)
                if k ≠ i && !is_ghost_vertex(k)
                    push!(new_neighbours, k)
                end
            end
        end
        union!(neighbours, new_neighbours)
    end
    return neighbours
end

"""
    get_area(tri::Triangulation) -> Number 

Returns the area of `tri`.
"""
function get_area(tri::Triangulation)
    F = number_type(tri)
    A = zero(F)
    for T in each_solid_triangle(tri)
        A += triangle_area(tri, T)
    end
    return A
end

"""
    norm(p) -> Number 

Assuming `p` is two-dimensional, computes the Euclidean norm of `p`.
"""
function norm(p)
    x, y = getxy(p)
    return hypot(x, y)
end

"""
    norm_sqr(p) -> Number

Assuming `p` is two-dimensional, computes the square of the Euclidean norm of `p`.
"""
function norm_sqr(p)
    x, y = getxy(p)
    return x^2 + y^2
end

"""
    dist(p, q) -> Number

Assuming `p` and `q` are two-dimensional, computes the Euclidean distance between `p` and `q`.
"""
function dist(p, q)
    px, py = getxy(p)
    qx, qy = getxy(q)
    vec = (qx - px, qy - py)
    return norm(vec)
end

"""
    dist_sqr(p, q) -> Number

Assuming `p` and `q` are two-dimensional, computes the square of the Euclidean distance between `p` and `q`.
"""
function dist_sqr(p, q)
    px, py = getxy(p)
    qx, qy = getxy(q)
    vec = (qx - px, qy - py)
    return norm_sqr(vec)
end

@doc """
    edge_length(tri::Triangulation, u, v) -> Number
    edge_length(tri::Triangulation, e) -> Number

Computes the length of the edge `e = (u, v)`.
"""
edge_length
function edge_length(tri::Triangulation, u, v)
    p, q = get_point(tri, u, v)
    return dist(p, q)
end
edge_length(tri::Triangulation, e) = edge_length(tri, initial(e), terminal(e))

@doc """
    edge_length_sqr(tri::Triangulation, u, v) -> Number
    edge_length_sqr(tri::Triangulation, e) -> Number

Computes the square of the length of the edge `e = (u, v)`.
"""
function edge_length_sqr(tri::Triangulation, u, v)
    p, q = get_point(tri, u, v)
    return dist_sqr(p, q)
end
edge_length_sqr(tri::Triangulation, e) = edge_length_sqr(tri, initial(e), terminal(e))

"""
    midpoint(p, q) -> Number or NTuple{2, Number}

Assuming `p` and `q` are either both numbers are both 2-vectors, computes their average.
"""
function midpoint(x::Number, y::Number)
    if max(abs(x), abs(y)) ≥ 1
        return x / 2 + y / 2
    else
        return (x + y) / 2
    end
end
function midpoint(p, q)
    px, py = getxy(p)
    qx, qy = getxy(q)
    return (midpoint(px, qx), midpoint(py, qy))
end

"""
    midpoint(tri::Triangulation, u, v) -> NTuple{2, Number}
    midpoint(tri::Triangulation, e) -> NTuple{2, Number}

Computes the midpoint of `e = (u, v)`.
"""
function midpoint(tri::Triangulation, u, v)
    p, q = get_point(tri, u, v)
    return midpoint(p, q)
end
midpoint(tri::Triangulation, e) = midpoint(tri, initial(e), terminal(e))

"""
    check_precision(x) -> Bool 

Returns `true` if `abs(x)` is less than or equal to `sqrt(eps(number_type(eps)))`.
"""
check_precision(x) = abs(x) ≤ ε(x)

"""
    check_absolute_precision(x, y) -> Bool

Returns `true` if `abs(x - y)` is less than or equal to `sqrt(eps(Float64))`.
"""
check_absolute_precision(x, y) = check_precision(x - y)

"""
    check_relative_precision(x, y) -> Bool

Returns `true` if `abs(x - y)/max(abs(x), abs(y))` is less than or equal to `sqrt(eps(Float64))`.
"""
function check_relative_precision(x, y)
    x, y = abs(x), abs(y)
    if x < y
        x, y = y, x
    end
    return !iszero(x) && check_precision(abs(x - y) / x)
end

"""
    check_ratio_precision(x, y) -> Bool 

Returns `true` if `abs(x/y)` is bounded between `0.99` and `1.01`.
"""
check_ratio_precision(x, y) = !iszero(y) && 0.99 < abs(x / y) < 1.01

"""
    get_boundary_chain(tri::Triangulation, i, j) -> Edges 

Given two boundary vertices `i` and `j` on a boundary with ghost vertex `ghost_vertex`, 
walks counter-clockwise from `i` to `j` along the boundary and returns the collection of all vertices encountered in 
counter-clockwise order.
"""
function get_boundary_chain(tri::Triangulation, i, j, ghost_vertex)
    I = integer_type(tri)
    chain = Vector{I}()
    push!(chain, i)
    w = I(∅)
    while w ≠ j
        w = get_adjacent(tri, i, ghost_vertex)
        add_edge!(chain, w)
        i = w
    end
    return chain
end

"""
    dist(tri::Triangulation, p) -> Number

Given a point `p`, returns the distance from `p` to the triangulation, using the 
conventions from [`distance_to_polygon`](@ref):

- `δ > 0`: If the returned distance is positive, then `p` is inside the triangulation.
- `δ < 0`: If the returned distance is negative, then `p` is outside the triangulation.
- `δ = 0`: If the returned distance is zero, then `p` is on the boundary of the triangulation.

Where we say distance, we are referring to the distance from `p` to the boundary of the triangulation.
"""
function dist(tri::Triangulation, p)
    points = get_points(tri)
    if has_boundary_nodes(tri)
        return distance_to_polygon(p, points, get_boundary_nodes(tri))
    else
        return distance_to_polygon(p, points, get_convex_hull_vertices(tri))
    end
end

"""
    adjust_θ(θ₁, θ₂, positive) -> (Number, Number)

Given two angles `θ₁` and `θ₂` in radians, adjusts the angles to new angles `θ₁′`, `θ₂′` so that 
`θ₁′ ≤ θ₂′` if `positive` is `true`, and `θ₁′ ≥ θ₂′` if `positive` is `false`.
"""
function adjust_θ(θ₁, θ₂, positive)
    if positive
        if θ₁ < θ₂
            return (θ₁, θ₂)
        elseif θ₁ > θ₂
            return (θ₁, θ₂ + 2π)
        else
            return (θ₁, θ₂ + 2π)
        end
    else
        if θ₁ < θ₂
            return (θ₁, θ₂ - 2π)
        elseif θ₁ > θ₂
            return (θ₁, θ₂)
        else
            return (θ₁, θ₂ - 2π)
        end
    end
end

"""
    uniquetol(A::Vector{Float64}; tol=1e-12) -> Vector{Float64}

Returns the unique elements of `A` up to tolerance `tol`. We say that two values `x` and `y` are within tolerance if `abs(u - v) ≤ M*tol`, where `M = maximum(abs.(A))`. It 
is assumed that `A` is sorted - this is NOT checked.
"""
function uniquetol(A::Vector{Float64}; tol = 1.0e-12)
    isempty(A) && return Float64[]
    M = max(abs(A[begin]), abs(A[end])) # assuming A is sorted
    intol = (x, y) -> abs(x - y) ≤ M * tol
    Auniq = [A[1]]
    for i in 2:lastindex(A)
        a = A[i]
        if !intol(a, Auniq[end])
            push!(Auniq, a)
        end
    end
    return Auniq
end

"""
    eval_fnc_at_het_tuple_element(f, tup, idx) 

Evaluates `f(tup[idx])` in a type-stable way. If `idx > length(tup)`, then `f` is evaluated on the last element of `tup`.
"""
@inline function eval_fnc_at_het_tuple_element(f::F, tup::T, idx) where {F, T}
    return _eval_fnc_at_het_tuple_element(f, idx, tup...)
end
@inline function _eval_fnc_at_het_tuple_element(f::F, idx, el::E, tup...) where {F, E}
    idx == 1 && return _eval_fnc_at_het_tuple_element(f, 1, el)
    return _eval_fnc_at_het_tuple_element(f, idx - 1, tup...)
end
@inline function _eval_fnc_at_het_tuple_element(f::F, idx, el::E) where {F, E}
    return f(el)
end

"""
    eval_fnc_at_het_tuple_two_elements(f, tup, idx1, idx2)

Evaluates `f(tup[idx1], tup[idx2])` in a type-stable way. 
"""
@inline function eval_fnc_at_het_tuple_two_elements(f::F, tup::T, idx1, idx2) where {F, T <: Tuple}
    return _eval_fnc_at_het_tuple_two_elements(f, idx2, tup, idx1, tup...)
end
@inline function _eval_fnc_at_het_tuple_two_elements(f::F, idx2, next_tup::T, idx1, el::E, tup...) where {F, E, T <: Tuple}
    idx1 == 1 && return _eval_fnc_at_het_tuple_two_elements(f, idx2, next_tup, 1, el)
    return _eval_fnc_at_het_tuple_two_elements(f, idx2, next_tup, idx1 - 1, tup...)
end
@inline function _eval_fnc_at_het_tuple_two_elements(f::F, idx2, next_tup::T, idx1, el::E) where {F, E, T <: Tuple}
    return _eval_fnc_at_het_tuple_two_elements(f, idx2, el, next_tup...)
end
@inline function _eval_fnc_at_het_tuple_two_elements(f::F, idx2, el::E, el2::V, tup...) where {F, E, V}
    idx2 == 1 && return _eval_fnc_at_het_tuple_two_elements(f, 1, el, el2)
    return _eval_fnc_at_het_tuple_two_elements(f, idx2 - 1, el, tup...)
end
@inline function _eval_fnc_at_het_tuple_two_elements(f::F, idx2, el::E, el2::V) where {F, E, V}
    return f(el, el2)
end

"""
    eval_fnc_at_het_tuple_element_with_arg(f, tup, arg, idx)

Evaluates `f(tup[idx], arg...)` in a type-stable way. If `idx > length(tup)`, then `f` is evaluated on the last element of `tup`.
"""
@inline function eval_fnc_at_het_tuple_element_with_arg(f::F, tup::T, arg, idx) where {F, T}
    return _eval_fnc_at_het_tuple_element_with_arg(f, idx, arg, tup...)
end
@inline function _eval_fnc_at_het_tuple_element_with_arg(f::F, idx, arg, el::E, tup...) where {F, E}
    idx == 1 && return _eval_fnc_at_het_tuple_element_with_arg(f, 1, arg, el)
    return _eval_fnc_at_het_tuple_element_with_arg(f, idx - 1, arg, tup...)
end
@inline function _eval_fnc_at_het_tuple_element_with_arg(f::F, idx, arg, el::E) where {F, E}
    return f(el, arg...)
end

"""
    eval_fnc_at_het_tuple_element_with_arg_and_prearg(f, tup, prearg, arg, idx)

Evaluates `f(prearg, tup[idx], arg...)` in a type-stable way. If `idx > length(tup)`, then `f` is evaluated on the last element of `tup`.
"""
@inline function eval_fnc_at_het_tuple_element_with_arg_and_prearg(f::F, tup::T, prearg, arg, idx) where {F, T}
    return _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f, idx, prearg, arg, tup...)
end
@inline function _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f::F, idx, prearg, arg, el::E, tup...) where {F, E}
    idx == 1 && return _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f, 1, prearg, arg, el)
    return _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f, idx - 1, prearg, arg, tup...)
end
@inline function _eval_fnc_at_het_tuple_element_with_arg_and_prearg(f::F, idx, prearg, arg, el::E) where {F, E}
    return f(prearg, el, arg...)
end

"""
    eval_fnc_in_het_tuple(tup, arg, idx)

Evaluates `tup[idx](arg...)` in a type-stable way. If `idx > length(tup)`, then `tup[end](arg...)` is evaluated.
"""
@inline function eval_fnc_in_het_tuple(tup::T, arg::A, idx) where {T, A}
    return eval_fnc_at_het_tuple_element_with_arg(self_eval, tup, arg, idx)
end

"""
    self_eval(f, args...)

Evaluates `f(args...)`.
"""
@inline self_eval(f, args...) = f(args...)

"""
    _to_val(v) -> Val 

Wraps `v` in a `Val`, or if `v isa Val` simply returns `v`.
"""
@inline _to_val(v::V) where {V} = Val(v)::Val{v}
@inline _to_val(v::Val{B}) where {B} = v
