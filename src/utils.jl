is_true(b::Bool) = b
is_true(b::Val{true}) = true
is_true(b::Val{false}) = false
is_true(b::Type{Val{true}}) = true
is_true(b::Type{Val{false}}) = false

"""
    number_type(x)

Given a container `x`, returns the number type used for storing coordinates.
"""
number_type(x) = number_type(typeof(x))
number_type(::Type{T}) where {T<:AbstractArray} = number_type(eltype(T))
number_type(::Type{NTuple{N,T}}) where {N,T} = T
number_type(::Type{T}) where {T<:Number} = T

function get_boundary_index(i, j, k)
    is_boundary_index(i) && return i
    is_boundary_index(j) && return j
    is_boundary_index(k) && return k
    throw(ArgumentError("No boundary indices provided."))
end
function get_boundary_index(i, j)
    is_boundary_index(i) && return i
    is_boundary_index(j) && return j
    throw(ArgumentError("No boundary indices provided."))
end

function rotate_ghost_triangle_to_standard_form(i, j, k)
    if is_boundary_index(i)
        return (j, k, i)
    elseif is_boundary_index(j)
        return (k, i, j)
    elseif is_boundary_index(k)
        return (i, j, k)
    end
end
function rotate_ghost_triangle_to_standard_form(T::V) where {V}
    i, j, k = indices(T)
    u, v, w = rotate_ghost_triangle_to_standard_form(i, j, k)
    return construct_triangle(V, u, v, w)
end

"""
    get_right_boundary_node(adj::Adjacent{I,E}, k, boundary_index, boundary_index_ranges, check_existence::C) where {I,E,C}

Given an [`Adjacent`](@ref) map, a boundary node index `k`, a `boundary_index` corresponding to the curve, 
and `boundary_index_ranges` from [`construct_boundary_index_ranges`](@ref),
returns the node on the boundary to the right of `k`. `check_existence` can be used if you need to check 
over all outer boundary indices, in case there are multiple segments and thus multiple possible boundary indices 
on the boundary.

See also [`get_left_boundary_node`](@ref).
"""
function get_right_boundary_node(adj::Adjacent{I,E}, k, boundary_index,
                                 boundary_index_ranges, check_existence::C) where {I,E,C}
    if is_true(check_existence)
        boundary_index_range = map_boundary_index(boundary_index_ranges, boundary_index)
        for index in boundary_index_range
            i = get_adjacent(adj, k, index)
            edge_exists(i) && return i
        end
        throw("No boundary node found.")
    else
        return get_adjacent(adj, k, boundary_index)
    end
end

"""
    get_left_boundary_node(adj::Adjacent{I,E}, k, boundary_index, boundary_index_ranges, check_existence::C) where {I,E,C}

Given an [`Adjacent`](@ref) map, a boundary node index `k`, a `boundary_index` corresponding to the curve, 
and `boundary_index_ranges` from [`construct_boundary_index_ranges`](@ref),
returns the node on the boundary to the left of `k`. `check_existence` can be used if you need to check 
over all outer boundary indices, in case there are multiple segments and thus multiple possible boundary indices 
on the boundary.

See also [`get_right_boundary_node`](@ref).
"""
function get_left_boundary_node(adj::Adjacent{I,E}, k, boundary_index,
                                boundary_index_ranges, check_existence::C) where {I,E,C}
    if is_true(check_existence)
        boundary_index_range = map_boundary_index(boundary_index_ranges, boundary_index)
        for index in boundary_index_range
            i = get_adjacent(adj, index, k)
            edge_exists(i) && return i
        end
        throw("No boundary node found.")
    else
        return get_adjacent(adj, boundary_index, k)
    end
end

"""
    find_edge(T, points, ℓ)

Given a triangle `T` and a set of points `points`, with 
the `ℓ`th point of `points` on an edge of `T`, returns 
the edge `(u, v)` that the point is on.
"""
function find_edge(T, points, ℓ)
    r = get_point(points, ℓ)
    for (u, v) in triangle_edges(T)
        p, q = get_point(points, u, v)
        cert = point_position_relative_to_line(p, q, r)
        is_collinear(cert) && return (u, v)
    end
    throw("The point $(get_point(points, ℓ)) is not on an edge of $T.")
end

"""
    choose_uvw(e1, e2, e3, u, v, w)

Choose values for `(u, v, w)` based on the Booleans `(e1, e2, e3)`, 
assuming only one is true. The three cases are: 

- If `e1`, returns `(u, v, w)`.
- If `e2`, returns `(v, w, u)`.
- If `e3`, returns `(w, u, v)`.
"""
function choose_uvw(e1, e2, e3, u, v, w)
    e1 && return (u, v, w)
    e2 && return (v, w, u)
    e3 && return (w, u, v)
    throw(ArgumentError("All the Booleans were false."))
end

"""
    is_circular(A)

Tests if `A[begin] == A[end]`.
"""
is_circular(A) = A[begin] == A[end]

"""
    circular_equality(A, B)

Tests if the arrays `A` and `B` are equal up to a circular shift.
"""
function circular_equality(A, B)
    @assert is_circular(A) && is_circular(B) "The input arrays must satisfy x[begin] == x[end]."
    length(A) ≠ length(B) && return false
    _A = @views A[begin:(end - 1)]
    _B = @views B[begin:(end - 1)]
    same_idx = findfirst(==(_A[begin]), _B)
    same_idx === nothing && return false
    n = length(_A)
    for (i, a) in pairs(_A)
        j = mod1(i + same_idx - 1, n)
        b = _B[j]
        a ≠ b && return false
    end
    return true
end
