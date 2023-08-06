"""
    is_true(b)

Returns `true` if `b` is `true`, `Val{true}`, or `Val(true)`. Returns `false` otherwise.
"""
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

"""
    get_boundary_index(i, j, k)
    get_boundary_index(i, j)

Given three indices `i`, `j`, and `k`, returns the index corresponding to a boundary index. If no boundary index is provided, an `ArgumentError` is thrown.
Similarly for the second method, which takes two indices.
"""
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

"""
    rotate_ghost_triangle_to_standard_form(i, j, k)
    rotate_ghost_triangle_to_standard_form(T::V) where {V}

Given a triangle `T = (i, j, k)`, rotates it to a new triangle `T′ = (u, v, w)`
such that `w` is a boundary index.
"""
function rotate_ghost_triangle_to_standard_form(i, j, k) # ghost index last
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
    rotate_triangle_to_standard_form(i, j, k)
    rotate_triangle_to_standard_form(T::V) where {V}

Given a triangle `T = (i, j, k)`, rotates it to a new triangle `T′ = (u, v, w)`
such that `w` is the minimum index.
"""
function rotate_triangle_to_standard_form(i, j, k) # minimum index last 
    min_ijk = min(i, j, k)
    if min_ijk == i
        return (j, k, i)
    elseif min_ijk == j
        return (k, i, j)
    else
        return (i, j, k)
    end
end
function rotate_triangle_to_standard_form(T::V) where {V}
    i, j, k = indices(T)
    u, v, w = rotate_triangle_to_standard_form(i, j, k)
    return construct_triangle(V, u, v, w)
end

"""
    get_right_boundary_node(adj::Adjacent{I,E}, k, boundary_index, boundary_index_ranges, check_existence::C) where {I,E,C}

Returns the node on the boundary that is to the right of `k`.

# Arguments 
- `adj::Adjacent`: The [`Adjacent`](@ref) map.
- `k`: The boundary node index.
- `boundary_index`: The boundary index corresponding to the curve.
- `boundary_index_ranges`: The boundary index ranges from [`construct_boundary_index_ranges`](@ref).
- `check_existence::C`: Whether to check over all boundary indices, in case there are multiple segments and thus multiple possible boundary indices on the boundary.

# Outputs 
- `i`: The node on the boundary to the right of `k`.

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

Returns the node on the boundary that is to the left of `k`.

# Arguments
- `adj::Adjacent`: The [`Adjacent`](@ref) map.
- `k`: The boundary node index.
- `boundary_index`: The boundary index corresponding to the curve.
- `boundary_index_ranges`: The boundary index ranges from [`construct_boundary_index_ranges`](@ref).
- `check_existence::C`: Whether to check over all boundary indices, in case there are multiple segments and thus multiple possible boundary indices on the boundary.

# Outputs
- `i`: The node on the boundary to the left of `k`.

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

Given a triangle `T` with indices corresponding to `points`, returns the edge of `T` that contains the point `ℓ`.
It is assumed that the point `ℓ` is on an edge of `T`. If this is not the case, an error is thrown.
"""
function find_edge(T, points, ℓ)
    r = get_point(points, ℓ)
    for (u, v) in triangle_edges(T)
        if !is_ghost_edge(u, v)
            p, q = get_point(points, u, v)
            cert = point_position_relative_to_line(p, q, r)
            is_collinear(cert) && return (u, v)
        end
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

Tests if `A[begin] == A[end]`. Also returns `true` if `A` is empty.
"""
is_circular(A) = isempty(A) || (A[begin] == A[end])

"""
    circular_equality(A, B, by = isequal)

Tests if the arrays `A` and `B` are equal up to a circular shift, assuming `A` 
and `B` are circular. The function `by` is used to test equality of elements.
"""
function circular_equality(A, B, by=isequal)
    @assert is_circular(A) && is_circular(B) "The input arrays must satisfy x[begin] == x[end]."
    length(A) ≠ length(B) && return false
    length(A) == length(B) == 0 && return true
    _A = @views A[begin:(end-1)]
    _B = @views B[begin:(end-1)]
    same_idx = findfirst(by(_A[begin]), _B)
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
    get_surrounding_polygon(adj::Adjacent{I,E}, graph::Graph, u, boundary_index_ranges, check_existence::C; skip_boundary_indices=false) where {I,E,C}

Given a point `u`, returns a vector `S` which gives a counter-clockwise sequence of the neighbours of `u`. 

# Arguments 
- `adj::Adjacent{I,E}`: The [`Adjacent`](@ref) map.
- `graph::Graph`: The [`Graph`](@ref).
- `u`: The vertex.
- `boundary_index_ranges`: The output of [`construct_boundary_index_ranges`](@ref).
- `check_existence::C`: Whether to check over all boundary indices, in case there are multiple segments and thus multiple possible boundary indices on the boundary.

# Keyword Arguments
- `skip_boundary_indices=false`: Whether to remove all boundary indices from the result at the end.

# Outputs 
- `S`: The surrounding polygon.

!!! note 

    - When `u` is an outer boundary index, the returned polygon is clockwise.
    - When `u` is a boundary vertex and you do not have ghost triangles, then this function may return an invalid polygon.
"""
function get_surrounding_polygon(adj::Adjacent{I,E}, graph::Graph, u, boundary_index_ranges, check_existence::C; skip_boundary_indices=false) where {I,E,C}
    neighbouring_vertices = get_neighbours(graph, u)
    v = first(neighbouring_vertices)
    if !is_boundary_index(u)
        k = num_neighbours(graph, u)
        nbnd = count(is_boundary_index, neighbouring_vertices)
        if nbnd > 0
            # If we are on a boundary that has multiple boundary indices, k will not 
            # be counted correctly. This step makes the adjustment.
            k = k - nbnd + 1 # + 1 because we do still want one boundary index
        end
    else
        # Need to be careful when there are multiple boundary indices - get_neighbours would only 
        # return the segment for a given boundary index, so we need to check all. The nodes do not uniquely 
        # appear in a single segment, so we cannot just keep adding to d above. To keep things 
        # a bit simple, we just take a union. This could be improved.
        boundary_index_range = map_boundary_index(boundary_index_ranges, u)
        for index in boundary_index_range
            neighbouring_vertices = neighbouring_vertices ∪ get_neighbours(graph, index)
        end
        k = length(neighbouring_vertices)
    end
    S = zeros(I, k)
    S[1] = v
    for i in 2:k
        v = get_adjacent(adj, u, v; check_existence, boundary_index_ranges)
        S[i] = v
    end
    skip_boundary_indices && filter!(!is_boundary_index, S)
    return S
end

"""
    sort_edge_by_degree(e::E, graph::Graph)

Given an edge `e` of a `graph`, say `e = (u, v)`,
returns:

- If `deg(u) ≤ deg(v)`, returns `e`;
- If `deg(u) > deg(v)`, returns `(v, u)`.

In particular, `e` is sorted so that `initial(e)` is the vertex of `e` 
with the smallest degree.
"""
function sort_edge_by_degree(e::E, graph::Graph) where {E}
    u = initial(e)
    v = terminal(e)
    d₁ = num_neighbours(graph, u)
    d₂ = num_neighbours(graph, v)
    if d₁ ≤ d₂
        return e
    else
        return construct_edge(E, v, u)
    end
end

"""
    split_constrained_edge!(constrained_edges, constrained_edge::E, collinear_segments) where {E}

Splits the `constrained_edge` at the segments in `collinear_segments`, updating `constrained_edges` accordingly.
"""
function split_constrained_edge!(constrained_edges, constrained_edge::E, collinear_segments) where {E}
    if num_edges(collinear_segments) > 0
        flipped_constrained_edge = reverse_edge(constrained_edge)
        delete_edge!(constrained_edges, constrained_edge)
        delete_edge!(constrained_edges, flipped_constrained_edge) # Delete both, in case we have reverse(e) rather than e
        for e in each_edge(collinear_segments)
            reverse_e = reverse_edge(e)
            if !(contains_edge(e, constrained_edges) || contains_edge(reverse_e, constrained_edges))
                add_edge!(constrained_edges, e)
            end
        end
    end
    return nothing
end

"""
    fix_segments!(segments::AbstractVector{E}, bad_indices) where {E}

Fixes the overlapping segments in `segments`, referred to via `bad_indices`.

## Example 

For example, if we had 

```julia-repl
julia> c = [(2, 15), (2, 28), (2, 41)]
```

then these edges come from connecting the start of a constrained segment with a point that 
it goes through, but they are not actual segments in the triangulation (because they 
all start with 2). So, using `bad_indices = [1, 2, 3]`, the function mutates `c` 
to give 

```julia-repl
julia> bad_indices = [1, 2, 3]
julia> fix_segments!(c, bad_indices)
julia> c
3-element Vector{Tuple{Int, Int}}:
 (2, 15)
 (15, 28)
 (28, 41)
```
"""
function fix_segments!(segments::AbstractVector{E}, bad_indices) where {E}
    for i in bad_indices
        if i == firstindex(segments) # If it starts with a bad index, then no problem, it connects two valid indices.
            continue
        else
            prev = segments[i-1]
            cur = segments[i]
            segments[i] = construct_edge(E, terminal(prev), terminal(cur))
        end
    end
    return nothing
end

"""
    connect_segments!(segments::AbstractVector{E}) where {E}

Given an ordered vector of `segments`, mutates so that the endpoints connect, preserving order.

## Example

```julia-repl
julia> C = [(7, 12), (12, 17), (17, 22), (32, 37), (37, 42), (42, 47)];

julia> DelaunayTriangulation.connect_segments!(C);

julia> C
7-element Vector{Tuple{Int, Int}}:
 (7, 12)
 (12, 17)
 (17, 22)
 (32, 37)
 (37, 42)
 (42, 47)
 (22, 32)
```
"""
function connect_segments!(segments::AbstractVector{E}) where {E}
    f = firstindex(segments)
    ℓ = lastindex(segments)
    I = integer_type(E)
    insert_idxs = Tuple{Int,I,I}[]
    shift_idx = 0
    for i in f:(ℓ-1)
        eᵢ = segments[i]
        eᵢ₊₁ = segments[i+1]
        if terminal(eᵢ) ≠ initial(eᵢ₊₁)
            push!(insert_idxs, (i + 1 + shift_idx, terminal(eᵢ), initial(eᵢ₊₁)))
            shift_idx += 1
        end
    end
    for (i, u, v) in insert_idxs
        insert!(segments, i, construct_edge(E, u, v)) # insert might be a bit slow, but segments should never be large so it's fine
    end
    return nothing
end

"""
    extend_segments!(segments::AbstractVector{E}, constrained_edge) where {E}

Given an ordered vector of `segments`, ensures that they also represent the 
replacement of `constrained_edge`.

## Example 

```julia-repl
julia> segments = [(2, 7), (7, 12), (12, 49)];
julia> constrained_edge = (1, 68);
julia> extend_segments!(segments, constrained_edge);
julia> segments
5-element Vector{Tuple{Int, Int}}:
 (1, 2)
 (2, 7)
 (7, 12)
 (12, 49)
 (49, 68)
```
"""
function extend_segments!(segments::AbstractVector{E}, constrained_edge) where {E}
    constrained_i, constrained_j = edge_indices(constrained_edge)
    first_i = initial(segments[begin])
    last_j = terminal(segments[end])
    if constrained_i ≠ first_i
        pushfirst!(segments, construct_edge(E, constrained_i, first_i))
    end
    if constrained_j ≠ last_j
        push!(segments, construct_edge(E, last_j, constrained_j))
    end
    return nothing
end

"""
    convert_boundary_points_to_indices(x, y; existing_points = NTuple{2, Float64}[], check_args=true, adjust=true)
    convert_boundary_points_to_indices(xy; existing_points = NTuple{2, Float64}[], check_args=true, adjust=true)

Given some points `(x, y)` representing a boundary, or `xy` giving the points combined rather than separated, converts their representation into a set of 
indices corresponding to each boundary. The points should match the specification of a boundary 
defined in the documentation. These points also get appended onto the set of points given by the 
`existing_points` keyword argument, which should be used if you have a pre-existing set of points.

The returned value is `(nodes, points)`, with `nodes` the indices and `points` the modified 
`existing_points` (which are mutated in-place also).
"""
function convert_boundary_points_to_indices(x::AAA, y::AAA; existing_points=NTuple{2,Float64}[], check_args=true, adjust=true) where {F<:Number,A<:AbstractVector{F},AA<:AbstractVector{A},AAA<:AbstractVector{AA}}
    check_args && @assert length(x) == length(y)
    nodes = [[Int[] for _ in eachindex(x[i])] for i in eachindex(x)]
    for i in eachindex(x)
        _nodes, _ = convert_boundary_points_to_indices(x[i], y[i]; existing_points=existing_points, check_args=true)
        copyto!(nodes[i], _nodes)
    end
    return nodes, existing_points
end
function convert_boundary_points_to_indices(x::AA, y::AA; existing_points=NTuple{2,Float64}[], check_args=true, adjust=true) where {F<:Number,A<:AbstractVector{F},AA<:AbstractVector{A}}
    if check_args
        @assert length(x) == length(y)
        @assert all(i -> length(x[i]) == length(y[i]), eachindex(x, y))
        @assert x[begin][begin] ≈ x[end][end]
        @assert y[begin][begin] ≈ y[end][end]
        @assert all(i -> x[i][end] ≈ x[i+1][begin], firstindex(x):(lastindex(x)-1))
        @assert all(i -> y[i][end] ≈ y[i+1][begin], firstindex(y):(lastindex(y)-1))
    end
    nodes = [Int[] for _ in eachindex(x)]
    for i in eachindex(x)
        _nodes, _ = convert_boundary_points_to_indices(x[i], y[i]; existing_points=existing_points, check_args=false, adjust=false)
        resize!(nodes[i], length(_nodes))
        copyto!(nodes[i], _nodes)
    end
    if adjust # needed so that the different segments connect
        for i in firstindex(nodes):(lastindex(nodes)-1)
            push!(nodes[i], nodes[i+1][begin])
        end
        push!(nodes[end], nodes[begin][begin])
    end
    return nodes, existing_points
end
function convert_boundary_points_to_indices(x::A, y::A; existing_points=NTuple{2,Float64}[], check_args=true, adjust=true) where {F<:Number,A<:AbstractVector{F}}
    if check_args
        @assert length(x) == length(y)
        @assert x[begin] ≈ x[end]
        @assert y[begin] ≈ y[end]
    end
    nodes = Int[]
    init = num_points(existing_points) + 1
    for i in firstindex(x):(lastindex(x)-1)
        push!(nodes, init)
        push_point!(existing_points, x[i], y[i])
        init += 1
    end
    adjust && push!(nodes, nodes[begin])
    return nodes, existing_points
end
function convert_boundary_points_to_indices(xy::A; existing_points=NTuple{2,Float64}[], check_args=true, adjust=true) where {F,A<:AbstractVector{F}}
    x = [getx(xy[i]) for i in eachindex(xy)]
    y = [gety(xy[i]) for i in eachindex(xy)]
    return convert_boundary_points_to_indices(x, y; existing_points=existing_points, check_args=check_args, adjust=adjust)
end
function convert_boundary_points_to_indices(xy::AA; existing_points=NTuple{2,Float64}[], check_args=true, adjust=true) where {F,A<:AbstractVector{F},AA<:AbstractVector{A}}
    x = [[getx(xy[i][j]) for j in eachindex(xy[i])] for i in eachindex(xy)]
    y = [[gety(xy[i][j]) for j in eachindex(xy[i])] for i in eachindex(xy)]
    return convert_boundary_points_to_indices(x, y; existing_points=existing_points, check_args=check_args, adjust=adjust)
end
function convert_boundary_points_to_indices(xy::AAA; existing_points=NTuple{2,Float64}[], check_args=true, adjust=true) where {F,A<:AbstractVector{F},AA<:AbstractVector{A},AAA<:AbstractVector{AA}}
    x = [[[getx(xy[i][j][k]) for k in eachindex(xy[i][j])] for j in eachindex(xy[i])] for i in eachindex(xy)]
    y = [[[gety(xy[i][j][k]) for k in eachindex(xy[i][j])] for j in eachindex(xy[i])] for i in eachindex(xy)]
    return convert_boundary_points_to_indices(x, y; existing_points=existing_points, check_args=check_args, adjust=adjust)
end

"""
    get_ordinal_suffix(i)

Returns the ordinal suffix for the given integer `i`. 

## Example

```julia-repl
julia> get_ordinal_suffix(1)
"st"
julia> get_ordinal_suffix(2)
"nd"
julia> get_ordinal_suffix(3)
"rd"    
julia> get_ordinal_suffix(4)
"th"
julia> get_ordinal_suffix(11)
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
    check_args(points, boundary_nodes)

Checks the arguments `points` and `boundary_nodes` to make sure that they are valid. If they are
not, an error is thrown. This function is called by `triangulate` if the `check_args` keyword
argument is set to `true`. If you are sure that your arguments are valid, you can set this
keyword argument to `false` to speed up the triangulation process.
"""
function check_args(points, boundary_nodes)
    @assert points_are_unique(points) "Duplicate points are not allowed."
    if !isnothing(boundary_nodes)
        if has_multiple_curves(boundary_nodes)
            areas = [polygon_features(points, get_boundary_nodes(boundary_nodes, i))[1] for i in 1:num_curves(boundary_nodes)]
            @assert areas[1] ≥ 0.0 "The outer boundary curve is clockwise when it should be counter-clockwise. "
            for i in 2:num_curves(boundary_nodes)
                @assert areas[i] ≤ 0.0 "The $(i)$(get_ordinal_suffix(i)) boundary curve is counter-clockwise when it should be clockwise. If this is a mistake, e.g. if this curve is inside of another one in which case it should be counter-clockwise, recall triangulate with check_arguments = false."
            end
            for i in 1:num_curves(boundary_nodes)
                curve_nodes = get_boundary_nodes(boundary_nodes, i)
                ns = num_segments(curve_nodes)
                for j in 1:(ns-1)
                    segment_nodes_cur = get_boundary_nodes(curve_nodes, j)
                    segment_nodes_next = get_boundary_nodes(curve_nodes, j + 1)
                    nnodes_cur = num_boundary_edges(segment_nodes_cur) + 1
                    @assert get_boundary_nodes(segment_nodes_cur, nnodes_cur) == get_boundary_nodes(segment_nodes_next, 1) "The $(j)$(get_ordinal_suffix(j)) segment of the $(i)$(get_ordinal_suffix(i)) curve does not connect with the start of the $(j+1)$(get_ordinal_suffix(j+1)) segment."
                end
                segment_nodes_first = get_boundary_nodes(curve_nodes, 1)
                segment_nodes_last = get_boundary_nodes(curve_nodes, ns)
                nnodes_last = num_boundary_edges(segment_nodes_last) + 1
                @assert get_boundary_nodes(segment_nodes_first, 1) == get_boundary_nodes(segment_nodes_last, nnodes_last) "The first segment of the $(i)$(get_ordinal_suffix(i)) curve does not connect with the end of the last segment."
            end
        else
            area = polygon_features(points, boundary_nodes)[1]
            @assert area ≥ 0.0 "The boundary curve is clockwise when it should be counter-clockwise."
            if has_multiple_segments(boundary_nodes)
                ns = num_segments(boundary_nodes)
                for j in 1:(ns-1)
                    segment_nodes_cur = get_boundary_nodes(boundary_nodes, j)
                    segment_nodes_next = get_boundary_nodes(boundary_nodes, j + 1)
                    nnodes_cur = num_boundary_edges(segment_nodes_cur) + 1
                    @assert get_boundary_nodes(segment_nodes_cur, nnodes_cur) == get_boundary_nodes(segment_nodes_next, 1) "The $(j)$(get_ordinal_suffix(j)) segment does not connect with the start of the $(j+1)$(get_ordinal_suffix(j+1)) segment."
                end
                segment_nodes_first = get_boundary_nodes(boundary_nodes, 1)
                segment_nodes_last = get_boundary_nodes(boundary_nodes, ns)
                nnodes_last = num_boundary_edges(segment_nodes_last) + 1
                @assert get_boundary_nodes(segment_nodes_first, 1) == get_boundary_nodes(segment_nodes_last, nnodes_last) "The first segment of the $(i)$(get_ordinal_suffix(i)) curve does not connect with the end of the last segment."
            else
                nnodes = num_boundary_edges(boundary_nodes) + 1
                @assert get_boundary_nodes(boundary_nodes, 1) == get_boundary_nodes(boundary_nodes, nnodes) "The first boundary node does not equal the last boundary node."
            end
        end
    end
    return true
end

"""
    min_max(a, b)

Returns `(min(a, b), max(a, b))`.
"""
function min_max(a, b)
    if b < a
        return b, a
    else
        return a, b
    end
end

"""
    min_med_max(a, b, c)

Returns `(min(a, b, c), med(a, b, c), max(a, b, c)), where `med(a, b, c)`
is the value that is neither `min(a, b, c)` or `max(a, b, c)`.
"""
function min_med_max(a, b, c)
    b, c = min_max(b, c)
    a, c = min_max(a, c)
    a, b = min_max(a, b)
    return a, b, c
end

"""
    balanced_power_of_two_ternary_split(ℓ)

Compute the the power of two that is closest 
to `ℓ/3` or to `ℓ/1.5`.
"""
function balanced_power_of_two_ternary_split(ℓ)
    balanced_split = 1.0
    while ℓ > 3balanced_split
        balanced_split = 2balanced_split
    end
    while ℓ < 1.5balanced_split
        balanced_split = 0.5balanced_split
    end
    return balanced_split
end

"""
    balanced_power_of_two_quarternary_split(ℓ)

Compute the the power of two that is closest
to `ℓ/4` or to `ℓ/0.5`.
"""
function balanced_power_of_two_quarternary_split(ℓ)
    balanced_split = 1.0
    while ℓ > 4balanced_split
        balanced_split = 2balanced_split
    end
    while ℓ < 2balanced_split
        balanced_split = 0.5balanced_split
    end
    return balanced_split
end

function nearest_power_of_two(ℓ)
    ispow2(ℓ) && return ℓ
    z = log2(ℓ)
    z⁺ = ceil(z)
    z⁻ = floor(z)
    if abs(ℓ - 2^z⁺) < abs(ℓ - 2^z⁻)
        return 2^z⁺
    else
        return 2^z⁻
    end
end

"""
    segment_vertices_adjoin_other_segments(tri::Triangulation, e)

Test if the segment `e`'s vertices adjoin other segments. Returns:

- `0`: No vertex adjoins another segment.
- `1`: One vertex adjoins another segment.
- `2`: Both vertices adjoin another segment.
"""
function segment_vertices_adjoin_other_segments(tri::Triangulation, e)
    u, v = edge_indices(e)
    count = 0
    for w in get_neighbours(tri, u)
        if w ≠ v
            if contains_constrained_edge(tri, u, w)
                count += 1
                break
            end
        end
    end
    for w in get_neighbours(tri, v)
        if w ≠ u
            if contains_constrained_edge(tri, v, w)
                count += 1
                break
            end
        end
    end
    return count
end

"""
    edge_lies_on_two_distinct_segments(tri::Triangulation, e)

Tests if the edge `(i, j)` lies on two distinct segments. The returned value is:

- `(true, common_vertex)`: If `e` lies on two distinct segments, and the common vertex is `common_vertex`.
- `(false, 0)`: Otherwise.

If there are multiple common vertices. In this case, the function returns the vertex that is closest to `e`.
"""
function edge_lies_on_two_distinct_segments(tri::Triangulation, i, j)
    # Need to eventually make this non-allocating by being a bit more 
    # clever with the geometry
    I = integer_type(tri)
    i_segments = Set{I}()
    for u in get_neighbours(tri, i)
        if u ≠ j
            if contains_constrained_edge(tri, i, u)
                push!(i_segments, u)
            end
        end
    end
    isempty(i_segments) && return false, I(0)
    j_segments = Set{I}()
    for u in get_neighbours(tri, j)
        if u ≠ i
            if contains_constrained_edge(tri, j, u)
                push!(j_segments, u)
            end
        end
    end
    isempty(j_segments) && return false, I(0)
    intersect!(i_segments, j_segments)
    if length(i_segments) == 0
        return false, I(0)
    elseif length(i_segments) == 1
        return true, first(i_segments)
    else
        F = number_type(tri)
        p, q = get_point(tri, i, j)
        px, py = _getxy(p)
        qx, qy = _getxy(q)
        min_dist = typemax(F)
        min_idx = I(0)
        for k in i_segments
            r = get_point(tri, k)
            rx, ry = _getxy(r)
            δ = squared_distance_to_segment(px, py, qx, qy, rx, ry)
            if δ < min_dist
                min_dist = δ
                min_idx = k
            end
        end
        return true, min_idx
    end
end

"""
    nextindex_circular(C, i)

Returns the next index in the collection `C` after `i`, wrapping around to the first index if `i` is the last index.
"""
function nextindex_circular(C, i)
    return i == lastindex(C) ? firstindex(C) : i + 1
end

"""
    previndex_circular(C, i)

Returns the previous index in the collection `C` before `i`, wrapping around to the second-last index if `i` is the first 
index, assuming that `is_circular(C)` so that `C[begin] == C[end]`.
"""
function previndex_circular(C, i)
    return i == firstindex(C) ? lastindex(C) - 1 : i - 1
end

"""
    is_first_boundary_index(cell, i)

Returns `true` if the index `cell[i]` is the first boundary index in the cell `cell`, assuming they come as a chain.
"""
function is_first_boundary_index(cell, i)
    prev = previndex_circular(cell, i)
    return !is_boundary_index(cell[prev])
end

"""
    is_last_boundary_index(cell, i)

Returns `true` if the index `cell[i]` is the last boundary index in the cell `cell`, assuming they come as a chain.
"""
function is_last_boundary_index(cell, i)
    prev = previndex_circular(cell, i)
    return is_boundary_index(cell[prev])
end

"""
    get_neighbouring_boundary_edges(tri::Triangulation, e)

Given an edge `e` on the boundary, returns the edges to the left 
and to the right of `e`, oriented so that `get_adjacent(tri, v)` 
is a boundary index, where `v` are the edges returned.
"""
function get_neighbouring_boundary_edges(tri::Triangulation, e)
    e = convert_to_boundary_edge(tri, e)
    u, v = edge_indices(e)
    bnd_idx = get_adjacent(tri, u, v)
    right_bnd = get_right_boundary_node(tri, u, bnd_idx)
    left_bnd = get_left_boundary_node(tri, v, bnd_idx)
    E = edge_type(tri)
    left_e = construct_edge(E, v, left_bnd)
    right_e = construct_edge(E, right_bnd, u)
    return left_e, right_e
end

"""
    convert_to_boundary_edge(tri::Triangulation, e)

Given an edge `e`, returns the edge that is on the boundary, oriented so that `get_adjacent(tri, v)`
is a boundary index, where `v` is the edge returned.
"""
function convert_to_boundary_edge(tri::Triangulation, e)
    if !is_boundary_edge(tri, e)
        return reverse_edge(e)
    else
        return e
    end
end

"""
    get_shared_vertex(e, f)

Given two edges `e` and `f`, returns the vertex that they share, or `DefaultAdjacentValue` if they do not share a vertex.
"""
function get_shared_vertex(e, f)
    u, v = edge_indices(e)
    w, x = edge_indices(f)
    if u == w || u == x
        return u
    elseif v == w || v == x
        return v
    else
        I = typeof(u)
        return I(DefaultAdjacentValue)
    end
end

"""
    replace_boundary_triangle_with_ghost_triangle(tri, V)

Given a triangulation `tri` and a boundary triangle `V`, returns the ghost triangle associated with the boundary edge.
"""
function replace_boundary_triangle_with_ghost_triangle(tri, V)
    u, v, w = indices(V)
    T = triangle_type(tri)
    is_boundary_edge(tri, v, u) && return construct_triangle(T, v, u, get_adjacent(tri, v, u))
    is_boundary_edge(tri, w, v) && return construct_triangle(T, w, v, get_adjacent(tri, w, v))
    return construct_triangle(T, u, w, get_adjacent(tri, u, w))
end

"""
    iterated_neighbourhood(tri, i, d)

Computes the `d`-times iterated neighbourhood of `i` in the triangulation `tri`. In particular,
this returns all indices that are within `d` edges of `i`, excluding `i` itself.
"""
function iterated_neighbourhood(tri, i, d)
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
        if !is_boundary_index(j)
            push!(neighbours, j)
        end
    end
    for _ in 2:d
        new_neighbours = Set{I}() # don't want to mutate the iterator while iterating
        for j in neighbours
            for k in get_neighbours(tri, j)
                if k ≠ i && !is_boundary_index(k)
                    push!(new_neighbours, k)
                end
            end
        end
        union!(neighbours, new_neighbours)
    end
    return neighbours
end