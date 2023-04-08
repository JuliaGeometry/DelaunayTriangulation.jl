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
    get_right_boundary_node(adj::Adjacent{I,E}, k, boundary_index, boundary_index_ranges, check_existence::C) where {I,E,C}

Given an [`Adjacent`](@ref) map, a boundary node index `k`, a `boundary_index` corresponding to the curve, 
and `boundary_index_ranges` from [`construct_boundary_index_ranges`](@ref),
returns the node on the boundary to the right of `k`. `check_existence` can be used if you need to check 
over all boundary indices, in case there are multiple segments and thus multiple possible boundary indices 
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
over all boundary indices, in case there are multiple segments and thus multiple possible boundary indices 
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
is_circular(A) = isempty(A) || (A[begin] == A[end])

"""
    circular_equality(A, B)

Tests if the arrays `A` and `B` are equal up to a circular shift.
"""
function circular_equality(A, B)
    @assert is_circular(A) && is_circular(B) "The input arrays must satisfy x[begin] == x[end]."
    length(A) ≠ length(B) && return false
    length(A) == length(B) == 0 && return true
    _A = @views A[begin:(end-1)]
    _B = @views B[begin:(end-1)]
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

"""
    get_surrounding_polygon(adj::Adjacent{I,E}, graph::Graph, u, boundary_index_ranges, check_existence::C; skip_boundary_indices=false) where {I,E,C}

Given an [`Adjacent`](@ref) map, a [`Graph`](@ref), a vertex `u`, `boundary_index_ranges` from [`construct_boundary_index_ranges`](@ref) 
for handling the case where `u` is on the boundary, returns a vector `S` which gives a counter-clockwise sequence of the neighbours of `u`.
`check_existence` can be used if you need to check over all boundary indices when `u` is a boundary node, 
in case there are multiple segments and thus multiple possible boundary indices on the boundary.

When `u` is an outer boundary index, the returned polygon is clockwise.

When `u` is a boundary vertex and you do not have ghost triangles, then this function may return an invalid polygon.

If you want to remove all boundary indices from the result at the end, set `skip_boundary_indices=true`.
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

Given a set of `constrained_edges` and a `constrained_edge` in the set,
and a vector of segments `collinear_segments` that are collinear with `constrained_edge`,
replaces `constrained_edge` with those segments in `collinear_segments`. 
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

If we have edges in `segments` that are `bad`, meaning are not valid edges, found via 
`bad_indices`, this function fixes them.

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
3-element Vector{Tuple{Int64, Int64}}:
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

```julia-repl
julia> C = [(7, 12), (12, 17), (17, 22), (32, 37), (37, 42), (42, 47)];

julia> DelaunayTriangulation.connect_segments!(C);

julia> C
7-element Vector{Tuple{Int64, Int64}}:
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
    insert_idxs = Tuple{Int64,I,I}[]
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
    convert_boundary_points_to_indices(x, y; existing_points = NTuple{2, Float64}[])

Given some points `(x, y)` representing a boundary, converts their representation into a 
a set of indices corresponding to each boundary. The points should match the specification 
given in [`generate_mesh`](@ref). These points also get appended onto the `existing_points` keyword 
argument, which should be used if you have a pre-existing set of points so that the  
boundary indices get adjusted accordingly. Points get added into `existing_points` via 
[`push_point!`](@ref).

The returned value is `(nodes, points)`, with `nodes` the indices and `points` the modified 
`existing_points` (which are mutated in-place also).
"""
function convert_boundary_points_to_indices(x::AAA, y::AAA; existing_points=NTuple{2,Float64}[], check_args=true, adjust=true) where {F<:Number,A<:AbstractVector{F},AA<:AbstractVector{A},AAA<:AbstractVector{AA}}
    check_args && @assert length(x) == length(y)
    nodes = [[Int64[] for _ in eachindex(x[i])] for i in eachindex(x)]
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
    nodes = [Int64[] for _ in eachindex(x)]
    for i in eachindex(x)
        _nodes, _ = convert_boundary_points_to_indices(x[i], y[i]; existing_points=existing_points, check_args=false, adjust=false)
        resize!(nodes[i], length(_nodes))
        copyto!(nodes[i], _nodes)
    end
    if adjust
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
    nodes = Int64[]
    init = num_points(existing_points) + 1
    for i in firstindex(x):(lastindex(x)-1)
        push!(nodes, init)
        push_point!(existing_points, x[i], y[i])
        init += 1
    end
    adjust && push!(nodes, nodes[begin])
    return nodes, existing_points
end

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
end