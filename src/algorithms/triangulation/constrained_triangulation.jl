"""
    fix_segments!(segments, bad_indices)

Fixes the overlapping segments in `segments`, referred to via `bad_indices`, by connecting consecutive edges where needed.

# Arguments 
- `segments`: The segments to fix.
- `bad_indices`: The indices of the segments to fix.

# Outputs
There are no outputs as `segments` is updated in-place.

# Example 
```jldoctest
julia> using DelaunayTriangulation

julia> segments = [(2, 15), (2, 28), (2, 41)]; # the edges all start with 2, so they are not actual segments in the triangulation, and so must be fixed

julia> bad_indices = [1, 2, 3];

julia> DelaunayTriangulation.fix_segments!(segments, bad_indices)
3-element Vector{Tuple{Int64, Int64}}:
 (2, 15)
 (15, 28)
 (28, 41)

julia> segments = [(2, 7), (2, 12), (12, 17), (2, 22), (2, 27), (2, 32), (32, 37), (2, 42), (42, 47)];

julia> bad_indices = [2, 4, 5, 6, 8]
5-element Vector{Int64}:
 2
 4
 5
 6
 8

julia> DelaunayTriangulation.fix_segments!(segments, bad_indices)
9-element Vector{Tuple{Int64, Int64}}:
 (2, 7)
 (7, 12)
 (12, 17)
 (17, 22)
 (22, 27)
 (27, 32)
 (32, 37)
 (37, 42)
 (42, 47)
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
    return segments
end

"""
    connect_segments!(segments)

Connects the ordered vector of `segments` so that the endpoints all connect, preserving order.

# Example 
```jldoctest
julia> using DelaunayTriangulation

julia> segments = [(7, 12), (12, 17), (17, 22), (32, 37), (37, 42), (42, 47)];

julia> DelaunayTriangulation.connect_segments!(segments)
7-element Vector{Tuple{Int64, Int64}}:
 (7, 12)
 (12, 17)
 (17, 22)
 (22, 32)
 (32, 37)
 (37, 42)
 (42, 47)
```
"""
function connect_segments!(segments::AbstractVector{E}) where {E}
    f = firstindex(segments)
    ℓ = lastindex(segments)
    I = number_type(E)
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
    return segments
end

"""
    extend_segments!(segments, segment)

Given an ordered vector of `segments`, ensures that they also represent the replacement of `segment`. In particular, 
suppose `segments` represents the sequence of edges

        ---(i₁, i₂)---(i₂, i₃)---(⋯, ⋯)---(iₙ₋₁, iₙ)---

and `segment` is `(i₀, iₙ₊₁)`. Then the extended sequence becomes 

        ---(i₀, i₁)---(i₁, i₂)---(i₂, i₃)---(⋯, ⋯)---(iₙ₋₁, iₙ)---(iₙ, iₙ₊₁)---

# Example 
```jldoctest
julia> using DelaunayTriangulation

julia> segments = [(2, 7), (7, 12), (12, 49)];

julia> segment = (1, 68);

julia> DelaunayTriangulation.extend_segments!(segments, segment)
5-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (2, 7)
 (7, 12)
 (12, 49)
 (49, 68)
```
"""
function extend_segments!(segments::AbstractVector{E}, segment) where {E}
    constrained_i, constrained_j = edge_vertices(segment)
    first_i = initial(segments[begin])
    last_j = terminal(segments[end])
    if constrained_i ≠ first_i
        pushfirst!(segments, construct_edge(E, constrained_i, first_i))
    end
    if constrained_j ≠ last_j
        push!(segments, construct_edge(E, last_j, constrained_j))
    end
    return segments
end

"""
    split_segment!(tri::Triangulation, segment, collinear_segments)
    split_segment!(segments, segment, collinear_segments)

Splits `segment` at the segments in `collinear_segments`, which are assumed to be collinear with `segment`.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `segments`: The underlying set of segments. This is `get_interior_segments(tri)` if `tri` is a [`Triangulation`](@ref).
- `segment`: The segment to split.
- `collinear_segments`: The segments that are collinear with `segment`.

# Outputs 
There is no output, as `segments` is updated in-place.

# Example
```jldoctest
julia> using DelaunayTriangulation

julia> segments = Set(((2, 3), (3, 5), (10, 12)))
Set{Tuple{Int64, Int64}} with 3 elements:
  (2, 3)
  (3, 5)
  (10, 12)

julia> collinear_segments = [(2, 10), (11, 15), (2, 3)]
3-element Vector{Tuple{Int64, Int64}}:
 (2, 10)
 (11, 15)
 (2, 3)

julia> segment = (3, 5)
(3, 5)

julia> DelaunayTriangulation.split_segment!(segments, segment, collinear_segments)
Set{Tuple{Int64, Int64}} with 4 elements:
  (2, 10)
  (2, 3)
  (11, 15)
  (10, 12)
```
"""
function split_segment!(segments, segment, collinear_segments)
    if !isempty(collinear_segments)
        delete_unoriented_edge!(segments, segment)
        for e in each_edge(collinear_segments)
            if !contains_unoriented_edge(e, segments)
                add_edge!(segments, e)
            end
        end
    end
    return segments
end
function split_segment!(tri::Triangulation, segment, collinear_segments)
    return split_segment!(get_interior_segments(tri), segment, collinear_segments)
end

@doc """
    convert_boundary_points_to_indices(x, y; existing_points = NTuple{2, Float64}[], check_args=true) -> (Vector, Vector)
    convert_boundary_points_to_indices(xy; existing_points = NTuple{2, Float64}[], check_args=true) -> (Vector, Vector)

Converts a boundary represented by `(x, y)` or `xy`, where the points are combined rather than as separate sets of coordinates, 
into a set of boundary nodes for use in [`triangulate`](@ref).

# Arguments 
- `x`, `y`: The `x` and `y`-coordinates for the boundary points. The individual vectors must match the specification required for boundaries outlined in the documentation.
- `xy`: As above, except the coordinates are combined rather than given as separate vectors. 

# Keyword Arguments
- `existing_points`: The existing points to append the boundary points to. This is useful if you have a pre-existing set of points.
- `check_args`: Whether to check that the arguments match the specification in the documentation. 

# Outputs 
- `boundary_nodes`: The boundary nodes. 
- `points`: The point set, which is the same as `existing_points` but with the boundary points appended to it. 
"""
convert_boundary_points_to_indices
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

function convert_boundary_points_to_indices(xy; existing_points=NTuple{2,Float64}[], check_args=true, adjust=true)
    if is_point2(xy[1])
        x = [getx(xy[i]) for i in eachindex(xy)]
        y = [gety(xy[i]) for i in eachindex(xy)]
    elseif is_point2(xy[1][1])
        x = [[getx(xy[i][j]) for j in eachindex(xy[i])] for i in eachindex(xy)]
        y = [[gety(xy[i][j]) for j in eachindex(xy[i])] for i in eachindex(xy)]
    elseif is_point2(xy[1][1][1])
        x = [[[getx(xy[i][j][k]) for k in eachindex(xy[i][j])] for j in eachindex(xy[i])] for i in eachindex(xy)]
        y = [[[gety(xy[i][j][k]) for k in eachindex(xy[i][j])] for j in eachindex(xy[i])] for i in eachindex(xy)]
    else
        throw(ArgumentError("Invalid boundary coordinates provided. Please provide either a (1) vector of numbers, (2) vector of vector of numbers, or (3) a vector of vector of vector of numbers."))
    end
    return convert_boundary_points_to_indices(x, y; existing_points=existing_points, check_args=check_args, adjust=adjust)
end

"""
    remake_triangulation_with_constraints(tri::Triangulation, segments, boundary_nodes) -> (Dict, Dict, Triangulation)

Remakes the triangulation `tri` so that it contains `segments` and `boundary_nodes` in its fields.

See also [`replace_ghost_vertex_information`](@ref).

# Arguments
- `tri::Triangulation`: The triangulation to remake.
- `segments`: The segments to add to the triangulation.
- `boundary_nodes`: The boundary nodes to add to the triangulation.

# Outputs
- `new_ghost_vertex_map`: The new ghost vertex map. This will not yet be added to the triangulation.
- `new_ghost_vertex_ranges`: The new ghost vertex ranges. This will not yet be added to the triangulation.
- `new_tri::Triangulation`: The new triangulation, now containing `boundary_nodes` in the `boundary_nodes` field and `segments` in the `interior_segments` field.
"""
function remake_triangulation_with_constraints(tri::Triangulation, segments, boundary_nodes)
    points = get_points(tri)
    triangles = get_triangles(tri)
    boundary_nodes = isnothing(boundary_nodes) ? get_boundary_nodes(tri) : boundary_nodes
    interior_segments = isnothing(segments) ? get_interior_segments(tri) : segments
    all_segments = get_all_segments(tri)
    weights = get_weights(tri)
    adjacent = get_adjacent(tri)
    adjacent2vertex = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    boundary_curves = get_boundary_curves(tri)
    I = integer_type(tri)
    E = edge_type(tri)
    boundary_edge_map = construct_boundary_edge_map(boundary_nodes, I, E)
    ghost_vertex_map = get_ghost_vertex_map(tri)
    new_ghost_vertex_map = construct_ghost_vertex_map(boundary_nodes, I) # Delay putting these in until we are done with the triangulation
    ghost_vertex_ranges = get_ghost_vertex_ranges(tri)
    new_ghost_vertex_ranges = construct_ghost_vertex_ranges(boundary_nodes, I) # Delay putting these in until we are done with the triangulation
    ch = get_convex_hull(tri)
    representative_point_list = get_representative_point_list(tri)
    polygon_hierarchy = get_polygon_hierarchy(tri)
    boundary_enricher = get_boundary_enricher(tri)
    cache = get_cache(tri)
    return new_ghost_vertex_map, new_ghost_vertex_ranges, Triangulation(
        points,
        triangles,
        boundary_nodes,
        interior_segments,
        all_segments,
        weights,
        adjacent,
        adjacent2vertex,
        graph,
        boundary_curves,
        boundary_edge_map,
        ghost_vertex_map, # Delay putting these in until we are done with the triangulation
        ghost_vertex_ranges, # Delay putting these in until we are done with the triangulation
        ch,
        representative_point_list,
        polygon_hierarchy,
        boundary_enricher,
        cache)
end

"""
    replace_ghost_vertex_information(tri::Triangulation, ghost_vertex_map, ghost_vertex_ranges) -> Triangulation

Replaces the ghost vertex information in `tri` with `ghost_vertex_map` and `ghost_vertex_ranges`, using the results from 
[`remake_triangulation_with_constraints`](@ref).

# Arguments 
- `tri::Triangulation`: The triangulation to remake.
- `ghost_vertex_map`: The ghost vertex map to add to the triangulation.
- `ghost_vertex_ranges`: The ghost vertex ranges to add to the triangulation.

# Outputs
- `new_tri::Triangulation`: The new triangulation, now containing `ghost_vertex_map` in the `ghost_vertex_map` field and `ghost_vertex_ranges` in the `ghost_vertex_ranges` field.
"""
function replace_ghost_vertex_information(tri::Triangulation, ghost_vertex_map, ghost_vertex_ranges)
    points = get_points(tri)
    triangles = get_triangles(tri)
    boundary_nodes = get_boundary_nodes(tri)
    interior_segments = get_interior_segments(tri)
    all_segments = get_all_segments(tri)
    weights = get_weights(tri)
    adjacent = get_adjacent(tri)
    adjacent2vertex = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    boundary_curves = get_boundary_curves(tri)
    boundary_edge_map = get_boundary_edge_map(tri)
    ch = get_convex_hull(tri)
    representative_point_list = get_representative_point_list(tri)
    polygon_hierarchy = get_polygon_hierarchy(tri)
    boundary_enricher = get_boundary_enricher(tri)
    cache = get_cache(tri)
    return Triangulation(
        points,
        triangles,
        boundary_nodes,
        interior_segments,
        all_segments,
        weights,
        adjacent,
        adjacent2vertex,
        graph,
        boundary_curves,
        boundary_edge_map,
        ghost_vertex_map,
        ghost_vertex_ranges,
        ch,
        representative_point_list,
        polygon_hierarchy,
        boundary_enricher,
        cache)
end

"""
    is_vertex_closer_than_neighbours([predicates::AbstractPredicateKernel=AdaptiveKernel(),] tri::Triangulation, u, v, jᵢ, jᵢ₋₁, jᵢ₊₁) -> Bool
    is_vertex_closer_than_neighbours([predicates::AbstractPredicateKernel=AdaptiveKernel(),] tri::Triangulation, list::ShuffledPolygonLinkedList, u, v, j) -> Bool

Tests if the vertex `jᵢ` is closer to the line `(u, v)` than its neighbours `jᵢ₋₁` and `jᵢ₊₁`, assuming all these 
vertices are to the left of the line.

See also [`point_closest_to_line`](@ref).

# Arguments 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `u`, `v`: The vertices of the line. 
- `jᵢ`, `jᵢ₋₁`, `jᵢ₊₁`: The vertices to compare. 

The second method extracts these latter two vertices using the doubly-linked `list` of vertices.

# Outputs 
- `flag`: Whether `jᵢ` is closer to the line than `jᵢ₋₁` and `jᵢ₊₁`.
"""
function is_vertex_closer_than_neighbours(predicates::AbstractPredicateKernel, tri::Triangulation, u, v, jᵢ, jᵢ₋₁, jᵢ₊₁)
    prev_comp = point_closest_to_line(predicates, tri, u, v, jᵢ, jᵢ₋₁)
    is_further(prev_comp) && return false # If we just do is_closer, then we can run into an infinite loop when all neighbours are equidistant 
    next_comp = point_closest_to_line(predicates, tri, u, v, jᵢ, jᵢ₊₁)
    is_further(next_comp) && return false
    return true
end
function is_vertex_closer_than_neighbours(predicates::AbstractPredicateKernel, tri::Triangulation, list::ShuffledPolygonLinkedList, u, v, j)
    jᵢ, jᵢ₋₁, jᵢ₊₁ = get_triplet(list, j)
    return is_vertex_closer_than_neighbours(predicates, tri, u, v, jᵢ, jᵢ₋₁, jᵢ₊₁)
end
is_vertex_closer_than_neighbours(tri::Triangulation, u, v, jᵢ, jᵢ₋₁, jᵢ₊₁) = is_vertex_closer_than_neighbours(AdaptiveKernel(), tri, u, v, jᵢ, jᵢ₋₁, jᵢ₊₁)
is_vertex_closer_than_neighbours(tri::Triangulation, list::ShuffledPolygonLinkedList, u, v, j) = is_vertex_closer_than_neighbours(AdaptiveKernel(), tri, list, u, v, j)

"""
    select_random_vertex(tri::Triangulation, list::ShuffledPolygonLinkedList, u, v, range, rng) -> Vertex 

Selects a random vertex that is not closer to the line `(u, v)` than both of its neighbours.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `list::ShuffledPolygonLinkedList`: The linked list of polygon vertices.
- `u`, `v`: The vertices of the line.
- `range`: The range of indices of the vertices to select from.
- `rng::Random.AbstractRNG`: The random number generator to use.
- `predicates::AbstractPredicateKernel`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
- `j`: The selected vertex.
"""
function select_random_vertex(tri::Triangulation, list::ShuffledPolygonLinkedList, u, v, range, rng::Random.AbstractRNG, predicates::AbstractPredicateKernel)
    j = rand(rng, range)
    while is_vertex_closer_than_neighbours(predicates, tri, list, u, v, j)
        j = rand(rng, range)
    end
    return j
end

"""
    prepare_vertex_linked_list(V) -> ShuffledPolygonLinkedList

Given a list of polygon vertices `V`, returns the doubly-linked list of polygon vertices.

# Arguments 
- `V`: The list of polygon vertices.

# Outputs 
- `list::ShuffledPolygonLinkedList`: A [`ShuffledPolygonLinkedList`](@ref). In `list`, `prev[begin]`, `prev[end]`, `next[begin]`, and `next[end]` are all `0`
   as are `shuffled_indices[begin]` and `shuffled_indices[end]`. Moreover, `shuffled_indices` will not have been shuffled yet.
"""
function prepare_vertex_linked_list(V::AbstractArray{I}) where {I}
    Base.require_one_based_indexing(V)
    m = length(V)
    next = zeros(I, m)
    prev = zeros(I, m)
    shuffled_indices = zeros(I, m)
    for i = 2:(m-1)
        next[i] = i + 1
        prev[i] = i - 1
        shuffled_indices[i] = i
    end
    return ShuffledPolygonLinkedList(next, prev, shuffled_indices, m, V)
end

"""
    delete_polygon_vertices_in_random_order!(list::ShuffledPolygonLinkedList, tri::Triangulation, u, v, rng::Random.AbstractRNG, predicates::AbstractPredicateKernel)

Deletes vertices from the polygon defined by `list` in a random order.

# Arguments 
- `list::ShuffledPolygonLinkedList`: The linked list of polygon vertices.
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `u`, `v`: The vertices of the segment `(u, v)` that was inserted in order to define the polygon `V = list.S`.
- `rng::Random.AbstractRNG`: The random number generator to use.
- `predicates::AbstractPredicateKernel`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
There is no output, but `list` is updated in-place.
"""
function delete_polygon_vertices_in_random_order!(list::ShuffledPolygonLinkedList, tri::Triangulation, u, v, rng::Random.AbstractRNG, predicates::AbstractPredicateKernel)
    m = list.k
    for i in (m-1):-1:3
        j = select_random_vertex(tri, list, u, v, 2:i, rng, predicates)
        delete_vertex!(list, j)
        swap_permutation!(list, i, j)
    end
    return tri
end

"""
    setup_cavity_cdt(tri::Triangulation, V; rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel()) -> ShuffledPolygonLinkedList

Prepares the linked list required for triangulating a cavity excavated by segment insertion in a constrained triangulation.

See also [`prepare_vertex_linked_list`](@ref) and [`delete_polygon_vertices_in_random_order!`](@ref).

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `V`: The list of polygon vertices, given as a counter-clockwise list of vertices, defining the cavity. 

# Keyword Arguments
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator to use.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
- `list::ShuffledPolygonLinkedList`: The linked list of polygon vertices representing the cavity.
"""
function setup_cavity_cdt(tri::Triangulation, V; rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())
    v = V[begin]
    u = V[end]
    list = prepare_vertex_linked_list(V)
    delete_polygon_vertices_in_random_order!(list, tri, u, v, rng, predicates)
    return list
end

"""
    triangulate_cavity_cdt!(tri::Triangulation, V, marked_vertices; rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())

Triangulates the cavity `V` left behind when deleting triangles intersected in a triangulation by an edge, updating `tri` to do so.

# Arguments
- `tri::Triangulation`: The [`Triangulation`](@ref) to update. This should be an empty triangulation.
- `V`: The list of polygon vertices, given as a counter-clockwise list of vertices, defining the cavity.
- `tri_fan::Triangulation`: The [`Triangulation`](@ref) to use for the fan of triangles to be re-triangulated. This should be an empty triangulation.
- `marked_vertices`: Cache for marking vertices to re-triangulate during the triangulation.
- `fan_triangles`: A cache used for sorting and identifying triangles in a fan for retriangulation.

# Keyword Arguments
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator to use or [`setup_cavity_cdt`](@ref).
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs 
There is no output, but `tri` is updated in-place.
"""
function triangulate_cavity_cdt!(tri::Triangulation, V, tri_fan::Triangulation, marked_vertices, fan_triangles; rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())
    list = setup_cavity_cdt(tri, V; rng, predicates)
    add_triangle!(tri, V[begin], V[list.shuffled_indices[2]], V[end]; protect_boundary=true, update_ghost_edges=false)
    for i in 3:(list.k-1)
        a, b, c = get_triplet(list, i)
        add_point_cavity_cdt!(tri, a, b, c, marked_vertices, predicates)
        if !isempty(marked_vertices) && (last(marked_vertices) == a) || (a ∈ marked_vertices) # We try and insert a last, so the first check will sometimes be a nice boost
            split_marked_vertices!(fan_triangles, tri, marked_vertices)
            num_triangles(fan_triangles) ≤ 1 && continue
            empty!(marked_vertices)
            sort_fan!(marked_vertices, fan_triangles, tri) # reuse the marked_vertices cache for the fan
            retriangulate_fan!(tri, tri_fan, marked_vertices, fan_triangles; rng, predicates)
            empty!(marked_vertices)
            empty!(fan_triangles)
            empty_unconstrained_triangulation!(tri_fan)
        end
    end
    return tri
end

"""
    retriangulate_fan!(tri::Triangulation, tri_fan::Triangulation, fan, fan_triangles; predicates::AbstractPredicateKernel=AdaptiveKernel(), rng::Random.AbstractRNG=Random.default_rng())

Given a sorted set of vertices `fan` in a fan of triangles associated with `fan_triangles`, retriangulates the fan, updating `tri` to do so and 
using `tri_fan` as a temporary triangulation. (This implements Lines 17--19 and Line 28 of the algorithms in [this paper](http://dx.doi.org/10.1016/j.comgeo.2015.04.006).)

The `predicates` argument defines the method for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
"""
function retriangulate_fan!(tri::Triangulation, tri_fan::Triangulation, fan, fan_triangles; predicates::AbstractPredicateKernel=AdaptiveKernel(), rng::Random.AbstractRNG=Random.default_rng())
    for T in each_triangle(fan_triangles)
        u, v, w = triangle_vertices(T)
        delete_triangle!(tri, u, v, w; protect_boundary=true)
    end
    triangulate_convex!(tri_fan, fan; predicates, rng)
    for T in each_solid_triangle(tri_fan)
        u, v, w = triangle_vertices(T)
        add_triangle!(tri, u, v, w; protect_boundary=true, update_ghost_edges=false)
    end
    return tri
end

"""
    split_marked_vertices!(fan_triangles, tri::Triangulation, marked_vertices)

Given a set of `marked_vertices` indicating a crossed triangle (like in Figure 9 of [this paper](http://dx.doi.org/10.1016/j.comgeo.2015.04.006)),
finds all triangles whose three vertices are all in `marked_vertices` and places them into `fan_triangles`. 
"""
function split_marked_vertices!(fan_triangles, tri::Triangulation, marked_vertices)
    for u in marked_vertices, v in marked_vertices
        u == v && continue
        w = get_adjacent(tri, u, v)
        if edge_exists(w) && w ∈ marked_vertices
            i, j, k = sort_triangle(u, v, w)
            add_triangle!(fan_triangles, i, j, k)
        end
    end
    return fan_triangles
end

"""
    sort_fan!(fan, fan_triangles, tri::Triangulation)

Given a set of triangles in a fan, `fan_triangles`, associated with some triangulation `tri`, places all the triangle vertices 
and sorts them counter-clockwise, placing the results into `fan`.
"""
function sort_fan!(fan, fan_triangles, tri::Triangulation)
    for τ in each_triangle(fan_triangles)
        i, j, k = triangle_vertices(τ)
        push!(fan, i, j, k)
    end
    unique!(fan)
    sort_convex_polygon!(fan, tri)
    return fan
end

"""
    add_point_cavity_cdt!(tri::Triangulation, u, v, w, marked_vertices)

Adds a point to the cavity `V` left behind when deleting triangles intersected in a triangulation by an edge, updating `tri` to do so.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref) to update.
- `u`: The vertex to add.
- `v`: The vertex along the polygon that is next to `u`.
- `w`: The vertex along the polygon that is previous to `u`.
- `marked_vertices`: Cache for marking vertices to re-triangulate during the triangulation. This gets mutated.

# Outputs
There is no output, but `tri` is updated in-place, as is `marked_vertices` if necessary.
"""
function add_point_cavity_cdt!(tri::Triangulation, u, v, w, marked_vertices, predicates::AbstractPredicateKernel=AdaptiveKernel())
    (u == v || v == w || u == w) && return tri # For some pathological cases, the found cavities double back on itself, e.g. a right cavity [6, 22, 124, 96, 135, 96, 124, 26] is possible. (Take i = 2567; rng = StableRNG(i); points, edges, mat_edges = get_random_vertices_and_constrained_edges(40, 200, 20, rng); tri = triangulate(points; segments=edges, rng=StableRNG(i)))
    x = get_adjacent(tri, w, v)
    if !edge_exists(x)
        insert_flag = true
    else
        p, q, r, s = get_point(tri, w, v, x, u) # Don't want to deal with boundary handling here 
        incircle_test = point_position_relative_to_circle(predicates, p, q, r, s)
        orient_test = triangle_orientation(predicates, tri, u, v, w)
        insert_flag = !is_inside(incircle_test) && is_positively_oriented(orient_test)
    end
    if insert_flag
        add_triangle!(tri, u, v, w; protect_boundary=true, update_ghost_edges=false)
    else
        delete_triangle!(tri, w, v, x; protect_boundary=true, update_ghost_edges=false)
        add_point_cavity_cdt!(tri, u, v, x, marked_vertices, predicates)
        add_point_cavity_cdt!(tri, u, x, w, marked_vertices, predicates)
        if !is_inside(incircle_test)
            push!(marked_vertices, x, w, v, u)
        end
    end
    return tri
end

"""
    add_new_triangles!(tri_original::Triangulation, tris)

Adds the triangles from `tris` to `tri_original`.
"""
function add_new_triangles!(tri_original::Triangulation, tris)
    for tri in each_triangle(tris)
        add_triangle!(tri_original, tri; protect_boundary=true, update_ghost_edges=false)
    end
    return tri_original
end

"""
    locate_intersecting_triangles(tri::Triangulation, e, rotate=Val(true), rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel()) -> (Vector, Vector, Vector, Vector)

Find all the triangles intersected by an edge `e`.

See also [`find_triangle`](@ref).

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `e`: The edge going through the triangulation.
- `rotate=Val(true)`: Whether to rotate the edge so that the minimum degree vertex of `e` is first.
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator to use.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
- `intersecting_triangles`: The intersected triangles. 
- `collinear_segments`: Segments that are collinear with `e`.
- `left_vertices`: The vertices of the intersected triangles that are left of `e`.
- `right_vertices`: The vertices of the intersected triangles that are right of `e`.
"""
function locate_intersecting_triangles(tri::Triangulation, e, rotate=Val(true), rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())
    V = triangle_type(tri)
    E = edge_type(tri)
    I = integer_type(tri)
    e = is_true(rotate) ? sort_edge_by_degree(tri, e) : e # faster to start at the minimum degree vertex of the edge
    history = PointLocationHistory{V,E,I}()
    add_left_vertex!(history, initial(e))
    add_right_vertex!(history, initial(e))
    find_triangle(tri, get_point(tri, terminal(e)); predicates,
        m=nothing, k=initial(e), store_history=Val(true), history, rng
    )
    add_left_vertex!(history, terminal(e))
    add_right_vertex!(history, terminal(e))
    intersecting_triangles = history.triangles
    collinear_segments = history.collinear_segments
    bad_indices = history.collinear_point_indices
    left_vertices = history.left_vertices
    right_vertices = history.right_vertices
    reverse!(left_vertices) # counter-clockwise
    if !isempty(collinear_segments)
        fix_segments!(collinear_segments, bad_indices)
        connect_segments!(collinear_segments)
        extend_segments!(collinear_segments, e)
    end
    return intersecting_triangles, collinear_segments, left_vertices, right_vertices
end

"""
    delete_intersected_triangles!(tri, triangles)

Deletes the triangles in `triangles` from `tri`.
"""
function delete_intersected_triangles!(tri, triangles) # don't really _need_ this method, but maybe it makes the code a bit clearer?
    for τ in each_triangle(triangles)
        delete_triangle!(tri, τ; protect_boundary=true)
    end
    return tri
end

"""
    process_intersecting_triangles!(tri::Triangulation, e, collinear_segments; predicates::AbstractPredicateKernel=AdaptiveKernel(), rng::Random.AbstractRNG=Random.default_rng()) -> Bool

Given segments in `collinear_segments` that are collinear with an edge `e`, updates `tri` so that this edge `e` is instead 
split so that it is instead represented by `collinear_segments`. These new segments will be placed into the triangulation using 
[`add_segment!`](@ref).

The `predicates::AbstractPredicateKernel` argument defines the method for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

See also [`connect_segments!`](@ref), [`extend_segments!`](@ref), [`split_segment!`](@ref) and [`split_boundary_edge_at_collinear_segments!`](@ref).
"""
function process_collinear_segments!(tri::Triangulation, e, collinear_segments; predicates::AbstractPredicateKernel=AdaptiveKernel(), rng::Random.AbstractRNG=Random.default_rng())
    isempty(collinear_segments) && return false
    all_segments = get_all_segments(tri) # the difference between all_segments and segments is important since we need to be careful about what segments are already in the triangulation 
    delete_edge!(all_segments, e)
    connect_segments!(collinear_segments)
    extend_segments!(collinear_segments, e)
    split_segment!(tri, e, collinear_segments)
    if contains_boundary_edge(tri, e)
        split_boundary_edge_at_collinear_segments!(tri, collinear_segments)
    end
    for η in each_edge(collinear_segments)
        add_segment!(tri, η; predicates, rng)
    end
    return true
end

"""
    merge_segments(tri::Triangulation, ghost_vertex_map) -> Edges 

Creates a set of edges that merge all the boundary nodes in `tri` as well as its interior segments into a single collection. 

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `ghost_vertex_map`: The ghost vertex map to use.

# Outputs 
- `all_segments`: The set of edges that merge all the boundary nodes in `tri` as well as its interior segments into a single collection, with type equal to that of `get_interior_segments(tri)`'s.
"""
function merge_segments(ghost_vertex_map, boundary_nodes, segments::Es) where {Es}
    all_segments = Es()
    E = edge_type(all_segments)
    for segment_index in values(ghost_vertex_map)
        bn_nodes = get_boundary_nodes(boundary_nodes, segment_index)
        nedges = num_boundary_edges(bn_nodes)
        for edge_idx in 1:nedges
            vᵢ = get_boundary_nodes(bn_nodes, edge_idx)
            vᵢ₊₁ = get_boundary_nodes(bn_nodes, edge_idx + 1)
            e = construct_edge(E, vᵢ, vᵢ₊₁)
            add_edge!(all_segments, e)
        end
    end
    for e in each_edge(segments)
        add_edge!(all_segments, e)
    end
    return all_segments
end
function merge_segments(tri::Triangulation, ghost_vertex_map)
    return merge_segments(ghost_vertex_map, get_boundary_nodes(tri), get_interior_segments(tri))
end

"""
    constrained_triangulation!(tri::Triangulation, segments, boundary_nodes, predicates::AbstractPredicateKernel, full_polygon_hierarchy; rng=Random.default_rng(), delete_holes=true) -> Triangulation

Creates a constrained triangulation from `tri` by adding `segments` and `boundary_nodes` to it. This will be in-place, but a new triangulation is returned 
to accommodate the changed types.

# Arguments 
- `tri::Triangulation`: The [`Triangulation`](@ref).
- `segments`: The interior segments to add to the triangulation.
- `boundary_nodes`: The boundary nodes to add to the triangulation.
- `predicates::AbstractPredicateKernel`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.
- `full_polygon_hierarchy`: The [`PolygonHierarchy`](@ref) defining the boundary. This will get copied into the existing polygon hierarchy.

# Keyword Arguments
- `rng=Random.default_rng()`: The random number generator to use.
- `delete_holes=true`: Whether to delete holes in the triangulation. See [`delete_holes!`](@ref).

# Outputs
- `new_tri`: The new triangulation, now containing `segments` in the `interior_segments` field and `boundary_nodes` in the `boundary_nodes` field, and with the updated [`PolygonHierarchy`](@ref). See also [`remake_triangulation_with_constraints`](@ref) and [`replace_ghost_vertex_information`](@ref).
"""
function constrained_triangulation!(tri::Triangulation, segments, boundary_nodes, predicates::AbstractPredicateKernel, full_polygon_hierarchy; rng=Random.default_rng(), delete_holes=true)
    ghost_vertex_map, ghost_vertex_ranges, new_tri = remake_triangulation_with_constraints(tri, segments, boundary_nodes)
    all_segments = merge_segments(new_tri, ghost_vertex_map)
    for e in each_edge(all_segments)
        add_segment!(new_tri, e; predicates, rng)
    end
    new_tri_2 = replace_ghost_vertex_information(new_tri, ghost_vertex_map, ghost_vertex_ranges)
    if !(isnothing(boundary_nodes) || !has_boundary_nodes(boundary_nodes)) && delete_holes
        delete_holes!(new_tri_2)
        add_boundary_information!(new_tri_2)
        add_ghost_triangles!(new_tri_2) # fix the ghost triangles 
        polygon_hierarchy = get_polygon_hierarchy(new_tri_2)
        copyto!(polygon_hierarchy, full_polygon_hierarchy)
    end
    return new_tri_2
end

