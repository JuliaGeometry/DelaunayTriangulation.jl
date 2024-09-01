"""
    SmallAngleComplexMember{I}

Struct for representing a member of a small-angle complex.

# Fields 
- `parent_curve::I`: The index of the parent curve in the boundary curves assoicated with the member. If this is `$∅`, then this is instead a member of a complex around an interior segment.
- `next_edge::I`: The next vertex after the apex in the boundary nodes associated with the member.
"""
struct SmallAngleComplexMember{I}
    parent_curve::I
    next_edge::I
end
Base.:(==)(member₁::SmallAngleComplexMember, member₂::SmallAngleComplexMember) = get_parent_curve(member₁) == get_parent_curve(member₂) && get_next_edge(member₁) == get_next_edge(member₂)
Base.show(io::IO, member::SmallAngleComplexMember) = print(io, "SmallAngleComplexMember with parent curve ", get_parent_curve(member), " and next edge ", get_next_edge(member))
Base.show(io::IO, ::MIME"text/plain", member::SmallAngleComplexMember) = print(io, "SmallAngleComplexMember with parent curve ", get_parent_curve(member), " and next edge ", get_next_edge(member), ".")

"""
    replace_next_edge(member::SmallAngleComplexMember{I}, next_edge) where {I} -> SmallAngleComplexMember{I}

Returns a new `SmallAngleComplexMember` with the same parent curve as `member` but with `next_edge` as the next edge.
"""
function replace_next_edge(member::SmallAngleComplexMember{I}, next_edge) where {I}
    return SmallAngleComplexMember{I}(get_parent_curve(member), next_edge)
end

"""
    SmallAngleComplex{I}

Struct for representing a small-angle complex.

# Fields
- `apex::I`: The apex vertex of the complex.
- `members::Vector{SmallAngleComplexMember{I}}`: The members of the complex.

## Extended help 
A small-angle complex is a set of curves that form a contiguous set of small angles, i.e. the angle between 
each consecutive pair of curves is less than 60°. The apex of the complex is the vertex that is shared by all of the curves.
"""
struct SmallAngleComplex{I}
    apex::I
    members::Vector{SmallAngleComplexMember{I}}
end
Base.:(==)(complex₁::SmallAngleComplex, complex₂::SmallAngleComplex) = get_apex(complex₁) == get_apex(complex₂) && get_members(complex₁) == get_members(complex₂)
Base.show(io::IO, complex::SmallAngleComplex) = print(io, "SmallAngleComplex with apex ", get_apex(complex), " and ", length(get_members(complex)), " members")
Base.show(io::IO, ::MIME"text/plain", complex::SmallAngleComplex) = print(io, "SmallAngleComplex with apex ", get_apex(complex), " and ", length(get_members(complex)), " members.")

"""
    replace_next_edge!(complex::SmallAngleComplex, member_id, next_edge)

Replaces the next edge of the `member_id`th member of `complex` with `next_edge`. 

See also [`replace_next_edge`](@ref).
"""
function replace_next_edge!(complex::SmallAngleComplex, member_id, next_edge)
    members = get_members(complex)
    members[member_id] = replace_next_edge(members[member_id], next_edge)
    return nothing
end

"""
    push!(complex::SmallAngleComplex, member::SmallAngleComplexMember) 

Pushes `member` onto the members of `complex`.
"""
Base.push!(complex::SmallAngleComplex, member::SmallAngleComplexMember) = push!(get_members(complex), member)

"""
    append!(complex::SmallAngleComplex, new_complex::SmallAngleComplex)

Appends the members of `new_complex` onto the members of `complex`.
"""
Base.append!(complex::SmallAngleComplex, new_complex::SmallAngleComplex) = append!(get_members(complex), get_members(new_complex))

"""
    get_parent_curve(member::SmallAngleComplexMember{I}) where {I} -> I 

Returns the parent curve of `member`.
"""
get_parent_curve(member::SmallAngleComplexMember) = member.parent_curve

"""
    get_next_edge(member::SmallAngleComplexMember{I}) where {I} -> I

Returns the next edge of `member`.
"""
get_next_edge(member::SmallAngleComplexMember) = member.next_edge

"""
    get_apex(complex::SmallAngleComplex{I}) where {I} -> I

Returns the apex of `complex`.
"""
get_apex(complex::SmallAngleComplex) = complex.apex

"""
    get_members(complex::SmallAngleComplex{I}) where {I} -> Vector{SmallAngleComplexMember{I}}

Returns the members of `complex`.
"""
get_members(complex::SmallAngleComplex) = complex.members

"""
    get_small_angle_complexes(points, boundary_nodes, boundary_curves, segments=nothing; IntegerType=Int) -> Dict{IntegerType,Vector{SmallAngleComplex{IntegerType}}}

Returns a map from an apex vertex to a list of all curves that define a small angle complex associated with that apex vertex.
"""
function get_small_angle_complexes(points, boundary_nodes, boundary_curves, segments = nothing; IntegerType = Int)
    d = Dict{IntegerType, Vector{SmallAngleComplex{IntegerType}}}()
    if has_multiple_curves(boundary_nodes)
        _get_small_angle_complexes_multiple_curves!(d, boundary_nodes, boundary_curves, IntegerType)
    elseif has_multiple_sections(boundary_nodes)
        _get_small_angle_complexes_multiple_sections!(d, boundary_nodes, boundary_curves, 0, IntegerType)
    else
        _get_small_angle_complexes_contiguous!(d, boundary_nodes, boundary_nodes, boundary_curves, 1, 1, IntegerType)
    end
    !isnothing(segments) && _get_small_angle_complexes_segments!(d, segments, points, IntegerType)
    for (apex, complexes) in d
        complex = complexes[1] # currently, there is only a single yet-to-be-partitioned complex 
        sort_members!(complex, points)
        new_complex = partition_members(complexes, points)
        d[apex] = new_complex
    end
    return d
end
function _get_small_angle_complexes_multiple_curves!(d, boundary_nodes, boundary_curves, ::Type{I}) where {I}
    nc = num_curves(boundary_nodes)
    ctr = 0
    for k in 1:nc
        curve_nodes = get_boundary_nodes(boundary_nodes, k)
        _get_small_angle_complexes_multiple_sections!(d, curve_nodes, boundary_curves, ctr, I)
        ctr += num_sections(curve_nodes)
    end
    return nothing
end
function _get_small_angle_complexes_multiple_sections!(d, boundary_nodes, boundary_curves, init_index₁, ::Type{I}) where {I}
    ns = num_sections(boundary_nodes)
    first_section = get_boundary_nodes(boundary_nodes, 1)
    index₁ = 1 + init_index₁
    for _index₂ in 2:ns
        index₂ = _index₂ + init_index₁
        next_section = get_boundary_nodes(boundary_nodes, _index₂)
        _get_small_angle_complexes_contiguous!(d, first_section, next_section, boundary_curves, index₁, index₂, I)
        index₁ = index₂
        first_section = next_section
    end
    next_section = get_boundary_nodes(boundary_nodes, 1)
    index₂ = 1 + init_index₁
    _get_small_angle_complexes_contiguous!(d, first_section, next_section, boundary_curves, index₁, index₂, I)
    return nothing
end
function _get_small_angle_complexes_contiguous!(d, first_section, next_section, boundary_curves, index₁, index₂, ::Type{I}) where {I}
    θ = angle_between(boundary_curves, index₁, index₂)
    if θ ≤ π / 3
        n = num_boundary_edges(first_section)
        apex = get_boundary_nodes(next_section, 1)
        next₁ = get_boundary_nodes(next_section, 2)
        next₂ = get_boundary_nodes(first_section, n)
        member₁ = SmallAngleComplexMember(I(index₂), I(next₁))
        member₂ = SmallAngleComplexMember(I(index₁), I(next₂))
        complex_vec = get!(Vector{SmallAngleComplex{I}}, d, apex)
        if isempty(complex_vec)
            complex = SmallAngleComplex(I(apex), [member₁, member₂])
            push!(complex_vec, complex)
        else
            push!(complex_vec[1], member₁, member₂)
        end
    end
    return nothing
end
function _get_small_angle_complexes_segments!(d, segments, points, ::Type{I}) where {I}
    segment_map = construct_segment_map(segments, points, I)
    for (vertex, vertices) in segment_map
        length(vertices) == 1 && continue
        p = get_point(points, vertex)
        px, py = getxy(p)
        for i in eachindex(vertices)
            j = i == lastindex(vertices) ? firstindex(vertices) : i + 1
            u, v = vertices[i], vertices[j]
            q, r = get_point(points, u, v)
            qx, qy = getxy(q)
            rx, ry = getxy(r)
            b = (qx - px, qy - py)
            a = (rx - px, ry - py)
            θ = angle_between(b, a)
            if θ ≤ π / 3
                member₁ = SmallAngleComplexMember(I(∅), u)
                member₂ = SmallAngleComplexMember(I(∅), v)
                complex_vec = get!(Vector{SmallAngleComplex{I}}, d, vertex)
                if isempty(complex_vec)
                    complex = SmallAngleComplex(I(vertex), [member₁, member₂])
                    push!(complex_vec, complex)
                else
                    push!(complex_vec[1], member₁, member₂)
                end
                unique!(get_members(complex_vec[1])) # typically, sets are small enough that this doesn't matter for a typical user.
            end
        end
    end
    return nothing
end

"""
    construct_segment_map(segments, points, IntegerType) -> Dict{IntegerType, Vector{IntegerType}}

Returns the segment map of `segments`. This is a map that maps a vertex to all vertices that share a segment with that vertex.
Segments are stored twice. The vertices associated with a vertex are sorted counter-clockwise, using the `points` argument 
to define the coordinates.
"""
function construct_segment_map(segments, points, ::Type{I}) where {I}
    segment_map = Dict{I, Vector{I}}()
    for e in each_edge(segments)
        i, j = edge_vertices(e)
        iset = get!(Vector{I}, segment_map, i)
        jset = get!(Vector{I}, segment_map, j)
        push!(iset, j)
        push!(jset, i)
    end
    for (vertex, vertices) in segment_map
        length(vertices) == 1 && continue
        p = get_point(points, vertex)
        first_vertex = first(vertices)
        q = get_point(points, first_vertex)
        px, py = getxy(p)
        qx, qy = getxy(q)
        base = (qx - px, qy - py)
        sort!(
            vertices, by = _vertex -> begin
                _q = get_point(points, _vertex)
                _qx, _qy = getxy(_q)
                next_base = (_qx - px, _qy - py)
                return angle_between(base, next_base)
            end, rev = false,
        )
    end
    return segment_map
end

"""
    sort_members!(complex::SmallAngleComplex, points)

Sorts the members of `complex` in a counter-clockwise order around the apex of `complex`.
"""
function sort_members!(complex::SmallAngleComplex, points)
    members = get_members(complex)
    apex = get_apex(complex)
    first_member = first(members)
    first_edge = get_next_edge(first_member)
    p, q = get_point(points, apex, first_edge)
    px, py = getxy(p)
    qx, qy = getxy(q)
    base = (qx - px, qy - py)
    sort!(
        members, by = member -> begin
            _q = get_point(points, get_next_edge(member))
            _qx, _qy = getxy(_q)
            next_base = (_qx - px, _qy - py)
            return angle_between(base, next_base)
        end, rev = false,
    )
    return complex
end

"""
    partition_members(complexes::Vector{SmallAngleComplex{I}}, points) where {I} -> Vector{SmallAngleComplex{I}}

Partitions the members of each complex in `complexes` into a new set of complexes. The complexes in `complexes` are assumed to be sorted in a counter-clockwise order around the apex of each complex.
The partitioning is done so that the original set of members are now correctly split into their own complexes, since the original complexes might not have formed a properly contiguous set of small angles.
"""
function partition_members(complexes::Vector{SmallAngleComplex{I}}, points) where {I}
    # Setup
    new_complexes = SmallAngleComplex{I}[]
    complex = first(complexes)
    members = get_members(complex)
    sizehint!(new_complexes, length(members))
    apex = get_apex(complex)
    init_complex = SmallAngleComplex{I}(apex, SmallAngleComplexMember{I}[])
    push!(new_complexes, init_complex)
    # Setup the loop 
    member = first(members)
    next_edge = get_next_edge(member)
    p, q = get_point(points, apex, next_edge)
    px, py = getxy(p)
    qx, qy = getxy(q)
    base = (qx - px, qy - py)
    n = length(members)
    push!(init_complex, member)
    # Now partition
    for i in 2:n
        base = _partition_members_itr!(new_complexes, members, apex, points, i, base, px, py)
    end
    # Decide what we need to do between the last and first members 
    _partition_members_itr!(new_complexes, members, apex, points, 1, base, px, py)
    return new_complexes
end
function _partition_members_itr!(new_complexes::Vector{SmallAngleComplex{I}}, members, apex, points, i, base, px, py) where {I}
    member = members[i]
    next_edge = get_next_edge(member)
    q = get_point(points, next_edge)
    qx, qy = getxy(q)
    next_base = (qx - px, qy - py)
    θ = angle_between(base, next_base)
    if θ ≤ π / 3
        current_complex = new_complexes[end]
        if i == 1 && length(new_complexes) > 1
            current_complex = pop!(new_complexes)
            first_complex = first(new_complexes)
            append!(current_complex, first_complex)
            new_complexes[1] = current_complex
        elseif i ≠ 1
            push!(current_complex, member)
        end
    elseif i ≠ 1
        new_complex = SmallAngleComplex{I}(apex, [member])
        push!(new_complexes, new_complex)
    end
    base = next_base
    return base
end

"""
    get_minimum_edge_length(complex::SmallAngleComplex, points) -> Float64

Returns the minimum edge length in `complex` with respect to `points`.
"""
function get_minimum_edge_length(complex::SmallAngleComplex, points)
    apex = get_apex(complex)
    members = get_members(complex)
    p = get_point(points, apex)
    len = Inf
    for member in members
        next_edge = get_next_edge(member)
        q = get_point(points, next_edge)
        len = min(len, dist(getxy(p), getxy(q)))
    end
    return len
end

"""
    BoundaryEnricher{P,B,C,I,BM,S,E}

Struct used for performing boundary enrichment on a curve-bounded boundary. 

See also [`enrich_boundary!`](@ref).

# Fields 
- `points::P`: The point set. 
- `boundary_nodes::B`: The boundary nodes.
- `segments::S`: The segments.
- `boundary_curves::C`: The boundary curves.
- `polygon_hierarchy::PolygonHierarchy{I}`: The polygon hierarchy.
- `parent_map::Dict{NTuple{2,I},I}`: A map from an edge, represented as a `Tuple`, to the index of the parent curve in `boundary_curves`.
- `curve_index_map::Dict{I,I}`: A map from a curve index to the index of the curve in `boundary_curves`.
- `boundary_edge_map::B`: A map from a boundary node to the index of the curve in `boundary_curves` that it belongs to. See [`construct_boundary_edge_map`](@ref).
- `spatial_tree::BoundaryRTree{P}`: The [`BoundaryRTree`](@ref) used for spatial indexing.
- `queue::Queue{I}`: A queue used for processing vertices during enrichment.
- `small_angle_complexes::Dict{I,Vector{SmallAngleComplex{I}}}`: A map from an apex vertex to a list of all curves that define a small angle complex associated with that apex vertex.

The first three fields should be those associated with [`convert_boundary_curves!`](@ref).

# Constructor 

    BoundaryEnricher(points, boundary_nodes; IntegerType=Int, n=4096, coarse_n=0)

This constructor will use [`convert_boundary_curves!`](@ref) to convert `points` and `boundary_nodes` into a set of boundary curves and modified boundary nodes suitable for enrichment. The boundary nodes 
field will no longer aliased with the input `boundary_nodes`, although `points` will be. The polygon hierarchy is computed using [`construct_polygon_hierarchy`](@ref).
The argument `n` is used in [`polygonise`](@ref) for filling out the boundary temporarily in order to construct the [`PolygonHierarchy`](@ref). The argument `coarse_n` defines the initial coarse discretisation 
through [`coarse_discretisation!`](@ref); the default `n=0` means that the coarse discretisation will be performed until the maximum total variation of a subcurve is less than π/2.
"""
struct BoundaryEnricher{P, B, C, I, BM, S, E}
    points::P
    boundary_nodes::B
    segments::S
    boundary_curves::C
    polygon_hierarchy::PolygonHierarchy{I}
    parent_map::Dict{E, I}
    curve_index_map::Dict{I, I}
    boundary_edge_map::BM
    spatial_tree::BoundaryRTree{P}
    queue::Queue{I}
    small_angle_complexes::Dict{I, Vector{SmallAngleComplex{I}}}
end
function BoundaryEnricher(
        points::P, boundary_nodes::B, segments = nothing; n = 4096, coarse_n = 0,
        IntegerType::Type{I} = Int,
        EdgeType::Type{E} = isnothing(segments) ? NTuple{2, IntegerType} : (edge_type ∘ typeof)(segments),
        EdgesType::Type{Es} = isnothing(segments) ? Set{EdgeType} : typeof(segments),
    ) where {P, B, I, E, Es}
    boundary_curves, new_boundary_nodes = convert_boundary_curves!(points, boundary_nodes, I)
    polygon_hierarchy = construct_polygon_hierarchy(points, new_boundary_nodes, boundary_curves; IntegerType, n)
    return _construct_boundary_enricher(points, new_boundary_nodes, boundary_curves, polygon_hierarchy, segments, n, coarse_n, I, E, Es)
end
function _construct_boundary_enricher(points, boundary_nodes, boundary_curves, polygon_hierarchy, segments, n, coarse_n, ::Type{I}, ::Type{E}, ::Type{Es}) where {I, E, Es}
    expand_bounds!(polygon_hierarchy, ε(number_type(points)))
    coarse_discretisation!(points, boundary_nodes, boundary_curves; n = coarse_n)
    boundary_edge_map = construct_boundary_edge_map(boundary_nodes, I)
    parent_map = Dict{E, I}()
    curve_index_map = Dict{I, I}()
    spatial_tree = BoundaryRTree(points)
    queue = Queue{I}()
    small_angle_complexes = get_small_angle_complexes(points, boundary_nodes, boundary_curves, segments; IntegerType = I)
    _segments = isnothing(segments) ? Es() : segments
    enricher = BoundaryEnricher(points, boundary_nodes, _segments, boundary_curves, polygon_hierarchy, parent_map, curve_index_map, boundary_edge_map, spatial_tree, queue, small_angle_complexes)
    construct_parent_map!(enricher)
    construct_curve_index_map!(enricher)
    construct_tree!(enricher)
    return enricher
end

function Base.:(==)(enricher1::BoundaryEnricher, enricher2::BoundaryEnricher)
    get_points(enricher1) ≠ get_points(enricher2) && return false
    get_boundary_nodes(enricher1) ≠ get_boundary_nodes(enricher2) && return false
    get_boundary_curves(enricher1) ≠ get_boundary_curves(enricher2) && return false
    get_polygon_hierarchy(enricher1) ≠ get_polygon_hierarchy(enricher2) && return false
    get_parent_map(enricher1) ≠ get_parent_map(enricher2) && return false
    get_curve_index_map(enricher1) ≠ get_curve_index_map(enricher2) && return false
    get_boundary_edge_map(enricher1) ≠ get_boundary_edge_map(enricher2) && return false
    get_spatial_tree(enricher1) ≠ get_spatial_tree(enricher2) && return false
    get_queue(enricher1) ≠ get_queue(enricher2) && return false
    get_small_angle_complexes(enricher1) ≠ get_small_angle_complexes(enricher2) && return false
    return true
end

function Base.show(io::IO, ::MIME"text/plain", enricher::BoundaryEnricher)
    print(io, "BoundaryEnricher with $(length(get_boundary_curves(enricher))) boundary curves and $(num_points(get_points(enricher))) points")
end
function check_args(enricher::BoundaryEnricher)
    points = get_points(enricher)
    boundary_nodes = get_boundary_nodes(enricher)
    hierarchy = get_polygon_hierarchy(enricher)
    boundary_curves = get_boundary_curves(enricher)
    return check_args(points, boundary_nodes, hierarchy, boundary_curves)
end

"""
    get_points(boundary_enricher::BoundaryEnricher{P}) -> P

Returns the point set associated with `boundary_enricher`.
"""
get_points(boundary_enricher::BoundaryEnricher) = boundary_enricher.points

"""
    get_boundary_nodes(boundary_enricher::BoundaryEnricher{P,B}) -> B

Returns the boundary nodes associated with `boundary_enricher`.
"""
get_boundary_nodes(boundary_enricher::BoundaryEnricher) = boundary_enricher.boundary_nodes

"""
    get_boundary_curves(boundary_enricher::BoundaryEnricher{P,B,C}) -> C

Returns the boundary curves associated with `boundary_enricher`.
"""
get_boundary_curves(boundary_enricher::BoundaryEnricher) = boundary_enricher.boundary_curves

"""
    get_segments(boundary_enricher::BoundaryEnricher{P,B,C,I,BM,S}) -> S

Returns the segments associated with `boundary_enricher`.
"""
get_segments(boundary_enricher::BoundaryEnricher) = boundary_enricher.segments

"""
    has_segments(boundary_enricher::BoundaryEnricher -> Bool

Returns `true` if `boundary_enricher` has interior segments, and `false` otherwise.
"""
function has_segments(boundary_enricher::BoundaryEnricher)
    segments = get_segments(boundary_enricher)
    isnothing(segments) && return false
    return !isempty(segments)
end

"""
    is_segment(enricher::BoundaryEnricher, i, j) -> Bool 

Returns `true` if `(i, j)` or `(j, i)` is an interior segment of `enricher`, and `false` otherwise.
"""
function is_segment(enricher::BoundaryEnricher, i, j)
    !has_segments(enricher) && return false
    segments = get_segments(enricher)
    E = edge_type(segments)
    e = construct_edge(E, i, j)
    return contains_unoriented_edge(e, segments)
end

"""
    get_boundary_curve(boundary_enricher::BoundaryEnricher, curve_index) -> AbstractParametricCurve

Returns the `curve_index`th curve from the boundary curves in `boundary_enricher`.
"""
get_boundary_curve(boundary_enricher::BoundaryEnricher, curve_index) = get_boundary_curves(boundary_enricher)[curve_index]

"""
    get_polygon_hierarchy(boundary_enricher::BoundaryEnricher{P,B,C,I}) -> PolygonHierarchy{I}

Returns the polygon hierarchy associated with `boundary_enricher`.
"""
get_polygon_hierarchy(boundary_enricher::BoundaryEnricher) = boundary_enricher.polygon_hierarchy

"""
    get_parent_map(boundary_enricher::BoundaryEnricher{P,B,C,I}) -> Dict{NTuple{2,I},I}

Returns the parent map associated with `boundary_enricher`.
"""
get_parent_map(boundary_enricher::BoundaryEnricher) = boundary_enricher.parent_map

"""
    get_curve_index_map(boundary_enricher::BoundaryEnricher{P,B,C,I}) -> Dict{I,I}

Returns the curve index map associated with `boundary_enricher`.
"""
get_curve_index_map(boundary_enricher::BoundaryEnricher) = boundary_enricher.curve_index_map

"""
    map_curve_index(boundary_enricher::BoundaryEnricher, curve_index) -> Integer 

Returns the curve index in `boundary_enricher` associated with `curve_index`.
"""
map_curve_index(boundary_enricher::BoundaryEnricher, curve_index) = get_curve_index_map(boundary_enricher)[curve_index]

"""
    get_boundary_edge_map(boundary_enricher::BoundaryEnricher{P,B,C,I,BM}) -> BM 

Returns the boundary edge map associated with `boundary_enricher`.
"""
get_boundary_edge_map(boundary_enricher::BoundaryEnricher) = boundary_enricher.boundary_edge_map

"""
    get_boundary_edge_map(boundary_enricher::BoundaryEnricher, i, j)

Returns the value from the key `(i, j)` in the boundary edge map of `boundary_enricher`. The returned value is a `Tuple` 
`(position, index)` so that `boundary_nodes = get_boundary_nodes(get_boundary_nodes(boundary_enricher), position)` are the boundary nodes associated 
with the section that `(i, j)` resides on, and `i = get_boundary_nodes(boundary_nodes, index)` and 
`j = get_boundary_nodes(boundary_nodes, index + 1)`.
"""
get_boundary_edge_map(boundary_enricher::BoundaryEnricher, i, j) = get_boundary_edge_map(boundary_enricher)[(i, j)]

"""
    get_spatial_tree(boundary_enricher::BoundaryEnricher{P,B,C,I}) -> RTree

Returns the spatial tree associated with `boundary_enricher`.
"""
get_spatial_tree(boundary_enricher::BoundaryEnricher) = boundary_enricher.spatial_tree

"""
    get_queue(boundary_enricher::BoundaryEnricher{P,B,C,I}) -> Queue{I}

Returns the queue associated with `boundary_enricher`.
"""
get_queue(boundary_enricher::BoundaryEnricher) = boundary_enricher.queue

"""
    get_small_angle_complexes(boundary_enricher::BoundaryEnricher{P,B,C,I}) -> Dict{I,Vector{SmallAngleComplex{I}}}

Returns the small angle complexes associated with `boundary_enricher`.
"""
get_small_angle_complexes(boundary_enricher::BoundaryEnricher) = boundary_enricher.small_angle_complexes

"""
    is_small_angle_complex_apex(boundary_enricher::BoundaryEnricher, apex) -> Bool

Returns `true` if `apex` is the apex of a small angle complex in `boundary_enricher`, and `false` otherwise.
"""
is_small_angle_complex_apex(boundary_enricher::BoundaryEnricher, apex) = haskey(get_small_angle_complexes(boundary_enricher), apex)

"""
    get_small_angle_complex(boundary_enricher::BoundaryEnricher, apex) -> Vector{SmallAngleComplex}

Returns the small angle complexes in `boundary_enricher` associated with `apex`.
"""
function get_small_angle_complexes(boundary_enricher::BoundaryEnricher, apex)
    complexes = get_small_angle_complexes(boundary_enricher)
    return complexes[apex]
end

"""
    get_parent(boundary_enricher::BoundaryEnricher{P,B,C,I}, i::I, j::I) -> I

Returns the parent of the edge `(i, j)` in `boundary_enricher`. If the edge is not in the parent map, then `$∅` is returned.
"""
get_parent(boundary_enricher::BoundaryEnricher{P, B, C, I}, i, j) where {P, B, C, I} = get(get_parent_map(boundary_enricher), (i, j), I(∅))

"""
    set_parent!(boundary_enricher::BoundaryEnricher, i, j, k)

Sets the parent of the edge `(i, j)` in `boundary_enricher` to `k`.
"""
function set_parent!(boundary_enricher::BoundaryEnricher, i, j, k)
    get_parent_map(boundary_enricher)[(i, j)] = k
    return nothing
end

"""
    delete_edge!(boundary_enricher::BoundaryEnricher, i, j)

Deletes the edge `(i, j)` in `boundary_enricher`.
"""
function delete_edge!(boundary_enricher::BoundaryEnricher, i, j)
    delete!(get_parent_map(boundary_enricher), (i, j))
    return nothing
end

"""
    update_parent_map!(boundary_enricher::BoundaryEnricher, i, j, k)

Replaces the edge `(i, j)` in `boundary_enricher` with the edges `(i, k)` and `(k, j)` in the parent map.
"""
function update_parent_map!(boundary_enricher::BoundaryEnricher, i, j, k)
    parent = get_parent(boundary_enricher, i, j)
    delete_edge!(boundary_enricher, i, j)
    set_parent!(boundary_enricher, i, k, parent)
    set_parent!(boundary_enricher, k, j, parent)
    return nothing
end

"""
    each_boundary_edge(enricher::BoundaryEnricher) -> KeySet 

Returns the set of keys in the parent map of `enricher`, i.e. each boundary edge in `enricher`.
"""
function each_boundary_edge(enricher::BoundaryEnricher)
    return keys(get_parent_map(enricher))
end

"""
    construct_parent_map!(enricher::BoundaryEnricher)

Constructs the parent map for `enricher`, modifying the parent map field in-place.
"""
function construct_parent_map!(enricher::BoundaryEnricher)
    parent_map = get_parent_map(enricher)
    empty!(parent_map)
    boundary_nodes = get_boundary_nodes(enricher)
    if has_multiple_curves(boundary_nodes)
        _construct_parent_map_multiple_curves!(enricher)
    elseif has_multiple_sections(boundary_nodes)
        _construct_parent_map_multiple_sections!(enricher)
    else
        _construct_parent_map_contiguous!(enricher)
    end
    return enricher
end
function _construct_parent_map_multiple_curves!(enricher::BoundaryEnricher)
    ctr = 1
    boundary_nodes = get_boundary_nodes(enricher)
    for curve_index in 1:num_curves(boundary_nodes)
        curve_nodes = get_boundary_nodes(boundary_nodes, curve_index)
        for section_index in 1:num_sections(curve_nodes)
            section_nodes = get_boundary_nodes(curve_nodes, section_index)
            _construct_parent_map_contiguous!(enricher, section_nodes, ctr)
            ctr += 1
        end
    end
    return enricher
end
function _construct_parent_map_multiple_sections!(enricher::BoundaryEnricher, boundary_nodes = get_boundary_nodes(enricher), ctr = 1)
    ns = num_sections(boundary_nodes)
    for i in 1:ns
        section_nodes = get_boundary_nodes(boundary_nodes, i)
        _construct_parent_map_contiguous!(enricher, section_nodes, ctr)
        ctr += 1
    end
    return enricher
end
function _construct_parent_map_contiguous!(enricher::BoundaryEnricher, boundary_nodes = get_boundary_nodes(enricher), ctr = 1)
    n = num_boundary_edges(boundary_nodes)
    for i in 1:n
        u = get_boundary_nodes(boundary_nodes, i)
        v = get_boundary_nodes(boundary_nodes, i + 1)
        set_parent!(enricher, u, v, ctr)
    end
    return enricher
end

"""
    construct_curve_index_map!(enricher::BoundaryEnricher)

Constructs the curve index map for `enricher`, modifying the curve index map field in-place. 
"""
function construct_curve_index_map!(enricher::BoundaryEnricher)
    boundary_nodes = get_boundary_nodes(enricher)
    if has_multiple_curves(boundary_nodes)
        _construct_curve_index_map_multiple_curves!(enricher)
    elseif has_multiple_sections(boundary_nodes)
        _construct_curve_index_map_multiple_sections!(enricher)
    else
        _construct_curve_index_map_contiguous!(enricher)
    end
    return enricher
end
function _construct_curve_index_map_multiple_curves!(enricher::BoundaryEnricher)
    boundary_nodes = get_boundary_nodes(enricher)
    curve_index_map = get_curve_index_map(enricher)
    ctr = 1
    for j in 1:num_curves(boundary_nodes)
        curve_nodes = get_boundary_nodes(boundary_nodes, j)
        for _ in 1:num_sections(curve_nodes)
            curve_index_map[ctr] = j
            ctr += 1
        end
    end
    return enricher
end
function _construct_curve_index_map_multiple_sections!(enricher::BoundaryEnricher)
    boundary_nodes = get_boundary_nodes(enricher)
    curve_index_map = get_curve_index_map(enricher)
    for i in 1:num_sections(boundary_nodes)
        curve_index_map[i] = 1
    end
    return enricher
end
function _construct_curve_index_map_contiguous!(enricher::BoundaryEnricher)
    curve_index_map = get_curve_index_map(enricher)
    curve_index_map[1] = 1
    return enricher
end

"""
    is_piecewise_linear(enricher::BoundaryEnricher, curve_index) -> Bool

Returns `true` if the `curve_index`th curve in `enricher` is piecewise linear, and `false` otherwise.
"""
@inline function is_piecewise_linear(enricher::BoundaryEnricher, curve_index)
    boundary_curves = get_boundary_curves(enricher)
    return is_piecewise_linear(boundary_curves, curve_index)
end
@inline function is_piecewise_linear(boundary_curves::C, curve_index) where {C <: Tuple}
    isempty(boundary_curves) && return true
    return eval_fnc_at_het_tuple_element(is_piecewise_linear, boundary_curves, curve_index)
end

"""
    get_inverse(enricher::BoundaryEnricher, curve_index, q) -> Float64

Returns the inverse of the `curve_index`th curve at `q`.
"""
@inline function get_inverse(enricher::BoundaryEnricher, curve_index, q)
    boundary_curves = get_boundary_curves(enricher)
    return get_inverse(boundary_curves, curve_index, q)
end
@inline function get_inverse(boundary_curves::C, curve_index, q) where {C <: Tuple}
    return eval_fnc_at_het_tuple_element_with_arg(get_inverse, boundary_curves, (q,), curve_index)
end

"""
    get_equivariation_split(enricher::BoundaryEnricher, curve_index, t₁, t₂) -> Float64, Float64

Returns the equivariation split of the `curve_index`th curve between `t₁` and `t₂`. Also returns the total variation of the two pieces.
"""
@inline function get_equivariation_split(enricher::BoundaryEnricher, curve_index, t₁, t₂)
    boundary_curves = get_boundary_curves(enricher)
    return get_equivariation_split(boundary_curves, curve_index, t₁, t₂)
end
@inline function get_equivariation_split(boundary_curves::C, curve_index, t₁, t₂) where {C <: Tuple}
    return eval_fnc_at_het_tuple_element_with_arg(get_equivariation_split, boundary_curves, (t₁, t₂), curve_index)
end

"""
    get_equidistant_split(enricher::BoundaryEnricher, curve_index, t₁, t₂) -> Float64 

Returns the equidistant split of the `curve_index`th curve between `t₁` and `t₂`.
"""
@inline function get_equidistant_split(enricher::BoundaryEnricher, curve_index, t₁, t₂)
    boundary_curves = get_boundary_curves(enricher)
    return get_equidistant_split(boundary_curves, curve_index, t₁, t₂)
end
@inline function get_equidistant_split(boundary_curves::C, curve_index, t₁, t₂) where {C <: Tuple}
    return eval_fnc_at_het_tuple_element_with_arg(get_equidistant_split, boundary_curves, (t₁, t₂), curve_index)
end


"""
    eval_boundary_curve(enricher::BoundaryEnricher, curve_index, t) -> NTuple{2,Float64}

Returns the `curve_index`th boundary curve at `t`.
"""
@inline function eval_boundary_curve(enricher::BoundaryEnricher, curve_index, t)
    boundary_curves = get_boundary_curves(enricher)
    return eval_boundary_curve(boundary_curves, curve_index, t)
end
@inline function eval_boundary_curve(boundary_curves::C, curve_index, t) where {C <: Tuple}
    return eval_fnc_in_het_tuple(boundary_curves, t, curve_index)
end

"""
    point_position_relative_to_curve([kernel::AbstractPredicateKernel=AdaptiveKernel(),] enricher::BoundaryEnricher, curve_index, p) -> Certificate

Returns a [`Certificate`](@ref) which is 

- `Left`: If `p` is to the left of the `curve_index`th curve.
- `Right`: If `p` is to the right of the `curve_index`th curve.
- `On`: If `p` is on the `curve_index`th curve.

The `kernel` argument determines how this result is computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.
"""
@inline function point_position_relative_to_curve(kernel::AbstractPredicateKernel, enricher::BoundaryEnricher, curve_index, p)
    boundary_curves = get_boundary_curves(enricher)
    return point_position_relative_to_curve(kernel, boundary_curves, curve_index, p)
end
@inline point_position_relative_to_curve(enricher::BoundaryEnricher, curve_index, p) = point_position_relative_to_curve(AdaptiveKernel(), enricher, curve_index, p)
@inline function point_position_relative_to_curve(kernel::AbstractPredicateKernel, boundary_curves::C, curve_index, p) where {C <: Tuple}
    return eval_fnc_at_het_tuple_element_with_arg_and_prearg(point_position_relative_to_curve, boundary_curves, kernel, (p,), curve_index)
end
@inline point_position_relative_to_curve(boundary_curves::C, curve_index, p) where {C <: Tuple} = point_position_relative_to_curve(AdaptiveKernel(), boundary_curves, curve_index, p)

"""
    angle_between(enricher::BoundaryEnricher, curve_index1, curve_index2) -> Float64 

Evaluates [`angle_between`](@ref) on the curves with indices `curve_index1` and `curve_index2` in `enricher`.
"""
function angle_between(enricher::BoundaryEnricher, curve_index1, curve_index2)
    boundary_curves = get_boundary_curves(enricher)
    return angle_between(boundary_curves, curve_index1, curve_index2)
end
function angle_between(boundary_curves::Tuple, curve_index1, curve_index2)
    return eval_fnc_at_het_tuple_two_elements(angle_between, boundary_curves, curve_index1, curve_index2)
end

"""
    get_circle_intersection(enricher::BoundaryEnricher, curve_index, t₁, t₂, r) -> (Float64, NTuple{2,Float64})

Finds the intersection of the `curve_index`th curve with the circle centered at the curve evaluated at `t₁` with radius `r`. The argument 
`t₂` defines the end of the subcurve to consider. The returned tuple is `(t, p)` where `t` is the parameter value of the intersection and `p` is the point of intersection.
"""
function get_circle_intersection(enricher::BoundaryEnricher, curve_index, t₁, t₂, r)
    boundary_curves = get_boundary_curves(enricher)
    return get_circle_intersection(boundary_curves, curve_index, t₁, t₂, r)
end
function get_circle_intersection(boundary_curves::Tuple, curve_index, t₁, t₂, r)
    return eval_fnc_at_het_tuple_element_with_arg(get_circle_intersection, boundary_curves, (t₁, t₂, r), curve_index)
end

"""
    polygonise(points, boundary_nodes, boundary_curves; n=4096)

Fills out a set of points for a curve-bounded domain for use with [`PolygonHierarchy`](@ref).

!!! warning 

    If the boundary curves are complicated so that they take a lot of points in order to be accurately resolved, then you should increase 
    `n`.

# Arguments 
- `points`: The point set.
- `boundary_nodes`: The boundary nodes.
- `boundary_curves`: The boundary curves.

# Keyword Arguments 
- `n=4096`: The number of points to use for filling in each boundary curves.

# Output
- `new_points`: The points defining the filled out boundaries.
- `new_boundary_nodes`: The boundary nodes associated with `new_points`.

!!! warning "Aliasing"

    If the boundary is not curve bounded, then `new_points` and `new_boundary_nodes` remain aliased 
    with the input `points` and `boundary_nodes`.
"""
function polygonise(points, boundary_nodes, boundary_curves; n = 4096)
    new_points = deepcopy(points)
    new_boundary_nodes = deepcopy(boundary_nodes)
    coarse_discretisation!(new_points, new_boundary_nodes, boundary_curves; n)
    return new_points, new_boundary_nodes
end

"""
    construct_tree!(enricher::BoundaryEnricher)

Constructs the spatial tree for `enricher`, modifying the spatial tree field in-place. The parent map 
must be correctly configured in order for this to be valid.
"""
function construct_tree!(enricher::BoundaryEnricher)
    tree = get_spatial_tree(enricher)
    parent_map = get_parent_map(enricher)
    for (i, j) in keys(parent_map)
        insert!(tree, i, j)
    end
    segments = get_segments(enricher)
    if !isnothing(segments)
        for e in segments
            i, j = edge_vertices(e)
            insert!(tree, i, j)
        end
    end
    return enricher
end

"""
    reorient_edge(enricher::BoundaryEnricher, i, j) -> NTuple{2,Integer}

Given an edge `(i, j)`, reorients it so that it is correctly oriented with the boundary. If `(i, j)`
is instead an interior segment rather than a boundary edge, then `(i, j)` is returned.
"""
function reorient_edge(enricher::BoundaryEnricher, i, j)
    boundary_edge_map = get_boundary_edge_map(enricher)
    if haskey(boundary_edge_map, (i, j))
        return (i, j)
    else
        return (j, i)
    end
end

"""
    split_boundary_edge!(enricher::BoundaryEnricher, i, j, r, update_boundary_nodes = Val(true))

Updates the fields of `enricher` after splitting a boundary edge `(i, j)` at the `r`th vertex. The `update_boundary_nodes` argument
can be used to avoid inserting an additional boundary node when `boundary_nodes` was already updated somewhere else (e.g., we need this for 
mesh refinement which already updates the `boundary_nodes` which is aliased with the same field in the enricher).
"""
function split_boundary_edge!(enricher::BoundaryEnricher, i, j, r, update_boundary_nodes = Val(true))
    boundary_nodes = get_boundary_nodes(enricher)
    boundary_edge_map = get_boundary_edge_map(enricher)
    spatial_tree = get_spatial_tree(enricher)
    pos = get_boundary_edge_map(enricher, i, j)
    new_pos = (pos[1], pos[2] + 1)
    is_true(update_boundary_nodes) && insert_boundary_node!(boundary_nodes, new_pos, r)
    split_boundary_edge_map!(boundary_edge_map, boundary_nodes, pos, i, j)
    split_edge!(spatial_tree, i, j, r)
    update_parent_map!(enricher, i, j, r)
    return enricher
end

"""
    split_interior_segment!(enricher::BoundaryEnricher, i, j, r, update_segments = Val(true))

Updates the fields of `enricher` after splitting an interior segment `(i, j)` at the `r`th vertex. 
The `update_segments` argument can be used to avoid inserting an additional segment when `segments` was already updated somewhere else (e.g., we need this for
mesh refinement which already updates the `interior_segments` which is aliased with the `segments` field in the enricher).
"""
function split_interior_segment!(enricher::BoundaryEnricher, i, j, r, update_segments = Val(true))
    segments = get_segments(enricher)
    spatial_tree = get_spatial_tree(enricher)
    E = edge_type(segments)
    e = construct_edge(E, i, j)
    e1 = construct_edge(E, i, r)
    e2 = construct_edge(E, r, j)
    if is_true(update_segments)
        delete_unoriented_edge!(segments, e)
        add_edge!(segments, e1, e2)
    end
    split_edge!(spatial_tree, i, j, r)
    return enricher
end

"""
    split_edge!(enricher::BoundaryEnricher, i, j, r, update_boundary_nodes = Val(true), update_segments = Val(true), is_interior = is_segment(enricher, i, j))

Updates the fields of `enricher` after splitting an edge `(i, j)` at the `r`th vertex. The `update_boundary_nodes` argument 
can be used to avoid inserting an additional boundary node when `boundary_nodes` was already updated somewhere else (e.g., we need this for
mesh refinement which already updates the `boundary_nodes` which is aliased with the same field in the enricher). The same 
point goes for `update_segments` which can be used to avoid inserting an additional segment when `segments` was already updated somewhere else.
The `is_interior` argument can be used to specify whether the edge is an interior segment or a boundary edge. 

See also [`split_boundary_edge!`](@ref) and [`split_interior_segment!`](@ref).
"""
function split_edge!(enricher::BoundaryEnricher, i, j, r, update_boundary_nodes = Val(true), update_segments = Val(true), is_interior = is_segment(enricher, i, j))
    if is_interior
        split_interior_segment!(enricher, i, j, r, update_segments)
    else
        split_boundary_edge!(enricher, i, j, r, update_boundary_nodes)
    end
    return enricher
end

"""
    is_small_angle_complex_member(enricher::BoundaryEnricher, i, j) -> Bool, I, IntegerType, IntegerType

Returns `true` if the edge `(i, j)` is a member of a small angle complex in `enricher`, and `false` otherwise. 

# Outputs 
- `flag`: `true` if the edge is a member of a small angle complex, and `false` otherwise.
- `apex`: If the edge is a member of a small angle complex, then `apex` is the apex of the complex. Otherwise, `apex` is `0`.
- `complex_id`: If the edge is a member of a small angle complex, then `complex_id` is the index of the complex in the list of complexes associated with `apex`. Otherwise, `complex_id` is `0`.
- `member_id`: If the edge is a member of a small angle complex, then `member_id` is the index of the member in the list of members associated with `complex_id`. Otherwise, `member_id` is `0`.
"""
function is_small_angle_complex_member(enricher::BoundaryEnricher, i, j)
    if !is_small_angle_complex_apex(enricher, i) && !is_small_angle_complex_apex(enricher, j)
        return false, oftype(i, 0), 0, 0
    end
    for r in (i, j)
        r′ = r == i ? j : i
        !is_small_angle_complex_apex(enricher, r) && continue
        r_complexes = get_small_angle_complexes(enricher, r)
        for (complex_id, complex) in enumerate(r_complexes)
            for (member_id, member) in enumerate(get_members(complex))
                if get_next_edge(member) == r′
                    return true, r, complex_id, member_id
                end
            end
        end
    end
    return false, oftype(i, 0), 0, 0
end

"""
    replace_next_edge!(enricher::BoundaryEnricher, apex, complex_id, member_id, next_edge)

Replaces the next edge of the `member_id`th member of the `complex_id`th complex associated with `apex` with `next_edge`.
"""
function replace_next_edge!(enricher::BoundaryEnricher, apex, complex_id, member_id, next_edge)
    complexes = get_small_angle_complexes(enricher, apex)
    complex = complexes[complex_id]
    replace_next_edge!(complex, member_id, next_edge)
    return enricher
end
