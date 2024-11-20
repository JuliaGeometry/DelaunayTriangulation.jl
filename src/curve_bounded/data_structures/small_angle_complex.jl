struct SmallAngleComplexMember{I}
    parent_curve::I
    next_edge::I
end
Base.:(==)(member₁::SmallAngleComplexMember, member₂::SmallAngleComplexMember) = get_parent_curve(member₁) == get_parent_curve(member₂) && get_next_edge(member₁) == get_next_edge(member₂)

Base.copy(member::SmallAngleComplexMember) = member

function replace_next_edge(member::SmallAngleComplexMember{I}, next_edge) where {I}
    return SmallAngleComplexMember{I}(get_parent_curve(member), next_edge)
end

struct SmallAngleComplex{I}
    apex::I
    members::Vector{SmallAngleComplexMember{I}}
end
Base.:(==)(complex₁::SmallAngleComplex, complex₂::SmallAngleComplex) = get_apex(complex₁) == get_apex(complex₂) && get_members(complex₁) == get_members(complex₂)

Base.copy(complex::SmallAngleComplex) = SmallAngleComplex(get_apex(complex), copy(get_members(complex)))

function replace_next_edge!(complex::SmallAngleComplex, member_id, next_edge)
    members = get_members(complex)
    members[member_id] = replace_next_edge(members[member_id], next_edge)
    return nothing
end

Base.push!(complex::SmallAngleComplex, member::SmallAngleComplexMember) = (push!(get_members(complex), member); complex)

Base.append!(complex::SmallAngleComplex, new_complex::SmallAngleComplex) = (append!(get_members(complex), get_members(new_complex)); complex)

get_parent_curve(member::SmallAngleComplexMember) = member.parent_curve

get_next_edge(member::SmallAngleComplexMember) = member.next_edge

get_apex(complex::SmallAngleComplex) = complex.apex

get_members(complex::SmallAngleComplex) = complex.members

function get_small_angle_complexes(points, boundary_nodes, boundary_curves, segments=nothing; IntegerType=Int)
    d = Dict{IntegerType,Vector{SmallAngleComplex{IntegerType}}}()
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

function construct_segment_map(segments, points, ::Type{I}) where {I}
    segment_map = Dict{I,Vector{I}}()
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
            vertices, by=_vertex -> begin
                _q = get_point(points, _vertex)
                _qx, _qy = getxy(_q)
                next_base = (_qx - px, _qy - py)
                return angle_between(base, next_base)
            end, rev=false,
        )
    end
    return segment_map
end

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
        members, by=member -> begin
            _q = get_point(points, get_next_edge(member))
            _qx, _qy = getxy(_q)
            next_base = (_qx - px, _qy - py)
            return angle_between(base, next_base)
        end, rev=false,
    )
    return complex
end

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