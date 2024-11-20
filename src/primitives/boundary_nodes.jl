const AV = AbstractVector

# Interface
@inline has_multiple_curves(::AAA) where {F<:Number,A<:AV{F},AA<:AV{A},AAA<:AV{AA}} = true
@inline has_multiple_curves(::AA) where {F<:Number,A<:AV{F},AA<:AV{A}} = false
@inline has_multiple_curves(::A) where {F<:Number,A<:AV{F}} = false
@inline function has_multiple_curves(boundary_nodes::AV)
    #=
    Need this since, e.g.
    [
    [
        [1, 2, 3, 4, 5, 6, 7, 1]
    ],
    [
        [CircularArc((0.5, 0.0), (0.5, 0.0), (0.0, 0.0), positive=false)]
    ],
    ]
    has type Vector{Vector}, but it has multiple curves.
    =#
    flag = !isempty(boundary_nodes) && all(x -> typeof(x) <: AV && eltype(x) <: AV, boundary_nodes)
    return flag
end
@inline has_multiple_curves(::NTuple{N,<:Integer}) where {N} = false

@inline has_multiple_sections(::AAA) where {F<:Number,A<:AV{F},AA<:AV{A},AAA<:AV{AA}} = true
@inline has_multiple_sections(::AA) where {F<:Number,A<:AV{F},AA<:AV{A}} = true
@inline has_multiple_sections(::A) where {F<:Number,A<:AV{F}} = false
@inline function has_multiple_sections(boundary_nodes::AV)
    #=
    Need this since, e.g.
        2-element Vector{Any}:
        [1, 2, 3]
        EllipticalArc[EllipticalArc((0.0, 0.0), 2.0, 0.5, (0.0, 1.0), 0.0, 3.141592653589793, (2.0, 0.0), (-2.0, 0.0))]
    =#
    flag = !isempty(boundary_nodes) && all(x -> typeof(x) <: AV, boundary_nodes)
    return flag
end
@inline has_multiple_sections(::NTuple{N,<:Integer}) where {N} = false

@inline num_curves(boundary_nodes) = has_multiple_curves(boundary_nodes) ? length(boundary_nodes) : 1

@inline num_sections(boundary_nodes) = has_multiple_sections(boundary_nodes) ? length(boundary_nodes) : 1

@inline num_boundary_edges(boundary_nodes::AV) = max(0, length(boundary_nodes) - 1)
@inline num_boundary_edges(::NTuple{N,<:Integer}) where {N} = max(0, N - 1)

@inline get_boundary_nodes(boundary_nodes, m::Integer) = boundary_nodes[m]
@inline get_boundary_nodes(boundary_nodes, m::Integer, n::Integer) = get_boundary_nodes(get_boundary_nodes(boundary_nodes, m), n)
@inline get_boundary_nodes(boundary_nodes, (m, n)::NTuple{2,<:Integer}) = get_boundary_nodes(boundary_nodes, m, n) # for indexing from a boundary map 
@inline get_boundary_nodes(boundary_nodes::A, ::A) where {A} = boundary_nodes # for indexing from a boundary map
@inline get_boundary_nodes(boundary_nodes::A, ::A) where {A<:Tuple{<:Integer,<:Integer}} = boundary_nodes # ambiguity
@inline get_boundary_nodes(boundary_nodes::A, ::A) where {A<:Integer} = boundary_nodes # ambiguity

@inline each_boundary_node(boundary_nodes) = boundary_nodes

@inline function insert_boundary_node!(boundary_nodes, pos, node)
    nodes = get_boundary_nodes(boundary_nodes, pos[1])
    insert!(nodes, pos[2], node)
    return boundary_nodes
end

@inline function delete_boundary_node!(boundary_nodes, pos)
    nodes = get_boundary_nodes(boundary_nodes, pos[1])
    deleteat!(nodes, pos[2])
    return boundary_nodes
end

@inline function get_skeleton(boundary_nodes, ::Type{I}) where {I}
    if has_multiple_curves(boundary_nodes)
        return _get_skeleton_multiple_curves(boundary_nodes, I)
    elseif has_multiple_sections(boundary_nodes)
        return _get_skeleton_multiple_sections(boundary_nodes, I)
    else
        return _get_skeleton_contiguous(boundary_nodes, I)
    end
end
@inline function _get_skeleton_multiple_curves(boundary_nodes, ::Type{I}) where {I}
    boundary_nodes′ = Vector{Vector{Vector{I}}}(undef, num_curves(boundary_nodes))
    for curve_index in each_curve_index(boundary_nodes)
        boundary_nodes_curve = get_boundary_nodes(boundary_nodes, curve_index)
        boundary_nodes′[curve_index] = _get_skeleton_multiple_sections(boundary_nodes_curve, I)
    end
    return boundary_nodes′
end
@inline function _get_skeleton_multiple_sections(boundary_nodes, ::Type{I}) where {I}
    boundary_nodes′ = Vector{Vector{I}}(undef, num_sections(boundary_nodes))
    for section_index in each_section_index(boundary_nodes)
        boundary_nodes_section = get_boundary_nodes(boundary_nodes, section_index)
        boundary_nodes′[section_index] = _get_skeleton_contiguous(boundary_nodes_section, I)
    end
    return boundary_nodes′
end
@inline function _get_skeleton_contiguous(boundary_nodes, ::Type{I}) where {I}
    return I[]
end

# Derived
@inline each_curve_index(boundary_nodes) = 1:num_curves(boundary_nodes)

@inline each_section_index(boundary_nodes) = 1:num_sections(boundary_nodes)

@inline each_edge_index(boundary_nodes) = 1:num_boundary_edges(boundary_nodes)

@inline get_boundary_edge(boundary_nodes, u) = (get_boundary_nodes(boundary_nodes, u), get_boundary_nodes(boundary_nodes, u + 1))
