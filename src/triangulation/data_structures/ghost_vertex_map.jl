struct GhostVertexMap{T}
    map::Vector{T}
end
@inline get_map(g::GhostVertexMap) = g.map
@inline Base.getindex(g::GhostVertexMap, i::Integer) = get_map(g)[𝒢+i+1]
@inline Base.copy(g::GhostVertexMap) = GhostVertexMap(copy(get_map(g)))

@inline function GhostVertexMap(boundary_nodes, IntegerType::Type{I}=number_type(boundary_nodes)) where {I}
    map = if has_multiple_curves(boundary_nodes)
        _gvm_curves(boundary_nodes, IntegerType)
    elseif has_multiple_sections(boundary_nodes)
        _gvm_sections(boundary_nodes, IntegerType)
    else
        _gvm_contiguous(boundary_nodes)
    end
    return GhostVertexMap(map)
end
function _gvm_curves(boundary_nodes, ::Type{I}) where {I}
    map = NTuple{2,I}[]
    for m in each_curve_index(boundary_nodes)
        bn_m = get_boundary_nodes(boundary_nodes, m)
        for n in each_section_index(bn_m)
            push!(map, (m, n))
        end
    end
    return map
end
function _gvm_sections(boundary_nodes, ::Type{I}) where {I}
    map = I[]
    for m in each_section_index(boundary_nodes)
        push!(map, m)
    end
    return map
end
_gvm_contiguous(boundary_nodes) = [boundary_nodes]

@inline get_curve_index(ghost_vertex::Tuple) = ghost_vertex[1]
@inline get_curve_index(ghost_vertex) = 1
@inline get_curve_index(gvm::GhostVertexMap, ghost_vertex) = get_curve_index(gvm[ghost_vertex])

@inline get_section_index(ghost_vertex::Tuple) = ghost_vertex[2]
@inline get_section_index(ghost_vertex::Integer) = ghost_vertex
@inline get_section_index(ghost_vertex) = 1
@inline get_section_index(gvm::GhostVertexMap, ghost_vertex) = get_section_index(gvm[ghost_vertex])
