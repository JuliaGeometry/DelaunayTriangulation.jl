struct GhostVertexRanges{I}
    ranges::Vector{UnitRange{I}}
end
@inline get_ranges(g::GhostVertexRanges) = g.ranges
@inline Base.getindex(g::GhostVertexRanges, i::Integer) = g.ranges[𝒢+i+1]
@inline Base.copy(g::GhostVertexRanges) = GhostVertexRanges(copy(get_ranges(g)))

function GhostVertexRanges(boundary_nodes, IntegerType::Type{I}=number_type(boundary_nodes)) where {I}
    gvr = UnitRange{I}[]
    ranges = if has_multiple_curves(boundary_nodes)
        _gvr_curves!(gvr, boundary_nodes)
    elseif has_multiple_sections(boundary_nodes)
        _gvr_sections!(gvr, boundary_nodes)
    else
        _gvr_contiguous!(gvr)
    end
    return GhostVertexRanges(ranges)
end
function _gvr_curves!(gvr::Vector{UnitRange{I}}, boundary_nodes) where {I}
    start = I(𝒢)
    for m in each_curve_index(boundary_nodes)
        bn_m = get_boundary_nodes(boundary_nodes, m)
        n = num_sections(bn_m)
        stop = I(start - n + 1)
        for _ in each_section_index(bn_m)
            push!(gvr, stop:start)
        end
        start -= n
    end
    return gvr
end
function _gvr_sections!(gvr::Vector{UnitRange{I}}, boundary_nodes) where {I}
    start = I(𝒢)
    n = num_sections(boundary_nodes)
    stop = I(start - n + 1)
    range = stop:start
    for _ in each_section_index(boundary_nodes)
        push!(gvr, range)
    end
    return gvr
end
function _gvr_contiguous(gvr::Vector{UnitRange{I}}) where {I}
    push!(gvr, I(𝒢):I(𝒢))
    return gvr
end
