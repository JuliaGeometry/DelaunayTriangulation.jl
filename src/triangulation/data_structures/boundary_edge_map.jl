struct BoundaryEdgeMap{E,T} <: AbstractDict{E,T}
    map::Dict{E,T}
end
@inline get_map(b::BoundaryEdgeMap) = b.map
@inline Base.getindex(b::BoundaryEdgeMap, e) = get_map(b)[e]
@inline Base.haskey(b::BoundaryEdgeMap, e) = haskey(get_map(b), e)
@inline Base.iterate(b::BoundaryEdgeMap, state...) = iterate(get_map(b), state...)
@inline Base.length(b::BoundaryEdgeMap) = length(get_map(b))

@inline is_contiguous(b::BoundaryEdgeMap{E,Tuple{NTuple{2,I},I}}) where {I} = false
@inline is_contiguous(b::BoundaryEdgeMap{E,NTuple{2,I}}) where {I} = false
@inline is_contiguous(b::BoundaryEdgeMap{E,Tuple{A,I}}) where {A,I} = true

@inline edge_type(::BoundaryEdgeMap{E,T}) where {E,T} = E
@inline integer_type(::BoundaryEdgeMap{E,Tuple{NTuple{2,I},I}}) where {I} = I
@inline integer_type(::BoundaryEdgeMap{E,NTuple{2,I}}) where {I} = I
@inline integer_type(::BoundaryEdgeMap{E,Tuple{A,I}}) where {E,A,I} = I

@inline Base.copy(b::BoundaryEdgeMap) = _bemcopy(b::BoundaryEdgeMap; boundary_nodes=_get_bem_boundary_nodes(b))
@inline function _bemcopy(b::BoundaryEdgeMap; boundary_nodes=_get_bem_boundary_nodes(b))
    if isnothing(boundary_map)
        return BoundaryEdgeMap(copy(get_map(b)))
    else
        I, E = integer_type(b), edge_type(b)
        return BoundaryEdgeMap(boundary_nodes, I, E)
    end
end

@inline _get_bem_boundary_nodes(::BoundaryEdgeMap{E,Tuple{NTuple{2,I},I}}) where {I} = nothing
@inline _get_bem_boundary_nodes(::BoundaryEdgeMap{E,NTuple{2,I}}) where {I} = nothing
@inline _get_bem_boundary_nodes(b::BoundaryEdgeMap{E,Tuple{A,I}}) where {A,I} = first(values(b))[1]

@inline function BoundaryEdgeMap(boundary_nodes, IntegerType::Type{I}=number_type(boundary_nodes), EdgeType::Type{E}=NTuple{2,IntegerType}) where {I,E}
    map = if has_multiple_curves(boundary_nodes)
        _bem_curves(boundary_nodes, I, E)
    elseif has_multiple_sections(boundary_nodes)
        _bem_sections(boundary_nodes, I, E)
    else
        _bem_contiguous(boundary_nodes, I, E)
    end
    return BoundaryEdgeMap(map)
end
@inline function _bem_curves(boundary_nodes, ::Type{I}, ::Type{E}) where {I,E}
    dict = Dict{E,Tuple{NTuple{2,I},I}}()
    for m in each_curve_index(boundary_nodes)
        bn_m = get_boundary_nodes(boundary_nodes, m)
        for n in each_section_index(bn_m)
            bn_n = get_boundary_nodes(bn_m, n)
            for ℓ in each_edge_index(bn_n)
                u, v = get_boundary_edge(bn_n, ℓ)
                e = construct_edge(E, u, v)
                dict[e] = ((m, n), ℓ)
            end
        end
    end
    return dict
end
@inline function _bem_sections(boundary_nodes, ::Type{I}, ::Type{E}) where {I,E}
    dict = Dict{E,NTuple{2,I}}()
    for n in each_section_index(boundary_nodes)
        bn_n = get_boundary_nodes(boundary_nodes, n)
        for ℓ in each_edge_index(bn_n)
            u, v = get_boundary_edge(bn_n, ℓ)
            e = construct_edge(E, u, v)
            dict[e] = (n, ℓ)
        end
    end
    return dict
end
@inline function _bem_contiguous(boundary_nodes::A, ::Type{I}, ::Type{E}) where {A,I,E}
    dict = Dict{E,Tuple{A,I}}()
    for ℓ in each_edge_index(boundary_nodes)
        u, v = get_boundary_edge(boundary_nodes, ℓ)
        e = construct_edge(E, u, v)
        dict[e] = (boundary_nodes, ℓ)
    end
    return dict
end 