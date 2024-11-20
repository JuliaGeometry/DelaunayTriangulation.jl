struct Adjacent2VertexCandidates{I} <: AbstractDict{I,I}
    candidates::Dict{I,I}
end
@inline get_candidates(adj2v::Adjacent2VertexCandidates) = adj2v.candidates
@inline Base.sizehint!(adj2v::Adjacent2VertexCandidates, n) = sizehint!(get_candidates(adj2v), n)
@inline Base.getindex(adj2v::Adjacent2VertexCandidates, u) = get_candidates(adj2v)[u]
@inline Base.haskey(adj2v::Adjacent2VertexCandidates, u) = haskey(get_candidates(adj2v), u)
@inline Base.iterate(adj2v::Adjacent2VertexCandidates, state...) = iterate(get_candidates(adj2v), state...)
@inline Base.length(adj2v::Adjacent2VertexCandidates) = length(get_candidates(adj2v))

@inline get_candidate(adj2v::Adjacent2VertexCandidates, u) = get_candidates(adj2v)[u]

@inline Base.copy(adj2v::Adjacent2VertexCandidates) = Adjacent2VertexCandidates(copy(get_candidates(adj2v)))

@inline function add_candidate!(adj2v::Adjacent2VertexCandidates, u, v)
    candidates = get_candidates(adj2v)
    candidates[u] = v
    return adj2v
end

@inline function delete_candidate!(adj2v::Adjacent2VertexCandidates, u)
    candidates = get_candidates(adj2v)
    delete!(candidates, u)
    return adj2v
end

@inline function add_candidates!(adj2v, u, v, w)
    add_candidate!(adj2v, u, v)
    add_candidate!(adj2v, v, w)
    add_candidate!(adj2v, w, u)
    return adj2v
end