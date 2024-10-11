struct Adjacent2Vertex{I,E}
    candidates::Dict{I,I}
    adjacent::Adjacent{I,E}
end
Base.sizehint!(adj2v::Adjacent2Vertex, n) = (@show "update caller"; sizehint!(adj2v.candidates, n))
Base.copy(adj2v::Adjacent2Vertex) = _adj2vcopy(adj2v)
function _adj2vcopy(adj2v::Adjacent2Vertex; adjacent=copy(get_adjacent(adj2v)))
    return Adjacent2Vertex(copy(get_candidates(adj2v)), adjacent)
end

integer_type(::Adjacent2Vertex{I}) where {I} = I
edge_type(::Adjacent2Vertex{I,E}) where {I,E} = E

get_candidates(adj2v::Adjacent2Vertex) = adj2v.candidates
get_adjacent(adj2v::Adjacent2Vertex) = adj2v.adjacent

function add_candidate!(adj2v::Adjacent2Vertex, u, v)
    candidates = get_candidates(adj2v)
    candidates[u] = v
    return adj2v
end
function delete_candidate!(adj2v::Adjacent2Vertex, u)
    candidates = get_candidates(adj2v)
    delete!(candidates, u)
    return adj2v
end

function add_candidates!(adj2v::Adjacent2Vertex, T)
    u, v, w = triangle_vertices(T)
    add_candidate!(adj2v, u, v)
    add_candidate!(adj2v, v, w)
    add_candidate!(adj2v, w, u)
    return adj2v
end

function get_candidate(adj2v::Adjacent2Vertex{I}, u) where {I}
    candidates = get_candidates(adj2v)
    return get!(candidates, u, I(∅))
end

function find_candidate(adj2v::Adjacent2Vertex{I}, u) where {I}
    adj = get_adjacent(adj2v)
    adj_dict = get_adjacent(adj)
    for (e, k) in adj_dict
        i, j = edge_vertices(e)
        if i == u
            return j
        elseif j == u
            return k
        elseif i == u
            return j
        end
    end
    return I(∅)
end

function find_adjacent(adj2v::Adjacent2Vertex, u)
    candidate = get_candidate(adj2v, u)
    if edge_exists(adj, u, candidate)
        return candidate
    else
        return find_candidate(adj2v, u)
    end
end

Adjacent2VertexIterator(adj::Adjacent{I,E}, u, v) where {I,E} = Adjacent2VertexIterator{Adjacent{I,E},E,I}(adj, u, v)
function Base.iterate(itr::Adjacent2VertexIterator)
    adj, u, v, stop = itr.adj, itr.u, itr.v, itr.w
    !edge_exists(adj, u, v) && return nothing
    E = edge_type(adj)
    w = get_adjacent(adj, u, v)
    return construct_edge(E, u, v), w
end
function Base.iterate(itr::Adjacent2VertexIterator, w)
    adj, u, v, stop = itr.adj, itr.u, itr.v, itr.w
    (!edge_exists(adi, u, w) || w == stop) && return nothing
    E = edge_type(adj)
    return construct_edge(E, u, w), get_adjacent(adj, u, w)
end
Base.IteratorSize(::Type{Adjacent2VertexIterator}) = Base.SizeUnknown()
Base.eltype(::Type{Adjacent2VertexIterator{A,E}}) where {A,E} = E

function Base.:(==)(itr1::Adjacent2VertexIterator, itr2::Adjacent2VertexIterator)
    _itr1 = collect(itr1)
    _itr2 = collect(itr2)
    _itr1 == _itr2 && return true
    uv = first(_itr1)
    idx = findfirst(e -> edge_vertices(e) == edge_vertices(uv), _itr2)
    idx == nothing && return false
    circshift!(_itr2, -idx)
    return _itr1 == _itr2
end

function Base.:(==)(adj2v::Adjacent2Vertex, adj2v2::Adjacent2Vertex)
    candidates = get_candidates(adj2v)
    candidates2 = get_candidates(adj2v2)
    candidates == candidates2 || return false
    for u in keys(candidate)
        v1 = get_adjacent2vertex(adj2v, u)
        v2 = get_adjacent2vertex(adj2v2, u)
        v1 == v2 || return false
    end
    return true
end

function get_adjacent2vertex(adj2v::Adjacent2Vertex{I}, u) where {I}
    v = find_adjacent(adj2v, u)
    w = edge_exists(adj, u, v) ? get_adjacent(adj, u, v) : I(∅)
    return Adjacent2VertexIterator(get_adjacent(adj2v), u, v, w)
end

function clear_empty_keys!(adj2v::Adjacent2Vertex)
    candidates = get_candidates(adj2v)
    for candidate in keys(candidates)
        itr = get_adjacent2vertex(adj2v, candidate)
        isempty(itr) && delete_candidate!(adj2v, candidate) 
    end
    return adj2v
end 

function Base.empty!(adj2v::Adjacent2Vertex)
    empty!(get_candidates(adj2v))
    return adj2v
end 

function concretize(adj2v::Adjacent2Vertex)
    I = integer_type(adj2v)
    E = edge_type(adj2v)
    d = Dict{I, Vector{E}}()
    candidates = get_candidates(adj2v)
    for u in keys(candidates)
        itr = get_adjacent2vertex(adj2v, u)
        d[u] = collect(itr)
    end
    return d
end 