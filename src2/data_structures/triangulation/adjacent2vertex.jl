struct Adjacent2Vertex{I,E}
    candidates::Dict{I,I}
    adjacent::Adjacent{I,E}
end
Base.sizehint!(adj2v::Adjacent2Vertex, n) = sizehint!(get_candidates(adj2v), n)
Base.copy(adj2v::Adjacent2Vertex) = _adj2vcopy(adj2v)
function _adj2vcopy(adj2v::Adjacent2Vertex; adjacent=copy(get_adjacent(adj2v)))
    return Adjacent2Vertex(copy(get_candidates(adj2v)), adjacent)
end

Base.length(adj2v::Adjacent2Vertex) = length(get_candidates(adj2v))

# Below is based on the show methods for AbstractDict in base/show.jl in Julia
function Base.show(io::IO, ::MIME"text/plain", itr::Adjacent2Vertex{I,E}) where {I,E}
    length(itr) == 0 && return show(io, itr)

    _keytype = integer_type(itr)
    _valtype = Adjacent2VertexIterator{Adjacent{I,E},E,I}

    recur_io = IOContext(io, :SHOWN_DICT => itr)
    limit = Base.get(io, :limit, false)::Bool
    if !haskey(io, :compact)
        recur_io = IOContext(recur_io, :compact => true)
    end
    recur_io_k = IOContext(recur_io, :typeinfo=>_keytype)
    recur_io_v = IOContext(recur_io, :typeinfo=>_valtype)

    summary(io, itr)
    length(itr) == 0 && return 
    print(io, ":")
    _show_circular(io, itr) && return
    if limit 
        sz = displaysize(io)
        rows, cols = sz[1] - 3, sz[2]
        rows < 2 && (print(io, " …"); return)
        cols < 12 && (cols = 12)
        cols -= 6 
        rows -= 1 

        ks = Vector{String}(undef, min(rows, length(itr)))
        vs = Vector{String}(undef, min(rows, length(itr)))
        keywidth = 0
        valwidth = 0
        for (i, k) in enumerate(keys(get_candidates(itr)))
            v = collect(get_adjacent2vertex(itr, k))
            i > rows && break
            ks[i] = sprint(show, k, context=recur_io_k, sizehint=0)
            vs[i] = sprint(show, v, context=recur_io_v, sizehint=0)
            keywidth = clamp(textwidth(ks[i]), keywidth, cols)
            valwidth = clamp(textwidth(vs[i]), valwidth, cols)
        end
        if keywidth > max(div(cols, 2), cols - valwidth)
            keywidth = max(cld(cols, 3), cols - valwidth)
        end
    else 
        rows = cols = typemax(Int)
    end

    for (i, k) in enumerate(keys(get_candidates(itr)))
        v = collect(get_adjacent2vertex(itr, k))
        print(io, "\n  ")
        if i == rows < length(itr)
            print(io, rpad("⋮", keywidth), " => ⋮")
            break
        end
        key = sprint(show, k, context=recur_io_k, sizehint=0)
        print(recur_io, key)
        print(io, " => ")

        show(recur_io_v, v)
    end
end

function Adjacent2Vertex(adj::Adjacent{I,E}) where {I,E}
    candidates = Dict{I,I}()
    sizehint!(candidates, length(adj) ÷ 3)
    adj2v = Adjacent2Vertex{I,E}(candidates, adj)
    for (e, w) in adj
        u, v = edge_vertices(e)
        add_candidates!(adj2v, u, v, w)
    end
    return adj2v
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

function add_candidates!(adj2v::Adjacent2Vertex, u, v, w)
    add_candidate!(adj2v, u, v)
    add_candidate!(adj2v, v, w)
    add_candidate!(adj2v, w, u)
    return adj2v
end

function get_candidate(adj2v::Adjacent2Vertex{I}, u) where {I}
    candidates = get_candidates(adj2v)
    return get!(candidates, u, I(∅))
end

function find_candidate(adj2v, u) 
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
    return integer_type(adj2v)(∅)
end

edge_exists(adj2v::Adjacent2Vertex, i, j) = edge_exists(get_adjacent(adj2v), i, j)

function find_adjacent(adj2v, u)
    candidate = get_candidate(adj2v, u)
    if edge_exists(adj2v, u, candidate)
        return candidate
    else
        return find_candidate(adj2v, u)
    end
end

struct Adjacent2VertexIterator{A,E,I}
    adj::A # parametric so that this can be used for both Adjacent and Triangulation
    u::I
    v::I
end
Adjacent2VertexIterator(adj::Adjacent{I,E}, u, v) where {I,E} = Adjacent2VertexIterator{Adjacent{I,E},E,I}(adj, u, v)
Adjacent2VertexIterator(adj2v::Adjacent2Vertex, u, v) = Adjacent2VertexIterator(get_adjacent(adj2v), u, v)
function Base.iterate(itr::Adjacent2VertexIterator)
    adj, u, v = itr.adj, itr.u, itr.v
    w = get_adjacent(adj, u, v)
    !edge_exists(w) && return nothing
    E = edge_type(adj)
    return construct_edge(E, v, w), w
end
function Base.iterate(itr::Adjacent2VertexIterator, w)
    adj, u, v = itr.adj, itr.u, itr.v
    x = get_adjacent(adj, u, w)
    (!edge_exists(x) || w == v) && return nothing
    E = edge_type(adj)
    return construct_edge(E, w, x), x
end
Base.IteratorSize(::Type{<:Adjacent2VertexIterator}) = Base.SizeUnknown()
Base.eltype(::Type{Adjacent2VertexIterator{A,E,I}}) where {A,E,I} = E
function Base.length(itr::Adjacent2VertexIterator)
    n = 0
    for _ in itr
        n += 1
    end
    return n
end
integer_type(::Adjacent2VertexIterator{A,E,I}) where {A,E,I} = I
edge_type(::Adjacent2VertexIterator{A,E,I}) where {A,E,I} = E

function Random.rand(rng::Random.AbstractRNG, v::Random.SamplerTrivial{<:Adjacent2VertexIterator})
    itr = v[]
    I = integer_type(itr)
    n = length(itr)
    idx = rand(rng, 1:n)
    ctr = 0
    state = iterate(itr)[2]
    for _ in 2:(idx-1)
        state = iterate(itr, state)[2]
    end
    E = edge_type(itr)
    return iterate(itr, state)[1]::E 
end

# Below is based on the show methods for AbstractSet in base/show.jl in Julia
# TODO: Define _show and use generic show
function Base.summary(io::IO, itr::Adjacent2VertexIterator)
    # Since the iterator has an UnknownLength, we need to get the length by iterating 
    n = length(itr)
    Base.showarg(io, itr, true)
    print(io, " with ", n, (n == 1 ? " edge" : " edges"))
end
function Base.show(io::IO, ::MIME"text/plain", itr::Adjacent2VertexIterator)
    n = length(itr)

    isempty(itr) && return show(io, itr)
    recur_io = IOContext(io, :SHOWN_SET => itr)
    limit = Base.get(io, :limit, false)::Bool

    summary(io, itr)
    isempty(itr) && return
    print(io, ":")
    _show_circular(io, itr) && return
    if limit
        sz = displaysize(io)
        rows = sz[1] - 3
        rows < 2 && (print(io, " …"); return)
        rows -= 1
    else
        rows = typemax(Int)
    end

    for (i, v) in enumerate(itr)
        print(io, "\n  ")
        if i == rows < n
            print(io, rpad("⋮", 2))
        end
        show(recur_io, v)
    end
end
_show_circular(io::IO, @nospecialize(x)) = false
function _show_circular(io::IOContext, @nospecialize(x))
    d = 1
    for (k, v) in io.dict
        if k === :SHOWN_SET
            if v == x
                print(io, "#= circular reference @-$d =#")
                return true
            end
            d += 1
        end
    end
    return false
end

function Base.:(==)(itr1::Adjacent2VertexIterator, itr2::Adjacent2VertexIterator)
    _itr1 = collect(itr1)
    _itr2 = collect(itr2)
    _itr1 == _itr2 && return true
    uv = first(_itr1)
    idx = findfirst(e -> edge_vertices(e) == edge_vertices(uv), _itr2)
    idx == nothing && return false
    circshift!(_itr2, -idx+1)
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

function get_adjacent2vertex(adj2v, u) 
    v = find_adjacent(adj2v, u)
    return Adjacent2VertexIterator(adj2v, u, v)
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

function concretize_adjacent2vertex(adj2v)
    I = integer_type(adj2v)
    E = edge_type(adj2v)
    d = Dict{I,Vector{E}}()
    candidates = get_candidates(adj2v)
    for u in keys(candidates)
        itr = get_adjacent2vertex(adj2v, u)
        d[u] = collect(itr)
    end
    return d
end