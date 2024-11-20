# Interface
@inline construct_edge(::Type{NTuple{2,I}}, i, j) where {I} = (I(i), I(j))
@inline construct_edge(::Type{Vector{I}}, i, j) where {I} = I[i, j]
@inline construct_edge(::Type{A}, i, j) where {A<:AbstractVector} = A((I(i), I(j)))

@inline initial(e) = e[1]
@inline terminal(e) = e[2]

@inline num_edges(E) = length(E)

@inline edge_type(E) = eltype(E)

@inline contains_edge(e, E) = e ∈ E

@inline add_to_edges!(E, e) = push!(E, e)

@inline delete_from_edges!(E, e) = delete!(E, e)

@inline each_edge(E) = E
@inline each_edge(E::AbstractMatrix) = eachcol(E)

@inline random_edge(E) = rand(E)
@inline random_edge(rng::Random.AbstractRNG, E) = rand(rng, E)

# Derived 
@inline edge_vertices(e) = (initial(e), terminal(e))

@inline reverse_edge(e::E) where {E} = construct_edge(E, terminal(e), initial(e))

@inline function compare_unoriented_edges(u, v)
    u₁, u₂ = edge_vertices(u)
    v₁, v₂ = edge_vertices(v)
    return (u₁, u₂) == (v₁, v₂) || (u₁, u₂) == (v₂, v₁)
end

@inline function contains_edge(i, j, E)
    et = edge_type(E)
    e = construct_edge(et, i, j)
    return contains_edge(e, E)
end

@inline function contains_unoriented_edge(e, E)
    e′ = reverse_edge(e)
    return contains_edge(e, E) || contains_edge(e′, E)
end

@inline function add_edge!(E, e...)
    foreach(e) do v
        add_to_edges!(E, v)
    end
    return E
end

@inline function delete_edge!(E, e...)
    foreach(e) do v
        delete_from_edges!(E, v)
    end
    return E
end

@inline function delete_unoriented_edge!(E, e)
    e′ = reverse_edge(e)
    delete_edge!(E, e, e′)
    return E
end

function compare_unoriented_edge_collections(E, F)
    num_edges(E) ≠ num_edges(F) && return false
    for e in each_edge(E)
        contains_unoriented_edge(e, F) || return false
    end
    return true
end

@inline function get_shared_vertex(e, f)
    u, v = edge_vertices(e)
    w, x = edge_vertices(f)
    if u == w || u == x
        return u
    elseif v == w || v == x
        return v
    else
        return oftype(u, ∅)
    end
end

@inline function edges_are_disjoint(e, e′)
    w = get_shared_vertex(e, e′)
    return w == oftype(w, ∅)
end