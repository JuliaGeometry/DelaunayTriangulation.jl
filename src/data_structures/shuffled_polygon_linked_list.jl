"""
    ShuffledPolygonLinkedList{I,T}

Data structure for representing a polygon as a doubly-linked list. In the 
descriptions below, `π` is used to denote the `shuffled_indices` vector.

# Fields 
- `next::Vector{I}`: The `next` vertices, so that `next[π[i]]` is the vertex after `S[π[i]]`.
- `prev::Vector{I}`: The `prev` vertices, so that `prev[π[i]]` is the vertex before `S[π[i]]`.
- `shuffled_indices::Vector{I}`: The shuffled indices of the vertices, so that `S[π[i]]` is the `i`th vertex.
- `k::I`: The number of vertices in the polygon.
- `S::T`: The vertices of the polygon. This should not be a circular vector, i.e. `S[begin] ≠ S[end]`, and must use one-based indexing. Additionally, 
   the vertices must be provided in counter-clockwise order - this is NOT checked.

# Constructor 
To construct this, use 

    ShuffledPolygonLinkedList(S::Vector; rng::AbstractRNG=Random.default_rng())

The argument `rng` is used for shuffling the `shuffled_indices` vector.
"""
struct ShuffledPolygonLinkedList{I,T}
    next::Vector{I}
    prev::Vector{I}
    shuffled_indices::Vector{I}
    k::I
    S::T
    function ShuffledPolygonLinkedList(next::Vector{I}, prev::Vector{I}, shuffled_indices::Vector{I}, k::I, S::T) where {I,T}
        Base.require_one_based_indexing(S)
        @assert !is_circular(S) "S must not be circular."
        @assert length(next) == length(prev) == length(shuffled_indices) == k == length(S) "The lengths of next, prev, shuffled_indices, and S must be equal."
        new{I, T}(next, prev, shuffled_indices, k, S)
    end
end
function ShuffledPolygonLinkedList(S::AbstractVector{I}; rng::AbstractRNG=Random.default_rng()) where {I} 
    k = I(length(S))
    shuffled_indices = collect(I, eachindex(S))
    next = zeros(I, k)
    prev = zeros(I, k)
    list = ShuffledPolygonLinkedList(next, prev, shuffled_indices, k, S)
    reset!(list; rng)
    return list
end

"""
    reset!(list::ShuffledPolygonLinkedList; rng::AbstractRNG=Random.default_rng())

Resets the linked `list`, so that `list.next[i] = mod1(i+1, list.k)` and `list.prev[i] = mod1(i-1, list.k)`,
and also reshuffles the `list.shuffled_indices` vector.
"""
function reset!(list::ShuffledPolygonLinkedList; rng::AbstractRNG=Random.default_rng())
    for i in 1:list.k 
        list.next[i] = mod1(i+1, list.k)
        list.prev[i] = mod1(i-1, list.k)
    end
    shuffle!(rng, list.shuffled_indices)
    return list
end

"""
    get_triplet(list::ShuffledPolygonLinkedList, i) -> (Vertex, Vertex, Vertex)

Returns `(u, v, w) = (S[πᵢ], S[next[πᵢ]], S[prev[πᵢ]])`, where `πᵢ = list.shuffled_indices[i]` and `S = list.S`.
"""
function get_triplet(list::ShuffledPolygonLinkedList, i)
    πᵢ = list.shuffled_indices[i]
    u = list.S[πᵢ]
    v = list.S[list.next[πᵢ]]
    w = list.S[list.prev[πᵢ]]
    return u, v, w
end

"""
    delete_vertex!(list::ShuffledPolygonLinkedList, i)

Deletes the vertex `S[πᵢ]` from the linked `list`, where `πᵢ = list.shuffled_indices[i]`.
This deletion is done symbolically rather than by mutating the vectors in `list`. In particular, 
we perform

- `next[prev[πᵢ]] = next[πᵢ]`
- `prev[next[πᵢ]] = prev[πᵢ]`

which is the same as removing `S[πᵢ]` from the linked `list`.
"""
function delete_vertex!(list::ShuffledPolygonLinkedList, i)
    πᵢ = list.shuffled_indices[i] 
    list.next[list.prev[πᵢ]] = list.next[πᵢ]
    list.prev[list.next[πᵢ]] = list.prev[πᵢ]
    return list
end

"""
    swap_permutation!(list::ShuffledPolygonLinkedList, i, j)

Reorders the permutation `list.shuffled_indices` of the linked `list`, swapping 
`πᵢ` and `πⱼ` where `πₖ = list.shuffled_indices[k]`.
"""
function swap_permutation!(list::ShuffledPolygonLinkedList, i, j)
    list.shuffled_indices[i], list.shuffled_indices[j] = list.shuffled_indices[j], list.shuffled_indices[i]
    return list
end