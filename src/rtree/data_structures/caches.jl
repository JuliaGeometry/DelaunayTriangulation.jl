struct NodeCache{Node, Child} # similar to why we use TriangulationCache. Think of it like implementing some CapacityVector that fails to push if it's full.
    cache::Vector{Node}
    size_limit::Int
    function NodeCache{Node, Child}(size_limit::Int) where {Node, Child}
        cache = Node[]
        sizehint!(cache, size_limit)
        return new{Node, Child}(cache, size_limit)
    end
end

Base.:(==)(cache1::NodeCache, cache2::NodeCache) = get_size_limit(cache1) == get_size_limit(cache2) && cache1.cache == cache2.cache

function Base.copy(cache::NodeCache{Node, Child}) where {Node, Child}
    # Just generates a new cache
    return NodeCache{Node, Child}(get_size_limit(cache))
end

const BranchCache = NodeCache{Branch, Branch}

const TwigCache = NodeCache{Branch, Leaf{Branch}}

const LeafCache = NodeCache{Leaf{Branch}, DiametralBoundingBox}

Base.length(cache::NodeCache) = length(cache.cache)

Base.isempty(cache::NodeCache) = isempty(cache.cache)

Base.pop!(cache::NodeCache) = pop!(cache.cache)

Base.push!(cache::NodeCache, node) = push!(cache.cache, node)

get_size_limit(cache::NodeCache) = cache.size_limit

spawn_node(::BranchCache) = Branch()
spawn_node(::TwigCache) = Branch(nothing, Leaf{Branch})
spawn_node(::LeafCache) = Leaf()

spawn_node!(cache::NodeCache) = isempty(cache) ? spawn_node(cache) : pop!(cache)

is_full(cache::NodeCache) = length(cache) ≥ get_size_limit(cache)

function cache_node!(cache::NodeCache, node)
    is_full(cache) || push!(cache, node)
    return cache
end

struct RTreeIntersectionCache
    node_indices::Vector{Int}
    need_tests::BitVector
end
RTreeIntersectionCache() = RTreeIntersectionCache(Vector{Int}(), BitVector())
Base.sizehint!(cache::RTreeIntersectionCache, n) = (sizehint!(get_node_indices(cache), n); sizehint!(get_need_tests(cache), n))

Base.:(==)(cache1::RTreeIntersectionCache, cache2::RTreeIntersectionCache) = get_node_indices(cache1) == get_node_indices(cache2) && get_need_tests(cache1) == get_need_tests(cache2)

Base.copy(::RTreeIntersectionCache) = RTreeIntersectionCache() # Just generates a new cache

get_node_indices(cache::RTreeIntersectionCache) = cache.node_indices

get_need_tests(cache::RTreeIntersectionCache) = cache.need_tests