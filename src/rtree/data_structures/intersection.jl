struct RTreeIntersectionIterator
    tree::RTree
    bounding_box::BoundingBox
    cache::RTreeIntersectionCache
end
function RTreeIntersectionIterator(tree::RTree, bounding_box::BoundingBox, cache_id::Integer)
    cache = get_intersection_cache(tree)
    return RTreeIntersectionIterator(tree, bounding_box, cache[cache_id])
end
Base.IteratorSize(::Type{RTreeIntersectionIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{RTreeIntersectionIterator}) = Base.HasEltype()
Base.eltype(::Type{RTreeIntersectionIterator}) = DiametralBoundingBox

struct RTreeIntersectionIteratorState
    leaf::Leaf{Branch}
    node_indices::Vector{Int}
    need_tests::BitVector
end

@enum QueryResult begin
    Contains 
    Intersects 
end

get_state(state::RTreeIntersectionIteratorState) = get_child(state.leaf, state.node_indices[begin])
