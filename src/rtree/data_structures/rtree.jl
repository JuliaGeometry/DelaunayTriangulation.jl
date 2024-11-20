mutable struct RTree # linear
    root::Union{Branch, Leaf{Branch}}
    num_elements::Int
    @const branch_cache::BranchCache
    @const twig_cache::TwigCache
    @const leaf_cache::LeafCache
    @const fill_factor::Float64
    @const free_cache::BitVector
    @const detached_cache::Vector{Union{Branch, Leaf{Branch}}}
    @const intersection_cache::NTuple{2, RTreeIntersectionCache}
end
function RTree(; size_limit = 100, fill_factor = 0.7) # https://en.wikipedia.org/wiki/R-tree: "however best performance has been experienced with a minimum fill of 30%–40%)
    branch_cache = BranchCache(size_limit)
    twig_cache = TwigCache(size_limit)
    leaf_cache = LeafCache(size_limit)
    root = spawn_node!(leaf_cache)
    num_elements = 0
    free_cache = BitVector()
    sizehint!(free_cache, size_limit)
    detached_cache = Vector{Union{Branch, Leaf{Branch}}}()
    sizehint!(detached_cache, size_limit)
    cache1, cache2 = RTreeIntersectionCache(), RTreeIntersectionCache()
    sizehint!(cache1, ceil(Int, log2(size_limit)))
    sizehint!(cache2, ceil(Int, log2(size_limit)))
    return RTree(
        root,
        num_elements,
        branch_cache,
        twig_cache,
        leaf_cache,
        fill_factor,
        free_cache,
        detached_cache,
        (cache1, cache2),
    )
end

function Base.copy(tree::RTree)
    root = copy(get_root(tree))
    branch_cache = copy(get_branch_cache(tree))
    twig_cache = copy(get_twig_cache(tree))
    leaf_cache = copy(get_leaf_cache(tree))
    fill_factor = get_fill_factor(tree)
    free_cache = copy(get_free_cache(tree))
    detached_cache = copy(get_detached_cache(tree))
    intersection_cache = copy.(get_intersection_cache(tree))
    return RTree(root, num_elements(tree), branch_cache, twig_cache, leaf_cache, fill_factor, free_cache, detached_cache, intersection_cache)
end

function Base.:(==)(tree1::RTree, tree2::RTree)
    num_elements(tree1) ≠ num_elements(tree2) && return false
    get_root(tree1) ≠ get_root(tree2) && return false
    return true
end

get_root(tree::RTree) = tree.root

get_branch_cache(tree::RTree) = tree.branch_cache

get_twig_cache(tree::RTree) = tree.twig_cache

get_leaf_cache(tree::RTree) = tree.leaf_cache

get_fill_factor(tree::RTree) = tree.fill_factor

get_height(tree::RTree) = get_level(get_root(tree))

get_bounding_box(tree::RTree) = get_bounding_box(get_root(tree))

get_size_limit(tree::RTree) = get_size_limit(get_branch_cache(tree)) # we assume that the caches all have the same size limit

is_full(node::AbstractNode, tree::RTree) = num_children(node) ≥ get_size_limit(tree)

increment_num_elements!(tree::RTree) = tree.num_elements += 1

decrement_num_elements!(tree::RTree) = tree.num_elements -= 1

is_root(node::AbstractNode, tree::RTree) = node === get_root(tree)

set_root!(tree::RTree, node::AbstractNode) = tree.root = node

num_elements(tree::RTree) = tree.num_elements

get_min_nodes(tree::RTree) = ceil(Int, get_fill_factor(tree) * get_size_limit(tree))

get_free_cache(tree::RTree) = tree.free_cache

get_detached_cache(tree::RTree) = tree.detached_cache

get_intersection_cache(tree::RTree) = tree.intersection_cache

function spawn_leaf!(tree::RTree, bounding_box::BoundingBox)
    leaf = spawn_node!(get_leaf_cache(tree))
    set_bounding_box!(leaf, bounding_box)
    return leaf
end

function spawn_branch!(tree::RTree, bounding_box::BoundingBox, level)
    cache = if level == 2
        get_twig_cache(tree)
    else
        get_branch_cache(tree)
    end
    branch = spawn_node!(cache)
    set_bounding_box!(branch, bounding_box)
    set_level!(branch, level)
    return branch
end 

function spawn_node!(tree::RTree, ::Type{N}, bounding_box::BoundingBox, level) where {N}
    if N <: Branch
        return spawn_branch!(tree, bounding_box, level)
    else
        return spawn_leaf!(tree, bounding_box)
    end
end
spawn_node!(tree::RTree, ::Type{N}, level) where {N} = spawn_node!(tree, N, InvalidBoundingBox, level)

function cache_node!(tree::RTree, leaf::Leaf)
    set_parent!(leaf, nothing)
    empty!(get_children(leaf))
    cache = get_leaf_cache(tree)
    cache_node!(cache, leaf)
    return tree
end
function cache_node!(tree::RTree, branch::Branch)
    set_parent!(branch, nothing)
    empty!(get_children(branch))
    cache = if get_child_type(branch) <: Leaf
        get_twig_cache(tree)
    else
        get_branch_cache(tree)
    end
    cache_node!(cache, branch)
    return tree
end