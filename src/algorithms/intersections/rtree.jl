#=
This file has been based heavily on https://github.com/alyst/SpatialIndexing.jl (v0.15). It is a simplified version that only aims to support 
bounding boxes in 2D, and only supports add and delete operations as well as basic intersection queries. Moreover, we only implement a linear 
R-tree without bulk loading. The license for SpatialIndexing.jl is below (https://github.com/alyst/SpatialIndexing.jl/blob/135a456c108503527491923b5c202c7ae85ce5ed/LICENSE).

MIT License

Copyright (c) 2018 Alexey Stukalov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

"""
    insert!(tree::RTree, bounding_box[, level = 1]) -> Bool 

Inserts `bounding_box` into `tree`. Returns `true` if the `tree`'s bounding boxes had to be adjusted and `false` otherwise.
"""
function Base.insert!(tree::RTree, bounding_box) # See e.g. https://rutgers-db.github.io/cs541-fall19/slides/notes4.pdf for discussions about overflow
    increment_num_elements!(tree)
    return insert!(tree::RTree, bounding_box, 1)
end
function Base.insert!(tree::RTree, bounding_box, level)
    subtree = find_subtree(tree, get_bounding_box(bounding_box), level)
    return insert!(subtree, bounding_box, tree)
end

"""
    insert!(node::AbstractNode, child, tree::RTree) -> Bool

Inserts `child` into `node` in `tree`. Returns `true` if the `tree`'s bounding boxes had to be adjusted and `false` otherwise.
"""
function Base.insert!(node::AbstractNode, child, tree::RTree)
    if !is_full(node, tree)
        original_bounding_box = get_bounding_box(node)
        append!(node, child)
        child_rect = get_bounding_box(child)
        if (child_rect ∉ original_bounding_box) && has_parent(node)
            update_bounding_box!(get_parent(node), find_position_in_parent(node), original_bounding_box, tree)
            return true
        else
            return false
        end
    else
        return overflow_insert!(node, child, tree)
    end
end

"""
    append!(node::AbstractNode, child)

Appends `child` to `node`'s children. Also updates `node`'s bounding box.
"""
function Base.append!(node::AbstractNode, child::C) where {C}
    original_bounding_box = get_bounding_box(node)
    flag = has_children(node)
    C <: AbstractNode && set_parent!(child, node)
    add_child!(node, child)
    child_rect = get_bounding_box(child)
    new_bounding_box = flag ? (original_bounding_box ∪ child_rect) : child_rect
    set_bounding_box!(node, new_bounding_box)
    return node
end

"""
    find_subtree(tree, bounding_box, level) -> Union{Branch,Leaf{Branch}}

Returns the subtree of `tree` at `level` that `bounding_box` should be inserted into.
"""
function find_subtree(tree, bounding_box, level)
    node = get_root(tree)::Union{Branch,Leaf{Branch}}
    while get_level(node) > level
        min_enlargement = minimise_enlargement(node, bounding_box)
        node = get_child(node, min_enlargement.idx)::Union{Branch,Leaf{Branch}}
    end
    return node
end

"""
    EnlargementValues

Type for representing the values used in the minimisation of enlargement.

# Fields
- `enlargement::Float64`: The enlargement of the bounding box of the child being compared with.
- `idx::Int`: The index of the child being compared with.
- `area::Float64`: The area of the child being compared with.
- `bounding_box::BoundingBox`: The bounding box being compared with for enlargement.

# Constructor 

    EnlargementValues(enlargement, idx, area, bounding_box)
    EnlargementValues(bounding_box)
"""
struct EnlargementValues
    enlargement::Float64
    idx::Int
    area::Float64
    bounding_box::BoundingBox # the bounding_box being compared with for enlargement
end
EnlargementValues(bounding_box) = EnlargementValues(Inf, 0, Inf, bounding_box)

"""
    update_enlargement(values::EnlargementValues, child, idx) -> EnlargementValues

Compare the enlargement state in `values` with `child` and return the updated `values` if necessary.
"""
function update_enlargement(values, child, idx)
    child_rect = get_bounding_box(child)
    child_area = get_area(child_rect)
    enlargement = get_area(child_rect ∪ values.bounding_box) - child_area
    if enlargement < values.enlargement || (enlargement == values.enlargement && child_area < values.area)
        return EnlargementValues(enlargement, idx, child_area, values.bounding_box)
    else
        return values
    end
end

"""
    minimise_enlargement(node, bounding_box) -> EnlargementValues

Returns the [`EnlargementValues`](@ref) associated with the child of `node` that minimises the enlargement, where enlargement 
is defined as the difference between the area of `bounding_box` and the area of the child's bounding box.
"""
function minimise_enlargement(node, bounding_box)
    values = EnlargementValues(bounding_box)
    for (idx, child) in enumerate(get_children(node))
        values = update_enlargement(values, child, idx)
    end
    return values
end

"""
    update_bounding_box!(node)

Updates the bounding box of `node` to be the union of its children's bounding boxes.
"""
function update_bounding_box!(node)
    first_child = get_child(node, 1)
    bounding_box = get_bounding_box(first_child)
    for child_idx in 2:num_children(node)
        child = get_child(node, child_idx)
        child_rect = get_bounding_box(child)
        bounding_box = bounding_box ∪ child_rect
    end
    set_bounding_box!(node, bounding_box)
    return node
end

"""
    update_bounding_box!(node, idx, original_bounding_box, tree)

Updates the bounding box of `node` to be the union of its children's bounding boxes. If `node` has a parent, updates the
bounding box of the parent as well if needed.
"""
function update_bounding_box!(node, idx, original_bounding_box, tree)
    child = get_child(node, idx)
    child_rect = get_bounding_box(child)
    node_rect = get_bounding_box(node)
    rect_is_old = (child_rect ∉ original_bounding_box) || is_touching(node_rect, original_bounding_box)
    rect_is_old && update_bounding_box!(node)
    if rect_is_old && has_parent(node)
        update_bounding_box!(get_parent(node), find_position_in_parent(node), original_bounding_box, tree)
    end
    return node
end

"""
    overflow_insert!(node, child, tree) -> Bool

Inserts `child` into `node` in `tree` when `node` is full. Returns `true`.
"""
function overflow_insert!(node, child, tree)
    original_bounding_box = get_bounding_box(node)
    append!(node, child)
    left, right = split!(node, tree)
    if is_root(node, tree)
        new_root = spawn_node!(tree, Branch, get_level(node) + 1)
        append!(new_root, left)
        append!(new_root, right)
        set_root!(tree, new_root)
    else
        replace!(node, left, right, original_bounding_box, tree)
    end
    return true
end

"""
    split!(node::AbstractNode, tree::RTree) -> Branch, Branch

Splits `node` into two other nodes using a linear splitting rule. Returns the two new nodes.
"""
function split!(node::N, tree) where {N}
    m = get_min_nodes(tree)
    level = get_level(node)
    i, j = split_seeds(node)
    i_spawn = spawn_node!(tree, N, level)
    j_spawn = spawn_node!(tree, N, level)
    childᵢ = get_child(node, i)
    childⱼ = get_child(node, j)
    left = append!(i_spawn, childᵢ)
    right = append!(j_spawn, childⱼ)
    free_cache = get_free_cache(tree)
    resize!(free_cache, num_children(node))
    fill!(free_cache, true)
    free_cache[i] = free_cache[j] = false
    nfree = num_children(node) - 2
    while nfree > 0
        left_at_min = num_children(left) + nfree == m
        right_at_min = num_children(right) + nfree == m
        if left_at_min || right_at_min
            sink = left_at_min ? left : right
            idx = findfirst(free_cache)
            while !isnothing(idx)
                append!(sink, get_child(node, idx))
                free_cache[idx] = false
                idx = findnext(free_cache, idx + 1)
            end
            nfree = 0
        else
            left_rect = get_bounding_box(left)
            left_area = get_area(left_rect)
            right_rect = get_bounding_box(right)
            right_area = get_area(right_rect)
            idx = findfirst(free_cache)
            child = get_child(node, idx)
            child_rect = get_bounding_box(child)
            left_enlargement = get_area(left_rect ∪ child_rect) - left_area
            right_enlargement = get_area(right_rect ∪ child_rect) - right_area
            select_left_enlargement = left_enlargement
            select_right_enlargement = right_enlargement
            select_child_idx = idx
            sink =
                select_left_enlargement < select_right_enlargement || (
                    select_left_enlargement == select_right_enlargement && (
                        left_area < right_area || (
                            left_area == right_area && num_children(left) ≤ num_children(right)
                        )
                    )
                ) ? left : right
            free_cache[select_child_idx] = false
            nfree -= 1
            append!(sink, get_child(node, select_child_idx))
        end
    end
    return left, right
end

"""
    split_seeds(node::AbstractNode) -> NTuple{2, Int}

Returns the indices of two children in `node` used to initiate the split in [`split!`](@ref).
"""
function split_seeds(node)
    i = j = 0
    max_separation = -Inf
    for d in 1:2 # ℝ²
        local argmax_low, argmin_high
        first_node = get_child(node, 1)
        first_rect = get_bounding_box(first_node)
        argmax_low = argmin_high = 1
        min_low = max_low = get_bl_corner(first_rect)[d]
        min_high = max_high = get_tr_corner(first_rect)[d]
        for i in 2:num_children(node)
            child = get_child(node, i)
            child_rect = get_bounding_box(child)
            low = get_bl_corner(child_rect)[d]
            high = get_tr_corner(child_rect)[d]
            if low > max_low
                max_low = low
                argmax_low = i
            elseif low < min_low
                min_low = low
            end
            if high < min_high
                min_high = high
                argmin_high = i
            elseif high > max_high
                max_high = high
            end
        end
        if max_high > min_low
            separation = (max_low - min_high) / (max_high - min_low)
        else
            separation = max_low - min_high
        end
        if separation > max_separation
            max_separation = separation
            i = argmin_high
            j = argmax_low
        end
    end
    if i == j
        j = j == 1 ? j + 1 : j - 1
    end
    return i, j
end

"""
    replace!(node::AbstractNode, left, right, original_bounding_box, tree::RTree)

Replaces the `node` in `tree` with `left` and `right`. Returns `true` if the `tree`'s bounding boxes had to be adjusted and `false` otherwise.
"""
function replace!(node, left, right, original_bounding_box, tree)
    parent = get_parent(node)
    node_idx = find_position_in_parent(node)
    old_parent_bounding_box = get_bounding_box(parent)
    left_rect = get_bounding_box(left)
    right_rect = get_bounding_box(right)
    rect_is_old = (left_rect ∉ old_parent_bounding_box) || (right_rect ∉ old_parent_bounding_box) || is_touching(old_parent_bounding_box, original_bounding_box)
    parent[node_idx] = left
    rect_is_old && update_bounding_box!(parent)
    flag = insert!(parent, right, tree)
    if !flag && rect_is_old && has_parent(parent)
        update_bounding_box!(get_parent(parent), find_position_in_parent(parent), old_parent_bounding_box, tree)
    end
    return flag || rect_is_old
end

"""
    delete!(tree::RTree, id_bounding_box::DiametralBoundingBox) 

Deletes `id_bounding_box` from `tree`.
"""
function Base.delete!(tree::RTree, id_bounding_box::DiametralBoundingBox)
    decrement_num_elements!(tree)
    leaf, idx = find_bounding_box(tree, id_bounding_box)
    detach!(leaf, idx)
    detached_cache = get_detached_cache(tree)
    collapse_after_deletion!(leaf, tree, detached_cache)
    insert_detached!(tree, detached_cache)
    return tree
end

"""
    find_bounding_box(tree::RTree, id_bounding_box::DiametralBoundingBox) -> Tuple{Leaf{Branch}, Int}

Returns the leaf node and the index in the leaf node's children that `id_bounding_box` is associated with.
"""
function find_bounding_box(tree::RTree, id_bounding_box::DiametralBoundingBox)
    idx_leaf = find_bounding_box(get_root(tree), id_bounding_box)::Tuple{Leaf{Branch},Int}
    return idx_leaf[1]::Leaf{Branch}, idx_leaf[2]::Int
end
function find_bounding_box(branch::Branch, id_bounding_box::DiametralBoundingBox)
    bounding_box = get_bounding_box(id_bounding_box)
    for child in get_children(branch)
        child_rect = get_bounding_box(child)
        if child_rect ∋ bounding_box
            res = find_bounding_box(child, id_bounding_box)
            !isnothing(res) && return res
        end
    end
    return nothing
end
function find_bounding_box(leaf::Leaf, id_bounding_box::DiametralBoundingBox)
    bounding_edge = get_edge(id_bounding_box)
    for (idx, child) in enumerate(get_children(leaf))
        child_edge = get_edge(child)
        compare_unoriented_edges(child_edge, bounding_edge) && return (leaf, idx)
    end
    return nothing
end

"""
    detach!(node::AbstractNode, idx) -> Bool

Detaches the `idx`th child of `node`. Returns `true` if the `node`'s bounding box had to be adjusted and `false` otherwise.
"""
function detach!(node::AbstractNode, idx)
    if num_children(node) > 1
        child = get_child(node, idx)
        child_rect = get_bounding_box(child)
        if idx < num_children(node)
            last_child = pop_child!(node)
            set_child!(node, last_child, idx)
        else
            pop_child!(node)
        end
        node_rect = get_bounding_box(node)
        touching = is_touching(node_rect, child_rect)
        touching && update_bounding_box!(node)
        return touching
    else
        pop_child!(node)
        set_bounding_box!(node, InvalidBoundingBox)
        return true
    end
    return false
end

"""
    collapse_after_deletion!(node::AbstractNode, tree::RTree, detached) 

Condenses `tree` after a deletion of one of `node`'s children. The `detached` argument will contain the nodes that were detached from `tree` during the condensing process.
"""
function collapse_after_deletion!(node::AbstractNode, tree, detached)
    m = get_min_nodes(tree)
    if is_root(node, tree)
        original_root = node
        min_root_level = 1
        for nd in detached
            min_root_level = max(min_root_level, get_level(nd))
        end
        while (get_level(node) > min_root_level) && (num_children(node) == 1)
            set_root!(tree, get_child(node, 1))
            set_parent!(get_root(tree), nothing)
            cache_node!(tree, node)
            node = get_root(tree)
        end
        if (get_level(node) ≥ min_root_level) && !has_children(node)
            new_root = if min_root_level == 1
                spawn_node!(tree, Leaf{Branch}, min_root_level)
            else
                spawn_node!(tree, Branch, min_root_level)
            end
            set_root!(tree, new_root)
        elseif node === original_root
            update_bounding_box!(node)
        end
        return node
    else
        idx = find_position_in_parent(node)
        if num_children(node) < m
            detach!(get_parent(node), idx)
            push!(detached, node)
        else
            update_bounding_box!(get_parent(node))
        end
        return collapse_after_deletion!(get_parent(node), tree, detached)
    end
end

"""
    insert_detached!(tree::RTree, detached)

Given the `detached` nodes from [`collapse_after_deletion!`](@ref), inserts them back into `tree`.
"""
function insert_detached!(tree::RTree, detached)
    isempty(detached) && return tree
    sort!(detached, by=get_level, rev=true)
    for node in detached
        for child in get_children(node)
            insert!(tree, child, get_level(node))
        end
        cache_node!(tree, node)
    end
    empty!(detached)
    return tree
end

"""
    get_intersections(tree::RTree, bounding_box::BoundingBox; cache_id=1) -> RTreeIntersectionIterator

Returns an [`RTreeIntersectionIterator`](@ref) over the elements in `tree` that intersect with `bounding_box`.
`cache_id` must be `1` or `2`, and determines what cache to use for the intersection query.
"""
function get_intersections(tree::RTree, bounding_box::BoundingBox; cache_id=1)
    return RTreeIntersectionIterator(tree, bounding_box, cache_id)
end

"""
    get_intersections(tree::RTree, point::NTuple{2,<:Number}; cache_id=1) -> RTreeIntersectionIterator

Returns an [`RTreeIntersectionIterator`](@ref) over the elements in `tree` that intersect with `point`, representing `point` 
as a [`BoundingBox`](@ref) with zero width and height centered at `point`.
`cache_id` must be `1` or `2`, and determines what cache to use for the intersection query.
"""
function get_intersections(tree::RTree, point::NTuple{2,<:Number}; cache_id=1)
    return RTreeIntersectionIterator(tree, BoundingBox(point), cache_id)
end

"""
    iterate(itr::RTreeIntersectionIterator, state...)

Iterate over the next state of `itr` to find more intersections with the bounding box in `RTreeIntersectionIterator`.
"""
function Base.iterate(itr::RTreeIntersectionIterator)
    tree = itr.tree
    root = get_root(tree)
    num_children(root) == 0 && return nothing
    root_match = test_intersection(root, itr)
    root_match == QR.Outside && return nothing
    cache = itr.cache
    node_indices = get_node_indices(cache)
    need_tests = get_need_tests(cache)
    height = get_height(tree)
    resize!(node_indices, height)
    resize!(need_tests, height)
    intersects = root_match == QR.Intersects
    for i in 1:height
        node_indices[i] = 1
        need_tests[i] = intersects
    end
    return _iterate(itr, root, node_indices, need_tests)
end
function Base.iterate(itr::RTreeIntersectionIterator, state::RTreeIntersectionIteratorState)
    start_idx = get_next_child(state.leaf, state.node_indices[begin] + 1, state.need_tests[begin], itr)[1]
    state.node_indices[begin] = start_idx
    if start_idx ≤ num_children(state.leaf)
        return get_state(state), state
    else
        return _iterate(itr, state.leaf, state.node_indices, state.need_tests)
    end
end

"""
    test_intersection(node::AbstractNode, itr::RTreeIntersectionIterator) -> QueryResult

Tests whether `node` intersects with the bounding box in `itr`, returning a [`QueryResult`](@ref).
"""
function test_intersection(node::AbstractNode, itr::RTreeIntersectionIterator)
    node_rect = get_bounding_box(node)
    if itr.bounding_box ∋ node_rect
        return QR.Contains
    elseif !isempty(node_rect ∩ itr.bounding_box)
        return QR.Intersects
    else
        return QR.Outside
    end
end

function test_intersection(element, itr::RTreeIntersectionIterator)
    element_rect = get_bounding_box(element)
    if !isempty(element_rect ∩ itr.bounding_box)
        return QR.Contains
    else
        return QR.Outside
    end
end

"""
    get_next_child(node::AbstractNode, start_idx, need_tests, itr::RTreeIntersectionIterator) -> Int, QueryResult

Returns the index of the next child of `node` that intersects with the bounding box in `itr` and the [`QueryResult`](@ref) of the intersection.
"""
function get_next_child(node::AbstractNode, start_idx, need_tests, itr)
    res = need_tests ? QR.Outside : QR.Contains
    !need_tests && return start_idx, res
    while start_idx ≤ num_children(node) && ((res = test_intersection(get_child(node, start_idx), itr)) == QR.Outside)
        start_idx += 1
    end
    return start_idx, res
end

function _iterate(itr::RTreeIntersectionIterator, node, node_indices, need_tests)
    while true
        local level, start_idx
        level = get_level(node)
        start_idx = node_indices[level]
        next_start_idx, query = get_next_child(node, start_idx, need_tests[level], itr)
        if next_start_idx > num_children(node)
            while next_start_idx > num_children(node)
                !has_parent(node) && return nothing
                node = get_parent(node)
                level = get_level(node)
                node_indices[level] += 1
                next_start_idx = node_indices[level]
            end
            start_idx = next_start_idx
        else
            if next_start_idx > start_idx
                node_indices[level] = next_start_idx
            end
            if typeof(node) <: Branch
                node = get_child(node, next_start_idx)
                level = get_level(node)
                start_idx = 1
                node_indices[level] = 1
                need_tests[level] = query ≠ QR.Contains
            else
                state = RTreeIntersectionIteratorState(node, node_indices, need_tests)
                return get_state(state), state
            end
        end
    end
end