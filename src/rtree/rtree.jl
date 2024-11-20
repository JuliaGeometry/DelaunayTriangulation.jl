function Base.insert!(tree::RTree, bounding_box) # See e.g. https://rutgers-db.github.io/cs541-fall19/slides/notes4.pdf for discussions about overflow
    increment_num_elements!(tree)
    return insert!(tree::RTree, bounding_box, 1)
end
function Base.insert!(tree::RTree, bounding_box, level)
    subtree = find_subtree(tree, get_bounding_box(bounding_box), level)
    return insert!(subtree, bounding_box, tree)
end

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

function find_subtree(tree, bounding_box, level)
    node = get_root(tree)::Union{Branch, Leaf{Branch}}
    while get_level(node) > level
        min_enlargement = minimise_enlargement(node, bounding_box)
        node = get_child(node, min_enlargement.idx)::Union{Branch, Leaf{Branch}}
    end
    return node
end

struct EnlargementValues
    enlargement::Float64
    idx::Int
    area::Float64
    bounding_box::BoundingBox # the bounding_box being compared with for enlargement
end
EnlargementValues(bounding_box) = EnlargementValues(Inf, 0, Inf, bounding_box)

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

function minimise_enlargement(node, bounding_box)
    values = EnlargementValues(bounding_box)
    for (idx, child) in enumerate(get_children(node))
        values = update_enlargement(values, child, idx)
    end
    return values
end

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

function Base.delete!(tree::RTree, id_bounding_box::DiametralBoundingBox)
    decrement_num_elements!(tree)
    leaf, idx = find_bounding_box(tree, id_bounding_box)
    detach!(leaf, idx)
    detached_cache = get_detached_cache(tree)
    collapse_after_deletion!(leaf, tree, detached_cache)
    insert_detached!(tree, detached_cache)
    return tree
end

function find_bounding_box(tree::RTree, id_bounding_box::DiametralBoundingBox)
    idx_leaf = find_bounding_box(get_root(tree), id_bounding_box)::Tuple{Leaf{Branch}, Int}
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
        return nothing # node
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
    return nothing
end

function insert_detached!(tree::RTree, detached)
    isempty(detached) && return tree
    sort!(detached, by = get_level, rev = true)
    for node in detached
        for child in get_children(node)
            insert!(tree, child, get_level(node))
        end
        cache_node!(tree, node)
    end
    empty!(detached)
    return tree
end

function get_intersections(tree::RTree, bounding_box::BoundingBox; cache_id = 1)
    return RTreeIntersectionIterator(tree, bounding_box, cache_id)
end

function get_intersections(tree::RTree, point::NTuple{2, <:Number}; cache_id = 1)
    return RTreeIntersectionIterator(tree, BoundingBox(point), cache_id)
end

function Base.iterate(itr::RTreeIntersectionIterator)
    tree = itr.tree
    root = get_root(tree)
    num_children(root) == 0 && return nothing
    root_match = test_intersection(root, itr)
    root_match == Outside && return nothing
    cache = itr.cache
    node_indices = get_node_indices(cache)
    need_tests = get_need_tests(cache)
    height = get_height(tree)
    resize!(node_indices, height)
    resize!(need_tests, height)
    intersects = root_match == Intersects
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

function test_intersection(node::AbstractNode, itr::RTreeIntersectionIterator)
    node_rect = get_bounding_box(node)
    if itr.bounding_box ∋ node_rect
        return Contains
    elseif !isempty(node_rect ∩ itr.bounding_box)
        return Intersects
    else
        return Outside
    end
end

function test_intersection(element, itr::RTreeIntersectionIterator)
    element_rect = get_bounding_box(element)
    if !isempty(element_rect ∩ itr.bounding_box)
        return Contains
    else
        return Outside
    end
end

function get_next_child(node::AbstractNode, start_idx, need_tests, itr)
    res = need_tests ? Outside : Contains
    !need_tests && return start_idx, res
    while start_idx ≤ num_children(node) && ((res = test_intersection(get_child(node, start_idx), itr)) == Outside)
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
                need_tests[level] = query ≠ Contains
            else
                state = RTreeIntersectionIteratorState(node, node_indices, need_tests)
                return get_state(state), state
            end
        end
    end
end
