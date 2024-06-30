# Originally wanted this for the Bentley-Ottmann algorithm, but I didn't get around to implementing it yet.
# Keeping this around until we do get around to implementing it.

"""
    mutable struct BalancedBSTNode{K}

Struct representing a node in a balanced binary search tree.

# Fields
- `key::K`: The key associated with the node.
- `height::Int8`: The height of the node.
- `count::Int32`: The number of nodes in the subtree rooted at this node, including this node.
- `parent::Union{Nothing,BalancedBSTNode{K}}`: The parent of the node.
- `left::Union{Nothing,BalancedBSTNode{K}}`: The left child of the node.
- `right::Union{Nothing,BalancedBSTNode{K}}`: The right child of the node.

# Constructor

    BalancedBSTNode(key::K) where {K}

Constructs a new node with key `key`, height `1`, and count `1`. The parent, left child, and right child are set to `nothing`.
"""
mutable struct BalancedBSTNode{K}
    key::K
    height::Int8
    count::Int32
    parent::Union{Nothing,BalancedBSTNode{K}}
    left::Union{Nothing,BalancedBSTNode{K}}
    right::Union{Nothing,BalancedBSTNode{K}}
    function BalancedBSTNode(key::K) where {K}
        new{K}(key, 1, 1, nothing, nothing, nothing)
    end
end
function Base.show(io::IO, ::MIME"text/plain", node::BalancedBSTNode)
    print(io, "BalancedBSTNode with key $(get_key(node)), height $(get_height(node)), and count $(get_count(node))")
end

function inorder!(keys, node)
    left = get_left(node)
    inorder!(keys, left)
    push!(keys, get_key(node))
    right = get_right(node)
    inorder!(keys, right)
    return keys
end
function inorder!(keys, node::Nothing)
    return keys
end

"""
    get_key(node::BalancedBSTNode{K}) -> K

Returns the key associated with `node`.
"""
get_key(node::BalancedBSTNode) = node.key

"""
    set_key!(node::BalancedBSTNode{K}, key::K)

Sets the key associated with `node` to `key`.
"""
set_key!(node::BalancedBSTNode, key) = node.key = key

"""
    get_height(node::Union{Nothing, BalancedBSTNode}) -> Int8

Returns the height of `node`. If the `node` is `nothing`, returns `0`.
"""
get_height(node::BalancedBSTNode) = node.height
get_height(node::Nothing) = Int8(0)

"""
    set_height!(node::BalancedBSTNode[, height = compute_height(node)])

Sets the height of `node` to `height`. If `height` is not provided, the height is computed using [`compute_height`](@ref).
"""
set_height!(node::BalancedBSTNode) = set_height!(node, compute_height(node))
set_height!(node::BalancedBSTNode, height) = node.height = height

"""
    get_count(node::Union{Nothing, BalancedBSTNode}) -> Int32

Returns the count of `node`. If the `node` is `nothing`, returns `0`.
"""
get_count(node::BalancedBSTNode) = node.count
get_count(node::Nothing) = Int32(0)

"""
    set_count!(node::BalancedBSTNode[, count = compute_count(node)])

Sets the count of `node` to `count`. If `count` is not provided, the count is computed using [`compute_count`](@ref).
"""
set_count!(node::BalancedBSTNode) = set_count!(node, compute_count(node))
set_count!(node::BalancedBSTNode, count) = node.count = count

"""
    get_parent(node::BalancedBSTNode) -> Union{Nothing, BalancedBSTNode}

Returns the parent of `node`. If the `node` is `nothing`, returns `nothing`.
"""
get_parent(node::BalancedBSTNode) = node.parent

"""
    set_parent!(node::BalancedBSTNode, parent::Union{Nothing, BalancedBSTNode})

Sets the parent of `node` to `parent`. 
"""
set_parent!(node::BalancedBSTNode, parent) = node.parent = parent

"""
    get_left(node::BalancedBSTNode) -> Union{Nothing, BalancedBSTNode}

Returns the left child of `node`. If the `node` is `nothing`, returns `nothing`.
"""
get_left(node::BalancedBSTNode) = node.left

"""
    set_left!(node::BalancedBSTNode, left::Union{Nothing, BalancedBSTNode})

Sets the left child of `node` to `left`.
"""
set_left!(node::BalancedBSTNode, left) = node.left = left

"""
    get_right(node::BalancedBSTNode) -> Union{Nothing, BalancedBSTNode}

Returns the right child of `node`. If the `node` is `nothing`, returns `nothing`.
"""
get_right(node::BalancedBSTNode) = node.right

"""
    set_right!(node::BalancedBSTNode, right::Union{Nothing, BalancedBSTNode})

Sets the right child of `node` to `right`.
"""
set_right!(node::BalancedBSTNode, right) = node.right = right

"""
    has_parent(node::BalancedBSTNode) -> Bool

Returns `true` if `node` has a parent, `false` otherwise.
"""
has_parent(node::BalancedBSTNode) = !isnothing(get_parent(node))

"""
    has_left(node::BalancedBSTNode) -> Bool

Returns `true` if `node` has a left child, `false` otherwise.
"""
has_left(node::BalancedBSTNode) = !isnothing(get_left(node))

"""
    has_right(node::BalancedBSTNode) -> Bool

Returns `true` if `node` has a right child, `false` otherwise.
"""
has_right(node::BalancedBSTNode) = !isnothing(get_right(node))

"""
    mutable struct BalancedBST{K}

Struct representing a balanced binary search tree.

# Fields
- `root::Union{Nothing,BalancedBSTNode{K}}`: The root of the tree.
- `count::Int32`: The number of nodes in the tree.

!!! warning "Duplicate keys"

    Nodes with duplicate keys are not supported. If a duplicate key is inserted, the tree will not be modified.
"""
mutable struct BalancedBST{K}
    root::Union{Nothing,BalancedBSTNode{K}}
    count::Int32
    BalancedBST(root::BalancedBSTNode{K}, count) where {K} = new{K}(root, count)
    BalancedBST{K}(root::BalancedBSTNode{K}, count) where {K} = new{K}(root, count)
    BalancedBST{K}(::Nothing, count) where {K} = new{K}(nothing, count)
    # need to separate out the constructors to avoid unbound type argument issues from Aqua
end
function Base.show(io::IO, ::MIME"text/plain", tree::BalancedBST)
    count = get_count(tree)
    if count == 1
        print(io, "BalancedBST with 1 node")
    else
        print(io, "BalancedBST with $count nodes")
    end
end

"""
    inorder(tree::BalancedBST{K}) -> Vector{K}

Returns the inorder traversal of `tree`.
"""
function inorder(tree::BalancedBST{K}) where {K}
    keys = K[]
    sizehint!(keys, get_count(tree))
    root = get_root(tree)
    inorder!(keys, root)
    return keys
end

"""
    BalancedBST{K}() where {K}

Constructs a new empty balanced binary search tree.
"""
BalancedBST{K}() where {K} = BalancedBST{K}(nothing, 0)

"""
    get_root(tree::BalancedBST{K}) -> Union{Nothing,BalancedBSTNode{K}}

Returns the root of `tree`. If `tree` is empty, returns `nothing`.
"""
get_root(tree::BalancedBST) = tree.root

"""
    set_root!(tree::BalancedBST{K}, root::Union{Nothing,BalancedBSTNode{K}})

Sets the root of `tree` to `root`.
"""
set_root!(tree::BalancedBST, root) = tree.root = root

"""
    get_count(tree::BalancedBST{K}) -> Int32
"""
get_count(tree::BalancedBST) = tree.count

"""
    set_count!(tree::BalancedBST{K}, count::Int32)

Sets the count of `tree` to `count`.
"""
set_count!(tree::BalancedBST, count) = tree.count = count

"""
    has_root(tree::BalancedBST{K}) -> Bool

Returns `true` if `tree` has a root, `false` otherwise.
"""
has_root(tree::BalancedBST) = !isnothing(get_root(tree))

"""
    push!(tree::BalancedBST{K}, key::K)

Inserts `key` into `tree` if it is not already present.
"""
function Base.push!(tree::BalancedBST, key)
    haskey(tree, key) && return tree
    root = get_root(tree)
    new_root = insert_node!(root, key)
    set_root!(tree, new_root)
    set_count!(tree, get_count(tree) + 1)
    return tree
end

"""
    insert_node!(node::Union{Nothing,BalancedBSTNode{K}}, key::K) -> Union{Nothing,BalancedBSTNode{K}}

Inserts `key` into the subtree rooted at `node` if it is not already present. Returns the new root of the subtree.
"""
insert_node!(node::Nothing, key) = BalancedBSTNode(key)
function insert_node!(node::BalancedBSTNode, key)
    if key < get_key(node)
        left = get_left(node)
        new_left = insert_node!(left, key)
        set_left!(node, new_left)
    else
        right = get_right(node)
        new_right = insert_node!(right, key)
        set_right!(node, new_right)
    end
    set_count!(node)
    set_height!(node)
    subtree_balance = compute_balance(node)
    if subtree_balance > 1
        if key < get_key(get_left(node))
            return rotate_right!(node)
        else
            left = get_left(node)
            rotated_left = rotate_left!(left)
            set_left!(node, rotated_left)
            return rotate_right!(node)
        end
    elseif subtree_balance < -1
        if key > get_key(get_right(node))
            return rotate_left!(node)
        else
            right = get_right(node)
            rotated_right = rotate_right!(right)
            set_right!(node, rotated_right)
            return rotate_left!(node)
        end
    end
    return node
end

"""
    compute_count(node::Union{Nothing,BalancedBSTNode{K}}) -> Int32

Computes the count of the subtree rooted at `node`, i.e. the number of nodes in the subtree rooted at `node`, including `node`.
"""
function compute_count(node::Union{Nothing,BalancedBSTNode})
    if isnothing(node)
        return Int32(0)
    else
        left = get_left(node)
        right = get_right(node)
        left_count = get_count(left)
        right_count = get_count(right)
        return Int32(1) + left_count + right_count
    end
end

"""
    compute_height(node::Union{Nothing,BalancedBSTNode{K}}) -> Int8

Computes the height of the subtree rooted at `node`.
"""
function compute_height(node::Union{Nothing,BalancedBSTNode})
    if isnothing(node)
        return Int8(0)
    else
        left = get_left(node)
        right = get_right(node)
        left_height = get_height(left)
        right_height = get_height(right)
        return Int8(1) + max(left_height, right_height)
    end
end

"""
    compute_balance(node::Union{Nothing,BalancedBSTNode{K}}) -> Int8

Computes the balance of the subtree rooted at `node`. This is the difference between the left and right heights.
"""
function compute_balance(node::Union{Nothing,BalancedBSTNode})
    if isnothing(node)
        return Int8(0)
    else
        left = get_left(node)
        right = get_right(node)
        left_height = get_height(left)
        right_height = get_height(right)
        return left_height - right_height
    end
end

"""
    rotate_left!(parent::BalancedBSTNode{K}) -> BalancedBSTNode{K}

Rotates a subtree rooted at `parent` to the left, returning the new root of the subtree.
This local operation is used to preserve the binary search tree property after inserting or deleting a node. 
"""
function rotate_left!(parent::BalancedBSTNode)
    right_child = get_right(parent)
    left_grandchild = get_left(right_child)
    set_left!(right_child, parent)
    set_right!(parent, left_grandchild)
    set_height!(parent)
    set_height!(right_child)
    set_count!(parent)
    set_count!(right_child)
    return right_child
end

"""
    rotate_right!(parent::BalancedBSTNode{K}) -> BalancedBSTNode{K}

Rotates a subtree rooted at `parent` to the right, returning the new root of the subtree.
This local operation is used to preserve the binary search tree property after inserting or deleting a node. 
"""
function rotate_right!(parent::BalancedBSTNode)
    left_child = get_left(parent)
    right_grandchild = get_right(left_child)
    set_right!(left_child, parent)
    set_left!(parent, right_grandchild)
    set_height!(parent)
    set_height!(left_child)
    set_count!(parent)
    set_count!(left_child)
    return left_child
end

"""
    findfirst(tree::BalancedBST{K}, key::K) -> Union{Nothing,BalancedBSTNode{K}}

Returns the node in `tree` with key `key`. If no such node exists, returns `nothing`.
"""
function Base.findfirst(tree::BalancedBST, key)
    prev = nothing
    node = get_root(tree)
    while !isnothing(node) && get_key(node) ≠ key
        prev = node
        if key < get_key(node)
            node = get_left(node)
        else
            node = get_right(node)
        end
    end
    if isnothing(node)
        return prev
    else
        return node
    end
end

"""
    haskey(tree::BalancedBST{K}, key::K) -> Bool

Returns `true` if `tree` has a node with key `key`, `false` otherwise.
"""
function Base.haskey(tree::BalancedBST, key)
    !has_root(tree) && return false
    node = findfirst(tree, key)
    return get_key(node) == key
end

"""
    delete!(tree::BalancedBST{K}, key::K) -> BalancedBST{K}

Deletes the node in `tree` with key `key` if it exists. Returns `tree`.
"""
function Base.delete!(tree::BalancedBST, key)
    !haskey(tree, key) && return tree
    root = get_root(tree)
    new_root = delete_node!(root, key)
    set_root!(tree, new_root)
    set_count!(tree, get_count(tree) - 1)
    return tree
end

"""
    delete_node!(node::Union{Nothing,BalancedBSTNode{K}}, key::K) -> Union{Nothing,BalancedBSTNode{K}}

Deletes the node with key `key` from the subtree rooted at `node` if it exists. Returns the new root of the subtree.
"""
function delete_node!(node::BalancedBSTNode, key)
    if key < get_key(node)
        left = get_left(node)
        new_left = delete_node!(left, key)
        set_left!(node, new_left)
    elseif key > get_key(node)
        right = get_right(node)
        new_right = delete_node!(right, key)
        set_right!(node, new_right)
    else
        if !has_left(node)
            return get_right(node)
        elseif !has_right(node)
            return get_left(node)
        else
            min = _minimum(get_right(node))
            set_key!(node, get_key(min))
            right = get_right(node)
            new_right = delete_node!(right, get_key(min))
            set_right!(node, new_right)
        end
    end
    set_count!(node)
    set_height!(node)
    subtree_balance = compute_balance(node)
    if subtree_balance > 1
        left_balance = compute_balance(get_left(node))
        if left_balance ≥ 0
            return rotate_right!(node)
        else
            left = get_left(node)
            rotated_left = rotate_left!(left)
            set_left!(node, rotated_left)
            return rotate_right!(node)
        end
    elseif subtree_balance < -1
        right_balance = compute_balance(get_right(node))
        if right_balance ≤ 0
            return rotate_left!(node)
        else
            right = get_right(node)
            rotated_right = rotate_right!(right)
            set_right!(node, rotated_right)
            return rotate_left!(node)
        end
    end
    return node
end

"""
    _minimum(node::Union{BalancedBSTNode,Nothing}) -> Union{BalancedBSTNode,Nothing}

Returns the node with the minimum key in the subtree rooted at `node`. If `node` is `nothing`, returns `nothing`.
"""
function _minimum(node::Union{BalancedBSTNode,Nothing})
    while !isnothing(node) && has_left(node)
        node = get_left(node)
    end
    return node
end