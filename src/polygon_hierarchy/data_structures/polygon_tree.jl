mutable struct PolygonTree{I}
    parent::Union{Nothing,PolygonTree{I}}
    @const children::Set{PolygonTree{I}}
    @const index::I
    height::Int
end
PolygonTree{I}(parent::Union{Nothing,PolygonTree{I}}, index, height) where {I} = PolygonTree{I}(parent, Set{PolygonTree{I}}(), index, height)
function hash_tree(tree::PolygonTree)
    height = get_height(tree)
    index = get_index(tree)
    parent_index = has_parent(tree) ? get_index(get_parent(tree)) : 0
    h = hash((parent_index, index, height))
    children = collect(get_children(tree))
    sort!(children, by=get_index)
    for child in children
        h = hash((h, hash_tree(child)))
    end
    return h
end
function Base.:(==)(tree1::PolygonTree, tree2::PolygonTree)
    children1 = get_children(tree1)
    children2 = get_children(tree2)
    length(children1) ≠ length(children2) && return false
    index1 = get_index(tree1)
    index2 = get_index(tree2)
    index1 ≠ index2 && return false
    height1 = get_height(tree1)
    height2 = get_height(tree2)
    height1 ≠ height2 && return false
    return hash_tree(tree1) == hash_tree(tree2)
end
function Base.deepcopy_internal(tree::PolygonTree{I}, dict::IdDict) where {I}
    haskey(dict, tree) && return dict[tree]
    copy_tree = PolygonTree{I}(nothing, get_index(tree), get_height(tree))
    dict[tree] = copy_tree
    new_parent = Base.deepcopy_internal(get_parent(tree), dict)
    if !isnothing(new_parent)
        set_parent!(copy_tree, new_parent)
    end
    for child in get_children(tree)
        add_child!(copy_tree, Base.deepcopy_internal(child, dict))
    end
    return copy_tree
end
function Base.copy(tree::PolygonTree{I}) where {I}
    parent = get_parent(tree)
    new_parent = isnothing(parent) ? nothing : copy(parent)
    return PolygonTree{I}(new_parent, copy(get_children(tree)), copy(get_index(tree)), copy(get_height(tree)))
end

get_parent(tree::PolygonTree) = tree.parent

get_children(tree::PolygonTree) = tree.children

get_index(tree::PolygonTree) = tree.index

get_height(tree::PolygonTree) = tree.height

has_parent(tree::PolygonTree) = !isnothing(get_parent(tree))

num_children(tree::PolygonTree) = length(get_children(tree))

has_children(tree::PolygonTree) = !isempty(get_children(tree))

add_child!(tree::PolygonTree, child::PolygonTree) = push!(get_children(tree), child)

set_parent!(tree::PolygonTree, parent::PolygonTree) = tree.parent = parent

set_height!(tree::PolygonTree, height) = tree.height = height

delete_child!(tree::PolygonTree, child) = delete!(get_children(tree), child)