abstract type AbstractNode end

set_child!(parent_node::AbstractNode, child_node, i::Integer) = parent_node.children[i] = child_node

set_parent!(child_node::AbstractNode, parent) = child_node.parent = parent

get_bounding_box(node::AbstractNode) = node.bounding_box

get_children(node::AbstractNode) = node.children

get_child(node::AbstractNode, i::Integer) = node.children[i]

get_parent(node::AbstractNode) = node.parent

has_parent(node::AbstractNode) = !isnothing(get_parent(node))

has_children(node::AbstractNode) = !isempty(get_children(node))

num_children(node::AbstractNode) = length(get_children(node))

get_child_type(node::AbstractNode) = eltype(get_children(node))

add_child!(node::AbstractNode, child) = push!(get_children(node), child)

set_bounding_box!(node::AbstractNode, bounding_box::BoundingBox) = node.bounding_box = bounding_box

pop_child!(node::AbstractNode) = pop!(get_children(node))

function find_position_in_parent(node::AbstractNode)
    if has_parent(node)
        return findfirst(==(node), get_children(get_parent(node)))::Int
    else
        return ∅
    end
end

mutable struct Leaf{Branch} <: AbstractNode
    parent::Union{Branch, Nothing}
    bounding_box::BoundingBox
    @const children::Vector{DiametralBoundingBox}
    Leaf(parent::Branch, bounding_box, children) where {Branch} = new{Branch}(parent, bounding_box, children)
    Leaf{Branch}(parent::Branch, bounding_box, children) where {Branch} = new{Branch}(parent, bounding_box, children)
    Leaf{Branch}(::Nothing, bounding_box, children) where {Branch} = new{Branch}(nothing, bounding_box, children)
    # need to separate out the constructors to avoid unbound type argument issues from Aqua
end
function Base.:(==)(leaf1::Leaf, leaf2::Leaf)
    xor(isnothing(leaf1), isnothing(leaf2)) && return false # if we test get_parent(leaf1) ≠ get_parent(leaf2), then we get a StackOverflowError
    get_bounding_box(leaf1) ≠ get_bounding_box(leaf2) && return false
    get_children(leaf1) ≠ get_children(leaf2) && return false
    return true
end

function Base.deepcopy_internal(leaf::Leaf{Branch}, dict::IdDict) where {Branch}
    haskey(dict, leaf) && return dict[leaf]
    new_parent = has_parent(leaf) ? Base.deepcopy_internal(get_parent(leaf), dict) : nothing
    bbox = get_bounding_box(leaf)
    cchildren = copy(get_children(leaf))
    new_leaf = Leaf{Branch}(new_parent, bbox, cchildren)
    dict[leaf] = new_leaf
    return new_leaf
end

function Base.copy(leaf::Leaf{B}) where {B}
    parent = get_parent(leaf)
    cparent = isnothing(parent) ? nothing : copy(parent)
    bbox = get_bounding_box(leaf)
    cchildren = copy(get_children(leaf))
    return Leaf{B}(cparent, bbox, cchildren)
end

mutable struct Branch <: AbstractNode
    parent::Union{Branch, Nothing}
    bounding_box::BoundingBox
    @const children::Union{Vector{Branch}, Vector{Leaf{Branch}}} # if we do e.g. Branch{C}, it makes resolving some of the other types a bit difficult, especially Leaf{Branch}.
    level::Int
end
Branch(parent::Union{Branch, Nothing} = nothing, ::Type{C} = Branch) where {C <: AbstractNode} = Branch(parent, InvalidBoundingBox, C[], 1)
Leaf(parent::Union{Branch, Nothing} = nothing) = Leaf{Branch}(parent, InvalidBoundingBox, DiametralBoundingBox[])

function Base.:(==)(branch1::Branch, branch2::Branch)
    xor(isnothing(branch1), isnothing(branch1)) && return false # if we test get_parent(branch1) ≠ get_parent(branch2), then we get a StackOverflowError
    get_bounding_box(branch1) ≠ get_bounding_box(branch2) && return false
    get_children(branch1) ≠ get_children(branch2) && return false
    get_level(branch1) ≠ get_level(branch2) && return false
    return true
end

function Base.deepcopy_internal(branch::Branch, dict::IdDict)
    haskey(dict, branch) && return dict[branch]
    new_parent = has_parent(branch) ? Base.deepcopy_internal(get_parent(branch), dict) : nothing
    bbox = get_bounding_box(branch)
    cchildren = copy(get_children(branch))
    new_branch = Branch(new_parent, bbox, cchildren, get_level(branch))
    dict[branch] = new_branch
    return new_branch
end

function Base.copy(branch::Branch)
    parent = get_parent(branch)
    cparent = isnothing(parent) ? nothing : copy(parent)
    bbox = get_bounding_box(branch)
    cchildren = copy(get_children(branch))
    return Branch(cparent, bbox, cchildren, get_level(branch))
end

get_level(leaf::Leaf) = 1
get_level(branch::Branch) = branch.level

set_level!(branch::Branch, level::Integer) = branch.level = level

function Base.setindex!(branch::Branch, child, i::Integer)
    set_child!(branch, child, i)
    set_parent!(child, branch)
    return child
end

