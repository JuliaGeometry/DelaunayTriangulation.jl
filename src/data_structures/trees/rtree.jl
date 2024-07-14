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
    abstract type AbstractBoundingShape

Abstract type for representing a bounding box.
"""
abstract type AbstractBoundingShape end

"""
    BoundingInterval <: AbstractBoundingShape

Type for representing a bounding interval `[a, b]`.

# Fields 
- `a::Float64`: The left endpoint.
- `b::Float64`: The right endpoint.
"""
struct BoundingInterval <: AbstractBoundingShape
    a::Float64
    b::Float64
end
Base.show(io::IO, I::BoundingInterval) = print(io, "[$(I.a) .. $(I.b)]")

"""
    InvalidBoundingInterval

A constant for representing an invalid interva[`BoundingInterval`](@ref), i.e. an interval with `NaN` endpoints.
"""
const InvalidBoundingInterval = BoundingInterval(NaN, NaN)

"""
    length(I::BoundingInterval) -> Float64 

Returns the length of the interval `I`.
"""
Base.length(I::BoundingInterval) = I.b - I.a

"""
    isempty(I::BoundingInterval) -> Bool

Returns `true` if `I` is empty, i.e. if `I.a` or `I.b` is `NaN` or if `length(I) < 0`.
"""
Base.isempty(I::BoundingInterval) = isnan(I.a) || isnan(I.b) || length(I) < 0

"""
    midpoint(I::BoundingInterval) -> Float64

Returns the midpoint of `I`.
"""
midpoint(I::BoundingInterval) = midpoint(I.a, I.b)

"""
    intersect(I::BoundingInterval, J::BoundingInterval) -> BoundingInterval
    I::BoundingInterval ∩ J::BoundingInterval -> BoundingInterval

Returns the intersection of `I` and `J`. If the intersection is empty, returns [`InvalidBoundingInterval`](@ref).
"""
function Base.:(∩)(I::BoundingInterval, J::BoundingInterval)
    a′ = max(I.a, J.a)
    b′ = min(I.b, J.b)
    if a′ ≤ b′
        return BoundingInterval(a′, b′)
    else
        return InvalidBoundingInterval
    end
end

"""
    union(I::BoundingInterval, J::BoundingInterval) -> BoundingInterval
    I::BoundingInterval ∪ J::BoundingInterval -> BoundingInterval

Returns the union of `I` and `J`, combining their bounds; i.e. the smallest interval that contains both `I` and `J`.
"""
function Base.:(∪)(I::BoundingInterval, J::BoundingInterval)
    a′ = min(I.a, J.a)
    b′ = max(I.b, J.b)
    return BoundingInterval(a′, b′)
end

"""
    in(a::Float64, I::BoundingInterval) -> Bool
    a::Float64 ∈ I::BoundingInterval -> Bool

Tests whether `a` is in `I`.
"""
Base.in(a::Float64, I::BoundingInterval) = I.a ≤ a ≤ I.b

"""
    in(I::BoundingInterval, J::BoundingInterval) -> Bool
    I::BoundingInterval ∈ J::BoundingInterval -> Bool

Tests whether the interval `I` is in the interval `J`.
"""
Base.in(I::BoundingInterval, J::BoundingInterval) = I.a ∈ J && I.b ∈ J

"""
    BoundingBox <: AbstractBoundingShape

Type for representing an axis-aligned bounding box, represented as a pair of interval `I` and `J` so that the bounding box is `I × J`.

# Fields 
- `x::BoundingInterval`: The interval for the x-axis.
- `y::BoundingInterval`: The interval for the y-axis.

# Constructors 
    
    BoundingBox(x::BoundingInterval, y::BoundingInterval)
    BoundingBox(a::Float64, b::Float64, c::Float64, d::Float64) = BoundingBox(BoundingInterval(a, b), BoundingInterval(c, d))
    BoundingBox(p::NTuple{2,<:Number}) = BoundingBox(p[1], p[1], p[2], p[2])
"""
struct BoundingBox <: AbstractBoundingShape
    x::BoundingInterval
    y::BoundingInterval
end
Base.show(io::IO, r::BoundingBox) = print(io, "[$(r.x.a), $(r.x.b)] × [$(r.y.a), $(r.y.b)]")

"""
    InvalidBoundingBox

A constant for representing an invalid rectangle, i.e. a rectangle with `NaN` endpoints.
"""
const InvalidBoundingBox = BoundingBox(InvalidBoundingInterval, InvalidBoundingInterval)
BoundingBox(a, b, c, d) = BoundingBox(BoundingInterval(a, b), BoundingInterval(c, d))
BoundingBox(p::NTuple{2,<:Number}) = BoundingBox(getx(p), getx(p), gety(p), gety(p))

"""
    hspan(r::BoundingBox) -> Float64

Returns the horizontal span of `r`, i.e. `length(r.x)`.
"""
hspan(r::BoundingBox) = length(r.x)

"""
    vspan(r::BoundingBox) -> Float64

Returns the vertical span of `r`, i.e. `length(r.y)`.
"""
vspan(r::BoundingBox) = length(r.y)

"""
    get_area(r::BoundingBox) -> Float64

Returns the area of `r`, i.e. `hspan(r) * vspan(r)`.
"""
get_area(r::BoundingBox) = hspan(r) * vspan(r)

"""
    isempty(r::BoundingBox) -> Bool

Returns `true` if `r` is empty, i.e. if `r.x` or `r.y` is empty.
"""
Base.isempty(r::BoundingBox) = isempty(r.x) || isempty(r.y)

"""
    midpoint(r::BoundingBox) -> NTuple{2,Float64}

Returns the center of `r`.
"""
midpoint(r::BoundingBox) = (midpoint(r.x), midpoint(r.y))

"""
    intersect(r1::BoundingBox, r2::BoundingBox) -> BoundingBox
    r1::BoundingBox ∩ r2::BoundingBox -> BoundingBox

Returns the intersection of `r1` and `r2`. If the intersection is empty, returns [`InvalidBoundingBox`](@ref).
"""
function Base.:(∩)(r1::BoundingBox, r2::BoundingBox)
    I = r1.x ∩ r2.x
    J = r1.y ∩ r2.y
    if isempty(I) || isempty(J)
        return InvalidBoundingBox
    else
        return BoundingBox(I, J)
    end
end

"""
    union(r1::BoundingBox, r2::BoundingBox) -> BoundingBox
    r1::BoundingBox ∪ r2::BoundingBox -> BoundingBox

Returns the union of `r1` and `r2`, i.e. the smallest bounding box that contains both `r1` and `r2`.
"""
function Base.:(∪)(r1::BoundingBox, r2::BoundingBox)
    I = r1.x ∪ r2.x
    J = r1.y ∪ r2.y
    return BoundingBox(I, J)
end

"""
    in(r1::BoundingBox, r2::BoundingBox) -> Bool
    r1::BoundingBox ∈ r2::BoundingBox -> Bool

Tests whether `r1` is in `r2`.
"""
Base.in(r1::BoundingBox, r2::BoundingBox) = (r1.x ∈ r2.x) && (r1.y ∈ r2.y)

"""
    in(p::NTuple{2,<:Number}, r::BoundingBox) -> Bool
    p::NTuple{2,<:Number} ∈ r::BoundingBox -> Bool

Tests whether `p` is in `r`.
"""
Base.in(p::NTuple{2,<:Number}, r::BoundingBox) = BoundingBox(p) ∈ r

"""
    is_touching(r1::BoundingBox, r2::BoundingBox) -> Bool

Tests whether `r1` and `r2` are touching, i.e. if they share a common boundary. This only considers interior touching.
"""
is_touching(r1::BoundingBox, r2::BoundingBox) = r1.x.a == r2.x.a || r1.x.b == r2.x.b || r1.y.a == r2.y.a || r1.y.b == r2.y.b # only considers interior touching

"""
    get_bl_corner(r::BoundingBox) -> NTuple{2,Float64}

Returns the bottom-left corner of `r`.
"""
get_bl_corner(r::BoundingBox) = (r.x.a, r.y.a)

"""
    get_tr_corner(r::BoundingBox) -> NTuple{2,Float64}

Returns the top-right corner of `r`.
"""
get_tr_corner(r::BoundingBox) = (r.x.b, r.y.b)

"""
    diametral_circle(p, q) -> (NTuple{2,Float64}, Float64)

Returns the circle with diameter `pq`.
"""
function diametral_circle(p, q)
    center = midpoint(p, q)
    radius = dist(p, q) / 2
    return center, radius
end

"""
    bounding_box(center, radius) -> BoundingBox

Returns the bounding box of the circle `(center, radius)`.
"""
function bounding_box(center::NTuple{2,<:Number}, radius::Number)
    cx, cy = getxy(center)
    return BoundingBox(cx - radius, cx + radius, cy - radius, cy + radius)
end

"""
    bounding_box(points) -> BoundingBox 

Gets the bounding box for a set of points.
"""
function bounding_box(points)
    xmin = Inf
    xmax = -Inf
    ymin = Inf
    ymax = -Inf
    for p in each_point(points)
        px, py = getxy(p)
        xmin = min(xmin, px)
        xmax = max(xmax, px)
        ymin = min(ymin, py)
        ymax = max(ymax, py)
    end
    return BoundingBox(xmin, xmax, ymin, ymax)
end

"""
    expand(box::BoundingBox, perc=0.10) -> BoundingBox

Expands the bounding box `box` by a factor `perc` in each direction.
"""
function expand(box::BoundingBox, perc=0.10)
    x = box.x
    y = box.y
    a, b = x.a, x.b
    c, d = y.a, y.b
    Δx = (b - a) * perc
    Δy = (d - c) * perc
    return BoundingBox(a - Δx, b + Δx, c - Δy, d + Δy)
end

"""
    bounding_box(p::NTuple, q::NTuple, r::NTuple) -> BoundingBox 

Returns the bounding box of the points `p`, `q` and `r`.
"""
function bounding_box(p::NTuple, q::NTuple, r::NTuple)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    xmin, _, xmax = min_med_max(px, qx, rx)
    ymin, _, ymax = min_med_max(py, qy, ry)
    return BoundingBox(xmin, xmax, ymin, ymax)
end

"""
    DiametralBoundingBox

Type for representing a bounding box generated from an edge's diametral circle.

# Fields
- `bounding_box::BoundingBox`: The bounding box.
- `edge::NTuple{2,Int}`: The generator edge.
"""
struct DiametralBoundingBox
    bounding_box::BoundingBox
    edge::NTuple{2,Int}
end

"""
    bounding_box(points, i, j) -> DiametralBoundingBox

Returns the bounding box of the diametral circle of the points `points[i]` and `points[j]` with generator edge `(i, j)`,
returned as an `DiametralBoundingBox`.
"""
function bounding_box(points, i, j)
    p, q = get_point(points, i, j)
    c, r = diametral_circle(p, q)
    bbox = bounding_box(c, r)
    edge = (Int(i), Int(j))
    return DiametralBoundingBox(bbox, edge)
end

Base.show(io::IO, id_bounding_box::DiametralBoundingBox) = print(io, "DiametralBoundingBox $(get_bounding_box(id_bounding_box)) with generator edge $(get_edge(id_bounding_box))")

"""
    get_bounding_box(id_bounding_box::DiametralBoundingBox) -> BoundingBox

Returns the bounding box of `id_bounding_box`.
"""
get_bounding_box(id_bounding_box::DiametralBoundingBox) = id_bounding_box.bounding_box

"""
    get_edge(id_bounding_box::DiametralBoundingBox) -> NTuple{2,Int}

Returns the generator edge of `id_bounding_box`.
"""
get_edge(id_bounding_box::DiametralBoundingBox) = id_bounding_box.edge

"""
    abstract type AbstractNode end

Abstract type for representing a node in an R-tree.
"""
abstract type AbstractNode end

"""
    set_child!(parent_node::AbstractNode, child_node, i::Integer)

Sets the `i`th child of `parent_node` to be `child_node`.
"""
set_child!(parent_node::AbstractNode, child_node, i::Integer) = parent_node.children[i] = child_node

"""
    set_parent!(child_node::AbstractNode, parent_node::AbstractNode)

Sets the parent of `child_node` to be `parent_node`.
"""
set_parent!(child_node::AbstractNode, parent) = child_node.parent = parent

"""
    get_bounding_box(node::AbstractNode) -> BoundingBox

Returns the bounding box of `node`.
"""
get_bounding_box(node::AbstractNode) = node.bounding_box

"""
    get_children(node::AbstractNode) -> Vector

Returns the children of `node`.
"""
get_children(node::AbstractNode) = node.children

"""
    get_child(node::AbstractNode, i::Integer) -> AbstractNode

Returns the `i`th child of `node`.
"""
get_child(node::AbstractNode, i::Integer) = node.children[i]

"""
    get_parent(node::AbstractNode) -> Union{Branch, Nothing}

Returns the parent of `node`.
"""
get_parent(node::AbstractNode) = node.parent

"""
    has_parent(node::AbstractNode) -> Bool

Returns `true` if `node` has a parent.
"""
has_parent(node::AbstractNode) = !isnothing(get_parent(node))

"""
    has_children(node::AbstractNode) -> Bool

Returns `true` if `node` has children.
"""
has_children(node::AbstractNode) = !isempty(get_children(node))

"""
    num_children(node::AbstractNode) -> Int

Returns the number of children of `node`.
"""
num_children(node::AbstractNode) = length(get_children(node))

"""
    get_child_type(node::AbstractNode) -> Union{Type{Leaf}, Type{Branch}}

Returns the type of the children of `node`.
"""
get_child_type(node::AbstractNode) = eltype(get_children(node))

"""
    add_child!(node::AbstractNode, child)

Adds `child` to `node`, i.e. appends `child` to the children of `node` via `push!`.
"""
add_child!(node::AbstractNode, child) = push!(get_children(node), child)

"""
    set_bounding_box!(node::AbstractNode, bounding_box::BoundingBox)

Sets the bounding box of `node` to be `bounding_box`.
"""
set_bounding_box!(node::AbstractNode, bounding_box::BoundingBox) = node.bounding_box = bounding_box

"""
    pop_child!(node::AbstractNode)

Removes the last child of `node` via `pop!`.
"""
pop_child!(node::AbstractNode) = pop!(get_children(node))

"""
    find_position_in_parent(node::AbstractNode) -> Int 

Returns the position of `node` in its parent's children. If `node` has no parent, returns `$∅`.
"""
function find_position_in_parent(node::AbstractNode)
    if has_parent(node)
        return findfirst(==(node), get_children(get_parent(node)))::Int
    else
        return ∅
    end
end

"""
    mutable struct Leaf <: AbstractNode

Type for representing a leaf node in an R-tree.

!!! danger "Type parametrisation"

    Technically, this type should be referred to by `Leaf{Branch}`. Due to a lack of support for mutually recursive types or 
    forward declarations, we have a parametric type in this struct's definition since Branch is not yet defined. In particular, 
    `Leaf` is not a concrete type, whereas `Leaf{Branch}` is.

# Fields
- `parent::Union{Branch, Nothing}`: The parent of the leaf node.
- `bounding_box::BoundingBox`: The bounding box of the leaf node.
- `children::Vector{DiametralBoundingBox}`: The children of the leaf node.

# Constructor 
    
    Leaf(parent::Union{Branch,Nothing}=nothing) = Leaf{Branch}(parent, InvalidBoundingBox, DiametralBoundingBox[])
"""
mutable struct Leaf{Branch} <: AbstractNode
    parent::Union{Branch,Nothing}
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

"""
    mutable struct Branch <: AbstractNode

Type for representing a branch node in an R-tree.

# Fields 
- `parent::Union{Branch, Nothing}`: The parent of the branch node.
- `bounding_box::BoundingBox`: The bounding box of the branch node.
- `children::Union{Vector{Branch},Vector{Leaf{Branch}}}`: The children of the branch node.
- `level::Int`: The level of the branch node.

# Constructor

    Branch(parent::Union{Branch,Nothing}=nothing, ::Type{C}=Branch) where {C<:AbstractNode} = new(parent, InvalidBoundingBox, C[], 1)
"""
mutable struct Branch <: AbstractNode
    parent::Union{Branch,Nothing}
    bounding_box::BoundingBox
    @const children::Union{Vector{Branch},Vector{Leaf{Branch}}} # if we do e.g. Branch{C}, it makes resolving some of the other types a bit difficult, especially Leaf{Branch}.
    level::Int
end
Branch(parent::Union{Branch,Nothing}=nothing, ::Type{C}=Branch) where {C<:AbstractNode} = Branch(parent, InvalidBoundingBox, C[], 1)
Leaf(parent::Union{Branch,Nothing}=nothing) = Leaf{Branch}(parent, InvalidBoundingBox, DiametralBoundingBox[])

function Base.:(==)(branch1::Branch, branch2::Branch)
    xor(isnothing(branch1), isnothing(branch1)) && return false # if we test get_parent(branch1) ≠ get_parent(branch2), then we get a StackOverflowError
    get_bounding_box(branch1) ≠ get_bounding_box(branch2) && return false
    get_children(branch1) ≠ get_children(branch2) && return false
    get_level(branch1) ≠ get_level(branch2) && return false
    return true
end

function Base.show(io::IO, ::MIME"text/plain", leaf::Leaf)
    parent = get_parent(leaf)
    nchildren = num_children(leaf)
    if has_parent(leaf)
        print(io, "Leaf with parent $(parent isa Leaf ? "Leaf" : "Branch") with $nchildren children")
    else
        print(io, "Leaf with $nchildren children")
    end
    return io
end
function Base.show(io::IO, ::MIME"text/plain", branch::Branch)
    parent = get_parent(branch)
    nchildren = num_children(branch)
    level = get_level(branch)
    if has_parent(branch)
        print(io, "Branch at level $level with parent $(parent isa Leaf ? "Leaf" : "Branch") with $nchildren children")
    else
        print(io, "Branch at level $level and $nchildren children")
    end
    return io
end
Base.show(io::IO, node::AbstractNode) = Base.show(io, MIME"text/plain"(), node)

@doc """
    get_level(node::AbstractNode) -> Int 

Returns the level of `node`. If `node` is a leaf, returns `1`.
"""
get_level

get_level(leaf::Leaf) = 1
get_level(branch::Branch) = branch.level

"""
    set_level!(branch::Branch, level::Integer)

Sets the level of `branch` to be `level`.
"""
set_level!(branch::Branch, level::Integer) = branch.level = level

"""
    setindex!(branch::Branch, child, i::Integer)
    branch[i] = child

Sets the `i`th child of `branch` to be `child`, also updating the parent of `child` to be `branch`.
"""
function Base.setindex!(branch::Branch, child, i::Integer)
    set_child!(branch, child, i)
    set_parent!(child, branch)
    return child
end

"""
    NodeCache{Node,Child}

Type for representing a cache of nodes whose children are of type `Child`. This is used for caching nodes that are detached from the R-tree, e.g. when a node is split.

# Fields 
- `cache::Vector{Node}`: The cache of nodes.
- `size_limit::Int`: The maximum number of nodes that can be cached.

# Constructor 
    
    NodeCache{Node,Child}(size_limit::Int) where {Node,Child} = new{Node,Child}(Node[], size_limit)
"""
struct NodeCache{Node,Child} # similar to why we use TriangulationCache. Think of it like implementing some CapacityVector that fails to push if it's full.
    cache::Vector{Node}
    size_limit::Int
    function NodeCache{Node,Child}(size_limit::Int) where {Node,Child}
        cache = Node[]
        sizehint!(cache, size_limit)
        return new{Node,Child}(cache, size_limit)
    end
end

"""
    BranchCache

Type for representing a cache of branch nodes.
"""
const BranchCache = NodeCache{Branch,Branch}

"""
    TwigCache

Type for representing a cache of twig nodes, i.e. branch nodes at level 2.
"""
const TwigCache = NodeCache{Branch,Leaf{Branch}}

"""
    LeafCache

Type for representing a cache of leaf nodes.
"""
const LeafCache = NodeCache{Leaf{Branch},DiametralBoundingBox}

"""
    length(cache::NodeCache) -> Int

Returns the number of nodes in `cache`.
"""
Base.length(cache::NodeCache) = length(cache.cache)

"""
    isempty(cache::NodeCache) -> Bool

Returns `true` if `cache` is empty.
"""
Base.isempty(cache::NodeCache) = isempty(cache.cache)

"""
    pop!(cache::NodeCache) -> Node

Removes and returns the last node in `cache`.
"""
Base.pop!(cache::NodeCache) = pop!(cache.cache)

"""
    push!(cache::NodeCache, node)
"""
Base.push!(cache::NodeCache, node) = push!(cache.cache, node)

"""
    get_size_limit(cache::NodeCache) -> Int

Returns the size limit of `cache`.
"""
get_size_limit(cache::NodeCache) = cache.size_limit

spawn_node(::BranchCache) = Branch()
spawn_node(::TwigCache) = Branch(nothing, Leaf{Branch})
spawn_node(::LeafCache) = Leaf()

"""
    spawn_node!(cache::NodeCache{Node}) where {Node} -> Node

Returns a node from `cache`. If `cache` is empty, returns a new node.
"""
spawn_node!(cache::NodeCache) = isempty(cache) ? spawn_node(cache) : pop!(cache)

"""
    is_full(cache::NodeCache) -> Bool

Returns `true` if `cache` is full, i.e. if `length(cache) ≥ get_size_limit(cache)`.
"""
is_full(cache::NodeCache) = length(cache) ≥ get_size_limit(cache)

"""
    cache_node!(cache::NodeCache, node)

Caches `node` in `cache` if `cache` is not full. Otherwise, does nothing.
"""
function cache_node!(cache::NodeCache, node)
    is_full(cache) || push!(cache, node)
    return cache
end

"""
    RTreeIntersectionCache

Type for representing a cache used for identifying intersections in an R-tree.

# Fields
- `node_indices::Vector{Int}`: A cache of indices used for identifying intersections.
- `need_tests::BitVector`: A `BitVector` cache for keeping track of which indices in `node_indices` need to be tested for intersections. 
"""
struct RTreeIntersectionCache
    node_indices::Vector{Int}
    need_tests::BitVector
end
RTreeIntersectionCache() = RTreeIntersectionCache(Vector{Int}(), BitVector())
Base.sizehint!(cache::RTreeIntersectionCache, n) = (sizehint!(get_node_indices(cache), n); sizehint!(get_need_tests(cache), n))

"""
    get_node_indices(cache::RTreeIntersectionCache) -> Vector{Int}

Returns the node indices of `cache`.
"""
get_node_indices(cache::RTreeIntersectionCache) = cache.node_indices

"""
    get_need_tests(cache::RTreeIntersectionCache) -> BitVector

Returns the `need_tests` cache of `tree`.
"""
get_need_tests(cache::RTreeIntersectionCache) = cache.need_tests

"""
    mutable struct RTree 

Type for representing an R-tree with linear splitting.

# Fields
- `root::Union{Branch,Leaf{Branch}}`: The root of the R-tree.
- `num_elements::Int`: The number of elements in the R-tree.
- `branch_cache::BranchCache`: The cache of branch nodes.
- `twig_cache::TwigCache`: The cache of twig nodes.
- `leaf_cache::LeafCache`: The cache of leaf nodes.
- `fill_factor::Float64`: The fill factor of the R-tree, i.e. the percentage of a node's capacity that should be filled after splitting.
- `free_cache::BitVector`: A bit vector for keeping track of which indices in `detached_cache` are free.
- `detached_cache::Vector{Union{Branch,Leaf{Branch}}}`: A cache of detached nodes, i.e. nodes that have been split from the R-tree. This is used for deleting nodes.
- `intersection_cache::NTuple{2,IntersectionCache}`: Cache used for identifying intersections. Each element of the `Tuple` is its own cache, allowing for up to 
    two intersection queries to be performed simultaneously. Note that this makes the R-tree non-thread-safe, and even non-safe when considering three or more 
    intersection queries simultaneously.

# Constructor 
    
        RTree(; size_limit=100, fill_factor=0.7)

The `size_limit` is the node capacity. All node types have the same capacity.
"""
mutable struct RTree # linear
    root::Union{Branch,Leaf{Branch}}
    num_elements::Int
    @const branch_cache::BranchCache 
    @const twig_cache::TwigCache 
    @const leaf_cache::LeafCache 
    @const fill_factor::Float64 
    @const free_cache::BitVector 
    @const detached_cache::Vector{Union{Branch,Leaf{Branch}}} 
    @const intersection_cache::NTuple{2,RTreeIntersectionCache}
    function RTree(; size_limit=100, fill_factor=0.7) # https://en.wikipedia.org/wiki/R-tree: "however best performance has been experienced with a minimum fill of 30%–40%)
        branch_cache = BranchCache(size_limit)
        twig_cache = TwigCache(size_limit)
        leaf_cache = LeafCache(size_limit)
        root = spawn_node!(leaf_cache)
        num_elements = 0
        free_cache = BitVector()
        sizehint!(free_cache, size_limit)
        detached_cache = Vector{Union{Branch,Leaf{Branch}}}()
        sizehint!(detached_cache, size_limit)
        cache1, cache2 = RTreeIntersectionCache(), RTreeIntersectionCache()
        sizehint!(cache1, ceil(Int, log2(size_limit)))
        sizehint!(cache2, ceil(Int, log2(size_limit)))
        return new(
            root,
            num_elements,
            branch_cache,
            twig_cache,
            leaf_cache,
            fill_factor,
            free_cache,
            detached_cache,
            (cache1, cache2)
        )
    end
end
function Base.:(==)(tree1::RTree, tree2::RTree)
    num_elements(tree1) ≠ num_elements(tree2) && return false
    get_root(tree1) ≠ get_root(tree2) && return false
    return true
end
function Base.show(io::IO, ::MIME"text/plain", tree::RTree)
    println(io, "RTree with height ", get_height(tree))
    println(io, "   Root: ", get_root(tree))
    return print(io, "   Num. elements: ", num_elements(tree))
end

"""
    get_root(tree::RTree) -> Union{Branch,Leaf{Branch}}

Returns the root of `tree`.
"""
get_root(tree::RTree) = tree.root

"""
    get_branch_cache(tree::RTree) -> BranchCache

Returns the branch cache of `tree`.
"""
get_branch_cache(tree::RTree) = tree.branch_cache

"""
    get_twig_cache(tree::RTree) -> TwigCache

Returns the twig cache of `tree`.
"""
get_twig_cache(tree::RTree) = tree.twig_cache

"""
    get_leaf_cache(tree::RTree) -> LeafCache

Returns the leaf cache of `tree`.
"""
get_leaf_cache(tree::RTree) = tree.leaf_cache

"""
    get_fill_factor(tree::RTree) -> Float64

Returns the fill factor of `tree`.
"""
get_fill_factor(tree::RTree) = tree.fill_factor

"""
    get_height(tree::RTree) -> Int

Returns the height of `tree`.
"""
get_height(tree::RTree) = get_level(get_root(tree))

"""
    get_bounding_box(tree::RTree) -> BoundingBox

Returns the bounding box of `tree`.
"""
get_bounding_box(tree::RTree) = get_bounding_box(get_root(tree))

"""
    get_size_limit(tree::RTree) -> Int

Returns the size limit of `tree`.
"""
get_size_limit(tree::RTree) = get_size_limit(get_branch_cache(tree)) # we assume that the caches all have the same size limit

"""
    is_full(node::AbstractNode, tree::RTree) -> Bool

Returns `true` if `node` is full, i.e. if `num_children(node) ≥ get_size_limit(tree)`.
"""
is_full(node::AbstractNode, tree::RTree) = num_children(node) ≥ get_size_limit(tree)

"""
    increment_num_elements!(tree::RTree)

Increments the number of elements in `tree` by 1.
"""
increment_num_elements!(tree::RTree) = tree.num_elements += 1

"""
    decrement_num_elements!(tree::RTree)

Decrements the number of elements in `tree` by 1.
"""
decrement_num_elements!(tree::RTree) = tree.num_elements -= 1

"""
    is_root(node::AbstractNode, tree::RTree) -> Bool

Returns `true` if `node` is the root of `tree`.
"""
is_root(node::AbstractNode, tree::RTree) = node === get_root(tree)

"""
    set_root!(tree::RTree, node::AbstractNode)

Sets the root of `tree` to be `node`.
"""
set_root!(tree::RTree, node::AbstractNode) = tree.root = node

"""
    num_elements(tree::RTree) -> Int

Returns the number of elements in `tree`.
"""
num_elements(tree::RTree) = tree.num_elements

"""
    get_min_nodes(tree::RTree) -> Int

Returns the minimum number of nodes that a node in `tree` can have.
"""
get_min_nodes(tree::RTree) = ceil(Int, get_fill_factor(tree) * get_size_limit(tree))

"""
    get_free_cache(tree::RTree) -> BitVector

Returns the free cache of `tree`.
"""
get_free_cache(tree::RTree) = tree.free_cache

"""
    get_detached_cache(tree::RTree) -> Vector{Union{Branch,Leaf{Branch}}}

Returns the detached cache of `tree`.
"""
get_detached_cache(tree::RTree) = tree.detached_cache

"""
    get_intersection_cache(tree::RTree) -> NTuple{2,RTreeIntersectionCache}

Returns the intersection cache of `tree`.
"""
get_intersection_cache(tree::RTree) = tree.intersection_cache

"""
    spawn_leaf!(tree::RTree, bounding_box::BoundingBox) -> Leaf{Branch}

Returns a new leaf node with bounding box `bounding_box` from `tree`.
"""
function spawn_leaf!(tree::RTree, bounding_box::BoundingBox)
    leaf = spawn_node!(get_leaf_cache(tree))
    set_bounding_box!(leaf, bounding_box)
    return leaf
end

"""
    spawn_branch!(tree::RTree, bounding_box::BoundingBox, level) -> Branch

Returns a new branch node with bounding box `bounding_box` and level `level` from `tree`.
"""
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

"""
    spawn_node!(tree::RTree, ::Type{N}, [bounding_box::BoundingBox], level) where {N} -> N

Returns a new node of type `N` with bounding box `bounding_box` and level `level` from `tree`. If 
`bounding_box` is not provided, it is replaced with [`InvalidBoundingBox`](@ref).
"""
function spawn_node!(tree::RTree, ::Type{N}, bounding_box::BoundingBox, level) where {N}
    if N <: Branch
        return spawn_branch!(tree, bounding_box, level)
    else
        return spawn_leaf!(tree, bounding_box)
    end
end
spawn_node!(tree::RTree, ::Type{N}, level) where {N} = spawn_node!(tree, N, InvalidBoundingBox, level)

@doc """
    cache_node!(tree::RTree, node::AbstractNode)

Caches `node` in into `tree`'s node caches.
"""
cache_node!

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

"""
    RTreeIntersectionIterator

Type for representing an iterator over the elements in an R-tree that intersect with a bounding box.

# Fields 
- `tree::RTree`: The R-tree.
- `bounding_box::BoundingBox`: The bounding box to test for intersections with.
- `cache::RTreeIntersectionCache`: The cache used for identifying intersections.    
"""
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

"""
    RTreeIntersectionIteratorState

The state of an [`RTreeIntersectionIterator`](@ref).

# Fields 
- `leaf::Leaf{Branch}`: The current leaf node.
- `node_indices::Vector{Int}`: The indices of the current node at each level.
- `need_tests::BitVector`: A `BitVector` cache for keeping track of which indices in `node_indices` need to be tested for intersections.
"""
struct RTreeIntersectionIteratorState
    leaf::Leaf{Branch}
    node_indices::Vector{Int}
    need_tests::BitVector
end

"""
    QueryResult 

An enum type for representing the result of an intersection query, with instances:

- `Contains`: The bounding box contains the element.
- `Intersects`: The bounding box intersects the element.
- `Outside`: The bounding box is outside the element.
"""
EnumX.@enumx QueryResult Contains Intersects Outside
const QR = QueryResult

"""
    get_state(state::RTreeIntersectionIteratorState) -> DiametralBoundingBox

Returns the current element in `state`.
"""
get_state(state::RTreeIntersectionIteratorState) = get_child(state.leaf, state.node_indices[begin])

"""
    BoundaryRTree{P}

Type for representing an R-tree of a boundary with an associated point set. The rectangular elements of the R-tree 
are the bounding box of the diametral circles of the boundary edges.

# Fields
- `tree::RTree`: The R-tree.
- `points::P`: The point set.

# Constructors 

    BoundaryRTree(points)
"""
struct BoundaryRTree{P}
    tree::RTree
    points::P
end
BoundaryRTree(points) = BoundaryRTree(RTree(), points)

function Base.:(==)(tree1::BoundaryRTree, tree2::BoundaryRTree)
    tree1.points ≠ tree2.points && return false
    tree1.tree ≠ tree2.tree && return false
    return true
end

"""
    bounding_box(tree::BoundaryRTree, i, j) -> DiametralBoundingBox

Returns the bounding box of the diametral circle of the edge between `i` and `j` in `tree`.
"""
function bounding_box(tree::BoundaryRTree, i, j)
    points = tree.points
    i′, j′ = minmax(i, j) # so that we can recompute safely and get the exact same box
    return bounding_box(points, i′, j′)
end

"""
    insert!(tree::BoundaryRTree, i, j) -> Bool

Inserts the diametral circle of the edge between `i` and `j` into `tree`. Returns `true` if the `tree`'s bounding boxes had to be adjusted and `false` otherwise.
"""
function Base.insert!(tree::BoundaryRTree, i, j)
    bbox = bounding_box(tree, i, j)
    return insert!(tree.tree, bbox)
end

"""
    delete!(tree::BoundaryRTree, i, j)

Deletes the bounding box of the diametral circle of the edge between `i` and `j` in `tree`.
"""
function Base.delete!(tree::BoundaryRTree, i, j)
    bbox = bounding_box(tree, i, j)
    return delete!(tree.tree, bbox)
end

"""
    split_edge!(tree::BoundaryRTree, i, j, r)

Splits the diametral bounding box associated with `(i, j)` into two new boxes associated 
with the diametral circles of `(i, r)` and `(j, r)`.
"""
function split_edge!(tree::BoundaryRTree, i, j, r)
    delete!(tree, i, j)
    insert!(tree, i, r)
    insert!(tree, r, j)
    return tree
end

"""
    get_intersections(tree::BoundaryRTree, i, j; cache_id=1) -> RTreeIntersectionIterator

Returns an [`RTreeIntersectionIterator`](@ref) over the elements in `tree` that intersect with the diametral circle of the edge between `i` and `j`.
`cache_id` must be `1` or `2`, and determines what cache to use for the intersection query.
"""
function get_intersections(tree::BoundaryRTree, i, j; cache_id=1)
    bbox = bounding_box(tree, i, j)
    return get_intersections(tree.tree, get_bounding_box(bbox); cache_id)
end

"""
    get_intersections(tree::BoundaryRTree, i, j, k; cache_id=1) -> RTreeIntersectionIterator

Returns an [`RTreeIntersectionIterator`](@ref) over the elements in `tree` that intersect with the bounding box of the triangle `(i, j, k)`.
`cache_id` must be `1` or `2`, and determines what cache to use for the intersection query.
"""
function get_intersections(tree::BoundaryRTree, i, j, k; cache_id=1)
    points = tree.points
    p, q, r = get_point(points, i, j, k)
    bbox = bounding_box(p, q, r)
    return get_intersections(tree.tree, bbox; cache_id)
end

"""
    get_intersections(tree::BoundaryRTree, i; cache_id=1) -> RTreeIntersectionIterator

Returns an [`RTreeIntersectionIterator`](@ref) over the elements in `tree` that intersect with the `i`th vertex.
`cache_id` must be `1` or `2`, and determines what cache to use for the intersection query.
"""
function get_intersections(tree::BoundaryRTree, i; cache_id=1)
    p = get_point(tree.points, i)
    return get_intersections(tree.tree, p; cache_id)
end

"""
    get_intersections(tree::BoundaryRTree, bbox::BoundingBox; cache_id=1) -> RTreeIntersectionIterator

Returns an [`RTreeIntersectionIterator`](@ref) over the elements in `tree` that intersect with `bbox`.
`cache_id` must be `1` or `2`, and determines what cache to use for the intersection query.
"""
function get_intersections(tree::BoundaryRTree, bbox::BoundingBox; cache_id=1)
    return get_intersections(tree.tree, bbox; cache_id)
end