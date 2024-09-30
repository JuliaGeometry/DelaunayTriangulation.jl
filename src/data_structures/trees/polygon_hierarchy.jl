"""
    mutable struct PolygonTree{I}

A tree structure used to define a polygon hierarchy.

# Fields
- `parent::Union{Nothing,PolygonTree{I}}`: The parent of the tree. If `nothing`, then the tree is the root.
- `children::Set{PolygonTree{I}}`: The children of the tree.
- `index::I`: The index of the tree. This is the index associated with the polygon.
- `height::Int`: The height of the tree. This is the number of polygons in the tree that `index` is inside of. The root has height `0`.

# Constructor

    PolygonTree{I}(parent::Union{Nothing,PolygonTree{I}}, index, height) where {I}

Constructs a [`PolygonTree`](@ref) with `parent`, `index`, and `height`, and no children.
"""
mutable struct PolygonTree{I}
    parent::Union{Nothing, PolygonTree{I}}
    @const children::Set{PolygonTree{I}}
    @const index::I
    height::Int
end
PolygonTree{I}(parent::Union{Nothing, PolygonTree{I}}, index, height) where {I} = PolygonTree{I}(parent, Set{PolygonTree{I}}(), index, height)
function hash_tree(tree::PolygonTree)
    height = get_height(tree)
    index = get_index(tree)
    parent_index = has_parent(tree) ? get_index(get_parent(tree)) : 0
    h = hash((parent_index, index, height))
    children = collect(get_children(tree))
    sort!(children, by = get_index)
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
@static if VERSION ≥ v"1.10"
    function Base.deepcopy(tree::PolygonTree) # without this definition, deepcopy would occassionally segfault 
        @warn "Deepcopy on PolygonTrees is not currently supported. Returning the tree without copying." maxlog = 1 # https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/129
        return tree
        #=
        parent = get_parent(tree)
        children = get_children(tree)
        index = get_index(tree)
        height = get_height(tree)
        new_parent = isnothing(parent) ? nothing : deepcopy(parent)
        new_children = deepcopy(children)
        new_index = deepcopy(index)
        return PolygonTree(new_parent, new_children, new_index, height)
        =# # Why does the above still segfault sometimes?
        # TODO: Fix
    end
end
function Base.show(io::IO, ::MIME"text/plain", tree::PolygonTree)
    print(io, "PolygonTree at height $(get_height(tree)) with index $(get_index(tree)) and $(length(get_children(tree))) children")
end
Base.show(io::IO, tree::PolygonTree) = Base.show(io, MIME"text/plain"(), tree)

"""
    get_parent(tree::PolygonTree) -> Union{Nothing,PolygonTree}

Returns the parent of `tree`.
"""
get_parent(tree::PolygonTree) = tree.parent

"""
    get_children(tree::PolygonTree) -> Set{PolygonTree}

Returns the children of `tree`.
"""
get_children(tree::PolygonTree) = tree.children

"""
    get_index(tree::PolygonTree{I}) -> I

Returns the index of `tree`.
"""
get_index(tree::PolygonTree) = tree.index

"""
    get_height(tree::PolygonTree) -> Int

Returns the height of `tree`.
"""
get_height(tree::PolygonTree) = tree.height

"""
    has_parent(tree::PolygonTree) -> Bool

Returns `true` if `tree` has a parent, and `false` otherwise.
"""
has_parent(tree::PolygonTree) = !isnothing(get_parent(tree))

"""
    has_children(tree::PolygonTree) -> Bool

Returns `true` if `tree` has children, and `false` otherwise.
"""
num_children(tree::PolygonTree) = length(get_children(tree))

"""
    has_children(tree::PolygonTree) -> Bool

Returns `true` if `tree` has children, and `false` otherwise.
"""
has_children(tree::PolygonTree) = !isempty(get_children(tree))

"""
    add_child!(tree::PolygonTree, child::PolygonTree)

Adds `child` to `tree`.
"""
add_child!(tree::PolygonTree, child::PolygonTree) = push!(get_children(tree), child)

"""
    set_parent!(tree::PolygonTree, parent::PolygonTree)

Sets the parent of `tree` to `parent`.
"""
set_parent!(tree::PolygonTree, parent::PolygonTree) = tree.parent = parent

"""
    set_height!(tree::PolygonTree, height::Int)

Sets the height of `tree` to `height`.
"""
set_height!(tree::PolygonTree, height) = tree.height = height

"""
    delete_child!(tree::PolygonTree, child::PolygonTree)

Deletes `child` from `tree`.
"""
delete_child!(tree::PolygonTree, child) = delete!(get_children(tree), child)

"""
    PolygonHierarchy{I}

Struct used to define a polygon hierarchy. The hierarchy is represented as a forest of [`PolygonTree`](@ref)s.

!!! danger "Overlapping polygons"

    The polygons must not intersect any other polygon's boundaries.

# Fields
- `polygon_orientations::BitVector`: A `BitVector` of length `n` where `n` is the number of polygons in the hierarchy. The `i`th entry is `true` if the `i`th polygon is positively oriented, and `false` otherwise.
- `bounding_boxes::Vector{BoundingBox}`: A `Vector` of [`BoundingBox`](@ref)s of length `n` where `n` is the number of polygons in the hierarchy. The `i`th entry is the [`BoundingBox`](@ref) of the `i`th polygon.
- `trees::Dict{I,PolygonTree{I}}`: A `Dict` mapping the index of a polygon to its [`PolygonTree`](@ref). The keys of `trees` are the roots of each individual tree, i.e. the outer-most polygons.
- `reorder_cache::Vector{PolygonTree{I}}`: A `Vector used for caching trees to be deleted in [`reorder_subtree!`](@ref).

!!! note "One-based indexing"

    Note that the vector definitions for `polygon_orientations` and `bounding_boxes` are treating the curves with the assumption that they are 
    enumerated in the order 1, 2, 3, ....

# Constructor

    PolygonHierarchy{I}() where {I}

Constructs a [`PolygonHierarchy`](@ref) with no polygons.
"""
struct PolygonHierarchy{I}
    polygon_orientations::BitVector
    bounding_boxes::Vector{BoundingBox}
    trees::Dict{I, PolygonTree{I}}
    reorder_cache::Vector{PolygonTree{I}}
end
PolygonHierarchy{I}() where {I} = PolygonHierarchy{I}(BitVector(), BoundingBox[], Dict{I, PolygonTree{I}}(), PolygonTree{I}[])
@static if VERSION ≥ v"1.10"
    function Base.deepcopy(hierarchy::PolygonHierarchy{I}) where {I} # without this definition, deepcopy would occassionally segfault 
        polygon_orientations = get_polygon_orientations(hierarchy)
        bounding_boxes = get_bounding_boxes(hierarchy)
        trees = get_trees(hierarchy)
        reorder_cache = get_reorder_cache(hierarchy)
        new_polygon_orientations = copy(polygon_orientations)
        new_bounding_boxes = copy(bounding_boxes)
        new_trees = Dict{I, PolygonTree{I}}()
        for (index, tree) in trees
            new_trees[index] = deepcopy(tree)
        end
        new_reorder_cache = similar(reorder_cache)
        for (index, tree) in enumerate(reorder_cache)
            new_reorder_cache[index] = deepcopy(tree)
        end
        return PolygonHierarchy{I}(new_polygon_orientations, new_bounding_boxes, new_trees, new_reorder_cache)
    end
end
function Base.empty!(hierarchy::PolygonHierarchy)
    polygon_orientations = get_polygon_orientations(hierarchy)
    bounding_boxes = get_bounding_boxes(hierarchy)
    trees = get_trees(hierarchy)
    reorder_cache = get_reorder_cache(hierarchy)
    empty!(polygon_orientations)
    empty!(bounding_boxes)
    empty!(trees)
    empty!(reorder_cache)
    return hierarchy
end
function Base.copyto!(hierarchy::PolygonHierarchy, full_hierarchy::PolygonHierarchy)
    polygon_orientations = get_polygon_orientations(hierarchy)
    full_polygon_orientations = get_polygon_orientations(full_hierarchy)
    bounding_boxes = get_bounding_boxes(hierarchy)
    full_bounding_boxes = get_bounding_boxes(full_hierarchy)
    trees = get_trees(hierarchy)
    full_trees = get_trees(full_hierarchy)
    reorder_cache = get_reorder_cache(hierarchy)
    full_reorder_cache = get_reorder_cache(full_hierarchy)
    resize!(polygon_orientations, length(full_polygon_orientations))
    resize!(bounding_boxes, length(full_bounding_boxes))
    empty!(trees)
    resize!(reorder_cache, length(full_reorder_cache))
    copyto!(polygon_orientations, full_polygon_orientations)
    copyto!(bounding_boxes, full_bounding_boxes)
    for (index, tree) in full_trees
        set_tree!(hierarchy, index, tree)
    end
    copyto!(reorder_cache, full_reorder_cache)
    return hierarchy
end
function Base.:(==)(hierarchy1::PolygonHierarchy, hierarchy2::PolygonHierarchy)
    polygon_orientations1 = get_polygon_orientations(hierarchy1)
    polygon_orientations2 = get_polygon_orientations(hierarchy2)
    polygon_orientations1 ≠ polygon_orientations2 && return false
    bounding_boxes1 = get_bounding_boxes(hierarchy1)
    bounding_boxes2 = get_bounding_boxes(hierarchy2)
    bounding_boxes1 ≠ bounding_boxes2 && return false
    trees1 = get_trees(hierarchy1)
    trees2 = get_trees(hierarchy2)
    length(trees1) ≠ length(trees2) && return false
    keys(trees1) ≠ keys(trees2) && return false
    for (index, tree1) in trees1
        tree2 = trees2[index]
        tree1 ≠ tree2 && return false
    end
    return true
end

"""
    get_polygon_orientations(hierarchy::PolygonHierarchy) -> BitVector

Returns the polygon orientations of `hierarchy`.
"""
get_polygon_orientations(hierarchy::PolygonHierarchy) = hierarchy.polygon_orientations

"""
    get_polygon_orientation(hierarchy::PolygonHierarchy, index) -> Bool

Returns the polygon orientation of the `index`th polygon in `hierarchy`.
"""
get_polygon_orientation(hierarchy::PolygonHierarchy, index) = get_polygon_orientations(hierarchy)[index]

"""
    get_bounding_boxes(hierarchy::PolygonHierarchy) -> Vector{BoundingBox}

Returns the bounding boxes of `hierarchy`.
"""
get_bounding_boxes(hierarchy::PolygonHierarchy) = hierarchy.bounding_boxes

"""
    get_bounding_box(hierarchy::PolygonHierarchy, index) -> BoundingBox

Returns the bounding box of the `index`th polygon in `hierarchy`.
"""
get_bounding_box(hierarchy::PolygonHierarchy, index) = get_bounding_boxes(hierarchy)[index]

"""
    get_exterior_curve_indices(hierarchy::PolygonHierarchy) -> KeySet

Returns the indices of the exterior curves of `hierarchy`.
"""
get_exterior_curve_indices(hierarchy::PolygonHierarchy) = keys(get_trees(hierarchy))

"""
    get_positive_curve_indices(hierarchy::PolygonHierarchy) -> Generator 

Returns the indices of the positively oriented curves of `hierarchy` as a generator, i.e. 
as a lazy result.
"""
function get_positive_curve_indices(hierarchy::PolygonHierarchy) 
    orientations = get_polygon_orientations(hierarchy)
    return (index for (index, orientation) in enumerate(orientations) if orientation)
end

"""
    get_trees(hierarchy::PolygonHierarchy) -> Dict{I,PolygonTree{I}}

Returns the trees of `hierarchy`, mapping the index of an exterior polygon to its [`PolygonTree`](@ref).
"""
get_trees(hierarchy::PolygonHierarchy) = hierarchy.trees

"""
    get_tree(hierarchy::PolygonHierarchy, index) -> PolygonTree{I}

Returns the [`PolygonTree`](@ref) of the `index`th polygon in `hierarchy`. The `index` must be associated with 
an exterior polygon.
"""
get_tree(hierarchy::PolygonHierarchy, index) = get_trees(hierarchy)[index]

"""
    set_tree!(hierarchy::PolygonHierarchy, index, tree)

Sets the [`PolygonTree`](@ref) of the `index`th polygon in `hierarchy` to `tree`, or adds it if it is not an existing key. 
The `index` must be associated with an exterior polygon.
"""
set_tree!(hierarchy::PolygonHierarchy, index, tree) = get_trees(hierarchy)[index] = tree

"""
    delete_tree!(hierarchy::PolygonHierarchy, index)

Deletes the [`PolygonTree`](@ref) of the `index`th polygon in `hierarchy`. The `index` must be associated with an exterior polygon.
"""
delete_tree!(hierarchy::PolygonHierarchy, index) = delete!(get_trees(hierarchy), index)

"""
    get_reorder_cache(hierarchy::PolygonHierarchy) -> Vector{PolygonTree{I}}

Returns the reorder cache of `hierarchy`.
"""
get_reorder_cache(hierarchy::PolygonHierarchy) = hierarchy.reorder_cache

"""
    set_orientation!(hierarchy::PolygonHierarchy, index, orientation)

Sets the polygon orientation of the `index`th polygon in `hierarchy` to `orientation`. If `index` is greater than the length of the
polygon orientations vector, the vector is resized.
"""
function set_orientation!(hierarchy::PolygonHierarchy, index, orientation)
    orientations = get_polygon_orientations(hierarchy)
    if length(orientations) < index
        resize!(orientations, index)
    end
    orientations[index] = orientation
    return hierarchy
end

"""
    set_bounding_box!(hierarchy::PolygonHierarchy, index, bounding_box)

Sets the bounding box of the `index`th polygon in `hierarchy` to `bounding_box`. If `index` is greater than the length of the
bounding boxes vector, the vector is resized.
"""
function set_bounding_box!(hierarchy::PolygonHierarchy, index, bounding_box)
    bounding_boxes = get_bounding_boxes(hierarchy)
    if length(bounding_boxes) < index
        resize!(bounding_boxes, index)
    end
    bounding_boxes[index] = bounding_box
    return hierarchy
end

"""
    construct_polygon_hierarchy(points; IntegerType=Int) -> PolygonHierarchy{IntegerType}

Returns a [`PolygonHierarchy`](@ref) defining the polygon hierarchy for a given set of `points`. This defines a hierarchy with a single polygon.
"""
function construct_polygon_hierarchy(points; IntegerType = Int)
    hierarchy = PolygonHierarchy{IntegerType}()
    return construct_polygon_hierarchy!(hierarchy, points)
end
function construct_polygon_hierarchy!(hierarchy::PolygonHierarchy{I}, points) where {I}
    empty!(hierarchy)
    set_orientation!(hierarchy, 1, true)
    tree = PolygonTree{I}(nothing, one(I), 0)
    set_tree!(hierarchy, one(I), tree)
    bbox = bounding_box(points)
    set_bounding_box!(hierarchy, one(I), bbox)
    return hierarchy
end

"""
    construct_polygon_hierarchy(points, boundary_nodes; IntegerType=Int) -> PolygonHierarchy{IntegerType}
    
Returns a [`PolygonHierarchy`](@ref) defining the polygon hierarchy for a given set of `boundary_nodes` that define a set of piecewise 
linear curves. 
"""
function construct_polygon_hierarchy(points, boundary_nodes; IntegerType = Int)
    hierarchy = PolygonHierarchy{IntegerType}()
    return construct_polygon_hierarchy!(hierarchy, points, boundary_nodes)
end
construct_polygon_hierarchy(points, ::Nothing; IntegerType = Int) = construct_polygon_hierarchy(points; IntegerType)
function construct_polygon_hierarchy!(hierarchy::PolygonHierarchy{I}, points, boundary_nodes) where {I}
    if !has_boundary_nodes(boundary_nodes)
        return construct_polygon_hierarchy!(hierarchy, points)
    end
    empty!(hierarchy)
    if has_multiple_curves(boundary_nodes)
        _construct_polygon_hierarchy_multiple_curves!(hierarchy, points, boundary_nodes)
    else
        _construct_polygon_hierarchy_single_curve!(hierarchy, points, boundary_nodes)
    end
    return hierarchy
end

function _construct_polygon_hierarchy_single_curve!(hierarchy::PolygonHierarchy{I}, points, boundary_nodes) where {I}
    area = polygon_features(points, boundary_nodes)[1]
    is_pos = area > 0.0
    xmin, xmax, ymin, ymax = polygon_bounds(points, boundary_nodes)
    bbox = BoundingBox(xmin, xmax, ymin, ymax)
    set_orientation!(hierarchy, one(I), is_pos)
    set_bounding_box!(hierarchy, one(I), bbox)
    tree = PolygonTree{I}(nothing, one(I), 0)
    set_tree!(hierarchy, one(I), tree)
    return nothing
end
function _construct_polygon_hierarchy_multiple_curves!(hierarchy::PolygonHierarchy{I}, points, boundary_nodes) where {I}
    curve_nodes = get_boundary_nodes(boundary_nodes, 1)
    _construct_polygon_hierarchy_single_curve!(hierarchy, points, curve_nodes)
    nc = num_curves(boundary_nodes)
    for curve_index in 2:nc
        curve_nodes = get_boundary_nodes(boundary_nodes, curve_index)
        area = polygon_features(points, curve_nodes)[1]
        is_pos = area > 0.0
        xmin, xmax, ymin, ymax = polygon_bounds(points, curve_nodes)
        bbox = BoundingBox(xmin, xmax, ymin, ymax)
        set_orientation!(hierarchy, curve_index, is_pos)
        set_bounding_box!(hierarchy, curve_index, bbox)
        section_nodes = get_boundary_nodes(curve_nodes, 1)
        u = get_boundary_nodes(section_nodes, 1)
        p = get_point(points, u)
        tree = find_tree(hierarchy, points, boundary_nodes, p)
        if !isnothing(tree)
            new_tree = PolygonTree{I}(tree, curve_index, get_height(tree) + 1)
            reorder_subtree!(hierarchy, points, boundary_nodes, tree, new_tree) # If the boundary is in the tree but not any of its children, then all of the tree's children now belong to this new boundary
        else
            new_tree = PolygonTree{I}(nothing, curve_index, 0)
            reorder_hierarchy!(hierarchy, points, boundary_nodes, new_tree)
        end
    end
    return nothing
end

"""
    reorder_subtree!(hierarchy::PolygonHierarchy, points, boundary_nodes, tree::PolygonTree, new_tree)

Given a `new_tree` contained inside `tree`, adds it into `hierarchy`. The children of `tree` are reordered if necessary, in case 
they are now contained inside `new_tree`.

# Arguments 
- `hierarchy::PolygonHierarchy`: The [`PolygonHierarchy`](@ref) to add `new_tree` to.
- `points`: The point set.
- `boundary_nodes`: The boundary nodes.
- `tree::PolygonTree`: The [`PolygonTree`](@ref) to add `new_tree` to.
- `new_tree::PolygonTree`: The [`PolygonTree`](@ref) to add to `hierarchy` and `tree`.
"""
function reorder_subtree!(hierarchy::PolygonHierarchy, points, boundary_nodes, tree::PolygonTree, new_tree)
    reorder_cache = get_reorder_cache(hierarchy)
    empty!(reorder_cache)
    for child in get_children(tree)
        index = get_index(child)
        curve_nodes = get_boundary_nodes(boundary_nodes, index)
        section_nodes = get_boundary_nodes(curve_nodes, 1)
        u = get_boundary_nodes(section_nodes, 1)
        p = get_point(points, u)
        is_in = is_in_tree(hierarchy, points, boundary_nodes, new_tree, p)
        if is_in
            set_parent!(child, new_tree)
            add_child!(new_tree, child)
            increase_depth!(child)
            push!(reorder_cache, child)
        end
    end
    for child in reorder_cache
        delete_child!(tree, child)
    end
    add_child!(tree, new_tree)
    return nothing
end

"""
    increase_depth!(tree::PolygonTree)

Increases the height of `tree` and all of its children by `1`.
"""
function increase_depth!(tree::PolygonTree) # increase height of tree and all of its children
    set_height!(tree, get_height(tree) + 1)
    for child in get_children(tree)
        increase_depth!(child)
    end
    return nothing
end

"""
    reorder_hierarchy!(hierarchy::PolygonHierarchy, points, boundary_nodes, new_tree::PolygonTree)

Given a `new_tree` that is not contained inside any other polygon in `hierarchy`, adds it into `hierarchy`. The existing trees are 
checked to see if they are contained inside `new_tree`, and if so, they are added as children of `new_tree` and removed from `hierarchy`.

# Arguments
- `hierarchy::PolygonHierarchy`: The [`PolygonHierarchy`](@ref) to add `new_tree` to.
- `points`: The point set.
- `boundary_nodes`: The boundary nodes.
- `new_tree::PolygonTree`: The [`PolygonTree`](@ref) to add to `hierarchy`.
"""
function reorder_hierarchy!(hierarchy::PolygonHierarchy{I}, points, boundary_nodes, new_tree::PolygonTree) where {I}
    orig_curve_index = get_index(new_tree)
    for (curve_index, tree) in get_trees(hierarchy) # curve_index should never be orig_curve_index
        curve_nodes = get_boundary_nodes(boundary_nodes, curve_index)
        section_nodes = get_boundary_nodes(curve_nodes, 1)
        u = get_boundary_nodes(section_nodes, 1)
        p = get_point(points, u)
        is_in = is_in_tree(hierarchy, points, boundary_nodes, new_tree, p)
        if is_in
            is_disjoint = false
            set_parent!(tree, new_tree)
            add_child!(new_tree, tree)
            increase_depth!(tree)
            index = get_index(tree)
            delete_tree!(hierarchy, index)
        end
    end
    set_tree!(hierarchy, orig_curve_index, new_tree)
    return nothing
end

"""
    is_in_tree(hierarchy::PolygonHierarchy, points, boundary_nodes, tree::PolygonTree, p) -> Bool

Tests if the point `p` is inside `tree`.

# Arguments
- `hierarchy::PolygonHierarchy`: The [`PolygonHierarchy`](@ref) containing `tree`.
- `points`: The point set.
- `boundary_nodes`: The boundary nodes.
- `tree::PolygonTree`: The [`PolygonTree`](@ref) to test `p` against.
- `p`: The point to test.

# Output
- `true` if `p` is inside `tree`, and `false` otherwise.
"""
function is_in_tree(hierarchy::PolygonHierarchy, points, boundary_nodes, tree::PolygonTree, p)
    index = get_index(tree)
    if has_multiple_sections(boundary_nodes)
        curve_nodes = get_boundary_nodes(boundary_nodes, index)
    else
        curve_nodes = boundary_nodes
    end
    bbox = get_bounding_box(hierarchy, index)
    if p ∈ bbox
        δ = distance_to_polygon(p, points, curve_nodes)
        return δ > 0
    else
        return false
    end
end

"""
    find_tree(hierarchy::PolygonHierarchy, points, boundary_nodes, p) -> Union{Nothing,PolygonTree}
    
Finds a tree in `hierarchy` containing `p`.

# Arguments
- `hierarchy::PolygonHierarchy`: The [`PolygonHierarchy`](@ref) to search.
- `points`: The point set.
- `boundary_nodes`: The boundary nodes.
- `p`: The point to test the trees of `hierarchy` against.

# Output
- `nothing` if `p` is not inside any tree in `hierarchy`, and the [`PolygonTree`](@ref) containing `p` otherwise.
"""
function find_tree(hierarchy::PolygonHierarchy, points, boundary_nodes, p)
    found = nothing
    for (_, tree) in get_trees(hierarchy)
        is_in = is_in_tree(hierarchy, points, boundary_nodes, tree, p)
        if is_in
            found = tree
            break
        end
    end
    isnothing(found) && return nothing
    return find_tree(hierarchy, points, boundary_nodes, found, p)
end

"""
    find_tree(hierarchy::PolygonHierarchy, points, boundary_nodes, tree::PolygonTree, p) -> PolygonTree

Finds a tree in `hierarchy` containing `p` that is a child of `tree`, assuming `p` is inside `tree`.

# Arguments
- `hierarchy::PolygonHierarchy`: The [`PolygonHierarchy`](@ref) to search.
- `points`: The point set.
- `boundary_nodes`: The boundary nodes.
- `tree::PolygonTree`: The [`PolygonTree`](@ref) to search, assuming `p` is inside `tree`.
- `p`: The point to test the children of `tree` against.

# Output
- `tree` if `p` is inside `tree` but none of its children, and a child containing `p` otherwise.
"""
function find_tree(hierarchy::PolygonHierarchy, points, boundary_nodes, tree::PolygonTree, p)
    found = tree
    for sub_tree in get_children(tree)
        is_in = is_in_tree(hierarchy, points, boundary_nodes, sub_tree, p)
        is_in && return find_tree(hierarchy, points, boundary_nodes, sub_tree, p)
    end
    return found
end

"""
    expand_bounds!(hierarchy::PolygonHierarchy, perc=0.10) -> PolygonHierarchy

Expands the bounding boxes of `hierarchy` by a factor of `perc` in each direction.
"""
function expand_bounds!(hierarchy::PolygonHierarchy, perc = 0.1)
    bboxes = get_bounding_boxes(hierarchy)
    for (i, bbox) in enumerate(bboxes)
        bboxes[i] = expand(bbox, perc)
    end
    return hierarchy
end

"""
    construct_polygon_hierarchy(points, boundary_nodes, boundary_curves; IntegerType=Int, n=4096) -> PolygonHierarchy{IntegerType}

Returns a [`PolygonHierarchy`](@ref) defining the polygon hierarchy for a given set of `boundary_nodes` that define a curve-bounded domain 
from the curves in `boundary_curves`. Uses [`polygonise`](@ref) to fill in the boundary curves.

# Arguments
- `points`: The point set.
- `boundary_nodes`: The boundary nodes. These should be the output from [`convert_boundary_curves!`](@ref).
- `boundary_curves`: The boundary curves. These should be the output from [`convert_boundary_curves!`](@ref).

# Keyword Arguments
- `IntegerType=Int`: The integer type to use for indexing the polygons.
- `n=4096`: The number of points to use for filling in the boundary curves in [`polygonise`](@ref).
"""
function construct_polygon_hierarchy(points, boundary_nodes, boundary_curves; IntegerType = Int, n = 4096)
    new_points, new_boundary_nodes = polygonise(points, boundary_nodes, boundary_curves; n)
    hierarchy = PolygonHierarchy{IntegerType}()
    return construct_polygon_hierarchy!(hierarchy, new_points, new_boundary_nodes)
end
construct_polygon_hierarchy(points, boundary_nodes, ::Tuple{}; IntegerType = Int, n = 4096) = construct_polygon_hierarchy(points, boundary_nodes; IntegerType)
