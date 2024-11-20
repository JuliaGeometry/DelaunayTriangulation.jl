struct PolygonHierarchy{I}
    polygon_orientations::BitVector
    bounding_boxes::Vector{BoundingBox}
    trees::Dict{I,PolygonTree{I}}
    reorder_cache::Vector{PolygonTree{I}}
end
PolygonHierarchy{I}() where {I} = PolygonHierarchy{I}(BitVector(), BoundingBox[], Dict{I,PolygonTree{I}}(), PolygonTree{I}[])
function Base.deepcopy_internal(hierarchy::PolygonHierarchy{I}, dict::IdDict) where {I}
    haskey(dict, hierarchy) && return dict[hierarchy]
    copy_polygon_orientations = Base.deepcopy_internal(get_polygon_orientations(hierarchy), dict)
    copy_bounding_boxes = Base.deepcopy_internal(get_bounding_boxes(hierarchy), dict)
    copy_trees = Dict{I, PolygonTree{I}}()
    for (key, tree) in get_trees(hierarchy)
        copy_trees[key] = Base.deepcopy_internal(tree, dict)
    end
    copy_reorder_cache = Vector{PolygonTree{I}}()
    for tree in get_reorder_cache(hierarchy)
        push!(copy_reorder_cache, Base.deepcopy_internal(tree, dict))
    end
    copy_hierarchy = PolygonHierarchy(
        copy_polygon_orientations,
        copy_bounding_boxes,
        copy_trees,
        copy_reorder_cache
    )
    dict[hierarchy] = copy_hierarchy
    return copy_hierarchy
end
function Base.copy(hierarchy::PolygonHierarchy{I}) where {I}
    return PolygonHierarchy(
        copy(get_polygon_orientations(hierarchy)),
        copy(get_bounding_boxes(hierarchy)),
        copy(get_trees(hierarchy)),
        copy(get_reorder_cache(hierarchy))
    )
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

get_polygon_orientations(hierarchy::PolygonHierarchy) = hierarchy.polygon_orientations

get_polygon_orientation(hierarchy::PolygonHierarchy, index) = get_polygon_orientations(hierarchy)[index]

get_bounding_boxes(hierarchy::PolygonHierarchy) = hierarchy.bounding_boxes

get_bounding_box(hierarchy::PolygonHierarchy, index) = get_bounding_boxes(hierarchy)[index]

get_exterior_curve_indices(hierarchy::PolygonHierarchy) = keys(get_trees(hierarchy))

function get_positive_curve_indices(hierarchy::PolygonHierarchy)
    orientations = get_polygon_orientations(hierarchy)
    return (index for (index, orientation) in enumerate(orientations) if orientation)
end

get_trees(hierarchy::PolygonHierarchy) = hierarchy.trees

get_tree(hierarchy::PolygonHierarchy, index) = get_trees(hierarchy)[index]

set_tree!(hierarchy::PolygonHierarchy, index, tree) = get_trees(hierarchy)[index] = tree

delete_tree!(hierarchy::PolygonHierarchy, index) = delete!(get_trees(hierarchy), index)

get_reorder_cache(hierarchy::PolygonHierarchy) = hierarchy.reorder_cache

function set_orientation!(hierarchy::PolygonHierarchy, index, orientation)
    orientations = get_polygon_orientations(hierarchy)
    if length(orientations) < index
        resize!(orientations, index)
    end
    orientations[index] = orientation
    return hierarchy
end

function set_bounding_box!(hierarchy::PolygonHierarchy, index, bounding_box)
    bounding_boxes = get_bounding_boxes(hierarchy)
    if length(bounding_boxes) < index
        resize!(bounding_boxes, index)
    end
    bounding_boxes[index] = bounding_box
    return hierarchy
end