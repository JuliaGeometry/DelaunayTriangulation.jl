function construct_polygon_hierarchy(points; IntegerType=Int)
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

function construct_polygon_hierarchy(points, boundary_nodes; IntegerType=Int)
    hierarchy = PolygonHierarchy{IntegerType}()
    return construct_polygon_hierarchy!(hierarchy, points, boundary_nodes)
end
construct_polygon_hierarchy(points, ::Nothing; IntegerType=Int) = construct_polygon_hierarchy(points; IntegerType)
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

function increase_depth!(tree::PolygonTree) # increase height of tree and all of its children
    set_height!(tree, get_height(tree) + 1)
    for child in get_children(tree)
        increase_depth!(child)
    end
    return nothing
end

function reorder_hierarchy!(hierarchy::PolygonHierarchy{I}, points, boundary_nodes, new_tree::PolygonTree) where {I}
    orig_curve_index = get_index(new_tree)
    for (curve_index, tree) in get_trees(hierarchy) # curve_index should never be orig_curve_index
        curve_nodes = get_boundary_nodes(boundary_nodes, curve_index)
        section_nodes = get_boundary_nodes(curve_nodes, 1)
        u = get_boundary_nodes(section_nodes, 1)
        p = get_point(points, u)
        is_in = is_in_tree(hierarchy, points, boundary_nodes, new_tree, p)
        if is_in
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

function find_tree(hierarchy::PolygonHierarchy, points, boundary_nodes, tree::PolygonTree, p)
    found = tree
    for sub_tree in get_children(tree)
        is_in = is_in_tree(hierarchy, points, boundary_nodes, sub_tree, p)
        is_in && return find_tree(hierarchy, points, boundary_nodes, sub_tree, p)
    end
    return found
end

function expand_bounds!(hierarchy::PolygonHierarchy, perc=0.1)
    bboxes = get_bounding_boxes(hierarchy)
    for (i, bbox) in enumerate(bboxes)
        bboxes[i] = expand(bbox, perc)
    end
    return hierarchy
end

function construct_polygon_hierarchy(points, boundary_nodes, boundary_curves; IntegerType=Int, n=4096)
    new_points, new_boundary_nodes = polygonise(points, boundary_nodes, boundary_curves; n)
    hierarchy = PolygonHierarchy{IntegerType}()
    return construct_polygon_hierarchy!(hierarchy, new_points, new_boundary_nodes)
end
construct_polygon_hierarchy(points, boundary_nodes, ::Tuple{}; IntegerType=Int, n=4096) = construct_polygon_hierarchy(points, boundary_nodes; IntegerType)
