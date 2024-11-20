struct BoundaryRTree{P}
    tree::RTree
    points::P
end
BoundaryRTree(points) = BoundaryRTree(RTree(), points)

function Base.copy(tree::BoundaryRTree)
    ctree = copy(tree.tree)
    cpoints = copy(tree.points)
    return BoundaryRTree(ctree, cpoints)
end

function Base.:(==)(tree1::BoundaryRTree, tree2::BoundaryRTree)
    tree1.points ≠ tree2.points && return false
    tree1.tree ≠ tree2.tree && return false
    return true
end

function bounding_box(tree::BoundaryRTree, i, j)
    points = tree.points
    i′, j′ = minmax(i, j) # so that we can recompute safely and get the exact same box
    return bounding_box(points, i′, j′)
end

function Base.insert!(tree::BoundaryRTree, i, j)
    bbox = bounding_box(tree, i, j)
    return insert!(tree.tree, bbox)
end

function Base.delete!(tree::BoundaryRTree, i, j)
    bbox = bounding_box(tree, i, j)
    return delete!(tree.tree, bbox)
end

function split_edge!(tree::BoundaryRTree, i, j, r)
    delete!(tree, i, j)
    insert!(tree, i, r)
    insert!(tree, r, j)
    return tree
end

function get_intersections(tree::BoundaryRTree, i, j; cache_id = 1)
    bbox = bounding_box(tree, i, j)
    return get_intersections(tree.tree, get_bounding_box(bbox); cache_id)
end

function get_intersections(tree::BoundaryRTree, i, j, k; cache_id = 1)
    points = tree.points
    p, q, r = get_point(points, i, j, k)
    bbox = bounding_box(p, q, r)
    return get_intersections(tree.tree, bbox; cache_id)
end

function get_intersections(tree::BoundaryRTree, i; cache_id = 1)
    p = get_point(tree.points, i)
    return get_intersections(tree.tree, p; cache_id)
end

function get_intersections(tree::BoundaryRTree, bbox::BoundingBox; cache_id = 1)
    return get_intersections(tree.tree, bbox; cache_id)
end
