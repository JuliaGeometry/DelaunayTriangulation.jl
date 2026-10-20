# Comparisons of the package's R-tree with SpatialIndexing.jl.

export get_dt_rectangles, si_diametral_bounding_box

function Base.:(==)(rect1::SI.Rect, rect2::DT.BoundingBox)
    ab = rect1.low
    cd = rect1.high
    ab′ = DT.get_bl_corner(rect2)
    cd′ = DT.get_tr_corner(rect2)
    return (ab === ab′) && (cd === cd′)
end

Base.:(==)(rect2::DT.BoundingBox, rect1::SI.Rect) = rect1 == rect2

function Base.:(==)(rect1::SI.SpatialElem, rect2::DT.DiametralBoundingBox)
    id = rect1.val
    id′ = DT.get_edge(rect2)
    return (id == id′) && rect1.mbr == DT.get_bounding_box(rect2)
end

Base.:(==)(rect2::DT.DiametralBoundingBox, rect1::SI.SpatialElem) = rect1 == rect2

function Base.:(==)(tree1::SI.RTree, tree2::DT.RTree)
    si_rects, si_els = get_si_rectangles(tree1)
    dt_rects, dt_els = get_dt_rectangles(tree2)
    return (si_rects == dt_rects) && (si_els == dt_els)
end

Base.:(==)(tree2::DT.RTree, tree1::SI.RTree) = tree1 == tree2

Base.:(==)(int::DT.RTreeIntersectionIterator, int_si::SI.RTreeRegionQueryIterator) = (isempty(int) && isempty(int_si)) || collect(int) == collect(int_si)

Base.:(==)(int_si::SI.RTreeRegionQueryIterator, int::DT.RTreeIntersectionIterator) = int == int_si

function get_si_rectangles(tree::SI.RTree)
    rects = Any[]
    els = Any[]
    return get_si_rectangles!(rects, els, tree.root)
end

function get_dt_rectangles(tree::DT.RTree)
    rects = DT.BoundingBox[]
    els = DT.DiametralBoundingBox[]
    return get_dt_rectangles!(rects, els, DT.get_root(tree))
end

function get_si_rectangles!(rects, els, node::SI.Branch)
    rect = node.mbr
    push!(rects, rect)
    for child in node.children
        get_si_rectangles!(rects, els, child)
    end
    return identity.(rects), identity.(els)
end

function get_dt_rectangles!(rects, els, node::DT.Branch)
    rect = DT.get_bounding_box(node)
    push!(rects, rect)
    for child in DT.get_children(node)
        get_dt_rectangles!(rects, els, child)
    end
    return identity.(rects), identity.(els)
end

function get_si_rectangles!(rects, els, node::SI.Leaf)
    rect = node.mbr
    push!(rects, rect)
    for child in node.children
        push!(els, child)
    end
    return identity.(rects), identity.(els)
end

function get_dt_rectangles!(rects, els, node::DT.Leaf)
    rect = DT.get_bounding_box(node)
    push!(rects, rect)
    for child in DT.get_children(node)
        push!(els, child)
    end
    return identity.(rects), identity.(els)
end

function si_diametral_bounding_box(points, i, j)
    bbox = DT.bounding_box(points, i, j)
    bbox = SI.Rect(DT.get_bl_corner(bbox.bounding_box), DT.get_tr_corner(bbox.bounding_box))
    return bbox
end
