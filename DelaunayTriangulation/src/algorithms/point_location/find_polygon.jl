"""
    find_polygon(tri::Triangulation, q) -> Integer 

Given a point `q`, finds the index of the polygon in the triangulation `tri` that contains `q`. If 
`q` is on the boundary of the triangulation or outside the triangulation, the function returns `0`.

See also [`dist`](@ref) and [`distance_to_polygon`](@ref).
"""
function find_polygon(tri::Triangulation, q)
    points = get_points(tri)
    boundary_nodes = if has_boundary_nodes(tri)
        get_boundary_nodes(tri)
    else
        get_convex_hull_vertices(tri)
    end
    hierarchy = get_polygon_hierarchy(tri)
    tree = find_tree(hierarchy, points, boundary_nodes, q)
    return isnothing(tree) ? 0 : get_index(tree)
end

function find_polygon(hierarchy::PolygonHierarchy, points, boundary_nodes, q)
    tree = find_tree(hierarchy, points, boundary_nodes, q)
    return isnothing(tree) ? 0 : get_index(tree)
end
