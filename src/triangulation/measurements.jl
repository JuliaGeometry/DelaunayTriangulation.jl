function get_area(tri::Triangulation)
    F = number_type(tri)
    A = zero(F)
    for T in each_solid_triangle(tri)
        u, v, w = triangle_vertices(T)
        p, q, r = get_point(tri, u, v, w)
        A += triangle_area(p, q, r)
    end
    return A
end

function edge_length(tri::Triangulation, u, v)
    p, q = get_point(tri, u, v)
    return dist(p, q)
end
edge_length(tri::Triangulation, e) = edge_length(tri, initial(e), terminal(e))

function edge_length_sqr(tri::Triangulation, u, v)
    p, q = get_point(tri, u, v)
    return dist_sqr(p, q)
end
edge_length_sqr(tri::Triangulation, e) = edge_length_sqr(tri, initial(e), terminal(e))

function midpoint(tri::Triangulation, u, v)
    p, q = get_point(tri, u, v)
    return midpoint(p, q)
end
midpoint(tri::Triangulation, e) = midpoint(tri, initial(e), terminal(e))

function dist(tri::Triangulation, p)
    points = get_points(tri)
    if has_boundary_nodes(tri)
        return distance_to_polygon(p, points, get_boundary_nodes(tri))
    else
        return distance_to_polygon(p, points, get_convex_hull_vertices(tri))
    end
end
