"""
    is_boundary_edge(tri::Triangulation, ij)
    is_boundary_edge(tri::Triangulation, i, j)

Returns `is_boundary_edge(ij, get_adjacent(tri))` or `is_boundary_edge(i, j, get_adjacent(tri))`, respectively,
testing if `ij = (i, j)` belongs to the boundary of the triangulation `tri`, i.e. `(i, j)` adjoins a ghost vertex.
"""
@inline is_boundary_edge(tri::Triangulation, ij) = is_boundary_edge(ij, get_adjacent(tri))
@inline is_boundary_edge(tri::Triangulation, i, j) = is_boundary_edge(i, j, get_adjacent(tri))

"""
    is_boundary_triangle(tri::Triangulation, i, j, k)
    is_boundary_triangle(tri::Triangulation, T)

Returns `is_boundary_triangle(i, j, k, get_adjacent(tri))` or `is_boundary_triangle(T, get_adjacent(tri))`,
testing if at least one of the edges `(u, v)` of `T = (i, j, k)` satisfies `is_boundary_edge(tri, u, v)`.
"""
@inline is_boundary_triangle(tri::Triangulation, i, j, k) = is_boundary_triangle(i, j, k, get_adjacent(tri))
@inline is_boundary_triangle(tri::Triangulation, T) = is_boundary_triangle(tri, geti(T), getj(T), getk(T))

"""
    triangle_orientation(tri::Triangulation, i, j, k)
    triangle_orientation(tri::Triangulation, T)

Computes the orientation of the triangle `T = (i, j, k)` in `tri`.
"""
@inline function triangle_orientation(tri::Triangulation, i, j, k)
    points = get_points(tri)
    representative_point_list = get_representative_point_list(tri)
    boundary_map = get_boundary_map(tri)
    return triangle_orientation(i, j, k, points, representative_point_list, boundary_map)
end
@inline triangle_orientation(tri::Triangulation, T) = triangle_orientation(tri, geti(T), getj(T), getk(T))

"""
    point_position_relative_to_circumcircle(tri::Triangulation, i, j, k, ℓ)
    point_position_relative_to_circumcircle(tri::Triangulation, T, ℓ)

Computes the position of the `ℓ`th point of `tri` relative to the circumcircle of the triangle `T = (i, j, k)` in `tri`.
"""
@inline function point_position_relative_to_circumcircle(tri::Triangulation, i, j, k, ℓ)
    points = get_points(tri)
    representative_point_list = get_representative_point_list(tri)
    boundary_map = get_boundary_map(tri)
    return point_position_relative_to_circumcircle(i, j, k, ℓ, points, representative_point_list, boundary_map)
end
@inline point_position_relative_to_circumcircle(tri::Triangulation, T, ℓ) = point_position_relative_to_circumcircle(tri, geti(T), getj(T), getk(T), ℓ)

"""
    point_position_relative_to_line(tri::Triangulation, i, j, u)

Computes the position of the `u`th point of `tri` relative to the line segment with indices `(i, j)` in `tri`.
"""
@inline function point_position_relative_to_line(tri::Triangulation, i, j, u)
    points = get_points(tri)
    representative_point_list = get_representative_point_list(tri)
    boundary_map = get_boundary_map(tri)
    return point_position_relative_to_line(i, j, u, points, representative_point_list, boundary_map)
end

"""
    point_closest_to_line(tri::Triangulation, i, j, u, v) 

Tests which of the points `u` or `v` is closest to the line segment with indices `(i, j)` in `tri`.
"""
@inline point_closest_to_line(tri::Triangulation, i, j, u, v) = point_closest_to_line(i, j, u, v, get_points(tri))

"""
    point_position_on_line_segment(tri::Triangulation, i, j, u)

Given vertices `i`, `j`, and `u` that are collinear, computes the position of `u` 
on the line segment `(i, j)`.
"""
@inline point_position_on_line_segment(tri::Triangulation, i, j, u) = point_position_on_line_segment(i, j, u, get_points(tri))

"""
    line_segment_intersection_type(tri::Triangulation, i, j, u, v)

Given two lines with indices `(i, j)` and `(u, v)`, respectively, in `tri`, computes the number of 
intersections.
"""
@inline line_segment_intersection_type(tri::Triangulation, u, v, i, j) = line_segment_intersection_type(u, v, i, j, get_points(tri))

"""
    point_position_relative_to_triangle(tri::Triangulation, i, j, k, u)
    point_position_relative_to_triangle(tri::Triangulation, T, u)

Computes the position of the `u`th point of `tri` relative to the triangle `T = (i, j, k)` in `tri`.
"""
@inline function point_position_relative_to_triangle(tri::Triangulation, i, j, k, u)
    points = get_points(tri)
    representative_point_list = get_representative_point_list(tri)
    boundary_map = get_boundary_map(tri)
    return point_position_relative_to_triangle(i, j, k, u, points, representative_point_list, boundary_map)
end
@inline point_position_relative_to_triangle(tri::Triangulation, T, u) = point_position_relative_to_triangle(tri, geti(T), getj(T), getk(T), u)

"""
    triangle_line_segment_intersection(tri::Triangulation, i, j, k, u, v)
    triangle_line_segment_intersection(tri::Triangulation, T, e)

Computes the type of the intersection of the triangle `T = (i, j, k)` in `tri` with the line segment `e = (u, v)`.
"""
@inline triangle_line_segment_intersection(tri::Triangulation, i, j, k, u, v) = triangle_line_segment_intersection(i, j, k, u, v, get_points(tri))
@inline triangle_line_segment_intersection(tri::Triangulation, T, e) = triangle_line_segment_intersection(tri, geti(T), getj(T), getk(T), initial(e), terminal(e))

"""
    is_outer_boundary_index(tri::Triangulation, i)

Returns `true` if `is_boundary_index(i)` and this boundary index `i` corresponds to the outermost boundary of `tri`. 
Returns `false` otherwise.
"""
@inline is_outer_boundary_index(tri::Triangulation, i) = is_outer_boundary_index(i, get_boundary_map(tri))

"""
    is_outer_boundary_node(tri::Triangulation, i)

Returns `true` if `i` is a boundary node belonging to the outermost boundary of `tri`.
Returns `false` otherwise.
"""
@inline is_outer_boundary_node(tri::Triangulation, i) = is_outer_boundary_node(i, get_graph(tri), get_boundary_index_ranges(tri)) # not the same as the above function, since i ≤ BoundaryIndex in the above

"""
    is_boundary_node(tri::Triangulation, i) 

Returns `true` if `i` is a node belonging to any boundary of `tri`.
"""
@inline is_boundary_node(tri::Triangulation, i) = is_boundary_node(i, get_graph(tri), get_boundary_index_ranges(tri))

"""
    edge_exists(tri::Triangulation, i, j)
    edge_exists(tri::Triangulation, ij)

Returns `true` if the edge `(i, j)` exists in `tri`, and `false` otherwise.
"""
@inline edge_exists(tri::Triangulation, i, j) = edge_exists(get_adjacent(tri, i, j))
@inline edge_exists(tri::Triangulation, ij) = edge_exists(get_adjacent(tri, ij))

"""
    has_ghost_triangles(tri::Triangulation)

Returns `true` if `tri` has ghost triangles, and `false` otherwise.
"""
@inline has_ghost_triangles(tri::Triangulation) = has_ghost_triangles(get_adjacent(tri), get_adjacent2vertex(tri))

"""
    has_boundary_nodes(tri::Triangulation)

Returns `true` if `tri` has boundary nodes, and `false` otherwise. 
"""
@inline has_boundary_nodes(tri::Triangulation) = has_multiple_segments(tri) || num_boundary_edges(get_boundary_nodes(tri)) ≠ 0

"""
    is_legal(tri::Triangulation, i, j)

Returns `true` if the edge `(i, j)` is legal in `tri`, and `false` otherwise. We also define 
constrained edges to be legal, as are ghost edges.
"""
function is_legal(tri::Triangulation, i, j)
    (contains_constrained_edge(tri, i, j) ||
     is_boundary_edge(tri, i, j) ||
     is_boundary_edge(tri, j, i) ||
     !edge_exists(tri, i, j) ||
     !edge_exists(tri, j, i) ||
     is_ghost_edge(i, j)) && return Cert.Legal
    k = get_adjacent(tri, i, j)
    ℓ = get_adjacent(tri, j, i)
    p, q, r, s = get_point(tri, i, j, k, ℓ)
    cert = is_legal(p, q, r, s)
    return cert
end

"""
    find_edge(tri::Triangulation, T, ℓ)

Given a point `ℓ` that is on an edge of the triangle `T` in `tri`, returns the 
edge that `ℓ` is on. 
"""
@inline find_edge(tri::Triangulation, T, ℓ) = find_edge(T, get_points(tri), ℓ)

"""
    is_constrained(tri::Triangulation, i, j)

Returns `true` if `tri` has any constrained edges, and `false` otherwise.
"""
is_constrained(tri::Triangulation) = !is_empty(get_all_constrained_edges(tri))
