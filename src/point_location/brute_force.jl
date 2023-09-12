"""
    brute_force_search(tri::Triangulation, q; itr = each_triangle(tri)

Returns the triangle in `tri` containing `q` using brute force search.
"""
function brute_force_search(tri::Triangulation, q; itr = each_triangle(tri))
    for V in itr
        cert = point_position_relative_to_triangle(tri, V, q)
        !is_outside(cert) && return V 
    end
    throw("Failed to find the point $(get_point(tri, q)).")
end