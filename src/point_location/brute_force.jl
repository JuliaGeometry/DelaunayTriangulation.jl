"""
    brute_force_search(T, r, pts, representative_point_list, boundary_map)

Searches for the triangle in `T` containing the point `r` in `T` using brute force, simply 
searching over all triangles. Ghost triangles are handled via the `representative_point_list`
and `boundary_map` from [`construct_boudary_map`](@ref).
"""
function brute_force_search(T, r, pts, representative_point_list, boundary_map)
    for V in each_triangle(T)
        cert = point_position_relative_to_triangle(V, r, pts, representative_point_list, boundary_map)
        if !is_outside(cert)
            return V
        end
    end
    throw("Failed to find the point $(get_point(pts, r)).")
end
