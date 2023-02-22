"""
    brute_force_search(T, r, pts, boundary_map)

Given a collection of triangles `T`, a point `r`, a collection of points 
`pts`, and a boundary map handling the mapping of boundary indices from 
[`construct_boundary_map`](@ref), returns the triangle in `T` containing `r` 
by searching over all triangles.
"""
function brute_force_search(T, r, pts, boundary_map)
    for V in each_triangle(T)
        cert = point_position_relative_to_triangle(V, r, pts, boundary_map)
        if !is_outside(cert)
            return V
        end
    end
    throw("Failed to find the point $(get_point(pts, r)).")
end
