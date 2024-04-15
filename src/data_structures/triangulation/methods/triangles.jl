"""
    num_triangles(tri::Triangulation) -> Integer

Returns the number of triangles in `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of these triangles will be ghost triangles.
"""
num_triangles(tri::Triangulation) = num_triangles(get_triangles(tri))

"""
    contains_triangle(tri::Triangulation, T) -> (Triangle, Bool)
    contains_triangle(tri::Triangulation, i, j, k) -> (Triangle, Bool)

Tests whether `tri` contains `T = (i, j, k)` up to rotation, returning 

- `V`: The rotated form of `T` that is in `tri`, or simply `T` if `T` is not in `tri`.
- `flag`: `true` if `T` is in `tri`, and `false` otherwise.
"""
contains_triangle(tri::Triangulation, T) = contains_triangle(T, get_triangles(tri))
contains_triangle(tri::Triangulation, i, j, k) = contains_triangle(i, j, k, get_triangles(tri))

"""
    construct_triangle(tri::Triangulation, i, j, k) -> Triangle

Returns a triangle in `tri` from the vertices `i`, `j`, and `k` such that the triangle is positively oriented.
"""
construct_positively_oriented_triangle(tri::Triangulation, i, j, k) = construct_positively_oriented_triangle(triangle_type(tri), i, j, k, get_points(tri))

"""
    num_ghost_triangles(tri::Triangulation) -> Integer

Returns the number of ghost triangles in `tri`.
"""
function num_ghost_triangles(tri::Triangulation)
    if has_ghost_triangles(tri)
        num_ghosts = 0
        all_ghosts = all_ghost_vertices(tri)
        for i in all_ghosts
            !has_vertex(tri, i) && continue # in case delete_holes=false when triangulate was called, there could be issues here
            S = get_adjacent2vertex(tri, i)
            num_ghosts += num_edges(S)
        end
    else
        num_ghosts = 0
    end
    return num_ghosts
end

"""
    num_solid_triangles(tri::Triangulation) -> Integer

Returns the number of solid triangles in `tri`.
"""
num_solid_triangles(tri::Triangulation) = num_triangles(tri) - num_ghost_triangles(tri)