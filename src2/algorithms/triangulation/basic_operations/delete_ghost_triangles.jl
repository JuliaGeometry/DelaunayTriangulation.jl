"""
    delete_ghost_triangles!(tri::Triangulation)

Deletes all the ghost triangles from `tri`.

!!! warning "Ghost vertices"

    Ghost vertices are still used in the `keys` of the [`Adjacent2Vertex`](@ref) 
    of `tri`, and are still present in the [`Graph`](@ref). If you want to delete the 
    ghost vertex `keys` from the [`Adjacent2Vertex`](@ref), you need to use 
    [`delete_adjacent2vertex!`](@ref). For deleting the ghost vertices from the 
    [`Graph`](@ref), you need [`delete_ghost_vertices_from_graph!`](@ref). Additionally, 
    edges in [`Adjacent`](@ref) can still map to ghost vertices. If you also want to delete 
    those, you need to filter through the `values` of the [`Adjacent`](@ref) map 
    that are ghost vertices, and use [`delete_adjacent!`](@ref).
"""
function delete_ghost_triangles!(tri::Triangulation)
    for g in each_ghost_vertex(tri)
        for uv in each_edge(get_adjacent2vertex(tri, g))
            @show 1
            u, v = edge_vertices(uv)
            delete_adjacent!(tri, v, g)
            delete_adjacent!(tri, g, u)
        end
    end
    return tri
end
