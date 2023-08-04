"""
    get_graph(tri::Triangulation)

Returns `get_edges(get_graph(tri))`.
"""
@inline get_edges(tri::Triangulation) = get_edges(get_graph(tri))

"""
    get_vertices(tri::Triangulation)

Returns `get_vertices(get_graph(tri))`, giving all vertices present in the triangulation. 
This will include ghost vertices. 

See also [`each_vertex`](@ref), [`each_solid_vertex`](@ref) and [`each_ghost_vertex`](@ref).
"""
@inline get_vertices(tri::Triangulation) = get_vertices(get_graph(tri))

"""
    get_neighbours(tri::Triangulation)

Returns `get_neighbours(get_graph(tri))`, the set of neighbourhoods in the triangulation `tri`.
"""
@inline get_neighbours(tri::Triangulation) = get_neighbours(get_graph(tri))

"""
    get_neighbours(tri::Triangulation, u)

Returns `get_neighbours(get_graph(tri), u)`, the set of neighbours of `u`.
"""
@inline get_neighbours(tri::Triangulation, u) = get_neighbours(get_graph(tri), u)

"""
    add_vertex!(tri::Triangulation, u...)

Calls `add_vertex!(get_graph(tri), u...)`, adding the vertex `u` into `tri`.
"""
@inline add_vertex!(tri::Triangulation, u...) = add_vertex!(get_graph(tri), u...)

"""
    num_neighbours(tri::Triangulation, u)

Returns `num_neighbours(get_graph(tri), u)`, the number of neighbours of `u`.
"""
@inline num_neighbours(tri::Triangulation, u) = num_neighbours(get_graph(tri), u)

"""
    add_neighbour!(tri::Triangulation, u, v...)

Calls `add_neighbour!(get_graph(tri), u, v...)`, adding each `v` into the neighbourhood of each `u`.
"""
@inline add_neighbour!(tri::Triangulation, u, v...) = add_neighbour!(get_graph(tri), u, v...)

"""
    delete_neighbour!(tri::Triangulation, u, v...)

Calls `delete_neighbour!(get_graph(tri), u, v...)`, deleting each `v` from the neighbourhood of `u`.
"""
@inline delete_neighbour!(tri::Triangulation, u, v...) = delete_neighbour!(get_graph(tri), u, v...)

"""
    delete_vertex!(tri::Triangulation, u...)

Calls `delete_vertex!(get_graph(tri), u...)`, deleting each vertex `u` from `tri`.
"""
@inline delete_vertex!(tri::Triangulation, u...) = delete_vertex!(get_graph(tri), u...)

"""
    delete_boundary_vertices_from_graph!(tri::Triangulation)

Calls `delete_boundary_vertices_from_graph!(get_graph(tri))`, deleting all ghost vertices 
from the graph.
"""
@inline delete_boundary_vertices_from_graph!(tri::Triangulation) = delete_boundary_vertices_from_graph!(get_graph(tri))

"""
    each_vertex(tri::Triangulation)

Returns `each_vertex(get_graph(tri))`, an iterator over all vertices present 
in the triangulation.

!!! warning

    This iterator will include ghost vertices. If you want to exclude these, 
    see [`each_solid_vertex`](@ref). Alternatively, if you only want ghost vertices, 
    see [`each_ghost_vertex`](@ref).
"""
@inline each_vertex(tri::Triangulation) = each_vertex(get_graph(tri))

"""
    num_vertices(tri::Triangulation)

Returns `num_vertices(get_graph(tri))`, the number of vertices in the triangulation.
"""
@inline num_vertices(tri::Triangulation) = num_vertices(get_graph(tri))

"""
    has_vertex(tri::Triangulation, u)

Returns `true` if the vertex `u` is in the triangulation,
and `false` otherwise.
"""
@inline has_vertex(tri::Triangulation, u) = has_vertex(get_graph(tri), u)

"""
    has_boundary_vertices(tri::Triangulation)

Returns `true` if the triangulation has any ghost vertices, and `false` otherwise.
"""
@inline has_boundary_vertices(tri::Triangulation) = has_boundary_vertices(get_graph(tri))