"""
    get_edges(tri::Triangulation) -> Set{NTuple{2,Vertex}}

Returns the set of all edges in `tri`. Orientation is ignored, so that 
only one of `(i, j)` and `(j, i)` will appear in the result. Note that, 
if `has_ghost_triangles(tri)`, then some of these edges will be ghost edges.

See also [`each_edge`](@ref), [`each_solid_edge`](@ref), and [`each_ghost_edge`](@ref).
"""
get_edges(tri::Triangulation) = get_edges(get_graph(tri))


"""
    get_vertices(tri::Triangulation) -> Set{Vertex}

Returns the set of all vertices in `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of these vertices will be ghost vertices.

See also [`each_vertex`](@ref), [`each_solid_vertex`](@ref), and [`each_ghost_vertex`](@ref).
"""
get_vertices(tri::Triangulation) = get_vertices(get_graph(tri))


"""
    get_neighbours(tri::Triangulation) -> Dict{Vertex, Set{Vertex}}

Returns the `neighbours` map of `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of the neighbours and vertices will be ghost vertices.
"""
get_neighbours(tri::Triangulation) = get_neighbours(get_graph(tri))


"""
    get_neighbours(tri::Triangulation, u) -> Set{Vertex}

Returns the set of neighbours of `u` in `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of the neighbours and vertices will be ghost vertices.
"""
get_neighbours(tri::Triangulation, u) = get_neighbours(get_graph(tri), u)


"""
    add_vertex!(tri::Triangulation, u...)

Adds the vertices `u...` into the graph of `tri`.
"""
add_vertex!(tri::Triangulation, u...) = add_vertex!(get_graph(tri), u...)


"""
    has_vertex(tri::Triangulation, u) -> Bool

Returns `true` if `u` is a vertex in `tri`, and `false` otherwise.
"""
has_vertex(tri::Triangulation, u) = has_vertex(get_graph(tri), u)


"""
    num_neighbours(tri::Triangulation, u) -> Integer

Returns the number of neighbours of `u` in `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of the neighbours counted might be ghost vertices if `u` is a boundary vertex.
"""
num_neighbours(tri::Triangulation, u) = num_neighbours(get_graph(tri), u)


"""
    add_neighbour!(tri::Triangulation, u, v...)

Adds the neighbours `v...` to `u` in the graph of `tri`.
"""
add_neighbour!(tri::Triangulation, u, v...) = add_neighbour!(get_graph(tri), u, v...)


"""
    delete_neighbour!(tri::Triangulation, u, v...)

Deletes the neighbours `v...` from `u` in the graph of `tri`.
"""
delete_neighbour!(tri::Triangulation, u, v...) = delete_neighbour!(get_graph(tri), u, v...)


"""
    delete_vertex!(tri::Triangulation, u...)

Deletes the vertices `u...` from the graph of `tri`.
"""
delete_vertex!(tri::Triangulation, u...) = delete_vertex!(get_graph(tri), u...)


"""
    delete_ghost_vertices_from_graph!(tri::Triangulation)

Deletes all ghost vertices from the graph of `tri`.
"""
delete_ghost_vertices_from_graph!(tri::Triangulation) = delete_ghost_vertices_from_graph!(get_graph(tri))


"""
    num_vertices(tri::Triangulation) -> Integer

Returns the number of vertices in `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of these vertices will be ghost vertices.

See also [`num_solid_vertices`](@ref) and [`num_ghost_vertices`](@ref).
"""
num_vertices(tri::Triangulation) = num_vertices(get_graph(tri))


"""
    has_ghost_vertices(tri::Triangulation) -> Bool

Returns `true` if `tri` has ghost vertices, and `false` otherwise.
"""
has_ghost_vertices(tri::Triangulation) = has_ghost_vertices(get_graph(tri))


"""
    num_ghost_vertices(tri::Triangulation) -> Integer

Returns the number of ghost vertices in `tri`. 

See also [`num_solid_vertices`](@ref) and [`num_vertices`](@ref).
"""
num_ghost_vertices(tri::Triangulation) = length(all_ghost_vertices(tri)) * has_ghost_vertices(tri)


"""
    num_solid_vertices(tri::Triangulation) -> Integer

Returns the number of solid vertices in `tri`.

See also [`num_ghost_vertices`](@ref) and [`num_vertices`](@ref).
"""
num_solid_vertices(tri::Triangulation) = num_vertices(tri) - num_ghost_vertices(tri)


"""
    num_edges(tri::Triangulation) -> Integer

Returns the number of edges in `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of these edges will be ghost edges.

See also [`num_solid_edges`](@ref) and [`num_ghost_edges`](@ref).
"""
num_edges(tri::Triangulation) = num_edges(get_graph(tri))


"""
    num_ghost_edges(tri::Triangulation) -> Integer

Returns the number of ghost edges in `tri`.

See also [`num_solid_edges`](@ref) and [`num_edges`](@ref).
"""
function num_ghost_edges(tri::Triangulation)
    I = integer_type(tri)
    !has_vertex(tri, I(𝒢)) && return 0 # if it doesn't have 𝒢, it doesn't have any ghosts
    all_ghosts = all_ghost_vertices(tri)
    num_ghosts = 0
    for i in all_ghosts
        !has_vertex(tri, i) && continue # in case delete_holes=false when triangulate was called, there could be issues here
        bnd_ngh = get_neighbours(tri, i)
        num_ghosts += num_edges(bnd_ngh)
    end
    return num_ghosts
end


"""
    num_solid_edges(tri::Triangulation) -> Integer

Returns the number of solid edges in `tri`.

See also [`num_ghost_edges`](@ref) and [`num_edges`](@ref).
"""
num_solid_edges(tri::Triangulation) = num_edges(tri) - num_ghost_edges(tri)


"""
    sort_edge_by_degree(tri::Triangulation, e) -> Edge

Returns the edge `e` sorted so that `initial(e)` has the smaller degree of the two vertices.
"""
function sort_edge_by_degree(tri::Triangulation, e)
    u, v = edge_vertices(e)
    d₁ = num_neighbours(tri, u)
    d₂ = num_neighbours(tri, v)
    if d₁ ≤ d₂
        return e
    else
        return reverse_edge(e)
    end
end
