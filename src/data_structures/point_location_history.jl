"""
    PointLocationHistory{T,E,I}

History from using [`find_triangle`](@ref).

# Fields 
- `triangles::Vector{T}`: The visited triangles. 
- `collinear_segments::Vector{E}`: Segments collinear with the original line `pq` using to jump.
- `collinear_point_indices::Vector{I}`: This field contains indices to segments in `collinear_segments` that refer to points that were on the original segment, but there is no valid segment for them. We use manually fix this after the fact. For example, we could add an edge `(1, 14)`, when really we mean something like `(7, 14)` which isn't a valid edge.
- `left_vertices::Vector{I}`: Vertices from the visited triangles to the left of `pq`.
- `right_verices::Vector{I}`: Vertices from the visited triangles to the right of `pq`.
"""
struct PointLocationHistory{T,E,I}
    triangles::Vector{T}
    collinear_segments::Vector{E}
    collinear_point_indices::Vector{I}
    left_vertices::Vector{I}
    right_vertices::Vector{I}
    PointLocationHistory{T,E,I}() where {T,E,I} = new{T,E,I}(T[], E[], I[], I[], I[])
end

"""
    add_triangle!(history::PointLocationHistory, i, j, k)

Adds the triangle `(i, j, k)` to the `triangles` field of `history`.
"""
add_triangle!(history::Ts, i::I, j::I, k::I) where {Ts<:PointLocationHistory,I<:Integer} = add_triangle!(history.triangles, i, j, k)

"""
    add_edge!(history::PointLocationHistory{T,E}, i, j)

Adds the edge `(i, j)` to the `collinear_segments` field of `history`.
"""
add_edge!(history::PointLocationHistory{T,E}, i, j) where {T,E} = add_edge!(history.collinear_segments, construct_edge(E, i, j))

"""
    add_left_vertex!(history::PointLocationHistory, i)

Adds the vertex `i` to the `left_vertices` field of `history`.
""" 
add_left_vertex!(history::PointLocationHistory, i) = push!(history.left_vertices, i)

"""
    add_right_vertex!(history::PointLocationHistory, j)

Adds the vertex `j` to the `right_vertices` field of `history`.
"""
add_right_vertex!(history::PointLocationHistory, j) = push!(history.right_vertices, j)

"""
    add_index!(history::PointLocationHistory, i)

Adds the index `i` to the `collinear_point_indices` field of `history`.
"""
add_index!(history::PointLocationHistory, i) = push!(history.collinear_point_indices, i)

"""
    num_edges(history::PointLocationHistory) -> Integer

Returns the number of edges in `history.collinear_segments`.
"""
num_edges(history::PointLocationHistory) = num_edges(history.collinear_segments)