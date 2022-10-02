############################################
##
## INDEXING AND UPDATING OF THE TRIANGULATION DATA STRUCTURES
##
############################################
## Adjacent  
""" 
    add_edge!(adj::Adjacent, uv::AbstractEdge, w)
    add_edge!(adj::Adjacent{I,E}, u, v, w) where {I,E<:AbstractEdge}

Adds the edge `E(u, v)` to the adjacent map `adj`.
"""
function add_edge!(adj::Adjacent, uv::AbstractEdge, w)
    adjacent(adj)[uv] = w
    return nothing
end
@doc (@doc add_edge!(::Adjacent::AbstractEdge, ::Any))
function add_edge!(adj::Adjacent{I,E}, u, v, w) where {I,E<:AbstractEdge}
    add_edge!(adj, E(u, v), w)
    return nothing
end

"""
    get_edge(adj::Adjacent, uv::AbstractEdge)
    get_edge(adj::Adjacent{I,E}, u, v) where {I,E<:AbstractEdge}

Accessor for the key `E(u, v)` in the adjacent map `adj`, so that `Triangle(u, v, get_edge(adj, u, v))` 
is positively oriented.
"""
function get_edge(adj::Adjacent, uv::AbstractEdge)
    return adjacent(adj)[uv]
end
function get_edge(adj::Adjacent{I,E}, u, v) where {I,E<:AbstractEdge}
    return get_edge(adj, E(u, v))
end

"""
    delete_edge!(adj::Adjacent{I,E}, i, j; protect_boundary=true) where {I,E<:AbstractEdge{I}}

Deletes the edge `E(i, j)` from the map `adj`. This function also deletes `E(j, i)` from `adj`. If 
the boundary needs to be protected, so that `E(i, j)` is not deleted if `E(i, j) = BoundaryIndex`
(and similarly for `E(j, i)`), set `protect_boundary=true`. 
"""
function delete_edge!(adj::Adjacent{I,E}, i, j; protect_boundary=true) where {I,E<:AbstractEdge{I}}
    (!protect_boundary || get_edge(adj, i, j) ≠ BoundaryIndex) && delete!(adjacent(adj), E(i, j))
    (!protect_boundary || get_edge(adj, j, i) ≠ BoundaryIndex) && delete!(adjacent(adj), E(j, i))
    return nothing
end

"""
    add_triangle!(adj::Adjacent, T::AbstractTriangle...)

Adds the edges of the triangle `T` to the adjacent map `adj`.
"""
function add_triangle!(adj::Adjacent, T::AbstractTriangle)
    i, j, k = T
    add_edge!(adj, i, j, k)
    add_edge!(adj, j, k, i)
    add_edge!(adj, k, i, j)
    return nothing
end
@doc (@doc add_triangle!(::Adjacent, ::AbstractTriangle))
function add_triangle!(adj::Adjacent, T::AbstractTriangle...)
    for T in T
        add_triangle!(adj, T)
    end
    return nothing
end

## Adjacent2Vertex
"""
    add_edge!(adj2v::Adjacent2Vertex{I,E}, w, uv::E) where {I,E<:AbstractEdge{I}}
    add_edge!(adj2v::Adjacent2Vertex{I,E}, w, u, v) where {I,E<:AbstractEdge{I}}

Adds the edge `E(u, v)` to the set of edges in `adjacent2vertex(adj2v)[w]`.
"""
function add_edge!(adj2v::Adjacent2Vertex{I,E}, w, uv::E) where {I,E<:AbstractEdge{I}}
    existing_segments = get!(Set{E}, adjacent2vertex(adj2v), w)
    push!(existing_segments, uv)
    return nothing
end
@doc (@doc add_edge!(::Adjacent2Vertex{I,E}, ::Any, ::E) where {I,E<:AbstractEdge{I}})
function add_edge!(adj2v::Adjacent2Vertex{I,E}, w, u, v) where {I,E<:AbstractEdge{I}}
    add_edge!(adj2v, w, E(u, v))
    return nothing
end

"""
    get_edge(adj2v::Adjacent2Vertex, w)

Accessor for the set of edges `(u, v)` such that `Triangle(u, v, w)` is positively oriented.
"""
function get_edge(adj2v::Adjacent2Vertex, w)
    return adjacent2vertex(adj2v)[w]
end

"""
    delete_edge!(adj2v::Adjacent2Vertex{I,E}, w, uv::E) where {I,E<:AbstractEdge{I}}
    delete_edge!(adj2v::Adjacent2Vertex{I,E}, w, u, v) where {I,E<:AbstractEdge{I}}

Deletes the edge `E(u, v)` from the set of edges around `w` that together define a positively 
oriented triangle.
"""
function delete_edge!(adj2v::Adjacent2Vertex{I,E}, w, uv::E) where {I,E<:AbstractEdge{I}}
    delete!(get_edge(adj2v, w), uv)
    return nothing
end
@doc (@doc delete_edge!(::Adjacent2Vertex{I,E}, ::Any, ::AbstractEdge) where {I,E<:AbstractEdge{I}})
function delete_edge!(adj2v::Adjacent2Vertex{I,E}, w, u, v) where {I,E<:AbstractEdge{I}}
    delete_edge!(adj2v, w, E(u, v))
    return nothing
end

"""
    delete_point!(adj2v::Adjacent2Vertex, w) 

Deletes the point `w` from the map `adj2v`.
"""
function delete_point!(adj2v::Adjacent2Vertex, w)
    delete!(adjacent2vertex(adj2v), w)
    return nothing
end

"""
    update_after_flip!(adj2v::Adjacent2Vertex, i, j, k, r)

Updates the map `adj2v` after the edge `(i, j)` is flipped to `(k, r)`, where 
the triangles `(i, j, r)` and `(k, j, i)` are positively oriented.
"""
function update_after_flip!(adj2v::Adjacent2Vertex, i, j, k, r)
    # Delete the necessary edges
    delete_edge!(adj2v, i, k, j)
    delete_edge!(adj2v, i, j, r)
    delete_edge!(adj2v, j, r, i)
    delete_edge!(adj2v, j, i, k)
    delete_edge!(adj2v, k, j, i)
    delete_edge!(adj2v, r, i, j)
    # Add the new edges 
    add_edge!(adj2v, i, k, r)
    add_edge!(adj2v, j, r, k)
    add_edge!(adj2v, k, j, r)
    add_edge!(adj2v, k, r, i)
    add_edge!(adj2v, r, k, j)
    add_edge!(adj2v, r, i, k)
    return nothing
end

"""
    update_after_insertion!(adj2v::Adjacent2Vertex, i, j, k, r)

Updates the map `adj2v` after the point `r` is put into the positively 
oriented triangle `(i, j, k)`, thus subdividing the triangle into the 
new triangles `(i, j, r)`, `(j, k, r)`, and `(k, i, r)`.`
"""
function update_after_insertion!(adj2v::Adjacent2Vertex, i, j, k, r)
    # Delete old edges 
    delete_edge!(adj2v, i, j, k)
    delete_edge!(adj2v, j, k, i)
    delete_edge!(adj2v, k, i, j)
    # Add the new edges 
    add_edge!(adj2v, i, j, r)
    add_edge!(adj2v, i, r, k)
    add_edge!(adj2v, j, k, r)
    add_edge!(adj2v, j, r, i)
    add_edge!(adj2v, k, i, r)
    add_edge!(adj2v, k, r, j)
    add_edge!(adj2v, r, i, j)
    add_edge!(adj2v, r, k, i)
    add_edge!(adj2v, r, j, k)
    return nothing
end

"""
    update_after_split!(adj2v::Adjacent2Vertex, i, j, k, ℓ, r)

Updates the map `adj2v` after the edge `(i, j)`, incident to the two positively 
oriented triangles `(i, j, k)` and `(j, i, ℓ)`, is split into two so that the triangles 
`(i, j, k)` and `(j, i, ℓ`) are split from a point `r` to the points `k` and `ℓ`, 
respectively, assuming `r` is on the edge `(i, j)`.
"""
function update_after_split!(adj2v::Adjacent2Vertex, i, j, k, ℓ, r)
    # Delete old edges 
    delete_edge!(adj2v, i, j, k)
    delete_edge!(adj2v, i, ℓ, j)
    delete_edge!(adj2v, j, k, i)
    delete_edge!(adj2v, j, i, ℓ)
    delete_edge!(adj2v, k, i, j)
    delete_edge!(adj2v, ℓ, j, i)
    # Add new edges 
    add_edge!(adj2v, i, r, k)
    add_edge!(adj2v, i, ℓ, r)
    add_edge!(adj2v, j, k, r)
    add_edge!(adj2v, j, r, ℓ)
    add_edge!(adj2v, k, i, r)
    add_edge!(adj2v, k, r, j)
    add_edge!(adj2v, ℓ, r, i)
    add_edge!(adj2v, ℓ, j, r)
    add_edge!(adj2v, r, k, i)
    add_edge!(adj2v, r, j, k)
    add_edge!(adj2v, r, ℓ, r)
    add_edge!(adj2v, r, ℓ, j)
    return nothing
end

## DelaunayGraph 
"""
    add_point!(DG::DelaunayGraph, u...)

Adds the point(s) `u` into the vertex set of `DG`.
"""
function add_point!(DG::DelaunayGraph, u)
    add!(graph(DG), u)
    return nothing
end
@doc (@doc add_point!(::DelaunayGraph, ::Any))
function add_point!(DG::DelaunayGraph, u...)
    for v in u
        add_point!(DG, v)
    end
    return nothing
end

"""
    add_neighbour!(DG::DelaunayGraph{I}, u, v...) 

Adds an edge from `u` to the points in `v` into the graph `DG`. Note that, as 
`DG` is an undirected graph, this also adds an edge from `v` to `u`.
"""
function add_neighbour!(DG::DelaunayGraph, u, v)
    add!(graph(DG), u, v)
    return nothing
end
@doc (@doc add_neighbour!(::DelaunayGraph, ::Any, ::Any))
function add_neighbour!(DG::DelaunayGraph, u, v...)
    for w in v
        add_neighbour!(DG, u, w)
    end
    return nothing
end

"""
    get_neighbour(DG::DelaunayGraph, u)

Extracts the neighbours of `u` from the graph `DG`.
"""
function get_neighbour(DG::DelaunayGraph, u)
    return graph(DG).N[u]
end

"""
    delete_neighbour!(DG::DelaunayGraph, u, v)

Deletes `v` from the neighbourhood of `u`. Note that, as `DG` is an 
undirected graph, this also deletes `u` from the neighbourhood of `v`.
"""
function delete_neighbour!(DG::DelaunayGraph, u, v)
    delete!(graph(DG), u, v)
    return nothing
end

"""
    delete_point!(DG::DelaunayGraph, u...)

Deletes the point(s) `u` from the vertex set of `DG`.
"""
function delete_point!(DG::DelaunayGraph, u)
    delete!(graph(DG), u)
    return nothing
end
@doc (@doc delete_point!(::DelaunayGraph, ::Any))
function delete_point!(DG::DelaunayGraph, u...)
    for v in u
        delete_point!(DG, v)
    end
    return nothing
end

## HistoryGraph 
"""
    SimpleGraphs.out_neighbors(G::HistoryGraph, T)

Returns the set of triangles that `T` was replaced by in the triangulation.
"""
SimpleGraphs.out_neighbors(G::HistoryGraph, T) = graph(G).N[T] # The version in SimpleGraphs is allocating...
SimpleGraphs.in_neighbors(G::HistoryGraph, T) = SimpleGraphs.in_neighbors(graph(G), T) # Need the allocating version for find_root. Only place it's used anyway.

SimpleGraphs.in_deg(G::HistoryGraph, T) = SimpleGraphs.in_deg(graph(G), T)
SimpleGraphs.out_deg(G::HistoryGraph, T) = SimpleGraphs.out_deg(graph(G), T)

"""
    add_triangle!(D::HistoryGraph, T::AbstractTriangle...)

Adds the triangles in `T` into the graph `D`. Checks are made for 
circular shifts of the vertices of `T` to avoid duplicates.
"""
function add_triangle!(D::HistoryGraph, T::AbstractTriangle)
    has(graph(D), T) && return nothing
    has(graph(D), shift_triangle_1(T)) && return nothing
    has(graph(D), shift_triangle_2(T)) && return nothing
    add!(graph(D), T)
    return nothing
end
@doc (@doc add_triangle!(::HistoryGraph, ::AbstractTriangle))
function add_triangle!(D::HistoryGraph, T::AbstractTriangle...)
    for V in T
        add_triangle!(D, V)
    end
end

"""
    add_edge!(G::HistoryGraph, T::AbstractTriangle, V::AbstractTriangle...)

Adds a directed edge from `T` to the triangles in `V`, making checks for circular shifts 
in the vertices of `T` and `V`.
"""
function add_edge!(G::HistoryGraph, T::AbstractTriangle, V::AbstractTriangle)
    for i in 0:2
        Tshift = shift_triangle(T, i)
        if has(graph(G), Tshift)
            for j in 0:2
                Vshift = shift_triangle(V, j)
                if has(graph(G), Vshift)
                    SimpleGraphs.add!(graph(G), Tshift, Vshift)
                    return nothing
                end
            end
        end
    end
    return nothing
end
@doc (@doc add_edge!(::HistoryGraph, ::AbstractTriangle, ::AbstractTriangle))
function add_edge!(G::HistoryGraph, T::AbstractTriangle, V::AbstractTriangle...)
    for V in V
        add_edge!(G, T, V)
    end
    return nothing
end