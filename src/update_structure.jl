"""
    add_edge!(𝒜::Adjacent, u, v, w)

This adds the edge `(u, v)` to the adjacency list `𝒜`, mapping 
`(u, v)` to `w` so that `(u, v, w)` is a positively oriented triangle.
"""
function add_edge!(𝒜::Adjacent, u, v, w)
    𝒜[(u, v)] = w
    return nothing
end
"""
    add_edge!(𝒜⁻¹::Adjacent2Vertex, w, u, v)

This adds the edge `(u, v)` to `𝒜⁻¹[w]`.
"""
function add_edge!(𝒜⁻¹::Adjacent2Vertex, w, u, v)
    uv = get!(Adjacent2VertexVector, 𝒜⁻¹, w) # https://discourse.julialang.org/t/how-do-i-append-add-data-in-dictionary-with-same-key/15891/5?u=legola18
    push!(uv, (u, v))
    return nothing
end
"""
    add_edge!(G::HistoryDAG, u, v)

Adds an edge between triangle `u ∈ G` and `v ∈ G` to the graph `G`. Checks 
are made for circular permutations of `u` and `v`'s indices.
"""
function add_edge!(𝒟::HistoryDAG, 𝒯::TriangleType, 𝒯new::TriangleType)
    for i in 1:3
        for j in 1:3
            i₁, i₂, i₃ = i, (i % 3) + 1, ((i + 1) % 3) + 1
            j₁, j₂, j₃ = j, (j % 3) + 1, ((j + 1) % 3) + 1
            T = TriangleType((𝒯[i₁], 𝒯[i₂], 𝒯[i₃]))
            V = TriangleType((𝒯new[j₁], 𝒯new[j₂], 𝒯new[j₃]))
            if has(graph(𝒟), T) && has(graph(𝒟), V)
                SimpleGraphs.add!(graph(𝒟), T, V)
                return nothing
            end
        end
    end
    return nothing
end
function add_edge!(𝒟::HistoryDAG, 𝒯::TriangleType, 𝒯new::TriangleType...)
    for T in 𝒯new
        add_edge!(𝒟, 𝒯, T)
    end
end
"""
    add_edge!(𝒟𝒢::DelaunayGraph, u, v...)
"""
function add_edge!(𝒟𝒢::DelaunayGraph, u, v...)
    for w in v
        SimpleGraphs.add!(graph(𝒟𝒢), u, w)
    end
    return nothing
end
add_neighbour!(𝒟𝒢::DelaunayGraph, u, v...) = add_edge!(𝒟𝒢, u, v...)#alias

"""
    update_adjacent!(𝒜::Adjacent,, 𝒯::TriangleType...)

Updates the adjacency list `𝒜`, and its inverse `𝒜⁻¹`, with the edges from `𝒯`, keying on each edge. 
`𝒯` must be positively oriented.
"""
function update_adjacent!(𝒜::Adjacent, 𝒯::TriangleType)
    i, j, k = 𝒯
    add_edge!(𝒜, i, j, k)
    add_edge!(𝒜, j, k, i)
    add_edge!(𝒜, k, i, j)
    return nothing
end
function update_adjacent!(𝒜::Adjacent, 𝒯::TriangleType...)
    for T in 𝒯
        update_adjacent!(𝒜, T)
    end
    return nothing
end

"""
    update_adjacent2vertex_flip!(𝒜⁻¹::Adjacent2Vertex, Tᵣₖⱼ::TriangleType, Tᵣᵢₖ::TriangleType)

Updates `𝒜⁻¹` after an edge flip. See [`flip_edge!`](@ref).
"""
function update_adjacent2vertex_flip!(𝒜⁻¹::Adjacent2Vertex, i, j, k, r)
    # Delete the necessary edges
    delete_edge!(𝒜⁻¹, i, k, j)
    delete_edge!(𝒜⁻¹, i, j, r)
    delete_edge!(𝒜⁻¹, j, r, i)
    delete_edge!(𝒜⁻¹, j, i, k)
    delete_edge!(𝒜⁻¹, k, j, i)
    delete_edge!(𝒜⁻¹, r, i, j)
    # Add the new edges 
    add_edge!(𝒜⁻¹, i, k, r)
    add_edge!(𝒜⁻¹, j, r, k)
    add_edge!(𝒜⁻¹, k, j, r)
    add_edge!(𝒜⁻¹, k, r, i)
    add_edge!(𝒜⁻¹, r, k, j)
    add_edge!(𝒜⁻¹, r, i, k)
    return nothing
end

"""
    update_adjacent2vertex_addition!(𝒜⁻¹::Adjacent2Vertex, 𝒯ᵢⱼᵣ, 𝒯ⱼₖᵣ, 𝒯ₖᵢᵣ)

Updates `𝒜⁻¹` after the insertion of a point into the triangle `(i, j, k)`.
"""
function update_adjacent2vertex_addition!(𝒜⁻¹::Adjacent2Vertex, i, j, k, r)
    # Delete old edges 
    delete_edge!(𝒜⁻¹, i, j, k)
    delete_edge!(𝒜⁻¹, j, k, i)
    delete_edge!(𝒜⁻¹, k, i, j)
    # Add the new edges 
    add_edge!(𝒜⁻¹, i, j, r)
    add_edge!(𝒜⁻¹, i, r, k)
    add_edge!(𝒜⁻¹, j, k, r)
    add_edge!(𝒜⁻¹, j, r, i)
    add_edge!(𝒜⁻¹, k, i, r)
    add_edge!(𝒜⁻¹, k, r, j)
    add_edge!(𝒜⁻¹, r, i, j)
    add_edge!(𝒜⁻¹, r, k, i)
    add_edge!(𝒜⁻¹, r, j, k)
    return nothing
end
"""
    delete_triangle!(𝒯, T...)

Deletes the triangle(s) `T` from `𝒯`, and all circular shifts of their indices.
"""
function delete_triangle!(𝒯::Triangles, T::TriangleType)
    #deleteat!(𝒯, findall(x -> x == T || x == T[[2, 3, 1]] || x == T[[3, 1, 2]], 𝒯))
    delete!(𝒯.triangles, T)
    delete!(𝒯.triangles, (T[2], T[3], T[1]))
    delete!(𝒯.triangles, (T[3], T[1], T[2]))
    return nothing
end
function delete_triangle!(𝒯::Triangles, T::TriangleType...)
    for T in T
        delete_triangle!(𝒯, T)
    end
    return nothing
end

"""
    add_triangle!(𝒯, T::TriangleType...)

Adds the triangle(s) `T` to `𝒯`.
"""
function add_triangle!(𝒯::Triangles, T::TriangleType...)
    push!(𝒯.triangles, T...)
    return nothing
end
"""
    add_triangle!(G::HistoryDAG, u)

Add the triangle `u` to the graph `G`. Checks that cyclic permutations of `u`'s indices are not 
in `G` already, returning `false` otherwise.
"""
function add_triangle!(𝒟::HistoryDAG, 𝒯::TriangleType)
    for i = 1:3
        i₁, i₂, i₃ = i, (i % 3) + 1, ((i + 1) % 3) + 1
        has(graph(𝒟), TriangleType((𝒯[i₁], 𝒯[i₂], 𝒯[i₃]))) && return false
    end
    SimpleGraphs.add!(graph(𝒟), 𝒯)
    return nothing
end
function add_triangle!(𝒟::HistoryDAG, 𝒯::TriangleType...)
    for T in 𝒯
        add_triangle!(𝒟, T)
    end
end

"""
    delete_edge!(𝒜, 𝒜⁻¹::Adjacent2Vertex i, j; protect_boundary = true)

Deletes the keys `(i, j)` and `(j, i)` from `𝒜`. Use `protect_boundary=true` to avoid 
removing edges on the boundary.
"""
function delete_edge!(𝒜::Adjacent, i, j; protect_boundary=true)
    (!protect_boundary || 𝒜(i, j) ≠ BoundaryIdx) && delete!(𝒜.adjacent, (i, j))
    (!protect_boundary || 𝒜(j, i) ≠ BoundaryIdx) && delete!(𝒜.adjacent, (j, i))
    return nothing
end
"""
    delete_edge!(𝒜⁻¹::Adjacent2Vertex, w, u, v)

Deletes the edge `(u, v)` from `𝒜⁻¹[w]`.
"""
function delete_edge!(𝒜⁻¹::Adjacent2Vertex, w, u, v)
    delete!(𝒜⁻¹.adjacent2vertex[w], (u, v))
    return nothing
end

"""
    delete_point_from_neighbour!(𝒟𝒢, u, v)

Removes the point `v` from the neighbourhood of `u`, updating `𝒟𝒢` in-place.
"""
function delete_neighbour!(𝒟𝒢::DelaunayGraph, u, v)
    SimpleGraphs.delete!(graph(𝒟𝒢), u, v)
    return nothing
end

"""
    add_point!(DG::DelaunayGraph, u::Integer)
    
Adds the point `u` into the graph `DG`.    
"""
add_point!(DG::DelaunayGraph, u::Integer) = (add!(graph(DG), u); return nothing)
function add_point!(DG::DelaunayGraph, u::Integer...)
    for v in u
        add_point!(DG, v)
    end
    return nothing
end

"""
    delete_point!(𝒟𝒢::DelaunayGraph, u...)

Deletes the point(s) `u` from `𝒟𝒢`. Note that this also deletes `u` from the neighbours, as 
`𝒟𝒢` is an undirected graph.
"""
function delete_point!(𝒟𝒢::DelaunayGraph, u)
    delete!(graph(𝒟𝒢), u)
    return nothing
end
function delete_point!(𝒟𝒢::DelaunayGraph, u...)
    for v in u
        delete_point!(𝒟𝒢, v)
    end
    return nothing
end
"""
    delete_point!(𝒜⁻¹, w)

Deletes the key `w` from `𝒜⁻¹`.
"""
function delete_point!(𝒜⁻¹, w)
    delete!(𝒜⁻¹.adjacent2vertex, w)
    return nothing
end