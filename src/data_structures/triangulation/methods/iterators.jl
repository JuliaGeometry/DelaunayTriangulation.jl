"""
    each_point(tri::Triangulation) -> Points 

Returns an iterator over all points in `tri`. 

!!! danger "Missing vertices"

    If `tri` has vertices that are not yet present in the triangulation, e.g. if you have deleted vertices or have 
    some submerged vertices in a weighted triangulation, then the corresponding points will still be present in this 
    iterator. It is recommended that you instead consider [`each_vertex`](@ref), [`each_solid_vertex`](@ref), or 
    [`each_ghost_vertex`](@ref) together with [`get_point`](@ref) to get the coordinates.
"""
each_point(tri::Triangulation) = each_point(get_points(tri))

"""
    each_point_index(tri::Triangulation) -> Indices 

Returns an iterator over all point indices in `tri`.

!!! danger "Missing vertices"

    If `tri` has vertices that are not yet present in the triangulation, e.g. if you have deleted vertices or have 
    some submerged vertices in a weighted triangulation, then the corresponding point indices will still be present in this 
    iterator. It is recommended that you instead consider [`each_vertex`](@ref), [`each_solid_vertex`](@ref), or 
    [`each_ghost_vertex`](@ref).
"""
each_point_index(tri::Triangulation) = each_point_index(get_points(tri))

"""
    each_vertex(tri::Triangulation) -> Set{Vertex}

Returns an iterator over all vertices in `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of these vertices will be ghost vertices.

See also [`each_solid_vertex`](@ref) and [`each_ghost_vertex`](@ref).
"""
each_vertex(tri::Triangulation) = get_vertices(get_graph(tri))

"""
    AbstractEachVertex{V}

An abstract type for an iterator over vertices in a triangulation.
"""
abstract type AbstractEachVertex{V} end

"""
    EachSolidVertex{V,T}

An iterator over all solid vertices in a triangulation.

# Fields
- `vertices::V`: The iterator over all vertices in the triangulation.
- `tri::T`: The triangulation.
"""
struct EachSolidVertex{V,T} <: AbstractEachVertex{V}
    vertices::V
    tri::T
end

"""
    EachGhostVertex{V,T}

An iterator over all ghost vertices in a triangulation.

# Fields
- `vertices::V`: The iterator over all vertices in the triangulation.
- `tri::T`: The triangulation.
"""
struct EachGhostVertex{V,T} <: AbstractEachVertex{V}
    vertices::V
    tri::T
end
Base.IteratorSize(::Type{<:AbstractEachVertex}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractEachVertex{V}}) where {V} = Base.IteratorEltype(V)
Base.eltype(::Type{<:AbstractEachVertex{V}}) where {V} = eltype(V)
each_vertex(itr::AbstractEachVertex) = itr
num_vertices(itr::AbstractEachVertex) = length(itr)
Base.length(verts::EachSolidVertex) = num_solid_vertices(verts.tri)
Base.length(verts::EachGhostVertex) = num_ghost_vertices(verts.tri)

Base.show(io::IO, ::MIME"text/plain", itr::EachSolidVertex) = print(io, "EachSolidVertex iterator over vertices in a triangulation.")
Base.show(io::IO, ::MIME"text/plain", itr::EachGhostVertex) = print(io, "EachGhostVertex iterator over vertices in a triangulation.")

"""
    each_solid_vertex(tri::Triangulation) -> Set{Vertex}

Returns an iterator over all solid vertices in `tri`. 

See also [`each_vertex`](@ref) and [`each_ghost_vertex`](@ref).
"""
each_solid_vertex(tri::Triangulation) = EachSolidVertex(each_vertex(tri), tri)

"""
    each_ghost_vertex(tri::Triangulation) -> Set{Vertex}

Returns an iterator over all ghost vertices in `tri`.

See also [`each_vertex`](@ref) and [`each_solid_vertex`](@ref).
"""
each_ghost_vertex(tri::Triangulation) = EachGhostVertex(all_ghost_vertices(tri), tri)

function Base.iterate(itr::EachSolidVertex, state...)
    vertices_state = iterate(itr.vertices, state...)
    vertices_state === nothing && return nothing
    vertices, state = vertices_state
    while is_ghost_vertex(vertices)
        vertices_state = iterate(itr.vertices, state...)
        vertices_state === nothing && return nothing
        vertices, state = vertices_state
    end
    return vertices, state
end
Base.iterate(itr::EachGhostVertex, state...) = has_ghost_vertices(itr.tri) ? Base.iterate(itr.vertices, state...) : nothing

function Random.rand(rng::AbstractRNG, v::Random.SamplerTrivial{<:EachSolidVertex})
    itr = v[]
    verts = itr.vertices
    r = rand(rng, verts)
    while is_ghost_vertex(r)
        r = rand(rng, verts)
    end
    return r
end
function Random.rand(rng::AbstractRNG, v::Random.SamplerTrivial{<:EachGhostVertex})
    itr = v[]
    verts = itr.vertices
    r = rand(rng, verts)
    while !is_ghost_vertex(r)
        r = rand(rng, verts)
    end
    return r
end

"""
    each_triangle(tri::Triangulation) -> Triangles

Returns an iterator over all triangles in `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of these triangles will be ghost triangles.

See also [`each_solid_triangle`](@ref) and [`each_ghost_triangle`](@ref).
"""
each_triangle(tri::Triangulation) = each_triangle(get_triangles(tri))

"""
    AbstractEachTriangle{T}

An abstract type for an iterator over triangles in a triangulation.
"""
abstract type AbstractEachTriangle{T} end

"""
    EachSolidTriangle{V,T}

An iterator over all solid triangles in a triangulation.

# Fields
- `triangles::V`: The iterator over all triangles in the triangulation.
- `tri::T`: The triangulation.
"""
struct EachSolidTriangle{V,T} <: AbstractEachTriangle{V}
    triangles::V
    tri::T
end

"""
    EachGhostTriangle{V,T}

An iterator over all ghost triangles in a triangulation.

# Fields
- `triangles::V`: The iterator over all triangles in the triangulation.
- `tri::T`: The triangulation.
"""
struct EachGhostTriangle{V,T} <: AbstractEachTriangle{V}
    triangles::V
    tri::T
end
triangle_type(itr::AbstractEachTriangle) = triangle_type(itr.tri)
Base.IteratorSize(::Type{<:AbstractEachTriangle}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractEachTriangle{T}}) where {T} = Base.IteratorEltype(T)
Base.eltype(::Type{<:AbstractEachTriangle{T}}) where {T} = triangle_type(T)
each_triangle(itr::AbstractEachTriangle) = itr
Base.length(itr::EachSolidTriangle) = num_solid_triangles(itr.tri)
Base.length(itr::EachGhostTriangle) = num_ghost_triangles(itr.tri)

Base.show(io::IO, ::MIME"text/plain", itr::EachSolidTriangle) = print(io, "EachSolidTriangle iterator over triangles in a triangulation.")
Base.show(io::IO, ::MIME"text/plain", itr::EachGhostTriangle) = print(io, "EachGhostTriangle iterator over triangles in a triangulation.")

"""
    each_solid_triangle(tri::Triangulation) -> Triangles

Returns an iterator over all solid triangles in `tri`.

See also [`each_triangle`](@ref) and [`each_ghost_triangle`](@ref).
"""
each_solid_triangle(tri::Triangulation) = EachSolidTriangle(each_triangle(tri), tri)

"""
    each_ghost_triangle(tri::Triangulation) -> Triangles

Returns an iterator over all ghost triangles in `tri`.

See also [`each_triangle`](@ref) and [`each_solid_triangle`](@ref).
"""
each_ghost_triangle(tri::Triangulation) = EachGhostTriangle(each_triangle(tri), tri)

function Base.iterate(itr::EachSolidTriangle, state...)
    tri_state = iterate(itr.triangles, state...)
    tri_state === nothing && return nothing
    tri, state = tri_state
    while is_ghost_triangle(tri)
        tri_state = iterate(itr.triangles, state)
        tri_state === nothing && return nothing
        tri, state = tri_state
    end
    return tri, state
end
function Base.iterate(itr::EachGhostTriangle, state...)
    !has_ghost_triangles(itr.tri) && return nothing
    tri_state = iterate(itr.triangles, state...)
    tri_state === nothing && return nothing
    tri, state = tri_state
    while !is_ghost_triangle(tri)
        tri_state = iterate(itr.triangles, state)
        tri_state === nothing && return nothing
        tri, state = tri_state
    end
    return tri, state
end

function Random.rand(rng::AbstractRNG, v::Random.SamplerTrivial{<:EachSolidTriangle})
    itr = v[]
    tris = itr.triangles
    V = rand(rng, tris)
    while is_ghost_triangle(V)
        V = rand(rng, tris)
    end
    return V
end
function Random.rand(rng::AbstractRNG, v::Random.SamplerTrivial{<:EachGhostTriangle})
    itr = v[]
    tris = itr.triangles
    V = rand(rng, tris)
    while !is_ghost_triangle(V)
        V = rand(rng, tris)
    end
    return V
end

"""
    each_edge(tri::Triangulation) -> Edges

Returns an iterator over all edges in `tri`. Note that, if `has_ghost_triangles(tri)`,
then some of these edges will be ghost edges.

See also [`each_solid_edge`](@ref) and [`each_ghost_edge`](@ref).
"""
each_edge(tri::Triangulation) = get_edges(tri)

"""
    AbstractEachEdge{E}

An abstract type for an iterator over edges in a triangulation.
"""
abstract type AbstractEachEdge{E} end

"""
    EachSolidEdge{E,T}

An iterator over all solid edges in a triangulation.

# Fields
- `edges::E`: The iterator over all edges in the triangulation.
- `tri::T`: The triangulation.
"""
struct EachSolidEdge{E,T} <: AbstractEachEdge{E}
    edges::E
    tri::T
end

"""
    EachGhostEdge{E,T}

An iterator over all ghost edges in a triangulation.

# Fields
- `edges::E`: The iterator over all edges in the triangulation.
- `tri::T`: The triangulation.
"""
struct EachGhostEdge{E,T} <: AbstractEachEdge{E}
    edges::E
    tri::T
end
edge_type(itr::AbstractEachEdge) = edge_type(itr.tri)
Base.IteratorSize(::Type{<:AbstractEachEdge}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractEachEdge{E}}) where {E} = Base.IteratorEltype(E)
Base.eltype(::Type{<:AbstractEachEdge{E}}) where {E} = edge_type(E)
each_edge(itr::AbstractEachEdge) = itr
Base.length(itr::EachSolidEdge) = num_solid_edges(itr.tri)
Base.length(itr::EachGhostEdge) = num_ghost_edges(itr.tri)

Base.show(io::IO, ::MIME"text/plain", itr::EachSolidEdge) = print(io, "EachSolidEdge iterator over edges in a triangulation.")
Base.show(io::IO, ::MIME"text/plain", itr::EachGhostEdge) = print(io, "EachGhostEdge iterator over edges in a triangulation.")


"""
    each_solid_edge(tri::Triangulation) -> Edges

Returns an iterator over all solid edges in `tri`.

See also [`each_edge`](@ref) and [`each_ghost_edge`](@ref).
"""
each_solid_edge(tri::Triangulation) = EachSolidEdge(each_edge(tri), tri)

"""
    each_ghost_edge(tri::Triangulation) -> Edges

Returns an iterator over all ghost edges in `tri`.

See also [`each_edge`](@ref) and [`each_solid_edge`](@ref).
"""
each_ghost_edge(tri::Triangulation) = EachGhostEdge(each_edge(tri), tri)

function Base.iterate(itr::EachSolidEdge, state...)
    edges_state = iterate(itr.edges, state...)
    edges_state === nothing && return nothing
    edges, state = edges_state
    while is_ghost_edge(edges)
        edges_state = iterate(itr.edges, state)
        edges_state === nothing && return nothing
        edges, state = edges_state
    end
    return edges, state
end
function Base.iterate(itr::EachGhostEdge, state...)
    edges_state = iterate(itr.edges, state...)
    edges_state === nothing && return nothing
    edges, state = edges_state
    while !is_ghost_edge(edges)
        edges_state = iterate(itr.edges, state)
        edges_state === nothing && return nothing
        edges, state = edges_state
    end
    return edges, state
end

function Random.rand(rng::AbstractRNG, v::Random.SamplerTrivial{<:EachSolidEdge})
    itr = v[]
    edges = itr.edges
    e = rand(rng, edges)
    while is_ghost_edge(e)
        e = rand(rng, edges)
    end
    return e
end
function Random.rand(rng::AbstractRNG, v::Random.SamplerTrivial{<:EachGhostEdge})
    itr = v[]
    edges = itr.edges
    e = rand(rng, edges)
    while !is_ghost_edge(e)
        e = rand(rng, edges)
    end
    return e
end