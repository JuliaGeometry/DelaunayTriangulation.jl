"""
    VoronoiCell{I}

A structure to store the vertices of a Voronoi cell, stored 
in a counter-clockwise order that [`is_circular`](@ref). There is only 
one field, `vertices`, used for this storage.

!!! warning 

    No checks are made for orientation or circularity.
"""
struct VoronoiCell{I}
    vertices::Vector{I} # note that is_circular(vertices)
end
get_vertices(cell::VoronoiCell) = cell.vertices
get_vertices(cell::VoronoiCell, i) = get_vertices(cell)[i]
VoronoiCell{I}() where {I} = VoronoiCell{I}(I[])
add_vertex!(cell::VoronoiCell{I}, v::I) where {I} = push!(cell.vertices, v)
Base.empty!(cell::VoronoiCell) = empty!(cell.vertices)
Base.length(cell::VoronoiCell) = length(cell.vertices)
Base.iterate(cell::VoronoiCell, state...) = Base.iterate(get_vertices(cell), state...)
Base.eltype(::VoronoiCell{I}) where {I} = I
Base.getindex(cell::VoronoiCell, i) = get_vertices(cell, i)
Base.firstindex(cell::VoronoiCell) = firstindex(get_vertices(cell))
Base.lastindex(cell::VoronoiCell) = lastindex(get_vertices(cell))
Base.eachindex(cell::VoronoiCell) = eachindex(get_vertices(cell))

function Base.show(io::IO, vor::VoronoiCell)
    Base.show(io, get_vertices(vor))
end
