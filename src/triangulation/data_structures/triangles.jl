num_solid_triangles(tri::Triangulation) = num_solid_triangles(get_triangle_counts(tri))
num_ghost_triangles(tri::Triangulation) = num_ghost_triangles(get_triangle_counts(tri))
num_triangles(tri::Triangulation) = num_triangles(get_triangle_counts(tri))

struct Triangles{T,I} <: AbstractSet{NTuple{3,I}}
    triangulation::T
    Triangles(tri::T) where {T} = new{T, integer_type(tri)}(tri)
end
Base.length(triangles::Triangles) = num_triangles(get_triangulation(triangles))
function Base.iterate(triangles::Triangles, state...)
    tri = get_triangulation(triangles)
    adj = get_adjacent(tri)
    dict = get_adjacent(adj)
    
end

get_triangulation(triangles::Triangles) = triangles.triangulation