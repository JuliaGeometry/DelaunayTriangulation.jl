"""
    VoronoiTessellation{Tr<:Triangulation,P,I,T,S,E}

Struct for a Voronoi tessellation. 

# Fields 
- `triangulation::Tr`

The underlying triangulation, dual to the tessellation (unless there are constraints, in which case the tessellation is 
no longer dual, but it's still used).
- `generators::Dict{I, P}`

The generators for each polygon. These are simply the points present in the triangulation. The keys are the vertices, 
while the values are the coordinates; we need to use a `Dict` in case the triangulation is missing points.
- `polygon_points::Vector{P}`

The points defining the vertices of the polygons.
- `polygons::Dict{I, Vector{I}}`

A `Dict` mapping a polygon index (same as a generator index) to the vertices of the polygon. The polygons are given in counter-clockwise order,
and their first and last vertices are equal.
- `circumcenter_to_triangle::Dict{I, T}`

A `Dict` mapping a circumcenter index to the triangle that contains it. The triangles are sorted such that the minimum index is last.
- `triangle_to_circumcenter::Dict{T, I}`

A `Dict` mapping a triangle to its circumcenter index. The triangles are sorted such that the minimum index is last.
- `unbounded_polygons::Set{I}`

A `Set` of the indices of the unbounded polygons.
- `cocircular_circumcenters::S`

A `Set` of the indices of the circumcenters that come from triangles that are cocircular with another triangle's vertices, and adjoin said triangles. 
- `adjacent::Adjacent{I, E}`

An `Adjacent` struct that stores the adjacency information for the polygons, mapping an oriented edge to the polygon to which it belongs. 
- `boundary_polygons::Set{I}`

A `Set` of the indices of the polygons that are on the boundary of the tessellation. Only relevant for clipped tessellations, otherwise see 
`unbounded_polygons`.
"""
struct VoronoiTessellation{Tr<:Triangulation,P,I,T,S,E}
    triangulation::Tr
    generators::Dict{I,P}
    polygon_points::Vector{P} # not guaranteed to be unique if there a circumcenter appears on the boundary, but the vertices (below) do handle this automatically
    polygons::Dict{I,Vector{I}}
    circumcenter_to_triangle::Dict{I,T}
    triangle_to_circumcenter::Dict{T,I}
    unbounded_polygons::Set{I}
    cocircular_circumcenters::S
    adjacent::Adjacent{I,E}
    boundary_polygons::Set{I}
end
for n in fieldnames(VoronoiTessellation)
    name = String(n)
    @eval begin
        @doc """
        get_$($(name))(vor::VoronoiTessellation)

    Returns the $($name) field from the Voronoi tessellation `vor`.
    """ ($(Symbol("get_$n")))(vor::VoronoiTessellation) = vor.$n
    end
end
function Base.show(io::IO, ::MIME"text/plain", vor::VoronoiTessellation)
    println(io, "Voronoi Tessellation.")
    println(io, "    Number of generators: $(num_generators(vor))")
    println(io, "    Number of polygon vertices: $(num_polygon_vertices(vor))")
    print(io, "    Number of polygons: $(num_polygons(vor))")
end

## Type queries 
edge_type(::VoronoiTessellation{Tr,P,I,T,S,E}) where {Tr,P,I,T,S,E} = E
