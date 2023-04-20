"""
    VoronoiTessellation{P,I,E,Es,T,Tr<:Triangulation}

A structure to store the Voronoi tessellation of a set of points.

# Generators 
- `generators::Vector{P}`

The generators of the Voronoi tessellation.
- `cell_points::Vector{P}`

The points that define the boundaries of the Voronoi cells. This will include both the circumcenters and any 
points coming from the obstacles.
- `cells::Vector{VoronoiCell{I}}`

The Voronoi cells.
- `adjacent_cell::Adjacent{I, E}`

Map that takes an edge to its adjacent cell.
- `adjacent_vertex::Adjacent{I, E}`

Map that takes an edge to its adjacent vertex.
- `vertex_to_cells::Adjacent2Vertex{I, Set{I}, I}`

Map that takes a vertex to its adjacent cells.
- `vertex_to_vertices::Graph{I}`

Graph that takes a vertex to all other vertices that share an edge with it.
- `cell_to_cells::Graph{I}`

Graph that takes a cell to all other cells that share an edge with it.
- `obstacles::Es`

The obstacles that define the boundary of the Voronoi tessellation or any other 
blocked visibility.
- `boundary_cells::Set{I}`

The cells that are on the boundary of the Voronoi tessellation.
- `unbounded_cells::Set{I}`

The cells that are unbounded.
- `circumcenter_to_triangle::Dict{I, T}`

Map that takes a circumcenter to the triangle that it is the circumcenter of.
- `triangle_to_circumcenter::Dict{T, I}`

Map that takes a triangle to the circumcenter of that triangle.
- `triangulation::Tr`

The triangulation that the Voronoi tessellation is based on.

!!! note 

    While, in theory, the Voronoi tessellation could be accessed directly from a 
    triangulation and a set of circumcenters, we choose a design that keeps 
    the two as separate as possible to make working with the tessellation simpler. 
    Moreover, this allows for us to mutate the points as needed to work with 
    visibility constraints. We do still provide the `triangulation` field,
    though, as it can help with point location queries. (The main point 
    here being that, while `triangulation` is provided as a field,
    typically you shouldn't need to actually interact with it.)
"""
struct VoronoiTessellation{P,I,E,Es,T,Tr<:Triangulation}
    generators::Vector{P}
    cell_points::Vector{P}
    cells::Vector{VoronoiCell{I}}
    adjacent_cell::Adjacent{I,E}
    adjacent_vertex::Adjacent{I,E}
    vertex_to_cells::Adjacent2Vertex{I,Set{I},I}
    vertex_to_vertices::Graph{I}
    cell_to_cells::Graph{I}
    obstacles::Es
    boundary_cells::Set{I}
    unbounded_cells::Set{I}
    circumcenter_to_triangle::Dict{I,T}
    triangle_to_circumcenter::Dict{T,I}
    triangulation::Tr
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
    println(io, "    Constrained: $(is_constrained(vor))")
    println(io, "    Number of generators: $(num_generators(vor))")
    println(io, "    Number of polygon vertices: $(num_polygon_vertices(vor))")
    println(io, "    Number of cells: $(num_cells(vor))")
    println(io, "    Number of edges: $(num_edges(vor))")
    println(io, "    Number of obstacles: $(num_obstacles(vor))")
    println(io, "    Number of boundary cells: $(length(get_boundary_cells(vor)))")
    print(io, "    Number of unbounded cells: $(length(get_unbounded_cells(vor)))")
end

num_generators(vor::VoronoiTessellation) = length(get_generators(vor))
num_cells(vor::VoronoiTessellation) = length(get_cells(vor))
num_edges(vor::VoronoiTessellation) = num_edges(get_vertex_to_vertices(vor))
num_polygon_vertices(vor::VoronoiTessellation) = length(get_cell_points(vor))
num_obstacles(vor::VoronoiTessellation) = num_edges(get_obstacles(vor))
is_constrained(vor::VoronoiTessellation) = !is_empty(get_obstacles(vor))
num_boundary_cells(vor::VoronoiTessellation) = length(get_boundary_cells(vor))
num_unbounded_cells(vor::VoronoiTessellation) = length(get_unbounded_cells(vor))
integer_type(::VoronoiTessellation{P,I}) where {P,I} = I
number_type(::VoronoiTessellation{P}) where {P} = number_type(P)
triangle_type(::VoronoiTessellation{P,I,E,Es,T}) where {P,I,E,Es,T} = T
get_cell(vor::VoronoiTessellation, i) = get_cells(vor)[i]
get_circumcenter_to_triangle(vor::VoronoiTessellation, i) = get_circumcenter_to_triangle(vor)[i]
get_triangle_to_circumcenter(vor::VoronoiTessellation, i) = get_triangle_to_circumcenter(vor)[i]
each_generator(vor::VoronoiTessellation) = eachindex(get_generators(vor))
get_generator(vor::VoronoiTessellation, i) = getxy(get_generators(vor)[i])
get_generator(vor::VoronoiTessellation, i::Vararg{I,N}) where {I,N} = ntuple(j -> get_generator(vor, i[j]), Val(N))
get_cell_point(vor::VoronoiTessellation, i) = getxy(get_cell_points(vor)[i])
get_cell_point(vor::VoronoiTessellation, i::Vararg{I,N}) where {I,N} = ntuple(j -> get_cell_point(vor, i[j]), Val(N))
each_cell(vor::VoronoiTessellation) = eachindex(get_cells(vor))
each_obstacle(vor::VoronoiTessellation) = each_edge(get_obstacles(vor))
each_cell_point(vor::VoronoiTessellation) = eachindex(get_cell_points(vor))
function get_surrounding_polygon(vorn::VoronoiTessellation, i)
    tri = get_triangulation(vorn)
    S = get_surrounding_polygon(tri, i)
    push!(S, S[begin])
    return S
end
