"""
    VoronoiTessellation{Tr<:Triangulation,P,I,T,S,E}

Struct for representing a Voronoi tessellation.

See also [`voronoi`](@ref).

# Fields 
- `triangulation::Tr`: The underlying triangulation. The tessellation is dual to this triangulation, 
   although if the underlying triangulation is constrained then this is no longer the case (but it is 
   still used).
- `generators::Dict{I,P}`: A `Dict` that maps vertices of generators to coordinates. These are simply the 
   points present in the triangulation. A `Dict` is needed in case the triangulation is missing some points.
- `polygon_points::Vector{P}`: The points defining the coordinates of the polygons. The points are not guaranteed
   to be unique if a circumcenter appears on the boundary or you are considering a clipped tessellation. (See also [`get_polygon_coordinates`](@ref).)
- `polygons::Dict{I,Vector{I}}`: A `Dict` mapping polygon indices (which is the same as a generator vertex) to the vertices of a polygon. The polygons are given in 
   counter-clockwise order and the first and last vertices are equal.
- `circumcenter_to_triangle::Dict{I,T}`: A `Dict` mapping a circumcenter index to the triangle that contains it. The triangles are sorted such that the 
   minimum vertex is last. 
- `triangle_to_circumcenter::Dict{T,I}`: A `Dict` mapping a triangle to its circumcenter index. The triangles are sorted such that the minimum vertex is last. 
- `unbounded_polygons::Set{I}`: A `Set` of indices of the unbounded polygons. 
- `cocircular_circumcenters::S`: A `Set` of indices of circumcenters that come from triangles that are cocircular with another triangle's vertices, and adjoin said triangles. 
- `adjacent::Adjacent{I,E}`: The adjacent map. This maps an oriented edge to the polygon that it belongs to.
- `boundary_polygons::Set{I}`: A `Set` of indices of the polygons that are on the boundary of the tessellation. Only relevant for clipped 
   tessellations, otherwise see `unbounded_polygons`.
"""
struct VoronoiTessellation{Tr<:Triangulation,P,I,T,S,E}
    triangulation::Tr
    generators::Dict{I,P}
    polygon_points::Vector{P} # not guaranteed to be unique if a circumcenter appears on the boundary, but the vertices (below) do handle this automatically
    polygons::Dict{I,Vector{I}}
    circumcenter_to_triangle::Dict{I,T}
    triangle_to_circumcenter::Dict{T,I}
    unbounded_polygons::Set{I}
    cocircular_circumcenters::S
    adjacent::Adjacent{I,E}
    boundary_polygons::Set{I}
end

"""
    get_triangulation(vorn::VoronoiTessellation) -> Triangulation 

Gets the underlying triangulation of the Voronoi tessellation `vorn`.
"""
get_triangulation(vorn::VoronoiTessellation) = vorn.triangulation

"""
    get_generators(vorn::VoronoiTessellation) -> Dict{Vertex,Point}

Gets the generators of the Voronoi tessellation `vorn` as a `Dict`, mapping vertices to their coordinates. These 
coordinates are given as `Tuple`s of the form `(x, y)`.
"""
get_generators(vorn::VoronoiTessellation) = vorn.generators

"""
    get_polygon_points(vorn::VoronoiTessellation) -> Vector{Point}

Gets the polygon points of the Voronoi tessellation `vorn`. These are the vertices of the Voronoi polygons, and are
given as `Tuple`s of the form `(x, y)`.
"""
get_polygon_points(vorn::VoronoiTessellation) = vorn.polygon_points

"""
    get_polygons(vorn::VoronoiTessellation) -> Dict{Index,Vector{Vertex}}

Gets the polygons of the Voronoi tessellation `vorn` as a `Dict`, mapping polygon indices to their vertices, 
where the vertices refer to points in `get_polygon_points(vorn)`. The vertices are given in counter-clockwise order 
and the first and last vertices are equal. 
"""
get_polygons(vorn::VoronoiTessellation) = vorn.polygons

"""
    get_circumcenter_to_triangle(vorn::VoronoiTessellation) -> Dict{Index,Triangle}    

Gets the circumcenters of the Voronoi tessellation `vorn` as a `Dict`, mapping circumcenter indices to their
corresponding triangles. The triangles are sorted so that the minimum vertex is last.
"""
get_circumcenter_to_triangle(vorn::VoronoiTessellation) = vorn.circumcenter_to_triangle

"""
    get_triangle_to_circumcenter(vorn::VoronoiTessellation) -> Dict{Triangle,Index}

Gets the triangles of the Voronoi tessellation `vorn` as a `Dict`, mapping triangle indices to their corresponding
circumcenters. The circumcenters are given as their vertices, referring to points in `get_polygon_points(vorn)`.
"""
get_triangle_to_circumcenter(vorn::VoronoiTessellation) = vorn.triangle_to_circumcenter

"""
    get_unbounded_polygons(vorn::VoronoiTessellation) -> Set{Index}

Gets the unbounded polygons of the Voronoi tessellation `vorn` as a `Set` of polygon indices.
"""
get_unbounded_polygons(vorn::VoronoiTessellation) = vorn.unbounded_polygons

"""
    get_cocircular_circumcenters(vorn::VoronoiTessellation) -> Set

Gets the cocircular circumcenters of the Voronoi tessellation `vorn` as a `Set` of circumcenter indices. These are circumcenters 
that come from triangles that are cocircular with another adjoining triangle.
"""
get_cocircular_circumcenters(vorn::VoronoiTessellation) = vorn.cocircular_circumcenters

"""
    get_adjacent(vorn::VoronoiTessellation) -> Adjacent{Index,Edge}

Gets the adjacency information of the Voronoi tessellation `vorn` as an `Adjacent` object. This object maps oriented edges 
to the polygons that they belong to.
"""
get_adjacent(vorn::VoronoiTessellation) = vorn.adjacent

"""
    get_boundary_polygons(vorn::VoronoiTessellation) -> Set{Index}

Gets the boundary polygons of the Voronoi tessellation `vorn` as a `Set` of polygon indices.
"""
get_boundary_polygons(vorn::VoronoiTessellation) = vorn.boundary_polygons

function Base.show(io::IO, ::MIME"text/plain", vor::VoronoiTessellation)
    println(io, "Voronoi Tessellation.")
    println(io, "    Number of generators: $(num_generators(vor))")
    println(io, "    Number of polygon vertices: $(num_polygon_vertices(vor))")
    print(io, "    Number of polygons: $(num_polygons(vor))")
end

"""
    edge_type(vorn::VoronoiTessellation) -> DataType 

Type used for representing individual edges in the Voronoi tessellation.
"""
edge_type(::VoronoiTessellation{Tr,P,I,T,S,E}) where {Tr,P,I,T,S,E} = E

"""
    number_type(vorn::VoronoiTessellation) -> DataType

Type used for representing individual coordinates in the Voronoi tessellation.
"""
number_type(::VoronoiTessellation{Tr,P}) where {Tr,P} = number_type(P)

"""
    integer_type(vorn::VoronoiTessellation) -> DataType

Type used for representing indices in the Voronoi tessellation.
"""
integer_type(::VoronoiTessellation{Tr,P,I}) where {Tr,P,I} = I

"""
    triangle_type(vorn::VoronoiTessellation) -> DataType

Type used for representing individual triangles in the Voronoi tessellation.
"""
triangle_type(::VoronoiTessellation{Tr,P,I,T}) where {Tr,P,I,T} = T

"""
    get_generator(vor::VoronoiTessellation, i) -> NTuple{2, Number}
    get_generator(vor::VoronoiTessellation, i...) -> NTuple{length(i), NTuple{2, Number}}

Gets the coordinates for the generators `i...`, returned as `Tuple`s of the form `(x, y)` for each generator.
"""
get_generator(vor::VoronoiTessellation, i) = get_generators(vor)[i]
get_generator(vor::VoronoiTessellation, i, j::Vararg{I, N}) where {I, N} = (get_generator(vor, i), ntuple(k -> get_generator(vor, j[k]), Val(N))...)

"""
    get_polygon_point(vor::VoronoiTessellation, i) -> NTuple{2, Number}
    get_polygon_point(vor::VoronoiTessellation, i...) -> NTuple{length(i), NTuple{2, Number}}

Gets the coordinates corresponding to the vertices `i...` of the polygons, returned as `Tuple`s of the form `(x, y)` for each vertex.
"""
get_polygon_point(vor::VoronoiTessellation, i...) = get_point(get_polygon_points(vor), i...)

"""
    get_polygon(vor::VoronoiTessellation, i) -> Vector{Vertex}

Gets the vector of vertices corresponding to the `i`th polygon, given in counter-clockwise order and 
with the first and last vertices equal. To obtain the coordinates, see [`get_polygon_point`](@ref).
"""
get_polygon(vor::VoronoiTessellation, i) = get_polygons(vor)[i]

"""
    get_circumcenter_to_triangle(vor::VoronoiTessellation, i) -> Triangle

Gets the triangle associated with the `i`th circumcenter. The triangle is sorted so that the minimum vertex is last.
"""
get_circumcenter_to_triangle(vor::VoronoiTessellation, i) = get_circumcenter_to_triangle(vor)[i]

"""
    get_triangle_to_circumcenter(vor::VoronoiTessellation, T) -> Index

Gets the circumcenter index associated with the triangle `T`. The triangle should be sorted so that the minimum vertex is last.
"""
function get_triangle_to_circumcenter(vor::VoronoiTessellation, T)
    if !is_ghost_triangle(T)
        return get_triangle_to_circumcenter(vor)[T]
    else # check over all boundary indices, incase we are using a segmented boundary
        tri = get_triangulation(vor)
        u, v, w = triangle_vertices(T)
        range = get_ghost_vertex_range(tri, w) # w = ghost vertex
        V = triangle_type(tri)
        dict = get_triangle_to_circumcenter(vor)
        for j in range
            T = construct_triangle(V, u, v, j)
            haskey(dict, T) && return dict[T]
        end
        throw(KeyError(T))
    end
end

"""
    get_neighbouring_boundary_edges(vorn::VoronoiTessellation, e) -> (Edge, Edge)

Given a boundary edge `e`, returns the edges left and right of `e`. 
"""
get_neighbouring_boundary_edges(vorn::VoronoiTessellation, e) = get_neighbouring_boundary_edges(get_triangulation(vorn), e)

"""
    num_polygons(vor::VoronoiTessellation) -> Integer

Returns the number of polygons in the Voronoi tessellation `vor`.
"""
num_polygons(vor::VoronoiTessellation) = length(get_polygons(vor))

"""
    num_polygon_vertices(vor::VoronoiTessellation) -> Integer

Returns the number of polygon vertices in the Voronoi tessellation `vor`. This might include duplicate vertices if `get_polygon_points(vor)` has duplicates.
"""
num_polygon_vertices(vor::VoronoiTessellation) = num_points(get_polygon_points(vor))

"""
    num_generators(vor::VoronoiTessellation) -> Integer

Returns the number of generators in the Voronoi tessellation `vor`.
"""
num_generators(vor::VoronoiTessellation) = length(get_generators(vor))

"""
    add_polygon!(vor::VoronoiTessellation, B, i)

Adds, or replaces, the polygon associated with the index `i` with `B`. `B` should be a counter-clockwise 
sequence of vertices, with `B[begin] == B[end]`.
"""
add_polygon!(vor::VoronoiTessellation, B, i) = get_polygons(vor)[i] = B

"""
    push_polygon_point!(vor::VoronoiTessellation, p)
    push_polygon_point!(vor::VoronoiTessellation, x, y)

Adds the point `p = (x, y)` into the vector of polygon points of `vor`.
"""
push_polygon_point!(vor::VoronoiTessellation, p) = push_point!(get_polygon_points(vor), p)
push_polygon_point!(vor::VoronoiTessellation, x, y) = push_point!(get_polygon_points(vor), x, y)

"""
    add_unbounded_polygon!(vor::VoronoiTessellation, i)

Adds the index `i` to the set of unbounded polygons of `vor`.
"""
add_unbounded_polygon!(vor::VoronoiTessellation, i) = push!(get_unbounded_polygons(vor), i)

"""
    delete_unbounded_polygon!(vor::VoronoiTessellation, i)

Deletes the index `i` from the set of unbounded polygons of `vor`.
"""
delete_unbounded_polygon!(vor::VoronoiTessellation, i) = delete!(get_unbounded_polygons(vor), i)

"""
    add_boundary_polygon!(vor::VoronoiTessellation, i)

Adds the index `i` to the set of boundary polygons of `vor`.
"""
add_boundary_polygon!(vorn::VoronoiTessellation, i) = push!(get_boundary_polygons(vorn), i)

"""
    each_boundary_polygon(vorn::VoronoiTessellation) -> KeySet 

Returns an iterator over the boundary polygon indices of `vorn`.
"""
each_generator(vor::VoronoiTessellation) = keys(get_generators(vor))

"""
    each_polygon_vertex(vor::VoronoiTessellation) -> UnitRange

Returns an iterator over each polygon point index of `vor`.
"""
each_polygon_vertex(vor::VoronoiTessellation) = eachindex(get_polygon_points(vor))

"""
    each_polygon(vor::VoronoiTessellation) -> Set 

Returns an iterator over the polygon indices of `vor`.
"""
each_unbounded_polygon(vor::VoronoiTessellation) = get_unbounded_polygons(vor)

"""
    each_polygon(vor::VoronoiTessellation) -> ValueIterator

Returns an iterator over each set of polygon vertices of `vor`.
"""
each_polygon(vor::VoronoiTessellation) = values(get_polygons(vor))

"""
    each_polygon_index(vor::VoronoiTessellation) -> KeySet 

Returns an iterator over the polygon indices of `vor`.
"""
each_polygon_index(vor::VoronoiTessellation) = keys(get_polygons(vor))

"""
    get_adjacent(vor::VoronoiTessellation, ij) -> Index 
    get_adjacent(vor::VoronoiTessellation, i, j) -> Index

Gets the polygon index associated with the oriented edge `ij` in the Voronoi tessellation `vor`.
"""
get_adjacent(vor::VoronoiTessellation, e) = get_adjacent(get_adjacent(vor), e)
get_adjacent(vor::VoronoiTessellation, i, j) = get_adjacent(get_adjacent(vor), i, j)

"""
    add_adjacent!(vor::VoronoiTessellation, ij, k)
    add_adjacent!(vor::VoronoiTessellation, i, j, k)

Adds the adjacency relationship `(i, j) ⟹ k` between the oriented edge `(i, j)` and polygon index `k` to the Voronoi tessellation `vor`.
"""
add_adjacent!(vor::VoronoiTessellation, e, i) = add_adjacent!(get_adjacent(vor), e, i)
add_adjacent!(vor::VoronoiTessellation, i, j, k) = add_adjacent!(get_adjacent(vor), i, j, k)

"""
    delete_adjacent!(vor::VoronoiTessellation, ij)
    delete_adjacent!(vor::VoronoiTessellation, i, j)

Deletes the edge `(i, j)` from the adjacency map of `vor`.
"""
delete_adjacent!(vor::VoronoiTessellation, e) = delete_adjacent!(get_adjacent(vor), e)
delete_adjacent!(vor::VoronoiTessellation, i, j) = delete_adjacent!(get_adjacent(vor), i, j)

"""
    delete_polygon!(vor::VoronoiTessellation, polygon)

Deletes the adjacent polygons of the boundary edges of the polygon `polygon` from the adjacency list.
"""
function delete_polygon_adjacent!(vorn::VoronoiTessellation, polygon)
    polygon_vertices = get_polygon(vorn, polygon)
    ne = num_boundary_edges(polygon_vertices)
    for ℓ in 1:ne
        u = get_boundary_nodes(polygon_vertices, ℓ)
        v = get_boundary_nodes(polygon_vertices, ℓ + 1)
        delete_adjacent!(vorn, u, v)
    end
    return vorn
end

"""
    add_polygon_adjacent!(vorn::VoronoiTessellation, polygon)

Adds the adjacent polygons of the boundary edges of the polygon `polygon` to the adjacency list.
"""
function add_polygon_adjacent!(vorn::VoronoiTessellation, polygon)
    polygon_vertices = get_polygon(vorn, polygon)
    ne = num_boundary_edges(polygon_vertices)
    for ℓ in 1:ne
        u = get_boundary_nodes(polygon_vertices, ℓ)
        v = get_boundary_nodes(polygon_vertices, ℓ + 1)
        add_adjacent!(vorn, u, v, polygon)
    end
    return vorn
end

"""
    polygon_features(vor::VoronoiTessellation, i) -> (Number, NTuple{2, Number})

Gets the area and centroid of the polygon with index `i` in `vor`.
"""
function polygon_features(vor::VoronoiTessellation, i)
    polygon = get_polygon(vor, i)
    if i ∈ get_unbounded_polygons(vor)
        F = number_type(vor)
        return (typemax(F), (typemax(F), typemax(F)))
    end
    return polygon_features(get_polygon_points(vor), polygon)
end

"""
    get_area(vor::VoronoiTessellation, i) -> Number

Gets the area of the `i`th Voronoi polygon.
"""
get_area(vor::VoronoiTessellation, i) = polygon_features(vor, i)[1]

"""
    get_centroid(vor::VoronoiTessellation, i) -> NTuple{2, Number}

Gets the centroid of the `i`th Voronoi polygon, given as a `Tuple` of the form `(x, y)`.
"""
get_centroid(vor::VoronoiTessellation, i) = polygon_features(vor, i)[2]

"""
    polygon_bounds(vorn::VoronoiTessellation, unbounded_extension_factor=0.0; include_polygon_vertices=true) -> (Number, Number, Number, Number)

Gets the bounding box for the Voronoi tessellation `vorn`.

# Arguments 
- `vorn::VoronoiTessellation`: The Voronoi tessellation.
- `unbounded_extension_factor=0.0`: The factor by which to extend the bounding box for unbounded polygons.

# Keyword Arguments 
- `include_polygon_vertices=true`: If `true`, then the bounding box will also include the polygon vertices. Otherwise, only the generators are included.

# Output 
- `xmin`: Given by `xmin′ - unbounded_extension_factor * (xmin′ - xmin′)`, where `xmin′` is the original minimum `x`-coordinate of the computed bounding box and similarly for `xmax′`.
- `xmax`: Given by `xmax′ + unbounded_extension_factor * (xmax′ - xmax′)`, where `xmax′` is the original maximum `x`-coordinate of the computed bounding box and similarly for `xmin′`.
- `ymin`: Given by `ymin′ - unbounded_extension_factor * (ymin′ - ymin′)`, where `ymin′` is the original minimum `y`-coordinate of the computed bounding box and similarly for `ymax′`.
- `ymax`: Given by `ymax′ + unbounded_extension_factor * (ymax′ - ymax′)`, where `ymax′` is the original maximum `y`-coordinate of the computed bounding box and similarly for `ymin′`.
"""
function polygon_bounds(vorn::VoronoiTessellation, unbounded_extension_factor=0.0; include_polygon_vertices=true)
    F = number_type(vorn)
    xmin = typemax(F)
    xmax = typemin(F)
    ymin = typemax(F)
    ymax = typemin(F)
    if include_polygon_vertices
        for i in each_polygon_vertex(vorn)
            x, y = getxy(get_polygon_point(vorn, i))
            xmin = min(xmin, x)
            xmax = max(xmax, x)
            ymin = min(ymin, y)
            ymax = max(ymax, y)
        end
    end
    for i in each_generator(vorn)
        x, y = getxy(get_generator(vorn, i))
        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
    end
    _xmin = xmin - unbounded_extension_factor * (xmax - xmin)
    _xmax = xmax + unbounded_extension_factor * (xmax - xmin)
    _ymin = ymin - unbounded_extension_factor * (ymax - ymin)
    _ymax = ymax + unbounded_extension_factor * (ymax - ymin)
    return _xmin, _xmax, _ymin, _ymax
end

"""
    get_surrounding_polygon(vor::VoronoiTessellation, i) -> Vector{Vertex}

Gets the polygon surrounding the generator with index `i` in `vor`. 

You shouldn't need to use this, see [`get_polygon`](@ref) instead.
"""
function get_surrounding_polygon(vor::VoronoiTessellation, i)
    tri = get_triangulation(vor)
    S = get_surrounding_polygon(tri, i)
    push!(S, S[begin])
    return S
end

"""
    convert_to_edge_adjoining_ghost_vertex(vorn::VoronoiTessellation, e) -> Edge

Returns the edge `e` if it is not a boundary edge, and the edge `reverse(e)` if it is a boundary edge. 

See also [`is_boundary_edge`](@ref).
"""
convert_to_edge_adjoining_ghost_vertex(vorn::VoronoiTessellation, e) = convert_to_edge_adjoining_ghost_vertex(get_triangulation(vorn), e)
