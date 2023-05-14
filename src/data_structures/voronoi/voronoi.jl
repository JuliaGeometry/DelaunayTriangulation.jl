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

The points defining the vertices of the polygons. The points are not guaranteed to be unique if a circumcenter 
appears on the boundary and you are considering a clipped tessellation.
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
"""
    edge_type(vor::VoronoiTessellation)

Returns the type of the edges in the Voronoi tessellation `vor`.
"""
edge_type(::VoronoiTessellation{Tr,P,I,T,S,E}) where {Tr,P,I,T,S,E} = E

"""
    number_type(vor::VoronoiTessellation)

Returns the type of the numbers used in the Voronoi tessellation `vor`.
"""
number_type(::VoronoiTessellation{Tr,P}) where {Tr,P} = number_type(P)

"""
    integer_type(vor::VoronoiTessellation)

Returns the type of the integers used in the Voronoi tessellation `vor`.
"""
integer_type(::VoronoiTessellation{Tr,P,I}) where {Tr,P,I} = I

"""
    triangle_type(vor::VoronoiTessellation)

Returns the type of the triangles used in the triangulation underlying the Voronoi tessellation `vor`.
"""
triangle_type(::VoronoiTessellation{Tr,P,I,T}) where {Tr,P,I,T} = T

## Getters 
"""
    get_generators(vor::VoronoiTessellation, i...)

Gets the coordinates for the `i`th generator(s) of `vor`.
"""
get_generator(vor::VoronoiTessellation, i) = get_generators(vor)[i]
get_generator(vor::VoronoiTessellation, i::Vararg{I,N}) where {I,N} = ntuple(j -> get_generator(vor, i[j]), Val(N))

"""
    get_polygon_points(vor::VoronoiTessellation, i...)

Gets the coordinates for the `i`th polygon vertex(vertices) of `vor`.
"""
get_polygon_point(vor::VoronoiTessellation, i...) = get_point(get_polygon_points(vor), i...)

"""
    get_polygons(vor::VoronoiTessellation, i)

Gets the vertices of the `i`th polygon of `vor`. The vertices are given in counter-clockwise order, and the first and last vertices are equal.
"""
get_polygon(vor::VoronoiTessellation, i) = get_polygons(vor)[i]

"""
    get_circumcenter_to_triangle(vor::VoronoiTessellation, i)

Gets the triangle that contains the `i`th circumcenter of `vor`.
"""
get_circumcenter_to_triangle(vor::VoronoiTessellation, i) = get_circumcenter_to_triangle(vor)[i]

"""
    get_triangle_to_circumcenter(vor::VoronoiTessellation, T)

Gets the index of the circumcenter of the triangle `T` in `vor`.
"""
function get_triangle_to_circumcenter(vor::VoronoiTessellation, T)
    if !is_ghost_triangle(T)
        return get_triangle_to_circumcenter(vor)[T]
    else # check over all boundary indices, incase we are using a segmented boundary
        tri = get_triangulation(vor)
        u, v, w = indices(T)
        range = get_boundary_index_range(tri, w) # w = ghost index
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
    get_polygon_coordinates(vor::VoronoiTessellation, j, bbox=nothing, bbox_order=nothing)

Gets the coordinates of the `j`th polygon of `vor`. If `bbox` is given, then the coordinates are clipped to the bounding box `bbox`, with the coordinates in `bbox`
given in counter-clockwise order according to `bbox_order`.

See also [`polygon_bounds`](@ref) for `bbox`.
"""
function get_polygon_coordinates(vorn::VoronoiTessellation, j, bbox=nothing, bbox_order=nothing)
    C = get_polygon(vorn, j)
    F = number_type(vorn)
    coords = Vector{NTuple{2,F}}(undef, length(C) - 1)
    for i in firstindex(C):(lastindex(C)-1)
        if !is_boundary_index(C[i])
            coords[i] = get_polygon_point(vorn, C[i])
        else
            ghost_tri = get_circumcenter_to_triangle(vorn, C[i])
            u, v, _ = indices(ghost_tri) # w is the ghost vertex
            p, q = get_generator(vorn, u, v)
            px, py = getxy(p)
            qx, qy = getxy(q)
            m = (px + qx) / 2, (py + qy) / 2
            is_first = is_first_boundary_index(C, i)
            if is_first
                prev_index = previndex_circular(C, i)
                r = get_polygon_point(vorn, C[prev_index])
            else
                next_index = nextindex_circular(C, i)
                r = get_polygon_point(vorn, C[next_index])
            end
            if r == m # It's possible for the circumcenter to lie on the edge and exactly at the midpoint (e.g. [(0.0,1.0),(-1.0,2.0),(-2.0,-1.0)]). In this case, just rotate 
                mx, my = getxy(m)
                dx, dy = qx - mx, qy - my
                rotated_dx, rotated_dy = -dy, dx
                r = mx + rotated_dx, my + rotated_dy
                if is_right(point_position_relative_to_line(p, q, r))
                    rotated_dx, rotated_dy = dy, -dx
                    r = mx + rotated_dx, my + rotated_dy
                end
            end
            r = getxy(r)
            if is_left(point_position_relative_to_line(p, q, r))
                intersection = intersection_of_ray_with_boundary(bbox, bbox_order, m, r)
            else
                intersection = intersection_of_ray_with_boundary(bbox, bbox_order, r, m)
            end
            coords[i] = intersection
        end
    end
    push!(coords, coords[begin])
    return coords
end

"""
    get_neighbouring_boundary_edges(vor::VoronoiTessellation, e)

Gets the boundary edges that are adjacent to the boundary edge `e` in `vor`.
"""
get_neighbouring_boundary_edges(vorn::VoronoiTessellation, e) = get_neighbouring_boundary_edges(get_triangulation(vorn), e)


## Nums 
"""
    num_polygons(vor::VoronoiTessellation)

Gets the number of polygons in `vor`.
"""
num_polygons(vor::VoronoiTessellation) = length(get_polygons(vor))

"""
    num_polygon_vertices(vor::VoronoiTessellation)

Gets the number of vertices in all the polygons in `vor`.
"""
num_polygon_vertices(vor::VoronoiTessellation) = num_points(get_polygon_points(vor))

"""
    num_generators(vor::VoronoiTessellation)

Gets the number of generators in `vor`.
"""
num_generators(vor::VoronoiTessellation) = length(get_generators(vor))

## Adders 
"""
    add_polygon!(vor::VoronoiTessellation, B, i)

Adds the polygon `B` to `vor` at index `i`.
"""
add_polygon!(vor::VoronoiTessellation, B, i) = get_polygons(vor)[i] = B

"""
    push_polygon_point!(vor::VoronoiTessellation, p)
    push_polygon_point!(vor::VoronoiTessellation, x, y)

Pushes the point `p = (x, y)`` to the list of polygon points in `vor`.
"""
push_polygon_point!(vor::VoronoiTessellation, p) = push_point!(get_polygon_points(vor), p)
push_polygon_point!(vor::VoronoiTessellation, x, y) = push_point!(get_polygon_points(vor), x, y)

"""
    add_unbounded_polygon!(vor::VoronoiTessellation, i)

Adds the index `i` to the list of unbounded polygons in `vor`.
"""
add_unbounded_polygon!(vor::VoronoiTessellation, i) = push!(get_unbounded_polygons(vor), i)

"""
    delete_unbounded_polygon!(vor::VoronoiTessellation, i)

Deletes the index `i` from the list of unbounded polygons in `vor`.
"""
delete_unbounded_polygon!(vor::VoronoiTessellation, i) = delete!(get_unbounded_polygons(vor), i)

"""
    add_boundary_polygon!(vorn::VoronoiTessellation, i)

Adds the index `i` to the list of boundary polygons in `vorn`.
"""
add_boundary_polygon!(vorn::VoronoiTessellation, i) = push!(get_boundary_polygons(vorn), i)

## Iterators 
"""
    each_generator(vor::VoronoiTessellation)

Gets an iterator over the indices of the generators in `vor`.
"""
each_generator(vor::VoronoiTessellation) = keys(get_generators(vor))

"""
    each_polygon_vertex(vor::VoronoiTessellation)

Gets an iterator over the indices of the polygon points in `vor`.
"""
each_polygon_vertex(vor::VoronoiTessellation) = eachindex(get_polygon_points(vor))

"""
    each_unbounded_polygon(vor::VoronoiTessellation)

Gets an iterator over the indices of the unbounded polygons in `vor`.
"""
each_unbounded_polygon(vor::VoronoiTessellation) = get_unbounded_polygons(vor)

"""
    each_polygon(vor::VoronoiTessellation)

Gets an iterator over the polygons in `vor`, giving their vertices.
"""
each_polygon(vor::VoronoiTessellation) = values(get_polygons(vor))

"""
    each_polygon_index(vor::VoronoiTessellation)

Gets an iterator over the indices of the polygons in `vor`.
"""
each_polygon_index(vor::VoronoiTessellation) = keys(get_polygons(vor))

## Adjacent 
"""
    get_adjacent(vor::VoronoiTessellation, e)
    get_adjacent(vor::VoronoiTessellation, i, j)

Gets the adjacent polygon to the edge `e` or the edge `(i, j)`.
"""
get_adjacent(vor::VoronoiTessellation, e) = get_adjacent(get_adjacent(vor), e)
get_adjacent(vor::VoronoiTessellation, i, j) = get_adjacent(get_adjacent(vor), i, j)

"""
    add_adjacent!(vor::VoronoiTessellation, e, i)
    add_adjacent!(vor::VoronoiTessellation, i, j, k)

Adds the adjacent polygon `i` to the edge `e` or the edge `(i, j)` to the polygon `k`.
"""
add_adjacent!(vor::VoronoiTessellation, e, i) = add_adjacent!(get_adjacent(vor), e, i)
add_adjacent!(vor::VoronoiTessellation, i, j, k) = add_adjacent!(get_adjacent(vor), i, j, k)

"""
    delete_adjacent!(vor::VoronoiTessellation, e, i)
    delete_adjacent!(vor::VoronoiTessellation, i, j, k)

Deletes the adjacent polygon `i` from the edge `e` or the edge `(i, j)` from the polygon `k`.
"""
delete_adjacent!(vor::VoronoiTessellation, e, i) = delete_adjacent!(get_adjacent(vor), e, i)
delete_adjacent!(vor::VoronoiTessellation, i, j, k) = delete_adjacent!(get_adjacent(vor), i, j, k)

"""
    delete_polygon_adjacent!(vorn::VoronoiTessellation, polygon)

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
    return nothing
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
    return nothing
end

## Features

"""
    polygon_features(vor::VoronoiTessellation, i)

Gets the area and centroid of the polygon with index `i` in `vor`.
"""
function polygon_features(vor::VoronoiTessellation, i)
    polygon = get_polygon(vor, i)
    if any(is_boundary_index, polygon)
        F = number_type(vor)
        return (typemax(F), (typemax(F), typemax(F)))
    end
    return polygon_features(get_polygon_points(vor), polygon)
end

"""
    get_area(vor::VoronoiTessellation, i)

Gets the area of the polygon with index `i` in `vor`.
"""
get_area(vor::VoronoiTessellation, i) = polygon_features(vor, i)[1]

"""
    get_centroid(vor::VoronoiTessellation, i)

Gets the centroid of the polygon with index `i` in `vor`.
"""
get_centroid(vor::VoronoiTessellation, i) = polygon_features(vor, i)[2]

"""
    polygon_bounds(vorn::VoronoiTessellation, unbounded_extension_factor=0.0)

Gets the bounding box of the polygons in `vorn`. If `unbounded_extension_factor` is positive, the bounding box is extended by this factor in each direction,
proportional to the width of each axis.
"""
function polygon_bounds(vorn::VoronoiTessellation, unbounded_extension_factor=0.0)
    F = number_type(vorn)
    xmin = typemax(F)
    xmax = typemin(F)
    ymin = typemax(F)
    ymax = typemin(F)
    for i in each_polygon_vertex(vorn)
        x, y = getxy(get_polygon_point(vorn, i))
        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
    end
    for i in each_generator(vorn)
        x, y = getxy(get_generator(vorn, i))
        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
    end
    xmin -= unbounded_extension_factor * (xmax - xmin)
    xmax += unbounded_extension_factor * (xmax - xmin)
    ymin -= unbounded_extension_factor * (ymax - ymin)
    ymax += unbounded_extension_factor * (ymax - ymin)
    return xmin, xmax, ymin, ymax
end

"""
    jump_and_march(vor::VoronoiTessellation, q; kwargs...)

Finds the polygon containing the point `q` in the Voronoi tessellation `vor`. The keyword arguments are passed to `jump_and_march`.
"""
function jump_and_march(vor::VoronoiTessellation, q; kwargs...)
    #=
    V = jump_and_march(get_triangulation(vor), q; kwargs...)
    qx, qy = getxy(q)
    V = rotate_triangle_to_standard_form(V)
    i, j, k = V
    a, b = get_generator(vor, i, j)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    daq = (qx - ax)^2 + (qy - ay)^2
    dbq = (qx - bx)^2 + (qy - by)^2
    if !is_boundary_index(k)
        c = get_generator(vor, k)
        cx, cy = getxy(c)
        dcq = (qx - cx)^2 + (qy - cy)^2
    else
        dcq = typemax(number_type(vor))
    end
    min_dist = min(daq, dbq, dcq)
    if min_dist == daq
        return i
    elseif min_dist == dbq
        return j
    else
        return k
    end
    =#
    return jump_to_voronoi_polygon(get_triangulation(vor), q; kwargs...)
end

function jump_to_voronoi_polygon(tri::Triangulation, q; kwargs...)
    V = jump_and_march(tri, q; kwargs...)
    qx, qy = getxy(q)
    V = rotate_triangle_to_standard_form(V)
    i, j, k = V
    a, b = get_point(tri, i, j)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    daq = (qx - ax)^2 + (qy - ay)^2
    dbq = (qx - bx)^2 + (qy - by)^2
    if !is_boundary_index(k)
        c = get_point(tri, k)
        cx, cy = getxy(c)
        dcq = (qx - cx)^2 + (qy - cy)^2
    else
        dcq = typemax(number_type(tri))
    end
    min_dist = min(daq, dbq, dcq)
    if min_dist == daq
        return i
    elseif min_dist == dbq
        return j
    else
        return k
    end
end

## Utils
"""
    get_surrounding_polygon(vor::VoronoiTessellation, i)

Gets the polygon surrounding the generator with index `i` in `vor`. You shouldn't need 
to use this, see [`get_polygon`](@ref) instead.
"""
function get_surrounding_polygon(vor::VoronoiTessellation, i)
    tri = get_triangulation(vor)
    S = get_surrounding_polygon(tri, i)
    push!(S, S[begin])
    return S
end

"""
    convert_to_boundary_edge(vorn::VoronoiTessellation, e)

Converts the edge `e` in the triangulation of `vorn` to a boundary edge so that `get_adjacent(vorn, e)` is a boundary index.
"""
convert_to_boundary_edge(vorn::VoronoiTessellation, e) = convert_to_boundary_edge(get_triangulation(vorn), e)
