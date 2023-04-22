## Definition

"""
    VoronoiTessellation{Tr<:Triangulation, P, I, T, S, E}

A Voronoi tessellation, dual to a triangulation.

!!! warning 

    If the triangulation is constrained, it is not guaranteed that 
    the Voronoi tessellation is dual to the triangulation. This duality 
    is only guaranteed for unconstrained triangulations, and only for point 
    sets with no four points that are cocircular.

# Fields
- `triangulation::Tr`: The triangulation that the Voronoi tessellation is based on.
- `generators::Dict{I, P}`: The points that define the generators of the Voronoi tessellation, with indices mapping to coordinates. These are the same as the points of the triangulation, but dealiased to allow for mutation (e.g. if using Lloyd's algorithm). We use a `Dict` since the point set in the triangulation might not have every point included in the triangulation, so a `Dict` helps preserve the generator indices (rather than a vector which could shift the indices).
- `polygon_points::Vector{P}`: The points that define points on the boundaries of the Voronoi polygons. This will include both the circumcenters and any points coming from intersections with the boundary (if the tessellation is clipped).
- `polygons::Dict{I, Vector{I}}`: The Voronoi polygons. Each polygon is a polygon, defined by a list of indices into `polygon_point`. Each vector is a circular vector (see [`is_circular`](@ref)). The keys of the `Dict` enumerate the polygon, with the key `i` corresponding to the `i`th point (the generator) of the triangulation.
- `circumcenter_to_triangle::Dict{I,T}`: Map that takes a circumcenter to the triangle that it is the circumcenter of.
- `triangle_to_circumcenter::Dict{T,I}`: Map that takes a triangle to the circumcenter of that triangle.
- `unbounded_polygons::Set{I}`: The polygons that are unbounded.
- `cocircular_circumcenters:S`: Indices for the circumcenters that are given by more than just one triangle.
- `adjacent::Adjacent{I,E}`: The adjacency structure for the Voronoi tessellation, mapping polygonal edges to the cells they belong to.

!!! warning "Degeneracy"

        When there are adjoining triangles that have the same circumcenter, the `triangle_to_circumcenter` and 
        `circumcenter_to_triangle` dictionaries will only store one of the triangles.

        Moreover, the Voronoi tessellation is only dual to the triangulation if there are no four points that are cocircular.
        While we do handle cocircular subsets automatically in the construction of the tessellation, we do not 
        adjust the triangulation accordingly so that triangles arising from adjoining triangles with the 
        same circumcenter are morphed into a polyhedron (returning a Delaunay subdivision). If we did make 
        this adjustment, it would indeed be dual in general. However, this would be a significant change to
        the triangulation code, and would require a lot of additional code to handle the morphing of the
        triangulation. Moreover, the data structures we use are not adequate for representing a combination of 
        triangles and polyhedrons. (For more detail, see the discussion around Fig 2.3 in the book by
        Cheng, Dey, and Shewchuk, and also the discussion on p. 156 right before Section 7.2.)
"""
struct VoronoiTessellation{Tr<:Triangulation,P,I,T,S,E}
    triangulation::Tr
    generators::Dict{I,P}
    polygon_points::Vector{P}
    polygons::Dict{I,Vector{I}}
    circumcenter_to_triangle::Dict{I,T}
    triangle_to_circumcenter::Dict{T,I}
    unbounded_polygons::Set{I}
    cocircular_circumcenters::S
    adjacent::Adjacent{I,E}
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

num_polygons(vor::VoronoiTessellation) = length(get_polygons(vor))
num_polygon_vertices(vor::VoronoiTessellation) = num_points(get_polygon_points(vor))
get_generator(vor::VoronoiTessellation, i) = get_generators(vor)[i]
get_generator(vor::VoronoiTessellation, i::Vararg{I,N}) where {I,N} = ntuple(j -> get_generator(vor, i[j]), Val(N))
get_polygon_point(vor::VoronoiTessellation, i...) = get_point(get_polygon_points(vor), i...)
get_polygon(vor::VoronoiTessellation, i) = get_polygons(vor)[i]
get_circumcenter_to_triangle(vor::VoronoiTessellation, i) = get_circumcenter_to_triangle(vor)[i]
get_triangle_to_circumcenter(vor::VoronoiTessellation, i) = get_triangle_to_circumcenter(vor)[i]
number_type(vor::VoronoiTessellation) = number_type(get_polygon_points(vor))
num_generators(vor::VoronoiTessellation) = length(get_generators(vor))
integer_type(::VoronoiTessellation{Tr,P,I}) where {Tr,P,I} = I
number_type(::VoronoiTessellation{Tr,P}) where {Tr,P} = number_type(P)
triangle_type(::VoronoiTessellation{Tr,P,I,T}) where {Tr,P,I,T} = T
each_generator(vor::VoronoiTessellation) = keys(get_generators(vor))
each_polygon_vertex(vor::VoronoiTessellation) = eachindex(get_polygon_points(vor))
each_unbounded_polygon(vor::VoronoiTessellation) = get_unbounded_polygons(vor)
each_polygon(vor::VoronoiTessellation) = values(get_polygons(vor))
contains_boundary_edge(vor::VoronoiTessellation, e) = contains_boundary_edge(get_triangulation(vor), e)
contains_boundary_edge(vor::VoronoiTessellation, i, j) = contains_boundary_edge(get_triangulation(vor), i, j)
is_encroached(vor::VoronoiTessellation, e) = is_encroached(get_triangulation(vor), e)
is_encroached(vor::VoronoiTessellation, i, j) =
    let tri = get_triangulation(vor)
        is_encroached(tri, construct_edge(edge_type(tri), i, j))
    end
push_polygon_point!(vor::VoronoiTessellation, p) = push_point!(get_polygon_points(vor), p)
push_polygon_point!(vor::VoronoiTessellation, x, y) = push_point!(get_polygon_points(vor), x, y)
each_polygon_index(vor::VoronoiTessellation) = keys(get_polygons(vor))
get_adjacent(vor::VoronoiTessellation, e) = get_adjacent(get_adjacent(vor), e)
get_adjacent(vor::VoronoiTessellation, i, j) = get_adjacent(get_adjacent(vor), i, j)
add_adjacent!(vor::VoronoiTessellation, e, i) = add_adjacent!(get_adjacent(vor), e, i)
add_adjacent!(vor::VoronoiTessellation, i, j, k) = add_adjacent!(get_adjacent(vor), i, j, k)
delete_adjacent!(vor::VoronoiTessellation, e, i) = delete_adjacent!(get_adjacent(vor), e, i)
delete_adjacent!(vor::VoronoiTessellation, i, j, k) = delete_adjacent!(get_adjacent(vor), i, j, k)

function polygon_features(vor::VoronoiTessellation, i)
    polygon = get_polygon(vor, i)
    if any(is_boundary_index, polygon)
        F = number_type(vor)
        return (typemax(F), (typemax(F), typemax(F)))
    end
    return polygon_features(get_polygon_points(vor), polygon)
end

get_area(vor::VoronoiTessellation, i) = polygon_features(vor, i)[1]
get_centroid(vor::VoronoiTessellation, i) = polygon_features(vor, i)[2]

function get_surrounding_polygon(vor::VoronoiTessellation, i)
    tri = get_triangulation(vor)
    S = get_surrounding_polygon(tri, i)
    push!(S, S[begin])
    return S
end

# not exactly accurate for unbounded cells due to the positioning of the 
# ghost triangles differing from solid triangles. But since we only typically need 
# this for points inside the domain or right off an edge (due to floating point 
# arithmetic), this is sufficient. If you do need point location queries for external 
# points, post an issue with some ideas.
function jump_and_march(vor::VoronoiTessellation, q; kwargs...)
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
end

## Plotting 
MakieCore.@recipe(Voronoiplot, vorn) do scene
    th = MakieCore.default_theme(scene, Mesh)
    return MakieCore.Attributes(;
        markersize=11,
        show_generators=true,
        generator_color=:black,
        strokecolor=:black,
        polygon_color=(:white, 0),
        strokewidth=1,
        unbounded_edge_extension_factor=2.0,
        colormap=th.colormap,
        colorrange=get(th.attributes, :colorrange, MakieCore.automatic),
        cycle=[:color]
    )
end

function MakieCore.plot!(p::Voronoiplot)
    ## Extract 
    vorn = p[:vorn]
    markersize = p[:markersize]
    show_generators = p[:show_generators]
    generator_color = p[:generator_color]
    strokecolor = p[:strokecolor]
    polygon_color = p[:polygon_color]
    strokewidth = p[:strokewidth]
    unbounded_edge_extension_factor = p[:unbounded_edge_extension_factor]
    colormap = p[:colormap]
    colorrange = p[:colorrange]
    cycle = p[:cycle]
    bbox = polygon_bounds(vorn[], unbounded_edge_extension_factor[])
    xmin, xmax, ymin, ymax = bbox
    bbox = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    bbox_order = [1, 2, 3, 4, 1]

    if !(typeof(polygon_color[]) <: AbstractVector)
        n = num_polygons(vorn[])
        polygon_color[] = [polygon_color[] for _ in 1:n]
    end

    ## Define all the necessary observables 
    generators_2f = MakieCore.Observable(NTuple{2,Float64}[])
    polygons = MakieCore.Observable([NTuple{2,Float64}[] for _ in each_polygon(vorn[])])

    ## Define the plotting function 
    function update_plot(vorn)
        empty!(generators_2f[])
        empty!(polygons[])
        for i in each_generator(vorn)
            push!(generators_2f[], get_generator(vorn, i))
            push!(polygons[], get_polygon_coordinates(vorn, i, bbox, bbox_order))
        end
    end

    ## Connect the plot so that it updates whenever we change a value 
    MakieCore.Observables.onany(vorn)

    ## Call it once to prepopulate with current values 
    update_plot(vorn[])

    ## Now plot 
    for i in eachindex(polygons[])
        poly!(p, polygons[][i], color=polygon_color[][i], strokecolor=strokecolor[],
            strokewidth=strokewidth[], colormap=colormap[],
            colorrange=colorrange[], cycle=cycle[])
    end
    if show_generators[]
        scatter!(p, generators_2f[], markersize=markersize[], color=generator_color[])
    end
    return p
end

## Utility functions 
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

function get_polygon_colors(vorn::VoronoiTessellation, cmap)
    F = number_type(vorn)
    gtr = [get_generator(vorn, i) for i in each_generator(vorn)]
    gtr_mat = reinterpret(reshape, F, gtr)
    colors = get(cmap, gtr_mat, :extrema)
    return [(a + b) / 2 for (a, b) in eachcol(colors)]
end

## Constructing a Voronoi tessellation 
function initialise_voronoi_tessellation(tri::Tr) where {Tr<:Triangulation}
    I = integer_type(tri)
    T = triangle_type(tri)
    F = number_type(tri)
    P = NTuple{2,F}
    polygon_points = Vector{P}()
    circumcenter_to_triangle = Dict{I,T}()
    triangle_to_circumcenter = Dict{T,I}()
    sizehint!(polygon_points, num_triangles(tri))
    sizehint!(circumcenter_to_triangle, num_triangles(tri))
    sizehint!(triangle_to_circumcenter, num_triangles(tri))
    cur_ghost_idx = I(0)
    cocircular_dict = Dict{P,I}()
    encountered_circumcenters = DefaultDict{P,I,I}(zero(I))
    cocircular_circumcenters = Set{I}()
    for V in each_triangle(tri)
        V = rotate_triangle_to_standard_form(V)
        if !is_ghost_triangle(V)
            cx, cy = triangle_circumcenter(tri, V)
            encountered_circumcenters[(cx, cy)] += 1
            if encountered_circumcenters[(cx, cy)] > 1 # If we've already encountered this circumcenter, don't push another
                idx = cocircular_dict[(cx, cy)]
                triangle_to_circumcenter[V] = idx
                push!(cocircular_circumcenters, idx)
            else
                push_point!(polygon_points, cx, cy)
                circumcenter_to_triangle[num_points(polygon_points)] = V
                triangle_to_circumcenter[V] = num_points(polygon_points)
                cocircular_dict[(cx, cy)] = num_points(polygon_points)
            end
        else
            circumcenter_to_triangle[I(BoundaryIndex)-cur_ghost_idx] = V
            triangle_to_circumcenter[V] = I(BoundaryIndex) - cur_ghost_idx
            cur_ghost_idx += I(1)
        end
    end
    polygons = Dict{I,Vector{I}}()
    sizehint!(polygons, num_solid_vertices(tri))
    unbounded_polygons = Set{I}()
    sizehint!(unbounded_polygons, num_ghost_edges(tri))
    generators = Dict{I,P}()
    sizehint!(generators, num_solid_vertices(tri))
    for i in each_solid_vertex(tri)
        generators[i] = get_point(tri, i)
    end
    E = edge_type(tri)
    adj = Adjacent{I,E}()
    return VoronoiTessellation{Tr,P,I,T,typeof(cocircular_circumcenters),E}(tri, generators, polygon_points, polygons, circumcenter_to_triangle, triangle_to_circumcenter, unbounded_polygons, cocircular_circumcenters, adj)
end

function voronoi(tri::Triangulation, clip=Val(has_boundary_nodes(tri)))
    has_ghost = has_ghost_triangles(tri)
    !has_ghost && add_ghost_triangles!(tri)
    vorn = initialise_voronoi_tessellation(tri)
    has_bnds = has_boundary_nodes(tri)
    if is_true(clip) && !has_bnds
        lock_convex_hull!(tri)
    end
    for i in each_generator(vorn)
        add_voronoi_polygon!(vorn, i)
    end
    is_true(clip) && clip_voronoi_tessellation!(vorn)
    if is_true(clip) && !has_bnds
        unlock_convex_hull!(tri)
    end
    !has_ghost && delete_ghost_triangles!(tri)
    return vorn
end

function prepare_add_voronoi_polygon(vorn::VoronoiTessellation, i)
    I = integer_type(vorn)
    S = get_surrounding_polygon(vorn, i)
    B = I[]
    sizehint!(B, length(S))
    unbounded_polygons = get_unbounded_polygons(vorn)
    return S, B, unbounded_polygons
end

function get_next_triangle_for_voronoi_polygon(vorn::VoronoiTessellation, i, k, S, m)
    T = triangle_type(vorn)
    j = k
    k = S[m]
    V = (rotate_triangle_to_standard_form ∘ construct_triangle)(T, i, j, k)
    ci = get_triangle_to_circumcenter(vorn, V)
    return ci, j, k
end

function connect_circumcenters!(B, ci)
    push!(B, ci)
    return nothing
end

function add_voronoi_polygon!(vorn::VoronoiTessellation, i)
    S, B, unbounded_polygons = prepare_add_voronoi_polygon(vorn, i)
    m = firstindex(S) + 1
    k = S[begin]
    encountered_duplicate_circumcenter = false
    ci, j, k = get_next_triangle_for_voronoi_polygon(vorn, i, k, S, m)
    if ci ∈ get_cocircular_circumcenters(vorn)
        encountered_duplicate_circumcenter = true
    end
    push!(B, ci)
    prev_ci = ci
    for m in (firstindex(S)+2):lastindex(S)
        ci, j, k = get_next_triangle_for_voronoi_polygon(vorn, i, k, S, m)
        if is_boundary_index(ci)
            push!(unbounded_polygons, i)
        end
        if ci ∈ get_cocircular_circumcenters(vorn)
            encountered_duplicate_circumcenter = true
        end
        connect_circumcenters!(B, ci)
        add_adjacent!(vorn, prev_ci, ci, i)
        prev_ci = ci
    end
    encountered_duplicate_circumcenter && unique!(B)
    push!(B, B[begin])
    add_adjacent!(vorn, prev_ci, B[begin], i)
    polygons = get_polygons(vorn)
    polygons[i] = B
    return B
end

function get_neighbouring_edges(vorn::VoronoiTessellation, e)
    tri = get_triangulation(vorn)
    i, j = edge_indices(e)
    if !is_boundary_edge(tri, e)
        i, j = j, i
    end
    bnd_idx = get_adjacent(tri, i, j)
    left_bnd = get_left_boundary_node(tri, j, bnd_idx)
    if left_bnd == i || left_bnd == j
        left_bnd = get_left_boundary_node(tri, i, bnd_idx)
    end
    right_bnd = get_right_boundary_node(tri, i, bnd_idx)
    if right_bnd == i || right_bnd == j
        right_bnd = get_right_boundary_node(tri, j, bnd_idx)
    end
    E = edge_type(tri)
    left_edge = construct_edge(E, j, left_bnd)
    right_edge = construct_edge(E, right_bnd, i)
    return left_edge, right_edge, left_bnd, right_bnd
end

"""
    segment_intersection_type(e, intersected_edge_cache)

Given edges intersected, `intersected_edge_cache`, with an edge `e` by some incident polygon's edge `(u, v)`, returns:

- `flag`: `true` if both of the edges in `intersected_edge_cache` are the same as `e`, and `false` otherwise.
- `index`: If `flag == true`, then `index == 0`, but otherwise this gives the index of the edge in `intersected_edge_cache` that is different from `e`.
- `shared_vertex`: If `flag == true`, then `shared_vertex == 0`, but otherwise this gives the vertex shared by `e` and the other edge in `intersected_edge_cache`.

For these edges, their orientation is ignored (so `(i, j)` is the same as `(j, i)`).
"""
function segment_intersection_type(e, intersected_edge_cache)
    e1, e2 = intersected_edge_cache
    i, j = edge_indices(e)
    i1, j1 = edge_indices(e1)
    i2, j2 = edge_indices(e2)
    flag_1 = (i == i1 && j == j1) || (i == j1 && j == i1)
    flag_2 = (i == i2 && j == j2) || (i == j2 && j == i2)
    if flag_1 && flag_2
        return true, 0, 0, 0
    elseif !flag_1 && flag_2
        return false, 1, i1 == i ? i1 : j1
    else # flag_1 && !flag_2
        return false, 2, i2 == i ? i2 : j2
    end
end

function arithmetic_average_unsorted(vorn::VoronoiTessellation, polygon_vertices) # polygons aren't sorted here, and they are not is_circular
    F = number_type(vorn)
    cx = zero(F)
    cy = zero(F)
    n = length(polygon_vertices)
    for i in polygon_vertices
        p = get_polygon_point(vorn, i)
        px, py = getxy(p)
        cx += px
        cy += py
    end
    cx /= n
    cy /= n
    return (cx, cy)
end

function add_segment_intersection!(segment_intersections, boundary_sites, p, i::I) where {I}
    intersection_indices = get!(Set{I}, boundary_sites, i)
    idx = findfirst(==(p), segment_intersections)
    if isnothing(idx)
        push!(segment_intersections, p)
        push!(intersection_indices, lastindex(segment_intersections))
    else
        push!(intersection_indices, idx)
    end
    return nothing
end

function process_ray_intersection!(vorn, u, v, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions)
    u_tri = get_circumcenter_to_triangle(vorn, u)
    a = geti(u_tri)
    b = getj(u_tri)
    p, q = get_generator(vorn, a, b)
    r = get_polygon_point(vorn, v)
    cert = point_position_relative_to_line(p, q, r)
    incident_polygon == 1 && @show cert, incident_polygon
    is_left(cert) && return num_intersections
    px, py = getxy(p)
    qx, qy = getxy(q)
    m = (px + qx) / 2, (py + qy) / 2
    add_segment_intersection!(segment_intersections, boundary_site_additions, m, incident_polygon)
    num_intersections += 1
    tri = get_triangulation(vorn)
    E = edge_type(tri)
    intersected_edge_cache[num_intersections] = construct_edge(E, a, b)
    return num_intersections
end

function process_segment_intersection!(vorn, u, v, e, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions, boundary_site_deletions)
    a, b = edge_indices(e)
    p, q = get_generator(vorn, a, b)
    tri = get_triangulation(vorn)
    r, s = get_polygon_point(vorn, u, v)
    intersects = line_segment_intersection_type(p, q, r, s)
    I = integer_type(tri)
    incident_polygon == 1 && @show intersects, incident_polygon
    if is_none(intersects) || is_touching(intersects)
        # Even if there are no intersections, we should still check if r and s 
        # appear outside of the boundary, since we need to declare these as being outside.
        #if is_left(point_position_on_line_segment(p, q, r)) && is_left(point_position_on_line_segment(p, q, s))
        #    # Can't have one be left and one right, since this would mean intersection. 
        #    bnd_deletions = get!(Set{I}, boundary_site_deletions, incident_polygon)
        #    push!(bnd_deletions, u, v)
        #end
        return num_intersections
    end
    m = segment_intersection_coordinates(p, q, r, s)
    if !is_boundary_edge(tri, e)
        p, q = q, p
        a, b = b, a
    end
    add_segment_intersection!(segment_intersections, boundary_site_additions, m, incident_polygon)
    num_intersections += 1
    intersected_edge_cache[num_intersections] = e
    incident_polygon == 1 && @show point_position_relative_to_line(p, q, r), point_position_relative_to_line(p, q, s), incident_polygon
    if is_left(point_position_relative_to_line(p, q, r))
        push!(boundary_site_deletions, u)
    else
        push!(boundary_site_deletions, v)
    end
    return num_intersections
end

function clip_voronoi_tessellation!(vorn::VoronoiTessellation)
    tri = get_triangulation(vorn)
    boundary_edges = (keys ∘ get_boundary_edge_map)(tri)
    I = integer_type(tri)
    boundary_sites = Set{I}()
    E = edge_type(tri)
    boundary_sites = Dict{I,E}()
    all_boundary_sites = Set{I}()
    for e in boundary_edges
        i, j = edge_indices(e)
        p, q = get_point(tri, i, j)
        px, py = getxy(p)
        qx, qy = getxy(q)

        m = (px + qx) / 2, (py + qy) / 2
        incident_polygon = jump_and_march(vorn, m; k=i)
        boundary_sites[i] = e
        boundary_sites[j] = e
        boundary_sites[incident_polygon] = e
        push!(all_boundary_sites, i, j, incident_polygon)
        # Turns out that just doing the above isn't good enough, since we might have an 
        # inner incident polygon that intersects the boundary edge at a point that 
        # is not the midpoint of the edge. So, rather than using jump_and_march(vorn, m; k = i),
        # let's just get the whole triangle. 
        #m = (px + qx) / 2, (py + qy) / 2
        #V = jump_and_march(get_triangulation(vorn), m; k=i)
        #for v in indices(V)
        #    if !is_boundary_index(v)
        #        boundary_sites[i] = e
        #        boundary_sites[j] = e
        #        boundary_sites[v] = e
        #    end
        #end
    end
    F = number_type(vorn)
    segment_intersections = NTuple{2,F}[]
    intersected_edge_cache = Vector{E}(undef, 2)
    boundary_site_additions = Dict{I,Set{I}}()
    boundary_site_deletions = Set{I}()
    incident_polygon_e_state = iterate(boundary_sites)
    while !isnothing(incident_polygon_e_state)
        (incident_polygon, e), state = incident_polygon_e_state
        left_edge, right_edge, left_bnd, right_bnd = get_neighbouring_edges(vorn, e)
        polygon = get_polygon(vorn, incident_polygon)
        nedges = num_boundary_edges(polygon)
        num_intersections = 0
        for ℓ in 1:nedges
            u = get_boundary_nodes(polygon, ℓ)
            v = get_boundary_nodes(polygon, ℓ + 1)
            if num_intersections < 2
                if is_boundary_index(u) && is_boundary_index(v)
                    continue
                elseif is_boundary_index(u) && !is_boundary_index(v)
                    _num_intersections = process_ray_intersection!(vorn, u, v, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions)
                elseif !is_boundary_index(u) && is_boundary_index(v)
                    _num_intersections = process_ray_intersection!(vorn, v, u, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions)
                else
                    for e in (e, left_edge, right_edge)
                        _num_intersections = process_segment_intersection!(vorn, u, v, e, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions, boundary_site_deletions)
                        if _num_intersections > num_intersections
                            num_intersections = _num_intersections
                            # break < -- Actually, don't break because a single line could go past multiple parts of the boundary (e.g. the fixed_shewchuk_example_constrained example).
                            # Need to also check the adjacent incident polygon here, to get the correct form of e. This won't duplicate the same code below, since we won't make it 
                            # through adjacent_incident_polygon ∉ all_boundary_sites.
                            adjacent_incident_polygon = get_adjacent(vorn, v, u)
                            if adjacent_incident_polygon ∉ all_boundary_sites
                                push!(all_boundary_sites, adjacent_incident_polygon)
                                boundary_sites[adjacent_incident_polygon] = e
                            end
                            num_intersections == 2 && break # This is the correct way to break.
                        end
                    end
                end
                if _num_intersections > num_intersections
                    num_intersections = _num_intersections
                else
                    adjacent_incident_polygon = get_adjacent(vorn, v, u)
                    if adjacent_incident_polygon ∉ all_boundary_sites
                        push!(all_boundary_sites, adjacent_incident_polygon)
                        boundary_sites[adjacent_incident_polygon] = e
                    end
                end
                for u in (u, v)
                    if !is_boundary_index(u) && u ∉ boundary_site_deletions
                        for e in (e, left_edge, right_edge)
                            ei, ej = edge_indices(e)
                            if !is_boundary_edge(tri, e)
                                ei, ej = ej, ei
                            end
                            p, q = get_generator(vorn, ei, ej)
                            incident_polygon == 1 && @show incident_polygon, point_position_relative_to_line(p, q, get_polygon_point(vorn, u))
                            if is_left(point_position_relative_to_line(p, q, get_polygon_point(vorn, u)))
                                push!(boundary_site_deletions, u)
                                break
                            end
                        end
                    end
                end
            else
                #Even if we have reached two intersections, we still need to clean up any vertices from the 
                #polygons that might be outside. To avoid doing extra processing in process_segment_intersection,
                #here we check both u and v (else we need to make sure we don't miss a u when we miss a segment 
                #intersection).)
                for u in (u, v)
                    if !is_boundary_index(u) && u ∉ boundary_site_deletions
                        for e in (e, left_edge, right_edge)
                            ei, ej = edge_indices(e)
                            if !is_boundary_edge(tri, e)
                                ei, ej = ej, ei
                            end
                            p, q = get_generator(vorn, ei, ej)
                            incident_polygon == 1 && @show incident_polygon, point_position_relative_to_line(p, q, get_polygon_point(vorn, u))
                            if is_left(point_position_relative_to_line(p, q, get_polygon_point(vorn, u)))
                                push!(boundary_site_deletions, u)
                                break
                            end
                        end
                    end
                end
            end
        end
        if num_intersections == 2
            interior_intersection, index, shared_vertex = segment_intersection_type(e, intersected_edge_cache)
            if !interior_intersection
                other_e = intersected_edge_cache[index]
                m = get_generator(vorn, shared_vertex)
                add_segment_intersection!(segment_intersections, boundary_site_additions, m, incident_polygon)
            end
        end
        delete!(boundary_sites, incident_polygon)
        incident_polygon_e_state = iterate(boundary_sites, state)
    end
    n = num_polygon_vertices(vorn)
    for p in each_point(segment_intersections)
        push_polygon_point!(vorn, p)
    end
    boundary_sites = all_boundary_sites
    for polygon in boundary_sites
        polygon_vertices = get_polygon(vorn, polygon)
        # First, reset the adjacencies
        ne = num_boundary_edges(polygon_vertices)
        for ℓ in 1:ne 
            u = get_boundary_nodes(polygon_vertices, ℓ)
            v = get_boundary_nodes(polygon_vertices, ℓ + 1)
            delete_adjacent!(vorn, u, v)
        end
        pop!(polygon_vertices) # remove the last vertex, which is a duplicate of the first
        filter!(!is_boundary_index, polygon_vertices)
        filter!(!∈(boundary_site_deletions), polygon_vertices)
        for new_vert in boundary_site_additions[polygon]
            push!(polygon_vertices, n + new_vert)
        end
        cx, cy = arithmetic_average_unsorted(vorn, polygon_vertices)
        θ = zeros(F, length(polygon_vertices))
        for (j, i) in pairs(polygon_vertices)
            p = get_polygon_point(vorn, i)
            px, py = getxy(p)
            θ[j] = atan(py - cy, px - cx)
        end
        idx = sortperm(θ)
        permute!(polygon_vertices, idx)
        push!(polygon_vertices, polygon_vertices[begin])
        # Now fix the adjacencies 
        ne = num_boundary_edges(polygon_vertices)
        for ℓ in 1:ne 
            u = get_boundary_nodes(polygon_vertices, ℓ)
            v = get_boundary_nodes(polygon_vertices, ℓ + 1)
            add_adjacent!(vorn, u, v, polygon)
        end
    end
    empty!(get_unbounded_polygons(vorn))
end
