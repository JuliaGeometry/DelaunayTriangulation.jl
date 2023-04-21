## Definition

"""
    VoronoiTessellation{Tr<:Triangulation, P, I, T}

A Voronoi tessellation, dual to a triangulation.

!!! warning 

    If the triangulation is constrained, it is not guaranteed that 
    the Voronoi tessellation is dual to the triangulation. This duality 
    is only guaranteed for unconstrained triangulations. 

# Fields
- `triangulation::Tr`: The triangulation that the Voronoi tessellation is based on.
- `generators::Dict{I, P}`: The points that define the generators of the Voronoi tessellation, with indices mapping to coordinates. These are the same as the points of the triangulation, but dealiased to allow for mutation (e.g. if using Lloyd's algorithm). We use a `Dict` since the point set in the triangulation might not have every point included in the triangulation, so a `Dict` helps preserve the generator indices (rather than a vector which could shift the indices).
- `polygon_points::Vector{P}`: The points that define points on the boundaries of the Voronoi polygons. This will include both the circumcenters and any points coming from intersections with the boundary (if the tessellation is clipped).
- `polygons::Dict{I, Vector{I}}`: The Voronoi polygons. Each polygon is a polygon, defined by a list of indices into `polygon_point`. Each vector is a circular vector (see [`is_circular`](@ref)). The keys of the `Dict` enumerate the polygon, with the key `i` corresponding to the `i`th point (the generator) of the triangulation.
- `circumcenter_to_triangle::Dict{I,T}`: Map that takes a circumcenter to the triangle that it is the circumcenter of.
- `triangle_to_circumcenter::Dict{T,I}`: Map that takes a triangle to the circumcenter of that triangle.
- `unbounded_polygons::Set{I}`: The polygons that are unbounded.
"""
struct VoronoiTessellation{Tr<:Triangulation,P,I,T}
    triangulation::Tr
    generators::Dict{I,P}
    polygon_points::Vector{P}
    polygons::Dict{I,Vector{I}}
    circumcenter_to_triangle::Dict{I,T}
    triangle_to_circumcenter::Dict{T,I}
    unbounded_polygons::Set{I}
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

function polygon_features(vor::VoronoiTessellation, i)
    polygon = get_polygon(vor, i)
    if any(is_boundary_index, polygon)
        F = number_type(vorn)
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
# points, post an issue.
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
    for V in each_triangle(tri)
        V = rotate_triangle_to_standard_form(V)
        if !is_ghost_triangle(V)
            cx, cy = triangle_circumcenter(tri, V)
            push_point!(polygon_points, cx, cy)
            circumcenter_to_triangle[num_points(polygon_points)] = V
            triangle_to_circumcenter[V] = num_points(polygon_points)
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
    return VoronoiTessellation{Tr,P,I,T}(tri, generators, polygon_points, polygons, circumcenter_to_triangle, triangle_to_circumcenter, unbounded_polygons)
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
    ci, j, k = get_next_triangle_for_voronoi_polygon(vorn, i, k, S, m)
    prev_ci = ci
    push!(B, ci)
    for m in (firstindex(S)+2):lastindex(S)
        ci, j, k = get_next_triangle_for_voronoi_polygon(vorn, i, k, S, m)
        if is_boundary_index(ci)
            push!(unbounded_polygons, i)
        end
        connect_circumcenters!(B, ci)
        prev_ci = ci
    end
    push!(B, B[begin])
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

Given edges intersected, `intersected_edge_cache`, with an edge `e` by some incident polygon, returns:

- `flag`: `true` if both of the edges in `intersected_edge_cache` are the same as `e`, and `false` otherwise.
- `index`: If `flag == true`, then `index == 0`, but otherwise this gives the index of the edge in `intersected_edge_cache` that is different from `e`.
- `vertex`: If `flag == true`, then `vertex == 0`, but otherwise this gives the other vertex on the other edge shared by `e`. 
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
        return false, 1, i1 == i ? j1 : i1, i1 == i ? i1 : j1
    else # flag_1 && !flag_2
        return false, 2, i2 == i ? j2 : i2, i2 == i ? i2 : j2
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
    r, s = get_polygon_point(vorn, u, v)
    intersects = line_segment_intersection_type(p, q, r, s)
    is_none(intersects) && return num_intersections
    m = segment_intersection_coordinates(p, q, r, s)
    add_segment_intersection!(segment_intersections, boundary_site_additions, m, incident_polygon)
    num_intersections += 1
    intersected_edge_cache[num_intersections] = e
    tri = get_triangulation(vorn)
    if !is_boundary_edge(tri, e)
        a, b = b, a
        p, q = q, p
    end
    if is_left(point_position_relative_to_line(p, q, r))
        # The finite segment intersects an edge, so we need to determine which 
        # part of the segment was outside of the boundary.
        boundary_site_deletions[incident_polygon] = u
    else
        boundary_site_deletions[incident_polygon] = v
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
    end
    F = number_type(vorn)
    segment_intersections = NTuple{2,F}[]
    intersected_edge_cache = Vector{E}(undef, 2)
    boundary_site_additions = Dict{I,Set{I}}()
    boundary_site_deletions = Dict{I,I}()
    for (incident_polygon, e) in boundary_sites
        left_edge, right_edge, left_bnd, right_bnd = get_neighbouring_edges(vorn, e)
        polygon = get_polygon(vorn, incident_polygon)
        nedges = num_boundary_edges(polygon)
        num_intersections = 0
        for ℓ in 1:nedges
            num_intersections == 2 && break
            u = get_boundary_nodes(polygon, ℓ)
            v = get_boundary_nodes(polygon, ℓ + 1)
            if is_boundary_index(u) && is_boundary_index(v)
                continue
            elseif is_boundary_index(u) && !is_boundary_index(v)
                num_intersections = process_ray_intersection!(vorn, u, v, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions)
            elseif !is_boundary_index(u) && is_boundary_index(v)
                num_intersections = process_ray_intersection!(vorn, v, u, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions)
            else
                for e in (e, left_edge, right_edge)
                    _num_intersections = process_segment_intersection!(vorn, u, v, e, incident_polygon, num_intersections, intersected_edge_cache, segment_intersections, boundary_site_additions, boundary_site_deletions)
                    if _num_intersections > num_intersections
                        num_intersections = _num_intersections
                        break
                    end
                end
            end
        end
        interior_intersection, index, adjacent_incident_polygon, shared_vertex = segment_intersection_type(e, intersected_edge_cache)
        if !interior_intersection
            other_e = intersected_edge_cache[index]
            m = get_generator(vorn, shared_vertex)
            add_segment_intersection!(segment_intersections, boundary_site_additions, m, incident_polygon)
        end
    end
    n = num_polygon_vertices(vorn)
    for p in each_point(segment_intersections)
        push_polygon_point!(vorn, p)
    end
    boundary_sites = keys(boundary_sites)
    for polygon in boundary_sites
        polygon_vertices = get_polygon(vorn, polygon)
        pop!(polygon_vertices) # remove the last vertex, which is a duplicate of the first
        filter!(!is_boundary_index, polygon_vertices)
        if haskey(boundary_site_deletions, polygon)
            deleteat!(polygon_vertices, findfirst(isequal(boundary_site_deletions[polygon]), polygon_vertices))
        end
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
    end
    empty!(get_unbounded_polygons(vorn))
end
