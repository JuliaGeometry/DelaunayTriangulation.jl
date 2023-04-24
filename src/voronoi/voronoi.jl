

num_polygons(vor::VoronoiTessellation) = length(get_polygons(vor))
num_polygon_vertices(vor::VoronoiTessellation) = num_points(get_polygon_points(vor))
get_generator(vor::VoronoiTessellation, i) = get_generators(vor)[i]
get_generator(vor::VoronoiTessellation, i::Vararg{I,N}) where {I,N} = ntuple(j -> get_generator(vor, i[j]), Val(N))
get_polygon_point(vor::VoronoiTessellation, i...) = get_point(get_polygon_points(vor), i...)
get_polygon(vor::VoronoiTessellation, i) = get_polygons(vor)[i]
get_circumcenter_to_triangle(vor::VoronoiTessellation, i) = get_circumcenter_to_triangle(vor)[i]
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
number_type(vor::VoronoiTessellation) = number_type(get_polygon_points(vor))
num_generators(vor::VoronoiTessellation) = length(get_generators(vor))
integer_type(::VoronoiTessellation{Tr,P,I}) where {Tr,P,I} = I
add_polygon!(vor::VoronoiTessellation, B, i) = get_polygons(vor)[i] = B
number_type(::VoronoiTessellation{Tr,P}) where {Tr,P} = number_type(P)
triangle_type(::VoronoiTessellation{Tr,P,I,T}) where {Tr,P,I,T} = T
each_generator(vor::VoronoiTessellation) = keys(get_generators(vor))
each_polygon_vertex(vor::VoronoiTessellation) = eachindex(get_polygon_points(vor))
each_unbounded_polygon(vor::VoronoiTessellation) = get_unbounded_polygons(vor)
each_polygon(vor::VoronoiTessellation) = values(get_polygons(vor))
contains_boundary_edge(vor::VoronoiTessellation, e) = contains_boundary_edge(get_triangulation(vor), e)
contains_boundary_edge(vor::VoronoiTessellation, i, j) = contains_boundary_edge(get_triangulation(vor), i, j)
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
add_unbounded_polygon!(vor::VoronoiTessellation, i) = push!(get_unbounded_polygons(vor), i)
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
delete_unbounded_polygon!(vor::VoronoiTessellation, i) = delete!(get_unbounded_polygons(vor), i)

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
        resize!(generators_2f[], num_generators(vorn))
        resize!(polygons[], num_generators(vorn))
        for i in each_generator(vorn)
            generators_2f[][i] = get_generator(vorn, i)
            polygons[][i] = get_polygon_coordinates(vorn, i, bbox, bbox_order)
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
            u, v, w = indices(V)
            p, q, r = get_point(tri, u, v, w)
            A = triangle_area(p, q, r)
            cx, cy = triangle_circumcenter(p, q, r, A)
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
    boundary_polygons = Set{I}()
    return VoronoiTessellation{Tr,P,I,T,typeof(cocircular_circumcenters),E}(tri, generators, polygon_points, polygons, circumcenter_to_triangle, triangle_to_circumcenter, unbounded_polygons, cocircular_circumcenters, adj, boundary_polygons)
end

function voronoi(tri::Triangulation, clip=Val(has_boundary_nodes(tri)))
    has_ghost = has_ghost_triangles(tri)
    !has_ghost && add_ghost_triangles!(tri)
    vorn = initialise_voronoi_tessellation(tri)
    has_bnds = has_boundary_nodes(tri)
    is_convex = has_bnds
    is_true(clip) && !has_bnds && lock_convex_hull!(tri)
    for i in each_generator(vorn)
        add_voronoi_polygon!(vorn, i)
    end
    is_true(clip) && clip_voronoi_tessellation!(vorn, is_convex)
    is_true(clip) && !has_bnds && unlock_convex_hull!(tri)
    !has_ghost && delete_ghost_triangles!(tri)
    return vorn
end

function prepare_add_voronoi_polygon(vorn::VoronoiTessellation, i)
    I = integer_type(vorn)
    S = get_surrounding_polygon(vorn, i)
    B = I[]
    sizehint!(B, length(S))
    return S, B
end

function get_next_triangle_for_voronoi_polygon(vorn::VoronoiTessellation, i, k, S, m)
    T = triangle_type(vorn)
    j = k
    k = S[m]
    V = (rotate_triangle_to_standard_form ∘ construct_triangle)(T, i, j, k)
    ci = get_triangle_to_circumcenter(vorn, V)
    return ci, k
end

function connect_circumcenters!(B, ci)
    push!(B, ci)
    return nothing
end

function add_edge_to_voronoi_polygon!(B, vorn::VoronoiTessellation, i, k, S, m, encountered_duplicate_circumcenter)
    ci, k = get_next_triangle_for_voronoi_polygon(vorn, i, k, S, m)
    is_boundary_index(ci) && add_unbounded_polygon!(vorn, i)
    (encountered_duplicate_circumcenter || ci ∈ get_cocircular_circumcenters(vorn)) && (encountered_duplicate_circumcenter = true)
    connect_circumcenters!(B, ci)
    return ci, encountered_duplicate_circumcenter, k
end

function close_voronoi_polygon!(vorn::VoronoiTessellation, B, i, encountered_duplicate_circumcenter, prev_ci)
    encountered_duplicate_circumcenter && unique!(B)
    connect_circumcenters!(B, B[begin])
    add_adjacent!(vorn, prev_ci, B[begin], i)
    add_polygon!(vorn, B, i)
    return nothing
end

function add_voronoi_polygon!(vorn::VoronoiTessellation, i)
    S, B = prepare_add_voronoi_polygon(vorn, i)
    m = firstindex(S) + 1
    k = S[begin]
    encountered_duplicate_circumcenter = false
    prev_ci, encountered_duplicate_circumcenter, k = add_edge_to_voronoi_polygon!(B, vorn, i, k, S, m, encountered_duplicate_circumcenter)
    for m in (firstindex(S)+2):lastindex(S)
        ci, encountered_duplicate_circumcenter, k = add_edge_to_voronoi_polygon!(B, vorn, i, k, S, m, encountered_duplicate_circumcenter)
        add_adjacent!(vorn, prev_ci, ci, i)
        prev_ci = ci
    end
    close_voronoi_polygon!(vorn, B, i, encountered_duplicate_circumcenter, prev_ci)
    return B
end

function add_segment_intersection!(segment_intersections, boundary_sites, intersection_point, incident_polygon::I) where {I}
    intersection_indices = get!(Set{I}, boundary_sites, incident_polygon)
    idx = findfirst(==(intersection_point), segment_intersections)
    if idx === nothing
        push_point!(segment_intersections, intersection_point)
        idx = num_points(segment_intersections)
    end
    push!(intersection_indices, idx)
    return idx
end

function add_to_intersected_edge_cache!(intersected_edge_cache::AbstractVector{V}, u, v, a, b) where {E,V<:Pair{E,E}}
    uv = construct_edge(E, u, v)
    ab = construct_edge(E, a, b)
    push!(intersected_edge_cache, uv => ab)
    return nothing
end

function process_ray_intersection!(
    vorn::VoronoiTessellation,
    u,
    v,
    incident_polygon,
    intersected_edge_cache,
    segment_intersections,
    boundary_sites,
    exterior_circumcenters,
    equal_circumcenter_mapping)
    u_tri = get_circumcenter_to_triangle(vorn, u)
    a = geti(u_tri)
    b = getj(u_tri)
    p, q = get_generator(vorn, a, b)
    r = get_polygon_point(vorn, v)
    _, intersection_coordinates = intersection_of_edge_and_bisector_ray(p, q, r)
    F = number_type(vorn)
    any(isnan, intersection_coordinates) && (push!(exterior_circumcenters, v); return (F(NaN), F(NaN)))
    idx = add_segment_intersection!(segment_intersections, boundary_sites, intersection_coordinates, incident_polygon)
    if intersection_coordinates == r
        equal_circumcenter_mapping[idx] = v
    end
    add_to_intersected_edge_cache!(intersected_edge_cache, u, v, a, b)
    return intersection_coordinates
end

function process_segment_intersection!(
    vorn,
    u,
    v,
    e,
    incident_polygon,
    intersected_edge_cache,
    segment_intersections,
    boundary_sites,
    exterior_circumcenters,
    equal_circumcenter_mapping)
    e = convert_to_boundary_edge(vorn, e)
    a, b = edge_indices(e)
    p, q = get_generator(vorn, a, b)
    r, s = get_polygon_point(vorn, u, v)
    intersection_cert, cert_u, cert_v, intersection_coordinates = classify_and_compute_segment_intersection(p, q, r, s)
    F = number_type(vorn)
    if is_none(intersection_cert) || is_touching(intersection_cert)
        if is_left(cert_u) && is_left(cert_v)
            push!(exterior_circumcenters, u, v)
        end
        return (F(NaN), F(NaN))
    end
    idx = add_segment_intersection!(segment_intersections, boundary_sites, intersection_coordinates, incident_polygon)
    if intersection_coordinates == r
        equal_circumcenter_mapping[idx] = u
    elseif intersection_coordinates == s
        equal_circumcenter_mapping[idx] = v
    end
    add_to_intersected_edge_cache!(intersected_edge_cache, u, v, a, b)
    return intersection_coordinates
end

function initialise_clipping_arrays(vorn::VoronoiTessellation)
    tri = get_triangulation(vorn)
    E = edge_type(vorn)
    I = integer_type(vorn)
    boundary_edges = (keys ∘ get_boundary_edge_map)(tri)
    edges_to_process = Set{E}()
    foreach(boundary_edges) do e
        push!(edges_to_process, e)
    end
    polygon_edge_queue = Queue{Tuple{E,I}}()
    boundary_sites = Dict{I,Set{I}}()
    F = number_type(vorn)
    segment_intersections = NTuple{2,F}[]
    processed_pairs = Set{Tuple{E,I}}()
    intersected_edge_cache = Pair{E,E}[]
    sizehint!(intersected_edge_cache, 2^3)
    exterior_circumcenters = Set{I}()
    left_edge_intersectors = Set{E}()
    right_edge_intersectors = Set{E}()
    current_edge_intersectors = Set{E}()
    equal_circumcenter_mapping = Dict{I,I}()
    return edges_to_process, polygon_edge_queue, boundary_sites, segment_intersections, processed_pairs, intersected_edge_cache, exterior_circumcenters,
    left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, equal_circumcenter_mapping
end

function enqueue_new_edge!(polygon_edge_queue, vorn::VoronoiTessellation, e)
    u, v = edge_indices(e)
    p, q = get_generator(vorn, u, v)
    px, py = getxy(p)
    qx, qy = getxy(q)
    m = (px + qx) / 2, (py + qy) / 2
    incident_polygon = jump_and_march(vorn, m; k=u)
    enqueue!(polygon_edge_queue, (e, incident_polygon))
    return nothing
end

is_segment_between_two_ghosts(u, v) = is_boundary_index(u) && is_boundary_index(v)
is_ray_going_in(u, v) = is_boundary_index(u) && !is_boundary_index(v)
is_ray_going_out(u, v) = !is_boundary_index(u) && is_boundary_index(v)
is_finite_segment(u, v) = !is_boundary_index(u) && !is_boundary_index(v)
get_neighbouring_boundary_edges(vorn::VoronoiTessellation, e) = get_neighbouring_boundary_edges(get_triangulation(vorn), e)
convert_to_boundary_edge(vorn::VoronoiTessellation, e) = convert_to_boundary_edge(get_triangulation(vorn), e)

function process_ray_intersection_with_other_edges!(vorn, u, v, e, left_edge, right_edge, r, segment_intersections,
    boundary_sites, incident_polygon, equal_circumcenter_mapping, intersected_edge_cache)
    if !any(isnan, r)
        E = edge_type(vorn)
        u_tri = get_circumcenter_to_triangle(vorn, u)
        a, b, _ = indices(u_tri)
        s = get_polygon_point(vorn, v)
        intersected_edge = construct_edge(E, a, b)
        for _e in (e, left_edge, right_edge)
            if !compare_unoriented_edge(intersected_edge, _e)
                i, j = edge_indices(_e)
                p, q = get_generator(vorn, i, j)
                intersection_cert, cert_u, cert_v, intersection_coordinates = classify_and_compute_segment_intersection(p, q, r, s)
                if !(is_none(intersection_cert) || is_touching(intersection_cert))
                    idx = add_segment_intersection!(segment_intersections, boundary_sites, intersection_coordinates, incident_polygon)
                    if intersection_coordinates == s # don't need to check r, since it would have been checked in process_ray_intersection! already
                        equal_circumcenter_mapping[idx] = v
                    end
                    add_to_intersected_edge_cache!(intersected_edge_cache, u, v, i, j)
                end
            end
        end
    end
    return nothing
end

function process_polygon!(vorn, e, incident_polygon, boundary_sites, segment_intersections, intersected_edge_cache, exterior_circumcenters, equal_circumcenter_mapping)
    left_edge, right_edge = get_neighbouring_boundary_edges(vorn, e)
    polygon_vertices = get_polygon(vorn, incident_polygon)
    nedges = num_boundary_edges(polygon_vertices)
    for ℓ in 1:nedges
        u = get_boundary_nodes(polygon_vertices, ℓ)
        v = get_boundary_nodes(polygon_vertices, ℓ + 1)
        if is_segment_between_two_ghosts(u, v)
            continue
        elseif is_ray_going_in(u, v)
            r = process_ray_intersection!(vorn, u, v, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping)
            # It's possible for an infinite ray to also intersect other boundary edges, e.g. look at 
            #   points = [0.290978 0.830755 0.0139574; 0.386411 0.630008 0.803881]
            # So, let's just look for intersections with other edges.
            process_ray_intersection_with_other_edges!(vorn, u, v, e, left_edge, right_edge, r, segment_intersections, boundary_sites, incident_polygon, equal_circumcenter_mapping, intersected_edge_cache)
        elseif is_ray_going_out(u, v)
            r = process_ray_intersection!(vorn, v, u, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping)
            process_ray_intersection_with_other_edges!(vorn, v, u, e, left_edge, right_edge, r, segment_intersections, boundary_sites, incident_polygon, equal_circumcenter_mapping, intersected_edge_cache)
        elseif is_finite_segment(u, v)
            for _e in (e, left_edge, right_edge)
                process_segment_intersection!(vorn, u, v, _e, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping)
            end
        end
    end
    return left_edge, right_edge, e
end

function classify_intersections!(intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, left_edge, right_edge, current_edge)
    for (uv, e) in intersected_edge_cache
        if compare_unoriented_edge(e, left_edge)
            push!(left_edge_intersectors, uv)
        elseif compare_unoriented_edge(e, right_edge)
            push!(right_edge_intersectors, uv)
        elseif compare_unoriented_edge(e, current_edge)
            push!(current_edge_intersectors, uv)
        end
    end
    return nothing
end

function process_intersection_points!(polygon_edge_queue, vorn, current_incident_polygon,
    left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
    left_edge, right_edge, current_edge, processed_pairs, segment_intersections, boundary_sites)
    all_indices = (initial(left_edge), terminal(left_edge),
        initial(right_edge), terminal(right_edge),
        initial(current_edge), terminal(current_edge))
    if num_polygon_vertices(vorn) > 1 # A single triangle is a special case that we add the corners into manually
        for (e, intersectors) in zip((left_edge, right_edge), (left_edge_intersectors, right_edge_intersectors))
            if (length(intersectors) > 0 && length(current_edge_intersectors) > 0) && ((e, current_incident_polygon) ∉ processed_pairs && (reverse_edge(e), current_incident_polygon) ∉ processed_pairs)
                i, j = edge_indices(e)
                enqueue!(polygon_edge_queue, (e, i))
                enqueue!(polygon_edge_queue, (e, j))
                if current_incident_polygon ∈ all_indices # only need to consider a corner point if the point we are considering is a point on the boundary
                    u = get_shared_vertex(e, current_edge)
                    if u == current_incident_polygon # The only way we can get a corner point like this if it corresponds to the same generator we already considering
                        p = get_generator(vorn, u)
                        add_segment_intersection!(segment_intersections, boundary_sites, p, current_incident_polygon)
                    end
                end
            end
        end
    end
    for (e, intersectors) in zip((left_edge, right_edge, current_edge), (left_edge_intersectors, right_edge_intersectors, current_edge_intersectors))
        for uv in intersectors
            u, v = edge_indices(uv)
            adjacent_incident_polygon = get_adjacent(vorn, v, u)
            if adjacent_incident_polygon == current_incident_polygon
                adjacent_incident_polygon = get_adjacent(vorn, u, v)
            end
            if (e, adjacent_incident_polygon) ∉ processed_pairs && (reverse_edge(e), adjacent_incident_polygon) ∉ processed_pairs
                enqueue!(polygon_edge_queue, (e, adjacent_incident_polygon))
            end
        end
    end
end

function dequeue_and_process!(vorn, polygon_edge_queue, edges_to_process,
    intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
    processed_pairs, boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping)
    if isempty(polygon_edge_queue)
        e = convert_to_boundary_edge(vorn, first(edges_to_process))
        enqueue_new_edge!(polygon_edge_queue, vorn, e)
    end
    e, incident_polygon = dequeue!(polygon_edge_queue)
    push!(processed_pairs, (e, incident_polygon))
    for cache in (intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors)
        empty!(cache)
    end
    left_edge, right_edge, e = process_polygon!(vorn, e, incident_polygon, boundary_sites, segment_intersections, intersected_edge_cache, exterior_circumcenters, equal_circumcenter_mapping)
    classify_intersections!(intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, left_edge, right_edge, e)
    process_intersection_points!(polygon_edge_queue, vorn, incident_polygon,
        left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
        left_edge, right_edge, e, processed_pairs, segment_intersections, boundary_sites)
    if contains_edge(e, edges_to_process)
        delete!(edges_to_process, e)
    elseif contains_edge(reverse_edge(e), edges_to_process)
        delete!(edges_to_process, reverse_edge(e))
    end
    return nothing
end

function find_all_intersections(vorn::VoronoiTessellation)
    edges_to_process,
    polygon_edge_queue,
    boundary_sites,
    segment_intersections,
    processed_pairs,
    intersected_edge_cache,
    exterior_circumcenters,
    left_edge_intersectors,
    right_edge_intersectors,
    current_edge_intersectors,
    equal_circumcenter_mapping = initialise_clipping_arrays(vorn)
    e = convert_to_boundary_edge(vorn, first(edges_to_process))
    enqueue_new_edge!(polygon_edge_queue, vorn, e)
    while !isempty(edges_to_process) || !isempty(polygon_edge_queue)
        dequeue_and_process!(vorn, polygon_edge_queue, edges_to_process,
            intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
            processed_pairs, boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping)
    end
    if num_polygon_vertices(vorn) == 1 # 1 triangle 
        for i in each_generator(vorn)
            p = get_generator(vorn, i)
            add_segment_intersection!(segment_intersections, boundary_sites, p, i)
        end
    end
    return boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping
end

function add_intersection_points!(vorn::VoronoiTessellation, segment_intersections)#, equal_circumcenter_mapping)
    n = num_polygon_vertices(vorn)
    for i in each_point_index(segment_intersections)
        p = get_point(segment_intersections, i)
        push_polygon_point!(vorn, p)
    end
    return n
end

function clip_polygon!(vorn::VoronoiTessellation, n, points, polygon, new_verts, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    delete_polygon_adjacent!(vorn, polygon)
    vertices = get_polygon(vorn, polygon)
    pop!(vertices) # vertices[begin] == vertices[end]
    filter!(v -> !is_boundary_index(v) && v ∉ exterior_circumcenters, vertices)
    for new_vert in new_verts
        if new_vert ∉ keys(equal_circumcenter_mapping) || equal_circumcenter_mapping[new_vert] ∉ vertices
            push!(vertices, n + new_vert)
        end
    end
    sort_convex_polygon!(vertices, points)
    push!(vertices, vertices[begin])
    add_polygon_adjacent!(vorn, polygon)
    delete_unbounded_polygon!(vorn, polygon)
end

function clip_all_polygons!(vorn::VoronoiTessellation, n, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    points = get_polygon_points(vorn)
    for (polygon, new_verts) in boundary_sites
        clip_polygon!(vorn, n, points, polygon, new_verts, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    end
end

function add_boundary_polygon!(vorn::VoronoiTessellation, i)
    push!(get_boundary_polygons(vorn), i)
end

function add_all_boundary_polygons!(vorn::VoronoiTessellation, boundary_sites)
    for i in keys(boundary_sites)
        add_boundary_polygon!(vorn, i)
    end
    return nothing
end

function clip_voronoi_tessellation!(vorn::VoronoiTessellation, is_convex=true)
    boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping = find_all_intersections(vorn)
    n = add_intersection_points!(vorn, segment_intersections)
    clip_all_polygons!(vorn, n, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    add_all_boundary_polygons!(vorn, boundary_sites)
    return nothing
end

function move_generator_to_centroid!(points, vorn::VoronoiTessellation, generator)
    c = get_centroid(vorn, generator)
    cx, cy = getxy(c)
    p = get_generator(vorn, generator)
    px, py = getxy(p)
    dist = sqrt((cx - px)^2 + (cy - py)^2)
    set_point!(points, generator, cx, cy)
    return dist
end

function default_displacement_tolerance(vorn::VoronoiTessellation)
    xmin, xmax, ymin, ymax = polygon_bounds(vorn)
    max_extent = max(xmax - xmin, ymax - ymin)
    return 1e-4max_extent
end

function _centroidal_smooth_itr(vorn::VoronoiTessellation, set_boundary_nodes, points, edges, boundary_nodes,
    I, E, V, Es, Ts, F, rng; kwargs...)
    max_dist = zero(F)
    for i in each_generator(vorn)
        if i ∉ set_boundary_nodes
            dist = move_generator_to_centroid!(points, vorn, i)
            max_dist = max(max_dist, dist)
        end
    end
    _tri = triangulate(points; edges, boundary_nodes,
        IntegerType=I, EdgeType=E, TriangleType=V, EdgesType=Es, TrianglesType=Ts, delete_ghosts=false,
        delete_empty_features=false, check_arguments=false, rng=rng, kwargs...)
    vorn = voronoi(_tri, true)
    return vorn, max_dist
end

function centroidal_smooth(vorn::VoronoiTessellation{Tr}; maxiters=1000, tol=default_displacement_tolerance(vorn), rng=Random.default_rng(), kwargs...) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL,Tr<:Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}}
    iter = 0
    F = number_type(vorn)
    max_dist = typemax(F)
    tri = get_triangulation(vorn)
    has_ghost = has_ghost_triangles(tri)
    !has_ghost && add_ghost_triangles!(tri)
    has_bnds = has_boundary_nodes(tri)
    !has_bnds && lock_convex_hull!(tri)
    set_boundary_nodes = get_all_boundary_nodes(tri)
    points = (deepcopy ∘ get_points)(tri)
    boundary_nodes = get_boundary_nodes(tri)
    edges = get_constrained_edges(tri)
    if isempty(edges)
        edges = nothing
    end
    V = triangle_type(Ts)
    while iter < maxiters && max_dist > tol
        vorn, max_dist = _centroidal_smooth_itr(vorn, set_boundary_nodes, points, edges, boundary_nodes,
            I, E, V, Es, Ts, F, rng; kwargs...)
        iter += 1
    end
    !has_bnds && unlock_convex_hull!(tri)
    !has_ghost && delete_ghost_triangles!(tri)
    return vorn
end