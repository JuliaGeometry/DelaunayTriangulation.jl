struct VoronoiTessellation{Tr<:Triangulation,P,I,T,S,E}
    triangulation::Tr
    generators::Dict{I,P}
    polygon_points::Vector{P}
    polygons::Dict{I,Vector{I}}
    circumcenters_to_triangle::Dict{I,T}
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

edge_type(::VoronoiTessellation{Tr,P,I,T,S,E}) where {Tr,P,I,T,S,E} = E
num_polygons(vor::VoronoiTessellration) = length(get_polygons(vor))
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
add_polygon!(vor::VoronoiTessellation, B, i) = get_polygons(vor)[i] = B
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
add_unbounded_polygon!(vor::VoronoiTessellation, i) = push!(get_unbounded_polygons(vor), i)

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

function add_edge_to_voronoi_polygon!(B, vorn::VoronoITessellation, i, k, S, m, encountered_duplicate_circumcenter)
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
        prev_ci, encountered_duplicate_circumcenter, k = add_edge_to_voronoi_polygon!(B, vorn, i, k, S, m, encountered_duplicate_circumcenter)
        add_adjacent!(vorn, prev_ci, ci, i)
    end
    close_voronoi_polygon!(vorn, B, i, encountered_duplicate_circumcenter, prev_ci)
    return B
end

function add_segment_intersection!(segment_intersections, boundary_sites, intersection_point, incident_polygon::I) where {I}
    intersection_indices = get!(Set{I}, boundary_sites, incident_polygon)
    idx = findfirst(==(intersection_point), segment_intersections)
    if idx === nothing
        push_point!(segment_intersections, intersection_point)
        idx = length(segment_intersections)
    end
    push!(intersection_indices, idx)
    return nothing
end

function process_ray_intersection!(
    vorn::VoronoiTessellation,
    u,
    v,
    incident_polygon,
    num_intersections,
    intersected_edge_cache,
    segment_intersections,
    boundary_site_additions)
    u_tri = get_circumcenter_to_triangle(vorn, u)
    a = geti(u_tri)
    b = getj(u_tri)
    p, q = get_generator(vorn, a, b)
    r = get_polygon_point(vorn, v)
    intersection_coordinates = intersection_of_edge_and_bisector_ray(p, q, r)
    any(isnan, intersection_coordinates) && return num_intersections
    add_segment_intersection!(segment_intersections, boundary_site_additions, intersection_coordinates, incident_polygon)
    num_intersections += 1
    E = edge_type(vorn)
    push!(intersected_edge_cache, construct_edge(E, a, b))
    return num_intersections
end

function process_segment_intersection!(
    vorn,
    u,
    v,
    e,
    incident_polygon,
    num_intersections,
    intersected_edge_cache,
    segment_intersections,
    boundary_site_additions)
    a, b = edge_indices(e)
    tri = get_triangulation(vorn)
    if !is_boundary_edge(tri, e)
        a, b = b, a
    end
    p, q = get_generator(vorn, a, b)
    r, s = get_polygon_point(vorn, u, v)
    intersection_cert, cert_u, cert_v, intersection_coordinates = classify_and_compute_segment_intersection(p, q, r, s)
    if is_none(intersection_cert) || is_touching(intersection_cert)
        return num_intersections
    end
    add_segment_intersection!(segment_intersections, boundary_site_additions, intersection_coordinates, incident_polygon)
    num_intersections += 1
    E = edge_type(vorn)
    push!(intersected_edge_cache, construct_edge(E, a, b))
    return num_intersections
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
    processed_pairs = Set{E}()
    intersected_edge_cache = Dict{E,E}()
    sizehint!(intersected_edge_cache, 2^3)
    exterior_circumcenters = Set{I}()
    return edges_to_process, polygon_edge_queue, boundary_sites, segment_intersections, processed_pairs, intersected_edge_cache, exterior_circumcenters
end

function enqueue_new_edge!(polygon_edge_queue, vorn::VoronoiTessellation, edges_to_process)
    e = first(edges_to_process)
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

function process_polygon!(vorn, e, incident_polygon, edges_to_process, polygon_edge_queue, boundary_sites, segment_intersections, processed_pairs, intersected_edge_cache)
    left_edge, right_edge = get_neighbouring_edges(vorn, e)
    polygon_vertices = get_polygon(vorn, incident_polygon)
    nedges = num_boundary_edges(polygon_vertices)
    for ℓ in 1:nedges
        u = get_boundary_nodes(polygon_vertices, ℓ)
        v = get_boundary_nodes(polygon_vertices, ℓ + 1)
        if is_segment_between_two_ghosts(u, v)
            continue
        elseif is_ray_going_in(u, v)
            process_ray_intersection!(vorn, u, v, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites)
        elseif is_ray_going_out(u, v)
            process_ray_intersection!(vorn, v, u, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites)
        elseif is_finite_segment(u, v)
            for _e in (e, left_edge, right_edge)
                process_segment_intersection!(vorn, u, v, _e, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites)
            end
        end
    end
    process_intersection_points!(polygon_edge_queue, vorn, intersected_edge_cache, e, incident_polygon)
    return nothing
end

function process_intersection_points!(polygon_edge_queue, vorn, intersected_edge_cache, e, incident_polygon)
    (uv₁, e₁), (uv₂, e₂) = intersected_edge_cache
    length(intersected_edge_cache) > 2 && throw("...")
    if compare_unoriented_edge(e₁, e₂)
        # Intersection points are on the same edge 
        u₁, v₁ = edge_indices(uv₁)
        u₂, v₂ = edge_indices(uv₂)
        adjacent_polygon_1 = get_adjacent(vorn, v₁, u₁)
        adjacent_polygon_2 = get_adjacent(vorn, v₂, u₂)
        enqueue!(polygon_edge_queue, (e₁, adjacent_polygon_1))
        enqueue!(polygon_edge_queue, (e₁, adjacent_polygon_2))
    else
        i, j = edge_indices(e)
        i₁, j₁ = edge_indices(e₁)
        i₂, j₂ = edge_indices(e₂)
        if i₁ == i
        end
    end
end

function clip_voronoi_tessellation!(vorn::VoronoiTessellation)
    edges_to_process,
    polygon_edge_queue,
    boundary_sites,
    segment_intersections,
    processed_pairs,
    intersected_edge_cache,
    exterior_circumcenters = initialise_clipping_arrays(vorn)
    enqueue_new_edge!(polygon_edge_queue, vorn, edges_to_process)
    while !isempty(edges_to_process)
        isempty(polygon_edge_queue) && enqueue_new_edge!(polygon_edge_queue, vorn, edges_to_process)
        empty!(intersected_edge_cache)
        e, incident_polygon = dequeue!(polygon_edge_queue)
        if e ∉ processed_pairs && reverse_edge(e) ∉ processed_pairs
            push!(processed_pairs, e)
            process_polygon!(vorn, e, incident_polygon, edges_to_process, polygon_edge_queue, boundary_sites, segment_intersections, processed_pairs, intersected_edge_cache)
            if contains_edge(e, edges_to_process)
                delete!(edges_to_process, e)
            elseif contains_edge(reverse_edge(e), edges_to_process)
                delete!(edges_to_process, reverse_edge(e))
            end
        end
    end
end