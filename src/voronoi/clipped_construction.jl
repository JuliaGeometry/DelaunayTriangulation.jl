"""
    add_segment_intersection!(segment_intersections, boundary_sites, intersection_point, incident_polygon::I) where {I}

Add the intersection point to the list of segment intersections and return the index of the intersection point in the list. 
If the intersection point already exists in the list, then the index of the existing point is returned and used instead.

# Arguments
- `segment_intersections`: The list of segment intersections.
- `boundary_sites`: A mapping from boundary sites to the indices of the segment intersections that are incident to the boundary site.
- `intersection_point`: The intersection point to add.
- `incident_polygon`: The index of the polygon that is incident to the intersection point.
"""
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

"""
    add_to_intersected_edge_cache!(intersected_edge_cache, u, v, a, b)

Add the edge `uv` to the list of intersected edges, where `uv` is the edge of the Voronoi polygon intersecting the edge of the boundary `ab`.
"""
function add_to_intersected_edge_cache!(intersected_edge_cache::AbstractVector{V}, u, v, a, b) where {E,V<:Pair{E,E}}
    uv = construct_edge(E, u, v)
    ab = construct_edge(E, a, b)
    push!(intersected_edge_cache, uv => ab)
    return nothing
end

"""
    process_ray_intersection!(
        vorn::VoronoiTessellation,
        u,
        v,
        incident_polygon,
        intersected_edge_cache,
        segment_intersections,
        boundary_sites,
        exterior_circumcenters,
        equal_circumcenter_mapping)

Process the intersection of the Voronoi polygon of the site `u` with the ray emanating from the circumcenter of the site `v`,
returning the coordinates of the intersection and updating via `add_segment_intersection`.

# Arguments
- `vorn`: The Voronoi tessellation.
- `u`: The index of the site `u`, given as a boundary index for the associated ghost triangle.
- `v`: The index of the site `v`.
- `incident_polygon`: The index of the Voronoi polygon of the site `u` that is incident to the ray emanating from the circumcenter of the site `v`.
- `intersected_edge_cache`: The list of intersected edges currently being considered.
- `segment_intersections`: The list of segment intersections.
- `boundary_sites`: A mapping from boundary sites to the indices of the segment intersections that are incident to the boundary site.
- `exterior_circumcenters`: The list of circumcenters of sites that are outside the boundary.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
"""
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

"""
    process_segment_intersection!(
        vorn::VoronoiTessellation,
        u,
        v,
        e,
        incident_polygon,
        intersected_edge_cache,
        segment_intersections,
        boundary_sites,
        exterior_circumcenters,
        equal_circumcenter_mapping)

Process the intersection of the Voronoi polygon's edge `(u, v)` with the edge `e` of the boundary, returning the coordinates of the intersection and updating via `add_segment_intersection`.

# Arguments
- `vorn`: The Voronoi tessellation.
- `u`: The index of the site `u`.
- `v`: The index of the site `v`.
- `e`: The edge `e` of the boundary.
- `incident_polygon`: The index of the Voronoi polygon currently being considered `(u, v)`.
- `intersected_edge_cache`: The list of intersected edges currently being considered.
- `segment_intersections`: The list of segment intersections.
- `boundary_sites`: A mapping from boundary sites to the indices of the segment intersections that are incident to the boundary site.
- `exterior_circumcenters`: The list of circumcenters of sites that are outside the boundary.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
"""
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

"""
    initialise_clipping_arrays(vorn::VoronoiTessellation)

Initialise the arrays used in the clipping algorithm.

# Outputs 
- `edges_to_process`: The set of edges that are to be processed.
- `polygon_edge_queue`: The queue of edges that are to be processed.
- `boundary_sites`: A mapping from boundary sites to the indices of the segment intersections that are incident to the boundary site.
- `segment_intersections`: The list of segment intersections.
- `processed_pairs`: The set of pairs of edges and polygons that have been processed.
- `intersected_edge_cache`: The list of intersected edges currently being considered.
- `exterior_circumcenters`: The list of circumcenters of sites that are outside the boundary.
- `left_edge_intersectors`: The set of sites that intersect the edge to the left of an edge currently being considered.
- `right_edge_intersectors`: The set of sites that intersect the edge to the right of an edge currently being considered.
- `current_edge_intersectors`: The set of sites that intersect the current edge being considered.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
"""
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

"""
    enqueue_new_edge!(polygon_edge_queue, vorn::VoronoiTessellation, e)

Enqueue the edge `e` of the boundary to be processed.
"""
function enqueue_new_edge!(polygon_edge_queue, vorn::VoronoiTessellation, e)
    u, v = edge_indices(e)
    p, q = get_generator(vorn, u, v)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    m = (px + qx) / 2, (py + qy) / 2
    incident_polygon = jump_and_march(vorn, m; k=u)
    enqueue!(polygon_edge_queue, (e, incident_polygon))
    return nothing
end

"""
    is_segment_between_two_ghosts(u, v)

Check if the segment between the sites `u` and `v` is between two ghost sites.
"""
is_segment_between_two_ghosts(u, v) = is_boundary_index(u) && is_boundary_index(v)

"""
    is_ray_going_in(u, v)

Check if the ray from the site `u` to the site `v` is going in.
"""
is_ray_going_in(u, v) = is_boundary_index(u) && !is_boundary_index(v)

"""
    is_ray_going_out(u, v)

Check if the ray from the site `u` to the site `v` is going out.
"""
is_ray_going_out(u, v) = !is_boundary_index(u) && is_boundary_index(v)

"""
    is_finite_segment(u, v)

Check if the segment between the sites `u` and `v` is finite.
"""
is_finite_segment(u, v) = !is_boundary_index(u) && !is_boundary_index(v)

"""
    process_ray_intersection_with_other_edges!(vorn::VoronoiTessellation, u, v, e, left_edge, right_edge, r, segment_intersections,
    boundary_sites, incident_polygon, equal_circumcenter_mapping, intersected_edge_cache)

Process the intersection of the ray from the ghost site `u` to the site `v` with the edges `e`, `left_edge` and `right_edge`.

# Arguments 
- `vorn::VoronoiTessellation`: The Voronoi tessellation.
- `u`: The index of the ghost site.
- `v`: The index of the site `u` is going to.
- `e`: The edge on the boundary being considered.
- `left_edge`: The edge to the left of `e` on the boundary.
- `right_edge`: The edge to the right of `e` on the boundary.
- `r`: The coordinates of the intersection of the ray from `u` to `v` with some edge. If `any(isnan, r)`, then the ray does not intersect any edge and we skip. 
- `segment_intersections`: The list of segment intersections.
- `boundary_sites`: The mapping from the indices of the sites on the boundary to the indices of the edges on the boundary that they intersect.
- `incident_polygon`: The index of the polygon that contains the intersection of the ray from `u` to `v` with the boundary.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
- `intersected_edge_cache`: A cache of the edges that have been intersected by the ray from `u` to `v`.
"""
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

"""
    process_polygon!(vorn::VoronoiTessellation, e, incident_polygon, boundary_sites, segment_intersections, intersected_edge_cache, exterior_circumcenters, equal_circumcenter_mapping)

Process the polygon `incident_polygon` for all of its intersections based on the boundary edge `e`, updating the caches in-place and returning `(left_edge, right_edge, e)`, where 
`left_edge` is the edge to the left of `e` on the boundary and `right_edge` is the edge to the right of `e` on the boundary.

# Arguments 
- `vorn::VoronoiTessellation`: The Voronoi tessellation.
- `e`: The edge on the boundary being considered.
- `incident_polygon`: The index of the polygon being considered.
- `boundary_sites`: The mapping from the indices of the sites on the boundary to the indices of the edges on the boundary that they intersect.
- `segment_intersections`: The list of segment intersections.
- `intersected_edge_cache`: A cache of the edges that have been intersected by the ray from `u` to `v`.
- `exterior_circumcenters`: A list of the circumcenters of the sites that are outside the convex hull of the sites on the boundary.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
"""
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

"""
    classify_intersections!(intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, left_edge, right_edge, current_edge)

Classify the intersections in `intersected_edge_cache` into `left_edge_intersectors`, `right_edge_intersectors`, and `current_edge_intersectors` based on whether they intersect `left_edge`, `right_edge`, or `current_edge`, respectively.
"""
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

"""
    process_intersection_points!(polygon_edge_queue, vorn, current_incident_polygon, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, left_edge, right_edge, current_edge, processed_pairs, segment_intersections, boundary_sites)

Process the intersection points in `left_edge_intersectors`, `right_edge_intersectors`, and `current_edge_intersectors` and add the new edges to `polygon_edge_queue` if necessary. Special 
care is taken to not miss any corner points. 

The rules are based on the paper "Efficient Computation of Clipped Voronoi Diagram for Mesh Generation" by Yan, Wang, Levy, and Liu. Namely, 
an edge that intersects a boundary edge and one to it has its shared vertex added to the queue together with the current polygon (`current_incident_polygon`) being 
considered, and any intersections have the adjacent polygon added to the queue together with the intersecting edge. These are not strictly 
the rules in the paper, but they are the rules that I was able to implement since they do not share their code.

# Arguments 
- `polygon_edge_queue`: The queue of edges that need to be processed.
- `vorn`: The Voronoi tessellation.
- `current_incident_polygon`: The index of the current polygon being processed.
- `left_edge_intersectors`: The intersection points of `left_edge` with other edges.
- `right_edge_intersectors`: The intersection points of `right_edge` with other edges.
- `current_edge_intersectors`: The intersection points of `current_edge` with other edges.
- `left_edge`: The left edge of the current polygon.
- `right_edge`: The right edge of the current polygon.
- `current_edge`: The current edge of the current polygon.
- `processed_pairs`: A set of pairs of edges and polygons that have already been processed.
- `segment_intersections`: A dictionary of segment intersections.
- `boundary_sites`: A dictionary of boundary sites.
"""
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

"""
    dequeue_and_process!(vorn, polygon_edge_queue, edges_to_process, intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, processed_pairs, boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping)

Dequeue an edge from `polygon_edge_queue` and process it. If `polygon_edge_queue` is empty, then we process the first edge in `edges_to_process`.

# Arguments
- `vorn`: The Voronoi tessellation.
- `polygon_edge_queue`: The queue of edges that need to be processed.
- `edges_to_process`: The edges that need to be processed.
- `intersected_edge_cache`: A cache of intersected edges.
- `left_edge_intersectors`: The intersection points of `left_edge` with other edges.
- `right_edge_intersectors`: The intersection points of `right_edge` with other edges.
- `current_edge_intersectors`: The intersection points of `current_edge` with other edges.
- `processed_pairs`: A set of pairs of edges and polygons that have already been processed.
- `boundary_sites`: A dictionary of boundary sites.
- `segment_intersections`: A dictionary of segment intersections.
- `exterior_circumcenters`: A dictionary of exterior circumcenters.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
"""
function dequeue_and_process!(vorn, polygon_edge_queue, edges_to_process,
    intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
    processed_pairs, boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping)
    if isempty(polygon_edge_queue)
        e = convert_to_boundary_edge(vorn, first(edges_to_process))
        enqueue_new_edge!(polygon_edge_queue, vorn, e)
    end
    e, incident_polygon = dequeue!(polygon_edge_queue)
    if (e, incident_polygon) ∈ processed_pairs || (reverse_edge(e), incident_polygon) ∈ processed_pairs
        return nothing
    end
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

"""
    find_all_intersections(vorn::VoronoiTessellation)

Find all intersections between the edges of the Voronoi tessellation and the boundary of the polygon.

# Outputs 
- `boundary_sites`: A dictionary of boundary sites.
- `segment_intersections`: The intersection points. 
- `exterior_circumcenters`: The circumcenters that are outside of the domain.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
"""
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

"""
    add_intersection_points!(vorn::VoronoiTessellation, segment_intersections)

Add the intersection points to the polygon vertices.
"""
function add_intersection_points!(vorn::VoronoiTessellation, segment_intersections)#, equal_circumcenter_mapping)
    n = num_polygon_vertices(vorn)
    for i in each_point_index(segment_intersections)
        p = get_point(segment_intersections, i)
        push_polygon_point!(vorn, p)
    end
    return n
end

"""
    clip_polygon!(vorn::VoronoiTessellation, n, points, polygon, new_verts, exterior_circumcenters, equal_circumcenter_mapping, is_convex)

Clip the polygon `polygon` by removing the vertices that are outside of the domain and adding the new vertices `new_verts` to the polygon.

# Arguments 
- `vorn`: The Voronoi tessellation.
- `n`: The number of vertices in the tessellation before clipping.
- `points`: The points of the tessellation. 
- `polygon`: The index of the polygon to be clipped.
- `new_verts`: The indices of the new vertices that are added to the polygon.
- `exterior_circumcenters`: Any exterior circumcenters to be filtered out. 
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
- `is_convex`: Whether the polygon is convex or not. Not currently used.
"""
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

"""
    clip_all_polygons!(vorn::VoronoiTessellation, n, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, is_convex)

Clip all of the polygons in the Voronoi tessellation.

# Arguments
- `vorn`: The Voronoi tessellation.
- `n`: The number of vertices in the tessellation before clipping.
- `boundary_sites`: A dictionary of boundary sites.
- `exterior_circumcenters`: Any exterior circumcenters to be filtered out.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
- `is_convex`: Whether the polygon is convex or not. Not currently used.
"""
function clip_all_polygons!(vorn::VoronoiTessellation, n, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    points = get_polygon_points(vorn)
    for (polygon, new_verts) in boundary_sites
        clip_polygon!(vorn, n, points, polygon, new_verts, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    end
end

"""
    add_all_boundary_polygons!(vorn::VoronoiTessellation, boundary_sites)

Add all of the boundary polygons to the Voronoi tessellation.
"""
function add_all_boundary_polygons!(vorn::VoronoiTessellation, boundary_sites)
    for i in keys(boundary_sites)
        add_boundary_polygon!(vorn, i)
    end
    return nothing
end