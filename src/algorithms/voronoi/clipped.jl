"""
    get_neighbouring_boundary_edges(tri::Triangulation, e) -> (Edge, Edge)

Returns the two boundary edges adjacent to the boundary edge `e` in the triangulation `tri`.

# Arguments 
- `tri::Triangulation`: a triangulation
- `e`: The boundary edge.

# Outputs 
- `left_e`: The left edge.
- `right_e`: The right edge.
"""
function get_neighbouring_boundary_edges(tri::Triangulation, e)
    e = convert_to_edge_adjoining_ghost_vertex(tri, e)
    u, v = edge_vertices(e)
    bnd_idx = get_adjacent(tri, u, v)
    right_bnd = get_right_boundary_node(tri, u, bnd_idx)
    left_bnd = get_left_boundary_node(tri, v, bnd_idx)
    E = edge_type(tri)
    left_e = construct_edge(E, v, left_bnd)
    right_e = construct_edge(E, right_bnd, u)
    return left_e, right_e
end

"""
    convert_to_edge_adjoining_ghost_vertex(tri::Triangulation, e) -> Edge

Returns the edge `e` if it is not a boundary edge, and the edge `reverse(e)` if it is a boundary edge. 

See also [`is_boundary_edge`](@ref).
"""
function convert_to_edge_adjoining_ghost_vertex(tri::Triangulation, e) # used to be convert_to_edge_adjoining_ghost_vertex, but after fixing the definition of is_boundary_edge this is more appropriate. We could fix the definition accordingly, but that makes some of the necessary changes to process_segment_intersection! a bit annoying
    if is_boundary_edge(tri, e)
        return reverse_edge(e)
    else
        return e
    end
end

"""
    get_shared_vertex(e, f) -> Vertex

Returns the vertex shared by the edges `e` and `f`, or `∅` if they do not share a vertex.

# Arguments 
- `e`: The first edge.
- `f`: The second edge.

# Outputs
- `u`: The shared vertex.

# Example
```jldoctest
julia> using DelaunayTriangulation

julia> DelaunayTriangulation.get_shared_vertex((1, 3), (5, 7))
0

julia> DelaunayTriangulation.get_shared_vertex((1, 3), (3, 7))
3

julia> DelaunayTriangulation.get_shared_vertex((10, 3), (10, 5))
10

julia> DelaunayTriangulation.get_shared_vertex((9, 4), (9, 5))
9
```
"""
function get_shared_vertex(e, f)
    u, v = edge_vertices(e)
    w, x = edge_vertices(f)
    if u == w || u == x
        return u
    elseif v == w || v == x
        return v
    else
        I = typeof(u)
        return I(∅)
    end
end

"""
    add_segment_intersection!(segment_intersections, boundary_sites, intersection_point, incident_polygon::I) where {I} -> Integer

Adds the `intersection_point` into the list of `segment_intersections`.

# Arguments 
- `segment_intersections`: The list of segment intersections.
- `boundary_sites`: A mapping from boundary sites to the indices of the segment intersections that are incident to the boundary site.
- `intersection_point`: The intersection point to add.
- `incident_polygon`: The index of the polygon that is incident to the intersection point.

# Outputs
- `idx`: The index of the intersection point in the list of segment intersections. iF the intersection point already exists in the list, then the index of the existing point is returned and used instead.
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

"""
    add_to_intersected_edge_cache!(intersected_edge_cache, u, v, a, b)

Add the edge `uv` to the list of intersected edges.

# Arguments 
- `intersected_edge_cache`: The list of intersected edges.
- `u`: The first vertex of the edge of the Voronoi polygon intersecting the edge `ab` of the boundary. 
- `v`: The second vertex of the edge of the Voronoi polygon intersecting the edge `ab` of the boundary.
- `a`: The first vertex of the edge of the boundary.
- `b`: The second vertex of the edge of the boundary.

# Outputs
There are no outputs, as `intersected_edge_cache` is modified in-place.
"""
function add_to_intersected_edge_cache!(intersected_edge_cache::AbstractVector{V}, u, v, a, b) where {E,V<:Pair{E,E}}
    uv = construct_edge(E, u, v)
    ab = construct_edge(E, a, b)
    push!(intersected_edge_cache, uv => ab)
    return intersected_edge_cache
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
        equal_circumcenter_mapping,
        predicates::AbstractPredicateKernel=AdaptiveKernel()) -> Point

Process the intersection of the Voronoi polygon of the site `u` with the ray emanating from the circumcenter of the site `v`.

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `u`: The index of the site `u`, given as a ghost vertex for the associated ghost triangle.
- `v`: The index of the site `v`.
- `incident_polygon`: The index of the Voronoi polygon of the site `u` that is incident to the ray emanating from the circumcenter of the site `v`.
- `intersected_edge_cache`: The list of intersected edges currently being considered.
- `segment_intersections`: The list of segment intersections.
- `boundary_sites`: A mapping from boundary sites to the indices of the segment intersections that are incident to the boundary site.
- `exterior_circumcenters`: The list of circumcenters of sites that are outside the boundary.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs 
- `p`: The coordinates of the intersection. 

In addition to the point `p`, [`add_segment_intersection!`](@ref) is also updated to incorporate the new intersection point, as is [`add_to_intersected_edge_cache!`](@ref).
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
    equal_circumcenter_mapping,
    predicates::AbstractPredicateKernel=AdaptiveKernel())
    u_tri = get_circumcenter_to_triangle(vorn, u)
    a = geti(u_tri)
    b = getj(u_tri)
    p, q = get_generator(vorn, a, b)
    r = get_polygon_point(vorn, v)
    _, intersection_coordinates = intersection_of_edge_and_bisector_ray(predicates, p, q, r)
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
        equal_circumcenter_mapping,
        predicates::AbstractPredicateKernel=AdaptiveKernel()) -> Point 

Process the intersection of the Voronoi polygon's edge `(u, v)` with the edge `e` of the boundary, returning the coordinates of the intersection and updating via [`add_segment_intersection!`](@ref).

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `u`: The index of the site `u`.
- `v`: The index of the site `v`.
- `e`: The edge `e` of the boundary.
- `incident_polygon`: The index of the Voronoi polygon currently being considered.
- `intersected_edge_cache`: The list of intersected edges currently being considered.
- `segment_intersections`: The list of segment intersections.
- `boundary_sites`: A mapping from boundary sites to the indices of the segment intersections that are incident to the boundary site.
- `exterior_circumcenters`: The list of circumcenters of sites that are outside the boundary.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
- `p`: The coordinates of the intersection. If there is no intersection, this is `(NaN, NaN)`.

In addition to the point `p`, [`add_segment_intersection!`](@ref) is also updated to incorporate the new intersection point, as is [`add_to_intersected_edge_cache!`](@ref).
"""
function process_segment_intersection!(
    vorn::VoronoiTessellation,
    u,
    v,
    e,
    incident_polygon,
    intersected_edge_cache,
    segment_intersections,
    boundary_sites,
    exterior_circumcenters,
    equal_circumcenter_mapping,
    predicates::AbstractPredicateKernel=AdaptiveKernel())
    e = convert_to_edge_adjoining_ghost_vertex(vorn, e)
    a, b = edge_vertices(e)
    p, q = get_generator(vorn, a, b)
    r, s = get_polygon_point(vorn, u, v)
    intersection_cert, cert_u, cert_v, intersection_coordinates = classify_and_compute_segment_intersection(predicates, p, q, r, s)
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
    initialise_clipping_arrays(vorn::VoronoiTessellation) -> (Set{E}, Queue{Tuple{E,I}}, Dict{I,Set{I}}, NTuple{2,F}[], Set{Tuple{E,I}}, Pair{E,E}[], Set{I}, Set{E}, Set{E}, Set{E}, Dict{I,I})

Initialise the arrays used in the clipping algorithm for the [`VoronoiTessellation`](@ref) `vorn`.

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).

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

# Arguments
- `polygon_edge_queue`: The queue of edges that are to be processed.
- `vorn`: The [`VoronoiTessellation`](@ref).
- `e`: The edge to be processed.
- `rng::AbstractRNG=Random.default_rng()`: Random number generator. Needed for [`get_nearest_neighbour`](@ref).
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods. Needed for [`get_nearest_neighbour`](@ref).

# Outputs
There are no outputs, as `polygon_edge_queue` is modified in-place.
"""
function enqueue_new_edge!(polygon_edge_queue, vorn::VoronoiTessellation, e, rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())
    u, v = edge_vertices(e)
    p, q = get_generator(vorn, u, v)
    m = midpoint(p, q)
    incident_polygon = get_nearest_neighbour(vorn, m; rng, predicates,k=u)
    push!(polygon_edge_queue, (e, incident_polygon))
    return polygon_edge_queue
end

"""
    is_segment_between_two_ghosts(u, v) -> Bool

Returns `true` if the segment `(u, v)` is between two ghost vertices, and `false` otherwise.
"""
is_segment_between_two_ghosts(u, v) = is_ghost_vertex(u) && is_ghost_vertex(v)

"""
    is_ray_going_in(u, v) -> Bool

Returns `true` if the ray `(u, v)` is coming in from infinity, and `false` otherwise.
"""
is_ray_going_in(u, v) = is_ghost_vertex(u) && !is_ghost_vertex(v)

"""
    is_ray_going_out(u, v) -> Bool

Returns `true` if the ray `(u, v)` is going out to infinity, and `false` otherwise.
"""
is_ray_going_out(u, v) = !is_ghost_vertex(u) && is_ghost_vertex(v)

"""
    is_finite_segment(u, v) -> Bool

Returns `true` if the segment `(u, v)` is finite, and `false` otherwise.
"""
is_finite_segment(u, v) = !is_ghost_vertex(u) && !is_ghost_vertex(v)

"""
    process_ray_intersection_with_other_edges!(vorn::VoronoiTessellation,
        u,
        v,
        e,
        left_edge,
        right_edge,
        r,
        segment_intersections,
        boundary_sites,
        incident_polygon,
        equal_circumcenter_mapping,
        intersected_edge_cache,
        predicates::AbstractPredicateKernel=AdaptiveKernel())

Process the intersection of the ray from the ghost site `u` to the site `v` with the edges `e`, `left_edge` and `right_edge`.

# Arguments
- `vorn::VoronoiTessellation`: The [`VoronoiTessellation`](@ref).
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
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
There are no outputs, but [`add_segment_intersection!`](@ref) and [`add_to_intersected_edge_cache!`](@ref) are used to update the intersection objects.
"""
function process_ray_intersection_with_other_edges!(vorn::VoronoiTessellation,
    u,
    v,
    e,
    left_edge,
    right_edge,
    r,
    segment_intersections,
    boundary_sites,
    incident_polygon,
    equal_circumcenter_mapping,
    intersected_edge_cache,
    predicates::AbstractPredicateKernel=AdaptiveKernel())
    if !any(isnan, r)
        E = edge_type(vorn)
        u_tri = get_circumcenter_to_triangle(vorn, u)
        a, b, _ = triangle_vertices(u_tri)
        s = get_polygon_point(vorn, v)
        intersected_edge = construct_edge(E, a, b)
        for _e in (e, left_edge, right_edge)
            if !compare_unoriented_edges(intersected_edge, _e)
                i, j = edge_vertices(_e)
                p, q = get_generator(vorn, i, j)
                intersection_cert, cert_u, cert_v, intersection_coordinates = classify_and_compute_segment_intersection(predicates, p, q, r, s)
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
    return vorn
end

"""
    process_polygon!(vorn::VoronoiTessellation, e, incident_polygon, boundary_sites, segment_intersections, intersected_edge_cache, exterior_circumcenters, equal_circumcenter_mapping, predicates::AbstractPredicateKernel=AdaptiveKernel()) -> (Edge, Edge, Edge)

Processes the polygon `incident_polygon` for all of its intersections based on the boundary edge `e`.

# Arguments
- `vorn::VoronoiTessellation`: The [`VoronoiTessellation`](@ref).
- `e`: The edge on the boundary being considered.
- `incident_polygon`: The index of the polygon being considered.
- `boundary_sites`: The mapping from the indices of the sites on the boundary to the indices of the edges on the boundary that they intersect.
- `segment_intersections`: The list of segment intersections.
- `intersected_edge_cache`: A cache of the edges that have been intersected by the ray from `u` to `v`.
- `exterior_circumcenters`: A list of the circumcenters of the sites that are outside the convex hull of the sites on the boundary.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
- `left_edge`: The edge to the left of `e` on the boundary.
- `right_edge`: The edge to the right of `e` on the boundary.
- `e`: The edge on the boundary being considered.

In addition to these outputs, the caches are also updated in-place.

# Extended help 
This function works as follows: 

1. First, for the current edge `e`, we get the edges `left_edge` and `right_edge` that neighbour it via [`get_neighbouring_boundary_edges`](@ref).
2. For each edge of the `incident_polygon`, we need to process it depending on whether the edge `(u, v)` is finite, between two ghosts, going out to infinity, or coming in from infinity. 
   If the edge is between two ghosts, we skip the edge. For rays that go out or in to infinity, we use [`process_ray_intersection!`](@ref) and [`process_ray_intersection_with_other_edges!`](@ref) 
   to process the intersection of the ray with the boundary edges. The function [`process_ray_intersection_with_other_edges!`](@ref) is needed since rays going out to infinity may have to go 
   through other boundary edges in order to do so, e.g. at a corner it may be that it crosses two boundary edges. For finite segments, [`process_segment_intersection!`](@ref) is used to process
   the intersection. We apply this function with each of `e`, `left_edge`, and `right_edge` to check for all intersections.
3. The function is done once each of the polygon edges has been considered.
"""
function process_polygon!(vorn, e, incident_polygon, boundary_sites, segment_intersections, intersected_edge_cache, exterior_circumcenters, equal_circumcenter_mapping, predicates::AbstractPredicateKernel=AdaptiveKernel())
    left_edge, right_edge = get_neighbouring_boundary_edges(vorn, e)
    polygon_vertices = get_polygon(vorn, incident_polygon)
    nedges = num_boundary_edges(polygon_vertices)
    for ℓ in 1:nedges
        u = get_boundary_nodes(polygon_vertices, ℓ)
        v = get_boundary_nodes(polygon_vertices, ℓ + 1)
        if is_segment_between_two_ghosts(u, v)
            continue
        elseif is_ray_going_in(u, v)
            r = process_ray_intersection!(vorn, u, v, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, predicates)
            # It's possible for an infinite ray to also intersect other boundary edges, e.g. look at 
            #   points = [0.290978 0.830755 0.0139574; 0.386411 0.630008 0.803881]
            # So, let's just look for intersections with other edges.
            process_ray_intersection_with_other_edges!(vorn, u, v, e, left_edge, right_edge, r, segment_intersections, boundary_sites, incident_polygon, equal_circumcenter_mapping, intersected_edge_cache, predicates)
        elseif is_ray_going_out(u, v)
            r = process_ray_intersection!(vorn, v, u, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, predicates)
            process_ray_intersection_with_other_edges!(vorn, v, u, e, left_edge, right_edge, r, segment_intersections, boundary_sites, incident_polygon, equal_circumcenter_mapping, intersected_edge_cache, predicates)
        elseif is_finite_segment(u, v)
            for _e in (e, left_edge, right_edge)
                process_segment_intersection!(vorn, u, v, _e, incident_polygon, intersected_edge_cache, segment_intersections, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, predicates)
            end
        end
    end
    return left_edge, right_edge, e
end

"""
    classify_intersections!(intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, left_edge, right_edge, current_edge)

Classify the intersections in `intersected_edge_cache` into `left_edge_intersectors`, `right_edge_intersectors`, and `current_edge_intersectors` based on whether they intersect `left_edge`, `right_edge`, or `current_edge`, respectively.

# Arguments
- `intersected_edge_cache`: The list of intersected edges currently being considered.
- `left_edge_intersectors`: The set of sites that intersect the edge to the left of an edge currently being considered.
- `right_edge_intersectors`: The set of sites that intersect the edge to the right of an edge currently being considered.
- `current_edge_intersectors`: The set of sites that intersect the current edge being considered.
- `left_edge`: The edge to the left of `e` on the boundary.
- `right_edge`: The edge to the right of `e` on the boundary.
- `current_edge`: The edge on the boundary being considered.

# Outputs
There are no outputs, but `left_edge_intersectors`, `right_edge_intersectors`, or `current_edge_intersectors` are updated all in-place depending on the type of intersection for each 
edge in `intersected_edge_cache`.
"""
function classify_intersections!(intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, left_edge, right_edge, current_edge)
    for (uv, e) in intersected_edge_cache
        if compare_unoriented_edges(e, left_edge)
            push!(left_edge_intersectors, uv)
        elseif compare_unoriented_edges(e, right_edge)
            push!(right_edge_intersectors, uv)
        elseif compare_unoriented_edges(e, current_edge)
            push!(current_edge_intersectors, uv)
        end
    end
    return intersected_edge_cache
end

"""
    process_intersection_points!(polygon_edge_queue, vorn, current_incident_polygon,
        left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
        left_edge, right_edge, current_edge, processed_pairs, segment_intersections, boundary_sites)

Process the intersection points in `left_edge_intersectors`, `right_edge_intersectors`, and `current_edge_intersectors` and add the new edges to `polygon_edge_queue` if necessary. Special
care is taken to not miss any corner points.

# Arguments 
- `polygon_edge_queue`: The queue of edges that need to be processed.
- `vorn`: The [`VoronoiTessellation`](@ref).
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

# Outputs
There are no outputs, but the caches and queues are updated in-place.

# Extended help 
The rules are based on the paper "Efficient Computation of Clipped Voronoi Diagram for Mesh Generation" by Yan, Wang, Levy, and Liu. Namely, 
an edge that intersects a boundary edge and one adjacent to it has its shared vertex added to the queue together with the current polygon (`current_incident_polygon`) being 
considered, and any intersections have the adjacent polygon added to the queue together with the intersecting edge. (These are not strictly 
the rules in the paper.)

This function works as follows:

1. First, assuming that there is more than one triangle in the underlying triangulation of `vorn`, we need to consider `left_edge` and `right_edge` individually. 
2. The procedure for each edge is the same, so here we just describe the `left_edge`. If there are any intersectors with the `left_edge`, and neither of 
   `(left_edge, current_incident_polygon)` or `(reverse_edge(left_edge), current_incident_polygon)` have already been processed (i.e., in `processed_pairs`), then we enqueue 
   `(left_edge, i)` and `(left_edge, j)` into `polygon_edge_queue`, where `i` and `j` are the vertices of `left_edge` which correspond to polygons. This will ensure that we can find intersections next to this polygon. 
3. After enqueueing these pairs, we also need to protect against corner points, which we check for by considering `current_incident_polygon ∈ all_indices`, where `all_indices` are the vertices of `left_edge`, 
   `right_edge`, and `current_edge`. If this is true, and if the shared vertex of `current_edge` and `left_edge` is equal to `current_incident_polygon`, then we need to add the point generator of `current_incident_polygon`
   as an intersection. This need comes from having to worry about corners, i.e. points where the two unbounded polygons meet and go directly left and right of a vertex so that that vertex is not considered an intersection;
   this point needs to be included.
4. Once the `left_edge` and `right_edge` have been processed as above, we need to then consider all of `left_edge`, `right_edge`, and `current_edge`, and each of the intersections through the respective edge. This step is done 
   regardless of whether there is a single triangle in the underlying triangulation. The procedure for each edge is the same, so let us just describe the `current_edge`. For each edge `uv` in the `current_edge_intersectors`,
   we need to get the polygon adjacent to that edge. Then, if `(current_edge, adjacent_incident_polygon)` or `(reverse_edge(current_edge), adjacent_incident_polygon)` have not been processed, we enqueue `(current_edge, adjacent_incident_polygon)`.
5. Once the edges have all been processed as above, we return.
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
                i, j = edge_vertices(e)
                push!(polygon_edge_queue, (e, i))
                push!(polygon_edge_queue, (e, j))
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
            u, v = edge_vertices(uv)
            adjacent_incident_polygon = get_adjacent(vorn, v, u)
            if adjacent_incident_polygon == current_incident_polygon
                adjacent_incident_polygon = get_adjacent(vorn, u, v)
            end
            if (e, adjacent_incident_polygon) ∉ processed_pairs && (reverse_edge(e), adjacent_incident_polygon) ∉ processed_pairs
                push!(polygon_edge_queue, (e, adjacent_incident_polygon))
            end
        end
    end
    return polygon_edge_queue
end

"""
    dequeue_and_process!(vorn, polygon_edge_queue, edges_to_process,
        intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
        processed_pairs, boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping,
        rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())

Dequeue an edge from `polygon_edge_queue` and process it. If `polygon_edge_queue` is empty, then we process the first edge in `edges_to_process`.

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).
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
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator. 
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
There are no outputs. Instead, the caches and queues are updated in-place.

# Extended help
This function works as follows:

1. Firstly, if there are no edges queued in `polygon_edge_queue`, then we enqueue the first in edge in `edges_to_process` using [`enqueue_new_edge!`](@ref).
2. We then dequeue the next edge to be processed. If the edge has already been processed, then we return early.
3. If we're still here, then we process the `(edge, polygon)` pair enqueued from `polygon_edge_queue` using [`process_polygon!`](@ref). This function checks for intersections of the `edge` with the `polygon`.
4. Once the polygon has been processed, we then needed to classify all of the intersections using [`classify_intersections!`](@ref), which determines, for each intersection, if the intersection is with `edge`, 
   or with the edge left of `edge`, or to the edge right of `edge`.
5. Then, [`process_intersection_points!`](@ref) is used to process the intersection points, enqueueing new edges when needed.
6. We then delete the edge from `edges_to_process` if it is in there and return.
"""
function dequeue_and_process!(vorn, polygon_edge_queue, edges_to_process,
    intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
    processed_pairs, boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping,
    rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())
    if isempty(polygon_edge_queue)
        e = convert_to_edge_adjoining_ghost_vertex(vorn, first(edges_to_process))
        enqueue_new_edge!(polygon_edge_queue, vorn, e, rng, predicates)
    end
    e, incident_polygon = popfirst!(polygon_edge_queue)
    if (e, incident_polygon) ∈ processed_pairs || (reverse_edge(e), incident_polygon) ∈ processed_pairs
        return vorn
    end
    push!(processed_pairs, (e, incident_polygon))
    for cache in (intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors)
        empty!(cache)
    end
    left_edge, right_edge, e = process_polygon!(vorn, e, incident_polygon, boundary_sites, segment_intersections, intersected_edge_cache, exterior_circumcenters, equal_circumcenter_mapping, predicates)
    classify_intersections!(intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors, left_edge, right_edge, e)
    process_intersection_points!(polygon_edge_queue, vorn, incident_polygon,
        left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
        left_edge, right_edge, e, processed_pairs, segment_intersections, boundary_sites)
    if contains_edge(e, edges_to_process)
        delete!(edges_to_process, e)
    elseif contains_edge(reverse_edge(e), edges_to_process)
        delete!(edges_to_process, reverse_edge(e))
    end
    return vorn
end

"""
    find_all_intersections(vorn::VoronoiTessellation; rng=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel()) -> (Dict, Vector, Set, Dict)

Find all intersections between the edges of the Voronoi tessellation and the boundary of the polygon.

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).

# Keyword Arguments 
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs
- `boundary_sites`: A dictionary of boundary sites.
- `segment_intersections`: The intersection points.
- `exterior_circumcenters`: The circumcenters that are outside of the domain.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.

# Extended help 
This algorithm works as follows:

1. First, using [`initialise_clipping_arrays`](@ref), we initialise the arrays that we will use to store the intersections, and queue up all boundary edges for processing.
2. Then, starting with the first edge in `edges_to_process`, we dequeue an edge from `polygon_edge_queue` and process it via [`dequeue_and_process!`](@ref).
3. We repeat step 2 until `polygon_edge_queue` and `edges_to_process` are both empty.
4. In the special case that there is just a single triangle in the underlying triangulation, we process the intersections using [`add_segment_intersection!`](@ref) directly.
5. We then return.
"""
function find_all_intersections(vorn::VoronoiTessellation; rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())
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
    e = convert_to_edge_adjoining_ghost_vertex(vorn, first(edges_to_process))
    enqueue_new_edge!(polygon_edge_queue, vorn, e, rng, predicates)
    while !isempty(edges_to_process) || !isempty(polygon_edge_queue)
        dequeue_and_process!(vorn, polygon_edge_queue, edges_to_process,
            intersected_edge_cache, left_edge_intersectors, right_edge_intersectors, current_edge_intersectors,
            processed_pairs, boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping,
            rng, predicates)
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
    add_intersection_points!(vorn::VoronoiTessellation, segment_intersections) -> Integer

Adds all of the `segment_intersections` into the polygon vertices of `vorn`.

# Arguments 
- `vorn`: The [`VoronoiTessellation`](@ref).
- `segment_intersections`: The intersection points from [`find_all_intersections`](@ref).

# Outputs
- `n`: The number of polygon vertices before the intersections were added.
"""
function add_intersection_points!(vorn::VoronoiTessellation, segment_intersections)
    n = num_polygon_vertices(vorn)
    for p in each_point(segment_intersections)
        push_polygon_point!(vorn, p)
    end
    return n
end

"""
    clip_polygon!(vorn::VoronoiTessellation, n, points, polygon, new_verts, exterior_circumcenters, equal_circumcenter_mapping, is_convex)

Clip the polygon `polygon` by removing the vertices that are outside of the domain and adding the new vertices `new_verts` to the polygon.

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).
- `n`: The number of vertices in the tessellation before clipping.
- `points`: The polygon points of the tessellation.
- `polygon`: The index of the polygon to be clipped.
- `new_verts`: The indices of the new vertices that are added to the polygon.
- `exterior_circumcenters`: Any exterior circumcenters to be filtered out.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
- `is_convex`: Whether the boundary is convex or not. Not currently used.

# Outputs
There are no outputs, but the polygon is clipped in-place.
"""
function clip_polygon!(vorn::VoronoiTessellation, n, points, polygon, new_verts, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    delete_polygon_adjacent!(vorn, polygon)
    vertices = get_polygon(vorn, polygon)
    pop!(vertices) # vertices[begin] == vertices[end]
    filter!(v -> !is_ghost_vertex(v) && v ∉ exterior_circumcenters, vertices)
    for new_vert in new_verts
        if new_vert ∉ keys(equal_circumcenter_mapping) || equal_circumcenter_mapping[new_vert] ∉ vertices
            push!(vertices, n + new_vert)
        end
    end
    sort_convex_polygon!(vertices, points)
    push!(vertices, vertices[begin])
    add_polygon_adjacent!(vorn, polygon)
    delete_unbounded_polygon!(vorn, polygon)
    return vorn
end

"""
    clip_all_polygons!(vorn::VoronoiTessellation, n, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, is_convex)

Clip all of the polygons in the Voronoi tessellation.

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).
- `n`: The number of vertices in the tessellation before clipping.
- `boundary_sites`: A dictionary of boundary sites.
- `exterior_circumcenters`: Any exterior circumcenters to be filtered out.
- `equal_circumcenter_mapping`: A mapping from the indices of the segment intersections that are equal to the circumcenter of a site to the index of the site.
- `is_convex`: Whether the boundary is convex or not. Not currently used.

# Outputs
There are no outputs, but the polygons are clipped in-place.
"""
function clip_all_polygons!(vorn::VoronoiTessellation, n, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    points = get_polygon_points(vorn)
    for (polygon, new_verts) in boundary_sites
        clip_polygon!(vorn, n, points, polygon, new_verts, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    end
    return vorn
end

"""
    add_all_boundary_polygons!(vorn::VoronoiTessellation, boundary_sites)

Add all of the boundary polygons to the Voronoi tessellation.

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).
- `boundary_sites`: A dictionary of boundary sites.

# Outputs
There are no outputs, but the boundary polygons are added in-place.
"""
function add_all_boundary_polygons!(vorn::VoronoiTessellation, boundary_sites)
    for i in keys(boundary_sites)
        add_boundary_polygon!(vorn, i)
    end
    return vorn
end

"""
    clip_voronoi_tessellation!(vorn::VoronoiTessellation, is_convex=true; rng=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())

Clip the Voronoi tessellation `vorn` to the convex hull of the generators in `vorn`. 

# Arguments
- `vorn`: The [`VoronoiTessellation`](@ref).
- `is_convex`: Whether the boundary is convex or not. Not currently used.

# Keyword Arguments
- `rng::Random.AbstractRNG=Random.default_rng()`: The random number generator.
- `predicates::AbstractPredicateKernel=AdaptiveKernel()`: Method to use for computing predicates. Can be one of [`FastKernel`](@ref), [`ExactKernel`](@ref), and [`AdaptiveKernel`](@ref). See the documentation for a further discussion of these methods.

# Outputs 
There are no outputs, but the Voronoi tessellation is clipped in-place.
"""
function clip_voronoi_tessellation!(vorn::VoronoiTessellation, is_convex=true; rng::Random.AbstractRNG=Random.default_rng(), predicates::AbstractPredicateKernel=AdaptiveKernel())
    boundary_sites, segment_intersections, exterior_circumcenters, equal_circumcenter_mapping = find_all_intersections(vorn; rng, predicates)
    n = add_intersection_points!(vorn, segment_intersections)
    clip_all_polygons!(vorn, n, boundary_sites, exterior_circumcenters, equal_circumcenter_mapping, is_convex)
    add_all_boundary_polygons!(vorn, boundary_sites)
    return vorn
end