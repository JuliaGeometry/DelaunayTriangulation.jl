"""
    initialise_voronoi_tessellation(tri::Triangulation)

Initialise a `VoronoiTessellation` from the triangulation `tri`.
"""
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

"""
    prepare_add_voronoi_polygon(vorn::VoronoiTessellation, i)

Prepare the arrays for adding the Voronoi polygon for the point `i` to the `VoronoiTessellation` `vorn`, returning:

- `S`: The vertices surrounding the generator `i`.
- `B`: An empty array for storing vertices.
"""
function prepare_add_voronoi_polygon(vorn::VoronoiTessellation, i)
    I = integer_type(vorn)
    S = get_surrounding_polygon(vorn, i)
    B = I[]
    sizehint!(B, length(S))
    return S, B
end

"""
    get_next_triangle_for_voronoi_polygon(vorn::VoronoiTessellation, i, k, S, m)

Get the next triangle for the Voronoi polygon for the point `i` in the `VoronoiTessellation` `vorn`, returning:

- `ci`: The index for the circumcenter of the next triangle.
- `k`: The index of the next vertex in `S`.
"""
function get_next_triangle_for_voronoi_polygon(vorn::VoronoiTessellation, i, k, S, m)
    T = triangle_type(vorn)
    j = k
    k = S[m]
    V = (rotate_triangle_to_standard_form ∘ construct_triangle)(T, i, j, k)
    ci = get_triangle_to_circumcenter(vorn, V)
    return ci, k
end

"""
    connect_circumcenters!(B, ci)

Add the circumcenter `ci` to the array `B`.
"""
function connect_circumcenters!(B, ci)
    push!(B, ci)
    return nothing
end

"""
    add_edge_to_voronoi_polygon!(B, vorn::VoronoiTessellation, i, k, S, m, encountered_duplicate_circumcenter)

Add the next edge to the Voronoi polygon for the point `i` in the `VoronoiTessellation` `vorn`, returning:

- `ci`: The index for the circumcenter of the triangle considered.
- `encountered_duplicate_circumcenter`: Whether or not a duplicate circumcenter has been encountered.
- `k`: The index of the latest vertex in `S`.
"""
function add_edge_to_voronoi_polygon!(B, vorn::VoronoiTessellation, i, k, S, m, encountered_duplicate_circumcenter)
    ci, k = get_next_triangle_for_voronoi_polygon(vorn, i, k, S, m)
    is_boundary_index(ci) && add_unbounded_polygon!(vorn, i)
    (encountered_duplicate_circumcenter || ci ∈ get_cocircular_circumcenters(vorn)) && (encountered_duplicate_circumcenter = true)
    connect_circumcenters!(B, ci)
    return ci, encountered_duplicate_circumcenter, k
end

"""
    close_voronoi_polygon!(vorn::VoronoiTessellation, B, i, encountered_duplicate_circumcenter, prev_ci)

Close the Voronoi polygon for the point `i` in the `VoronoiTessellation` `vorn`.
"""
function close_voronoi_polygon!(vorn::VoronoiTessellation, B, i, encountered_duplicate_circumcenter, prev_ci)
    encountered_duplicate_circumcenter && unique!(B)
    connect_circumcenters!(B, B[begin])
    add_adjacent!(vorn, prev_ci, B[begin], i)
    add_polygon!(vorn, B, i)
    return nothing
end

"""
    add_voronoi_polygon!(vorn::VoronoiTessellation, i)

Add the Voronoi polygon for the point `i` to the `VoronoiTessellation` `vorn`.
""" 
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