"""
    locate_intersecting_triangles(
        e,
        pts,
        adj,
        adj2v,
        graph::Graph{I},
        boundary_index_ranges,
        representative_point_list,
        boundary_map,
        TriangleType::Type{V},
        check_existence::C=Val(has_multiple_segments(boundary_map)),
        rng::AbstractRNG=Random.default_rng()) where {I,V,Vs,C}

Given an edge `e`, returns a set of 
triangles whose interior intersects the edge. 

# Arguments 
- `e`: The edge to find the intersection of.
- `pts`: The points of the triangulation.
- `adj`: The [`Adjacent`](@ref) data structure.
- `adj2v`: The [`Adjacent2Vertex`](@ref) data structure.
- `graph`: The [`Graph`](@ref) data structure.
- `boundary_index_ranges`: The boundary index ranges from [`construct_boundary_index_ranges`](@ref).
- `representative_point_list`: The representative point list.
- `boundary_map`: The boundary map from [`construct_boundary_node_map`](@ref).
- `TriangleType`: The type of triangle to use.
- `check_existence`: Whether to check for the existence of the edge in the triangulation when using [`get_adjacent`](@ref).
- `rng`: The random number generator to use.

# Outputs 
- `intersecting_triangles`: The set of triangles that intersects `e`.
- `collinear_segments`: The set of segments that are collinear with `e`.
- `left_vertices`: The vertices of the `intersecting_triangles` that are to the left of `e`.
- `right_vertices`: The vertices of the `intersecting_triangles` that are to the right of `e`.
"""
function locate_intersecting_triangles(
    e,
    pts,
    adj,
    adj2v::Adjacent2Vertex{I,Es,E},
    graph::Graph{I},
    boundary_index_ranges,
    representative_point_list,
    boundary_map,
    TriangleType::Type{V},
    check_existence::C=Val(has_multiple_segments(boundary_map)),
    rotate=Val(true),
    rng::AbstractRNG=Random.default_rng()) where {I,V,C,Es,E}
    e = is_true(rotate) ? sort_edge_by_degree(e, graph) : e # faster to start at the minimum degree vertex of the edge
    history = PointLocationHistory{V,E,I}()
    add_left_vertex!(history, initial(e))
    add_right_vertex!(history, initial(e))
    jump_and_march(
        pts,
        adj,
        adj2v,
        graph,
        boundary_index_ranges,
        representative_point_list,
        boundary_map,
        get_point(pts, terminal(e));
        m=nothing,
        k=initial(e),
        TriangleType,
        check_existence,
        store_history=Val(true),
        history,
        rng
    )
    add_left_vertex!(history, terminal(e))
    add_right_vertex!(history, terminal(e))
    intersecting_triangles = history.triangles
    collinear_segments = history.collinear_segments
    bad_indices = history.collinear_point_indices
    left_vertices = history.left_vertices
    right_vertices = history.right_vertices
    reverse!(left_vertices) # counter-clockwise
    if !isempty(collinear_segments)
        fix_segments!(collinear_segments, bad_indices)
        connect_segments!(collinear_segments)
        extend_segments!(collinear_segments, e)
    end
    return intersecting_triangles, collinear_segments, left_vertices, right_vertices
end

function delete_intersected_triangles!(tri, triangles) # don't really _need_ this method, but maybe it makes the code a bit clearer?
    for τ in each_triangle(triangles)
        delete_triangle!(tri, τ; protect_boundary=true)
    end
end