"""
    locate_intersecting_triangles(
        e,
        pts,
        adj,
        adj2v,
        graph::Graph{I},
        boundary_index_ranges,
        boundary_map,
        TriangleType::Type{V},
        check_existence::C=Val(has_multiple_segments(boundary_map)),
        rng::AbstractRNG=Random.default_rng()) where {I,V,Vs,C}

Given an edge `e` of a triangulation, returns a set of 
triangles whose interior intersects the edge. 
"""
function locate_intersecting_triangles(
    e,
    pts,
    adj,
    adj2v::Adjacent2Vertex{I,Es,E},
    graph::Graph{I},
    boundary_index_ranges,
    boundary_map,
    TriangleType::Type{V},
    check_existence::C=Val(has_multiple_segments(boundary_map)),
    rng::AbstractRNG=Random.default_rng()) where {I,V,C,Es,E}
    e = sort_edge_by_degree(e, graph) # faster to start at the minimum degree vertex of the edge
    history = PointLocationHistory{V,E,I}()
    add_left_vertex!(history, initial(e))
    add_right_vertex!(history, initial(e))
    jump_and_march(
        pts,
        adj,
        adj2v,
        graph,
        boundary_index_ranges,
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
    left_vertices = history.left_vertices
    right_vertices = history.right_vertices
    reverse!(left_vertices) # counter-clockwise
    if !isempty(collinear_segments)
        connect_segments!(collinear_segments)
        extend_segments!(collinear_segments, e)
    end
    return intersecting_triangles, collinear_segments, left_vertices, right_vertices
end