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
        TrianglesType::Type{Vs},
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
    TrianglesType::Type{Vs},
    check_existence::C=Val(has_multiple_segments(boundary_map)),
    rng::AbstractRNG=Random.default_rng()) where {I,V,Vs,C,Es,E}
    e = sort_edge_by_degree(e, graph) # faster to start at the minimum degree vertex of the edge
    intersecting_triangles = initialise_triangles(Vs)
    collinear_segments = Vector{E}() # use vector so that we get the ordering of the segments
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
        store_visited_triangles=Val(true),
        visited_triangles=intersecting_triangles,
        collinear_segments,
        rng
    )
    if !isempty(collinear_segments)
        connect_segments!(collinear_segments)
        extend_segments!(collinear_segments, e)
    end
    return intersecting_triangles, collinear_segments
end