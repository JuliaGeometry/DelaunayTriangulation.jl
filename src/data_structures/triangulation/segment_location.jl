"""
    locate_intersecting_triangles(tri::Triangulation, e, rotate=Val(true);
        check_existence::C=Val(has_multiple_segments(tri)),
        rng::AbstractRNG=Random.default_rng()) where {C}

Returns a list of triangles intersecting the segment `e` in `tri`. If `is_true(rotate)`, 
then `e` will be sorted such that `initial(e)` has the least degree in `tri`.

More precisely, the returned values are:
- `intersecting_triangles`: The triangles intersecting `e`.
- `collinear_segments`: Any segments collinear with `e`, giving in order of appearance.
- `left_vertices`: All vertices of `intersecting_triangles` appearing to the left of `e`.
- `right_vertices`: All vertices of `intersecting_triangles` appearing to the right of `e`.
"""
function locate_intersecting_triangles(tri::Triangulation, e, rotate=Val(true);
    check_existence::C=Val(has_multiple_segments(tri)),
    rng::AbstractRNG=Random.default_rng()) where {C}
    pts = get_points(tri)
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    boundary_index_ranges = get_boundary_index_ranges(tri)
    boundary_map = get_boundary_map(tri)
    T = triangle_type(tri)
    rep = get_representative_point_list(tri)
    return locate_intersecting_triangles(e,
        pts,
        adj,
        adj2v,
        graph,
        boundary_index_ranges,
        rep,
        boundary_map,
        T,
        check_existence,
        rotate,
        rng)
end