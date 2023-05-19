"""
    locate_intersecting_triangles(tri::Triangulation, e, rotate=Val(true); rng::AbstractRNG=Random.default_rng()) where {C}

Returns a list of triangles intersecting the segment `e` in `tri`. If `is_true(rotate)`, 
then `e` will be sorted such that `initial(e)` has smaller degree in `tri` than `terminal(e)`.

More precisely, the returned values are:
- `intersecting_triangles`: The triangles intersecting `e`.
- `collinear_segments`: Any segments collinear with `e`, giving in order of appearance.
- `left_vertices`: All vertices of `intersecting_triangles` appearing to the left of `e`.
- `right_vertices`: All vertices of `intersecting_triangles` appearing to the right of `e`.
"""
function locate_intersecting_triangles(
    tri::Triangulation{P,Ts,I,E,Es},
    e,
    rotate=Val(true),
    rng::AbstractRNG=Random.default_rng()) where {P,Ts,I,E,Es}
    V = triangle_type(tri)
    e = is_true(rotate) ? sort_edge_by_degree(tri, e) : e # faster to start at the minimum degree vertex of the edge
    history = PointLocationHistory{V,E,I}()
    add_left_vertex!(history, initial(e))
    add_right_vertex!(history, initial(e))
    jump_and_march(
        tri,
        get_point(tri, terminal(e));
        m=nothing,
        k=initial(e),
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