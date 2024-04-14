"""
    each_segment(tri::Triangulation) -> Edges

Returns an iterator over all segments in `tri`. This includes both interior and boundary segments. If you only want 
interior segments, then see [`get_interior_segments`](@ref),
"""
each_segment(tri::Triangulation) = (each_edge ∘ get_all_segments)(tri)

"""
    contains_segment(tri::Triangulation, ij) -> Bool 
    contains_segment(tri::Triangulation, i, j) -> Bool

Returns `true` if `(i, j)` is a segment in `tri`, and `false` otherwise. Both `(i, j)` and `(j, i)` are checked.
"""
function contains_segment(tri::Triangulation, e)
    segments = get_all_segments(tri)
    return contains_unoriented_edge(e, segments)
end
function contains_segment(tri::Triangulation, i, j)
    E = edge_type(tri)
    e = construct_edge(E, i, j)
    return contains_segment(tri, e)
end