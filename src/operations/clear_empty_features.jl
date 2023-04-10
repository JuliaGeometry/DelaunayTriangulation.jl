"""
    clear_empty_features!(tri::Triangulation)

Clears all empty features from the triangulation `tri`.
"""
function clear_empty_features!(tri)
    adj = get_adjacent(tri)
    adj2v = get_adjacent2vertex(tri)
    graph = get_graph(tri)
    clear_empty_keys!(adj)
    clear_empty_keys!(adj2v)
    clear_empty_points!(graph)
    return nothing
end



