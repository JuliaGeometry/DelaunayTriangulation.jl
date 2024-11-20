"""
    clear_empty_features!(tri::Triangulation)

Clears all empty features from the triangulation `tri`.
"""
function clear_empty_features!(tri::Triangulation)
    adj2v = get_adjacent2vertex(tri)
    clear_empty_keys!(adj2v)
    return tri
end
