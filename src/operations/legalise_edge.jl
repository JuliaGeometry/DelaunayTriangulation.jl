"""
    legalise_edge!(tri::Triangulation, i, j, r)

Given a triangulation `tri`, an edge `(i, j)`, and a 
point `r` that was added into a triangle that `(i, j)` 
belongs to, legalises the edge `(i, j)` and other neighbouring 
edges recursively.
"""
function legalise_edge!(tri::Triangulation, i, j, r)
    cert = is_legal(tri, i, j)
    if is_illegal(cert)
        e = get_adjacent(tri, j, i)
        flip_edge!(tri, i, j, e, r)
        legalise_edge!(tri, i, e, r)
        legalise_edge!(tri, e, j, r)
    end
    return nothing
end