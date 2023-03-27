"""
    sort_edge_by_degree(tri::Triangulation, e)

Given an edge `e` of a triangulation `tri`, say `e = (u, v)`,
returns:

- If `deg(u) < deg(v)`, returns `e`;
- If `deg(u) ≥ deg(v)`, returns `(v, u)`.

In particular, `e` is sorted so that `initial(e)` is the vertex of `e` 
with the smallest degree.
"""
function sort_edge_by_degree(tri::Triangulation, e)
    E = edge_type(tri)
    u = initial(e)
    v = terminal(v)
    d₁ = num_neighbours(tri, u)
    d₂ = num_neighbours(tri, v)
    if d₁ < d₂
        return e
    else
        return construct_edge(E, v, u)
    end
end