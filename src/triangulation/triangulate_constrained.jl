"""
    triangulate_constrained(tri::Triangulation; rng=Random.default_rng())

Inserts all constrained edges and boundary edges into `tri`.
"""
function triangulate_constrained!(tri::Triangulation, bn_map=get_boundary_map(tri); rng=Random.default_rng())
    all_constrained_edges = merge_constrained_edges(tri, bn_map)
    for e in each_edge(all_constrained_edges)
        add_edge!(tri, e; rng)
    end
    return nothing
end