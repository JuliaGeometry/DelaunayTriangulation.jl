"""
    get_ghost_vertex_range(tri::Triangulation, ℓ) -> UnitRange 

Given a ghost vertex `ℓ` of `tri`, returns the range of all 
ghost vertices corresponding to the same curve or section as `ℓ` does.
"""
get_ghost_vertex_range(tri::Triangulation, ℓ) = get_ghost_vertex_ranges(tri)[ℓ]

"""
    all_ghost_vertices(tri::Triangulation) -> KeySet 

Returns the set of all ghost vertices in `tri`.
"""
all_ghost_vertices(tri::Triangulation) = keys(get_ghost_vertex_ranges(tri))
