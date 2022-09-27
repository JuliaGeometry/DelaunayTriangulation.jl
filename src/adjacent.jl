"""
    Adjacent{A}

This is a data structure for the Delaunay triangulation. It stores a `Dict`, `adjacent`, so that 
`(u, v, adjacent[(u, v)])` is a positively oriented triangle. This struct is also callable, e.g. if 
`𝒜::Adjacent`, then `(u, v, 𝒜(u, v))` is a positively oriented triangle.
"""
struct Adjacent{A}
    adjacent::DirectedGraph{EdgeType}
    function Adjacent()
        A = DirectedGraph{EdgeType}()
        adj = new{typeof(A)}(A)
        return adj
    end
end
