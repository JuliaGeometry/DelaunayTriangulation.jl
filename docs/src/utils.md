```@meta
CurrentModule = DelaunayTriangulation
```

# Other Utilities 

Here are some docstrings for other utility functions.

```@docs 
get_right_boundary_node(::Adjacent{I, E}, ::Any, ::Any, ::Any, ::C) where {I,E,C}
get_left_boundary_node(::Adjacent{I, E}, ::Any, ::Any, ::Any, ::C) where {I,E,C}
find_edge(::Any, ::Any, ::Any)
choose_uvw(::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
is_circular 
circular_equality
get_surrounding_polygon(::Adjacent{I, E}, ::Graph, ::Any, ::Any, ::C) where {I,E,C}
sort_edge_by_degree(e::E, graph::Graph) where {E}
split_constrained_edge!(::Any, ::Any, ::Any)
```