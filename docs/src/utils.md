```@meta
CurrentModule = DelaunayTriangulation
```

# Other Utilities 

Here are some docstrings for other utility functions.

```@docs 
is_true
get_boundary_index
rotate_ghost_triangle_to_standard_form 
get_right_boundary_node(::Adjacent{I, E}, ::Any, ::Any, ::Any, ::C) where {I,E,C}
get_left_boundary_node(::Adjacent{I, E}, ::Any, ::Any, ::Any, ::C) where {I,E,C}
find_edge(::Any, ::Any, ::Any)
choose_uvw 
is_circular 
circular_equality 
get_surrounding_polygon(::Adjacent{I,E}, ::Graph, ::Any, ::Any, ::C) where {I,E,C}
sort_edge_by_degree(::E, ::Graph) where {E}
split_constrained_edge!(::Any, ::E, ::Any) where {E}
fix_segments!(::AbstractVector{E}, bad_indices) where {E}
connect_segments!(::AbstractVector{E}) where {E}
extend_segments!(::AbstractVector{E}, ::Any) where {E}
convert_boundary_points_to_indices(::AAA, ::AAA) where {F<:Number,A<:AbstractVector{F},AA<:AbstractVector{A},AAA<:AbstractVector{AA}}
get_ordinal_suffix 
check_args
min_max 
min_med_max 
balanced_power_of_two_ternary_split
balanced_power_of_two_quarternary_split 
segment_vertices_adjoin_other_segments
edge_lies_on_two_distinct_segments
get_shared_vertex
convert_to_boundary_index 
get_neighbouring_boundary_edges 
is_last_boundary_index 
is_first_boundary_index 
previndex_circular
nextindex_circular
```