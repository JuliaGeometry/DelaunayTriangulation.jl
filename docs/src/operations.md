```@meta
CurrentModule = DelaunayTriangulation
```

# Operations

We define some specific operations for acting on `Triangulation`s directly. These are listed below.

```@docs 
add_point!(::Triangulation, ::Any)
add_edge!(::Triangulation, ::Any)
delete_point!(::Triangulation, ::Any)
add_boundary_information!(::Triangulation)
delete_ghost_triangles!(::Triangulation)
add_ghost_triangles!(::Triangulation)
add_triangle!(::Ts, ::Integer, ::Integer, ::Integer) where {Ts<:Triangulation}
delete_triangle!(::Ts, ::Integer, ::Integer, ::Integer) where {Ts<:Triangulation}
split_edge!
legalise_split_edge!
complete_split_and_legalise!
split_triangle!
legalise_split_triangle!
flip_edge!
legalise_edge!
clear_empty_features!
```