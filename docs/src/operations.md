```@meta
CurrentModule = DelaunayTriangulation
```

# Operations

We define some specific operations for acting on `Triangulation`s directly. These are listed below.

```@docs 
add_boundary_information!(::Triangulation)
add_ghost_triangles!(::Triangulation)
add_point!(::Triangulation, ::Any)
add_triangle!(::Ts, ::Integer, ::Integer, ::Integer) where {Ts<:Triangulation}
delete_ghost_triangles!(::Triangulation)
delete_triangle!(::Ts, ::Integer, ::Integer, ::Integer) where {Ts<:Triangulation}
split_edge!
split_triangle!
flip_edge!
legalise_edge!
delete_point!
```