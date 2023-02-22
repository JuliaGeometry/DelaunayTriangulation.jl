```@meta
CurrentModule = DelaunayTriangulation
```

# Predicates 

The predicates that we use in this package are all built from ExactPredicates.jl, avoiding degeneracies from predicates owing to floating point arithmetic. The results from predicates are based on certificates, coming from a `Certificate` type defined with EnumX.jl. The definition of this is below.

```@docs 
Certificate
```

## General 

Below we list some general predicates.

```@docs 
orient_predicate
incircle_predicate 
sameside_predicate 
meet_predicate
triangle_orientation(::Any, ::Any, ::Any)
point_position_relative_to_circle 
point_position_relative_to_line(::Any, ::Any, ::Any) 
point_position_on_line_segment(::Any, ::Any, ::Any) 
line_segment_intersection_type(::Any, ::Any, ::Any, ::Any) 
point_position_relative_to_triangle(::Any, ::Any, ::Any, ::Any) 
point_position_relative_to_oriented_outer_halfplane
```

## Boundaries and Ghosts 

Below we list some predicates for working with boundaries and ghost triangles. 

```@docs 
is_boundary_index 
is_boundary_edge(::Any, ::Adjacent) 
is_boundary_triangle(::Any, ::Any, ::Any, ::Any) 
is_ghost_edge 
is_ghost_triangle 
is_interior_curve
is_outer_boundary_index(::Any, ::Any) 
is_outer_ghost_triangle 
is_outer_ghost_edge
is_outer_boundary_node(::Any, ::Graph{I}, ::Any) where {I} 
edge_exists(::I) where {I}
edge_exists(::Any, ::Adjacent{I,E}) where {I,E}
has_ghost_triangles(::Adjacent{I,E}, ::Any) where {I,E} 
```

## Index and Ghost Handling

Below we list methods for working with predicates that are used when we provide indices for points rather than points directly.

```@docs 
triangle_orientation(::Any, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_circumcircle(::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_line(::Any, ::Any, ::Any, ::Any, ::Any)
point_position_on_line_segment(::Any, ::Any, ::Any, ::Any)
line_segment_intersection_type(::Any, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_triangle(::Any, ::Any, ::Any, ::Any, ::Any, ::AbstractDict)
point_position_relative_to_triangle(::Any, ::Any, ::Any, ::AbstractDict)
```