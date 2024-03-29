```@meta
CurrentModule = DelaunayTriangulation
```

# Predicates 

The predicates that we use in this package are all built from ExactPredicates.jl, avoiding degeneracies from predicates owing to floating point arithmetic. The results from predicates are based on certificates, coming from a `Certificate` type defined with EnumX.jl. The definition of this is below.

```@docs 
Certificate
```

## General 

Below we list some general predicates. The core ones that all other predicates are based on are:

```@docs 
orient_predicate
incircle_predicate 
parallelorder_predicate
```

In code, these predicates could be defined by (the actual definition with ExactPredicates.jl is much more involved):

```julia
_det(a, b, c, d) = a * d - b * c
_det(a, b, c, d, e, f, g, h, i) = a * _det(e, f, h, i) - d * _det(b, c, h, i) + g * _det(b, c, e, f) # cofactor expansion 
function orient_predicate(a, b, c)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    o = _det(ax - cx, ay - cy, bx - cx, by - cy)
    return Int(sign(o)) # need Int for xor
end
function incircle_predicate(a, b, c, d)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    dx, dy = getxy(d)
    o = _det(ax - dx, ay - dy, (ax - dx)^2 + (ay - dy)^2,
        bx - dx, by - dy, (bx - dx)^2 + (by - dy)^2,
        cx - dx, cy - dy, (cx - dx)^2 + (cy - dy)^2)
    return Int(sign(o)) # need Int for xor
end
function parallelorder_predicate(a, b, p, q)
    return orient_predicate(b .- a, q .- p, (0.0, 0.0))
end
```

You could use this as a reference if you want to disconnect from using ExactPredicates.jl (or e.g. use the predicates also defined in GeometricalPredicates.jl). This could be useful if you are not too worried about robustness (although you should typically care about this, so be careful - proceed at your own peril) and just want fast code.

The other predicates are:

```@docs
sameside_predicate
opposite_signs
meet_predicate 
triangle_orientation(::Any, ::Any, ::Any)
point_position_relative_to_circle(::Any, ::Any, ::Any, ::Any)
point_position_relative_to_line(::Any, ::Any, ::Any)
point_closest_to_line(::Any, ::Any, ::Any, ::Any)
point_position_on_line_segment(::Any, ::Any, ::Any)
line_segment_intersection_type(::Any, ::Any, ::Any, ::Any)
point_position_relative_to_triangle(::Any, ::Any, ::Any, ::Any)
point_position_relative_to_oriented_outer_halfplane(::Any, ::Any, ::Any)
is_legal(::Any, ::Any, ::Any, ::Any)
triangle_line_segment_intersection(::Any, ::Any, ::Any, ::Any, ::Any)
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
is_outer_ghost_triangle(::Any, ::Any, ::Any, ::Any)
is_outer_ghost_edge(::Any, ::Any, ::Any)
is_outer_boundary_node(::Any, ::Graph{I}, ::Any) where {I}
is_boundary_node(::Any, ::Graph{I}, ::Any) where {I}
edge_exists(::I) where {I}
edge_exists(::Any, ::Adjacent{I, E}) where {I, E}
has_ghost_triangles(::Adjacent{I, E}, ::Any) where {I, E}
```

## Index and Ghost Handling

Below we list methods for working with predicates that are used when we provide indices for points rather than points directly.

```@docs 
triangle_orientation(::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_circumcircle(::Any, ::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_line(::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
point_closest_to_line(::Any, ::Any, ::Any, ::Any, ::Any)
point_position_on_line_segment(::Any, ::Any, ::Any, ::Any)
line_segment_intersection_type(::Any, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_triangle(::Any, ::Any, ::Any, ::Any, ::Any, ::Any, ::AbstractDict)
triangle_line_segment_intersection(::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
```