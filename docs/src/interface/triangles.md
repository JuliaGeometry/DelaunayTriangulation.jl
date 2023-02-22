```@meta
CurrentModule = DelaunayTriangulation
```

## Individual Triangles 

Triangles are assumed to be of the form `(i, j, k)`, with positive orientation, but we allow for customisation in the way these indices are represented. The following methods are used for working with triangles. First, we list the methods that must be defined, and then methods that extend these former methods. 

### Necessary Methods

```@docs 
construct_triangle 
geti 
getj 
getk 
integer_type 
```

### Generic Methods 

```@docs 
indices 
triangle_edges
rotate_triangle 
construct_positively_oriented_triangle(::Type{V}, ::Any, ::Any, ::Any, ::Any) where {V}
compare_triangles 
```

## Collection of Triangles 

A collection of triangles simply stores many triangles, and the collection itself must be mutable so that triangles can be added and deleted. The following methods are used for working with collections of triangles. First, we list the methods that must be defined, and then methods that extend these former methods.

### Necessary Methods 

```@docs 
initialise_triangles 
triangle_type 
num_triangles 
add_to_triangles! 
delete_from_triangles!
each_triangle 
```

You must also provide definitions for `Base.in` and `Base.sizehint!` for your type.
Note also that `Triangulation`s also define `each_solid_triangle` and `each_ghost_triangle`.

### Generic Methods 

```@docs 
contains_triangle 
add_triangle!(::Any, ::Vararg{F, N}) where {F, N}
add_triangle!(::Ts, ::Integer, ::Integer, ::Integer) where {Ts}
delete_triangle!(::Any, ::Vararg{F, N}) where {F, N}
delete_triangle!(::Ts, ::Integer, ::Integer, ::Integer) where {Ts}
compare_triangle_collections 
```