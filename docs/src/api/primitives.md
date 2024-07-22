```@meta 
CurrentModule = DelaunayTriangulation
```

# Primitive Interfaces 

## Points 

```@docs 
getx
gety
getxy
is_point2 
getpoint
get_point
push_point!
pop_point!
set_point!
```

## Edges 

```@docs 
construct_edge
initial
terminal
edge_vertices
reverse_edge
compare_unoriented_edges
edge_type
contains_edge
contains_unoriented_edge
add_to_edges!
delete_from_edges!
delete_unoriented_edge!
random_edge
edges_are_disjoint
```

## Triangles 

```@docs 
construct_triangle
geti
getj
getk
triangle_vertices
triangle_type
triangle_edges
compare_triangles
add_to_triangles!
delete_from_triangles!
```

## Boundary Nodes 

```@docs 
num_boundary_edges
each_boundary_node
```