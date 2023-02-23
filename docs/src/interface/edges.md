```@meta
CurrentModule = DelaunayTriangulation
```

## Individual Edges 

Edges are assumed to take the form `(i, j)`, with customisation available for how we represent `(i, j)`. The following methods are used, where we first list the methods that must be defined and then methods that extend these former methods. 

### Necessary Methods 

```@docs 
construct_edge 
initial
terminal
```

### Generic Methods 

```@docs 
edge_indices 
```

## Collection of Edges 

A collection of edges simply stores many edges, and this collection must be mutable so that edges can be deleted added and deleted. The following methods are used, where we first list the methods that must be defined and then methods that extend these former methods. 

### Necessary Methods 

```@docs 
initialise_edges 
edge_type 
num_edges 
contains_edge 
add_to_edges! 
delete_from_edges! 
each_edge 
```

### Generic Methods 

```@docs 
add_edge!
delete_edge!
random_edge 
is_empty 
```