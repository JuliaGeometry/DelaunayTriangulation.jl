```@meta
CurrentModule = DelaunayTriangulation
```

## Boundary Nodes

As mentioned at the start of this section, the interface for representing boundary nodes allows for support for a contiguous boundary, a segmented boundary, and multiple separate boundaries. This interface is customisable, and we define the following methods for this.

### Necessary Methods 

```@docs 
has_multiple_curves 
has_multiple_segments 
num_curves 
num_segments 
num_boundary_edges 
getboundarynodes 
each_boundary_node 
```

### Generic Methods 

```@docs 
get_boundary_nodes(::Any, ::Vararg{Any})
construct_boundary_map 
construct_boundary_index_ranges 
map_boundary_index(::Any, ::Any) 
get_curve_index 
get_segment_index 
num_outer_boundary_segments 
```
