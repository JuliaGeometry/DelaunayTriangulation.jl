```@meta
CurrentModule = DelaunayTriangulation
```

# Triangulation 

The most important structure of the package is the Triangulation data structure. Its complete definition is given below, and then afterwards we give all the docstrings for methods that are useful for working with triangulations.

```@docs 
Triangulation
```

Each field has its own accessor: 

```@docs 
get_points(::Triangulation) 
get_triangles(::Triangulation)
get_adjacent(::Triangulation)
get_adjacent2vertex(::Triangulation)
get_graph(::Triangulation)
get_boundary_nodes(::Triangulation)
get_boundary_edge_map(::Triangulation)
get_boundary_map(::Triangulation)
get_boundary_index_ranges(::Triangulation)
get_constrained_edges(::Triangulation)
get_all_constrained_edges(::Triangulation)
get_convex_hull(::Triangulation)
get_representative_point_list(::Triangulation)
```

There are several useful methods available for working with triangulations. We split these into several sections. 

## Adjacent Methods 

```@docs 
get_adjacent(::Triangulation, ::Any)
add_adjacent!(::Triangulation, ::Any, ::Any)
delete_adjacent!(::Triangulation, ::Any)
```

## Adjacent2Vertex Methods

```@docs 
get_adjacent2vertex(::Triangulation, ::Any)
add_adjacent2vertex!(::Triangulation, ::Any, ::Any)
delete_adjacent2vertex!(::Triangulation, ::Any, ::Any)
delete_adjacent2vertex!(::Triangulation, ::Any)
```

## Boundary Nodes Methods

```@docs 
has_multiple_curves(::Triangulation)
has_multiple_segments(::Triangulation)
num_curves(::Triangulation)
num_segments(::Triangulation)
get_boundary_nodes(::Triangulation, ::Any)
map_boundary_index(::Triangulation, ::Any)
get_curve_index(::Triangulation, ::Any)
get_segment_index(::Triangulation, ::Any)
num_outer_boundary_segments(::Triangulation)
get_right_boundary_node(::Triangulation, ::Any, ::Any)
get_left_boundary_node(::Triangulation, ::Any, ::Any)
get_boundary_index_range(::Triangulation, ::Any)
get_boundary_edge_map(::Triangulation, ::Any)
insert_boundary_node!(::Triangulation, ::Any, ::Any)
split_boundary_edge!(::Triangulation, ::Any, ::Any)
split_boundary_edge_at_collinear_segments!(::Triangulation, ::Any)
contains_boundary_edge(::Triangulation, ::Any)
merge_constrained_edges(::Any, ::Any, ::Es) where {Es}
get_all_boundary_nodes(::Triangulation)
all_boundary_indices(::Triangulation)
delete_boundary_node!(::Triangulation, ::Any)
merge_boundary_node!(::Triangulation, ::Any, ::Any)
```

## Convex Hull Methods

```@docs 
get_convex_hull_indices(::Triangulation)
convex_hull!(::Triangulation)
```

## Edges Methods

```@docs 
edge_type(::Triangulation{P,Ts,I,E}) where {P,Ts,I,E}
num_edges(::Triangulation)
each_edge(::Triangulation)
each_constrained_edge(::Triangulation)
contains_constrained_edge(::Triangulation, ::Any)
num_ghost_edges(::Triangulation)
num_solid_edges(::Triangulation)
each_solid_edge
each_ghost_edge 
sort_edge_by_degree(::Triangulation, ::Any)
split_constrained_edge!(::Triangulation, ::Any, ::Any)
```

## Graph Methods

```@docs 
get_edges(::Triangulation)
get_vertices(::Triangulation)
get_neighbours(::Triangulation)
get_neighbours(::Triangulation, ::Any)
add_vertex!(::Triangulation, ::Vararg{Any})
num_neighbours(::Triangulation, ::Any)
add_neighbour!(::Triangulation, ::Any, ::Vararg{Any})
delete_neighbour!(::Triangulation, ::Any, ::Vararg{Any})
delete_vertex!(::Triangulation, ::Vararg{Any})
delete_boundary_vertices_from_graph!(::Triangulation)
each_vertex(::Triangulation)
num_vertices(::Triangulation)
```

## Points Methods

```@docs 
get_point(::Triangulation, ::Any)
each_point_index(::Triangulation)
each_point(::Triangulation)
num_points(::Triangulation)
push_point!(::Triangulation, ::Any, ::Any)
each_solid_vertex 
each_ghost_vertex
num_solid_vertices 
num_ghost_vertices
```

## Triangles Methods

```@docs 
triangle_type(::Triangulation{P,Ts}) where {P,Ts}
num_triangles(::Triangulation)
each_triangle(::Triangulation)
contains_triangle(::Triangulation, ::Any)
construct_positively_oriented_triangle(::Triangulation, ::Any, ::Any, ::Any)
num_ghost_triangles(::Triangulation)
num_solid_triangles(::Triangulation)
each_solid_triangle 
each_ghost_triangle
```

## Predicates Methods

```@docs 
is_boundary_edge(::Triangulation, ::Any)
is_boundary_triangle(::Triangulation, ::Any, ::Any, ::Any)
triangle_orientation(::Triangulation, ::Any, ::Any, ::Any)
point_position_relative_to_circumcircle(::Triangulation, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_line(::Triangulation, ::Any, ::Any, ::Any)
point_closest_to_line(::Triangulation, ::Any, ::Any, ::Any, ::Any)
point_position_on_line_segment(::Triangulation, ::Any, ::Any, ::Any)
line_segment_intersection_type(::Triangulation, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_triangle(::Triangulation, ::Any, ::Any, ::Any, ::Any)
triangle_line_segment_intersection(::Triangulation, ::Any, ::Any, ::Any, ::Any, ::Any)
is_outer_boundary_index(::Triangulation, ::Any)
is_outer_boundary_node(::Triangulation, ::Any)
is_boundary_node(::Triangulation, ::Any)
edge_exists(::Triangulation, ::Any, ::Any)
has_ghost_triangles(::Triangulation)
has_boundary_nodes(::Triangulation)
is_legal(::Triangulation, ::Any, ::Any)
find_edge(::Triangulation, ::Any, ::Any)
is_constrained(::Triangulation)
```

## Representative Points Methods

```@docs 
compute_representative_points!(::Triangulation)
get_representative_point_coordinates(::Triangulation, ::Any)
reset_representative_points!(::Triangulation)
update_centroid_after_addition!(::Triangulation, ::Any, ::Any)
update_centroid_after_deletion!(::Triangulation, ::Any, ::Any)
new_representative_point!(::Triangulation, ::Any)
```