```@meta
CurrentModule = DelaunayTriangulation
```

# Voronoi Tessellation 

The data structure for a Voronoi tessellation is reasonably simple. The data structure here is probably not as exhaustive as it could be, but it is sufficient.

```@docs 
VoronoiTessellation 
```

Each field has its own accessor:

```@docs 
get_triangulation(::VoronoiTessellation) 
get_generators(::VoronoiTessellation) 
get_polygon_points(::VoronoiTessellation)
get_polygons(::VoronoiTessellation) 
get_circumcenter_to_triangle(::VoronoiTessellation)
get_triangle_to_circumcenter(::VoronoiTessellation) 
get_unbounded_polygons(::VoronoiTessellation) 
get_cocircular_circumcenters(::VoronoiTessellation) 
get_adjacent(::VoronoiTessellation) 
get_boundary_polygons(::VoronoiTessellation) 
```

There are several useful methods available for working with this data structure. We list some of these below; for functions that actually construct the tessellation, see the dedicated tessellation section in the sidebar.

## Type queries 

```@docs 
edge_type(::VoronoiTessellation{Tr,P,I,T,S,E}) where {Tr,P,I,T,S,E}
number_type(::VoronoiTessellation{Tr,P}) where {Tr,P}
integer_type(::VoronoiTessellation{Tr,P,I}) where {Tr,P,I}
triangle_type(::VoronoiTessellation{Tr,P,I,T})
```

## Getters 

```@docs 
get_generator(::VoronoiTessellation, ::Any)
get_polygon_point(::VoronoiTessellation, ::Any)
get_polygon(::VoronoiTessellation, ::Any)
get_circumcenter_to_triangle(::VoronoiTessellation, ::Any)
get_triangle_to_circumcenter(::VoronoiTessellation, ::Any)
get_polygon_coordinates 
get_neighbouring_boundary_edges 
```

## Nums 

```@docs 
num_polygons 
num_polygon_vertices 
num_generators 
``` 

## Adders 

```@docs
add_polygon! 
push_polygon_point! 
add_unbounded_polygon! 
delete_unbounded_polygon! 
add_boundary_polygon! 
```

## Iterators 

```@docs 
each_generator 
each_polygon_vertex 
each_unbounded_polygon 
each_polygon 
each_polygon_index 
```

## Adjacent 

```@docs 
get_adjacent(::VoronoiTessellation, ::Any)
add_adjacent!(::VoronoiTessellation, ::Any, ::Any)
delete_adjacent!(::VoronoiTessellation, ::Any, ::Any)
delete_polygon_adjacent!
add_polygon_adjacent!
```

## Features 

```@docs 
polygon_features(::VoronoiTessellation, ::Any)
get_area
get_centroid 
polygon_bounds(::VoronoiTessellation, ::Any)
jump_and_march(::VoronoiTessellation, ::Any)
```

## Utilities 

```@docs 
get_surrounding_polygon(::VoronoiTessellation, ::Any)
convert_to_boundary_edge(::VoronoiTessellation, :Any)
```