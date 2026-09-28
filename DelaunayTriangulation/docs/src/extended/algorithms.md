```@meta 
CurrentModule = DelaunayTriangulation
```

# Algorithm Internals 

Here we give some functions that are used for some of the algorithms in this package.

## Unconstrained Triangulations

Here are some of the internal functions used for the Bowyer-Watson algorithm, i.e. for the computation of an unconstrained triangulation.

```@docs 
postprocess_triangulate!
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/triangulation/unconstrained_triangulation.jl"]
```

## Triangulation Rectangular Domains

Here are some of the internal functions used for `triangulate_rectangle`.

```@docs 
get_lattice_triangles 
get_lattice_points 
get_lattice_boundary
```

## Triangulating Convex Polygons

Here are some functions used in `triangulate_convex` for triangulating a convex polygon.

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/triangulation/triangulate_convex.jl"]
```

## Constrained Triangulations

There are many functions to list for the computation of a constrained triangulation.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/triangulation/basic_operations/add_segment.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation] 
Pages = ["src/algorithms/triangulation/constrained_triangulation.jl"]
Filter = t -> !(t in (DelaunayTriangulation.convert_boundary_points_to_indices,))
```

```@autodocs
Modules = [DelaunayTriangulation] 
Pages = ["src/data_structures/triangulation/methods/boundary_edge_map.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation] 
Pages = ["src/data_structures/triangulation/methods/boundary_nodes.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation] 
Pages = ["src/data_structures/triangulation/methods/exterior_curve_indices.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation] 
Pages = ["src/data_structures/triangulation/methods/ghost_vertex_map.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation] 
Pages = ["src/data_structures/triangulation/methods/ghost_vertex_ranges.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation] 
Pages = ["src/data_structures/triangulation/methods/representative_point_list.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation] 
Pages = ["src/data_structures/triangulation/methods/segments.jl"]
```

## Weighted Triangulations

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/data_structures/triangulation/methods/weights.jl"]
Filter = t -> !(t in (DelaunayTriangulation.ZeroWeight,))
```

## Mesh Refinement 

Here are some functions involved with mesh refinement.

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/triangulation/mesh_refinement.jl"]
```

```@docs 
is_encroachmentfailure 
is_successfulinsertion
is_failedinsertion
is_precisionfailure
```

## RTree

The RTrees we work with during boundary enrichment have several functions associated with them.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/intersections/rtree.jl"]
```

## Point Location

Here are some functions related to point location.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/point_location/brute_force.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/point_location/jump_and_march.jl"]
```

## Triangulation Operations

Some of the triangulation operations have internal functions associated with them.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/triangulation/basic_operations/delete_holes.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/triangulation/basic_operations/delete_point.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/triangulation/basic_operations/split_edge.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/triangulation/basic_operations/split_triangle.jl"]
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/data_structures/triangulation/methods/adjacent.jl"]
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/data_structures/triangulation/methods/adjacent2vertex.jl"]
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/data_structures/triangulation/methods/graph.jl"]
Filter = t -> !(t in (DelaunayTriangulation.get_vertices, DelaunayTriangulation.get_edges))
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/data_structures/triangulation/methods/points.jl"]
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/data_structures/triangulation/methods/triangles.jl"]
```

## Voronoi Tessellations

Here are some functions related to the computation of unbounded Voronoi tessellations.

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/voronoi/unbounded.jl"]
```

## Clipped Voronoi Tessellations

Here are some functions related to the computation of clipped Voronoi tessellations.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/voronoi/clipped.jl", "src/algorithms/voronoi/clipped_coordinates.jl", "src/algorithms/polygon_clipping/liang_barsky.jl", "src/algorithms/polygon_clipping/sutherland_hodgman.jl"]
```

## Centroidal Voronoi Tessellations 

Here are some functions related to the computation of centroidal Voronoi tessellations.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/voronoi/centroidal.jl"]
```

## Triangulating Curve-Bounded Domains

We have several functions related to the triangulation of curve-bounded domains.

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/algorithms/triangulation/triangulate_curve_bounded.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/data_structures/triangulation/methods/boundary_curves.jl"]
```

```@docs 
is_visible
is_invisible
test_visibility
```