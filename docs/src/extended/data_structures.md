```@meta 
CurrentModule = DelaunayTriangulation
```

# Data Structures

Here we list the internal data structures used throughout the package, along with functions related to working with the respective data structures.

## MaxPriorityQueue 

The `MaxPriorityQueue` is our implementation of a maximum priority queue; see, for example, [the wikipedia page](https://en.wikipedia.org/wiki/Priority_queue). The `MaxPriorityQueue` works similarly to the priority queue in [DataStructures.jl](https://github.com/JuliaCollections/DataStructures.jl). This struct is used for determining which triangles and edges to prioritise during mesh refinement, as well as which cells to prioritise during the quadtree construction for the pole of inaccessibility. 

```@docs 
MaxPriorityQueue
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["data_structures/queue/max_priority_queue.jl"]
```

## Queue 

The `Queue` is our implementation of a basic FIFO queue; see, for example, [the wikipedia page](https://en.wikipedia.org/wiki/Queue_(abstract_data_type)). The implementation of a `Queue` in this package is very simple - it is just a vector. A block-based approach could also be implemented in the future, it just hasn't. The idea here is to avoid a dependency on DataStructures.jl. The `Queue` is used in a variety of places, such as during boundary enrichment in [`enrich_boundary!`](@ref) to determine which edges to process next.

```@docs 
Queue 
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/queue/queue.jl"]
```

## BalancedBST

The `BalancedBST` is our implementation of a balanced binary search tree; see, for example, [these notes](https://courses.csail.mit.edu/6.006/fall09/lecture_notes/lecture04.pdf). The `BalancedBST` is not currently used anywhere. Originally, it was intended for use with the Bentley-Ottman algorithm, but the complete algorithm is yet to be implemented; this algorithm was to be implemented to allow for users to provide intersecting segments that, through the algorithm, could be automatically split for the user. The `BalancedBST` implementation is fully functional for this purpose, the work in getting the actual algorithm up and going just has to be done.

```@docs 
BalancedBST
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/trees/bst.jl"]
```

## RTree

The `RTree` is our implementation of a simple `RTree` data structure, following heavily the original code in [SpatialIndexing.jl](https://github.com/alyst/SpatialIndexing.jl). This tree is needed for triangulating curve-bounded domains. In particular, since the enrichment phase of triangulating a curve-bounded domain has no triangulation data structure that can be used for finding intersections between an edge's diametral circle and another vertex, the `RTree` is used to find these intersections.

```@docs 
RTree
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/trees/rtree.jl"]
```

```@docs 
QueryResult
```

## PolygonHierarchy

The `PolygonHierarchy` is a data structure used for representing hierarchical information about a set of polygons, such as the piecewise linear approximation of a curve-bounded domain's boundary. This structure is implemented as a collection of [directed trees](https://en.wikipedia.org/wiki/Polytree), where the disjoint trees each represent the disjoint parts of a domain, and the root of each tree contains the boundary curves of all the boundaries that are inside this tree.

```@docs 
PolygonHierarchy
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["data_structures/trees/polygon_hierarchy.jl"]
```

## Adjacent 

The `Adjacent` data structure is used for mapping edges to vertices (for triangulations) or polygons (for tessellations). The `Adjacent` structure itself is in the public API, as is `get_adjacent`.

```@docs; canonical=false
Adjacent
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["data_structures/triangulation/adjacent.jl"]
```

## Adjacent2Vertex

The `Adjacent2Vertex` data structure is used for mapping vertices to edges for triangulations. The `Adjacent2Vertex` structure itself is in the public API, as is `get_adjacent2vertex`.

```@docs; canonical=false
Adjacent2Vertex
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/triangulation/adjacent2vertex.jl"]
```

## Graph 

The `Graph` data structure is used for storing connectivity information about vertices in a triangulation. The `Graph` structure itself is in the public API, as is `get_graph` and `get_neighbours`.

```@docs; canonical=false
Graph
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/triangulation/graph.jl"]
```

## Curves 

There are many data structures used to define the curves we provide in this package, all subtyping the `AbstractParametricCurve` type. This type, and its subtypes, are all in the public API with the exception of `PiecewiseLinear`.

```@docs; canonical=false 
AbstractParametricCurve
LineSegment
CircularArc
EllipticalArc
BezierCurve
BSpline
CatmullRomSpline
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/curves/abstract.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/curves/beziercurve.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/curves/bspline.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/curves/catmullromspline.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/curves/circulararc.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/curves/ellipticalarc.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/curves/linesegment.jl"]
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/curves/piecewiselinear.jl"]
```


## RepresentativeCoordinates

We use a `RepresentativeCoordinates` struct for storing the representative coordinates used for interpreting ghost vertices.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/representative_coordinates/representative_coordinates.jl"]
```

```@docs 
update_centroid_after_addition!
```

## Cell

We use a `Cell` struct for representing an individual square in a quadtree during the computation of the pole of inaccessibility.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/representative_coordinates/cell.jl"]
```

## CellQueue

We use a `CellQueue` struct for storing the `Cell`s in a quadtree during the computation of the pole of inaccessibility.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/representative_coordinates/cell_queue.jl"]
```

## ConvexHull

We use a `ConvexHull` struct to represent a [convex hull](https://en.wikipedia.org/wiki/Convex_hull). This struct is in the public API.

```@docs; canonical=false
ConvexHull
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/convex_hull.jl"]
```

## Triangulation

We use a `Triangulation` struct to represent a triangulation. This struct is in the public API.

```@docs; canonical=false
Triangulation
```

```@autodocs    
Modules = [DelaunayTriangulation]
Pages = ["data_structures/triangulation/triangulation.jl"]
```

## TriangulationCache 

The `TriangulationCache` is what we store in the `cache` field of a triangulation.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/triangulation/triangulation_cache.jl"]
```

## BoundaryEnricher 

We use a `BoundaryEnricher` struct for storing information using during the enrichment phase of the triangulation of a curve-bounded domain.

```@docs 
BoundaryEnricher 
```

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/boundary_enricher.jl"]
```

## AbstractEach(Vertex/Edge/Triangle) Iterators 

We use a variety of iterators for iterating over the vertices, edges, and triangles of a triangulation.

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["data_structures/triangulation/methods/iterators.jl"]
```

## PointLocationHistory 

We provide a means for storing the history of triangles encountered during point location, using a `PointLocationHistory` struct. The main motivation for this struct is for constrained triangulations.

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["data_structures/point_location_history.jl"]
```

## IndividualTriangleStatistics 

We provide an `IndividualTriangleStatistics` struct for storing statistics about individual triangles in a triangulation. This struct is in the public API, as listed in the [API](../api/overview.md).

```@docs; canonical=false 
IndividualTriangleStatistics 
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/statistics/individual_triangle_statistics.jl"]
```

## TriangulationStatistics 

We also provide a function for computing statistics about a triangulation, using the `TriangulationStatistics` struct. This struct is in the public API, as listed in the [API](../api/overview.md).

```@docs; canonical=false
TriangulationStatistics
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/statistics/triangulation_statistics.jl"]
```

## InsertionEventHistory 

For mesh refinement we need a way to identify what happens to a triangulation after a point is added, in case we need to reverse the insertion. For this, we use `InsertionEventHistory` internally.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/insertion_event_history.jl"]
```

## RefinementConstraints

For mesh refinement, internally we store the user-provided constraints inside their own struct `RefinementConstraints`.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/refinement_constraints.jl"]
```

## RefinementQueue 

The mesh refinement algorithm we use requires a method for prioritising which segments and triangles need to be split. Using a `MaxPriorityQueue`, so that large segments and large triangles are prioritised, we use a `RefinementQueue` internally to determine this ordering.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/refinement_queue.jl"]
```

## RefinementArguments 

There are many arguments that need to be passed around during mesh refinement. To simplify this, we use a `RefinementArguments` struct internally.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/mesh_refinement/refinement_arguments.jl"]
```

## VoronoiTessellation

We use a `VoronoiTessellation` struct to represent a Voronoi tessellation. This struct is in the public API.

```@docs; canonical=false
VoronoiTessellation
```

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/voronoi.jl"]
```

## Polygon

We sometimes find it useful to use a specific `Polygon` struct for representing polygons.

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["data_structures/polygon.jl"]
```

## ShuffledPolygonLinkedList

For computing constrained Delaunay triangulations, we use a `ShuffledPolygonLinkedList` struct to store the polygons involved.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["data_structures/shuffled_polygon_linked_list.jl"]
```

## Points (Primitive Interface)

Here are functions that are used for defining and working with points in the package.

```@autodocs 
Modules = [DelaunayTriangulation]
Pages = ["src/geometric_primitives/points.jl"]
```

## Edges (Primitive Interface)

Here are functions that are used for defining and working with edges in the package.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/geometric_primitives/edges.jl"]
```

## Triangles (Primitive Interface)

Here are functions that are used for defining and working with triangles in the package.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/geometric_primitives/triangles.jl"]
```

## Boundary Nodes (Primitive Interface)

Here are functions that are used for defining and working with boundary nodes in the package.

```@autodocs
Modules = [DelaunayTriangulation]
Pages = ["src/geometric_primitives/boundary_nodes.jl"]
```