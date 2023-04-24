```@meta
CurrentModule = DelaunayTriangulation
```

# Statistics 

We provide an interface for computing statistics of a triangulation using `statistics`. (If it is important to you, note that the methods used for computing these statistics do not use adaptive floating point arithmetic unfortunately.)

```@docs 
IndividualTriangleStatistics 
TriangulationStatistics
```

For example:

```julia
p1 = (0.0, 0.0)
p2 = (1.0, 0.0)
p3 = (1.0, 1.0)
p4 = (0.0, 1.0)
pts = [p1, p2, p3, p4]
tri = triangulate(pts)
```
```jula-repl
julia> stats = statistics(tri)
Delaunay Triangulation Statistics.
   Triangulation area: 1.0
   Number of vertices: 5
   Number of solid vertices: 4
   Number of ghost vertices: 1
   Number of edges: 9
   Number of solid edges: 5
   Number of ghost edges: 4
   Number of triangles: 2
   Number of solid triangles: 2
   Number of ghost triangles: 0
   Number of constrained boundary edges: 0
   Number of constrained interior edges: 0
   Number of constrained edges: 0
   Number of convex hull points: 4
   Smallest angle: 45.0°
   Largest angle: 90.0°
   Smallest area: 0.5
   Largest area: 0.5
   Smallest radius-edge ratio: 0.7071067811865476
   Largest radius-edge ratio: 0.7071067811865476
```

```julia
stats = refine!(tri, max_area = 1e-2) # refine! computes and returns stats = statistics(tri) post-refinement
```
```julia-repl
julia> stats
Delaunay Triangulation Statistics.
   Triangulation area: 0.9999999999999999
   Number of vertices: 98
   Number of solid vertices: 97
   Number of ghost vertices: 1
   Number of edges: 288
   Number of solid edges: 256
   Number of ghost edges: 32
   Number of triangles: 160
   Number of solid triangles: 160
   Number of ghost triangles: 0
   Number of constrained boundary edges: 0
   Number of constrained interior edges: 0
   Number of constrained edges: 0
   Number of convex hull points: 32
   Smallest angle: 30.231153476251915°
   Largest angle: 113.40705209438516°
   Smallest area: 0.0038884849847605116
   Largest area: 0.0098765990101033
   Smallest radius-edge ratio: 0.5799214659156856
   Largest radius-edge ratio: 0.9930687672569417
```

All the relevant docstrings for working with these structs are below.

```@docs 
get_all_stat
num_vertices(::TriangulationStatistics)
num_solid_vertices(::TriangulationStatistics)
num_ghost_vertices(::TriangulationStatistics)
num_edges(::TriangulationStatistics)
num_solid_edges(::TriangulationStatistics)
num_ghost_edges(::TriangulationStatistics)
num_triangles(::TriangulationStatistics)
num_solid_triangles(::TriangulationStatistics)
num_ghost_triangles(::TriangulationStatistics)
num_constrained_boundary_edges(::TriangulationStatistics)
num_constrained_interior_edges(::TriangulationStatistics)
num_constrained_edges(::TriangulationStatistics)
num_convex_hull_points(::TriangulationStatistics)
get_smallest_angle(::TriangulationStatistics)
get_largest_angle(::TriangulationStatistics)
get_smallest_area(::TriangulationStatistics)
get_largest_area(::TriangulationStatistics)
get_smallest_radius_edge_ratio(::TriangulationStatistics)
get_largest_radius_edge_ratio(::TriangulationStatistics)
get_total_area(::TriangulationStatistics)
get_individual_statistics(::TriangulationStatistics)
get_minimum_angle(::TriangulationStatistics, ::Any)
get_maximum_angle(::TriangulationStatistics, ::Any)
get_median_angle(::TriangulationStatistics, ::Any)
get_area(::TriangulationStatistics, ::Any)
get_lengths(::TriangulationStatistics, ::Any)
get_circumcenter(::TriangulationStatistics, ::Any)
get_circumradius(::TriangulationStatistics, ::Any)
get_angles(::TriangulationStatistics, ::Any)
get_radius_edge_ratio(::TriangulationStatistics, ::Any)
get_edge_midpoints(::TriangulationStatistics, ::Any)
get_aspect_ratio(::TriangulationStatistics, ::Any)
get_inradius(::TriangulationStatistics, ::Any)
get_perimeter(::TriangulationStatistics, ::Any)
get_centroid(::TriangulationStatistics, ::Any)
```