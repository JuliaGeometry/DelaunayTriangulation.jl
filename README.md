# DelaunayTriangulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/dev/)
[![Coverage](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl)
[![DOI](https://zenodo.org/badge/540660309.svg)](https://zenodo.org/badge/latestdoi/540660309)

This is a package for constructing Delaunay triangulations and Voronoi tessellations of planar point sets. Supports unconstrained and constrained triangulations, mesh refinement, Voronoi tessellations, and clipped and centroidal Voronoi tessellations. All geometric predicates are computed via [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl). Many features are available, some of these being:

- [Geometric predicates](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/predicates/) are implemented with [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl), and many predicates have been extended from ExactPredicates.jl.
- [Unconstrained](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/triangulations/unconstrained/) and [constrained](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/triangulations/constrained/) triangulations. Support is provided for many types of domains, as given in the docs.
- [Mesh refinement](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/triangulations/refinement/), with support for custom angle and area constraints.
- Dynamic point insertion, point deletion, and segment insertion, amongst many other [operations](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/operations/).
- Computation of convex hulls, either [from the triangulation itself](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/data_structures/convex_hull/) or using [the monotone chain algorithm](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/other_features/convex_hull/). 
- [Triangulation of convex polygons](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/tri_algs/convex/).
- [Efficient point location](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/other_features/point_location/) on convex geometries, even with interior holes. Partial support exists for non-convex geometries, although it is much slower and not perfect.
- [Computation of the pole of inaccessibility](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/other_features/pole_of_inaccessibility/), i.e. the point in a polygon that is furthest from a boundary (see e.g. [this blogpost](https://blog.mapbox.com/a-new-algorithm-for-finding-a-visual-center-of-a-polygon-7c77e6492fbc)).
- [Fully customisable interface](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/interface/interface/) for defining geometric primitives.
- [Simple iteration over mesh elements, including points, edges, or triangles](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/data_structures/triangulation/).
- Computation of [statistics](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/data_structures/statistics/) over individual triangular elements and over a complete triangulation.
- [Computation of Voronoi tessellations](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/tessellations/voronoi/), including [clipping of polygons to the convex hull](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/tessellations/clipped/). I hope to get this working for constrained triangulations, but it's difficult.
- Computation of [centroidal Voronoi tessellations](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/tessellations/lloyd/) using Lloyd's algorithm.

Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013). Feel free to use the issues tab for any suggestions, feedback, or if you have any questions about using the package, internals, etc.

## Quick Example 

See the docs for plenty of examples. For now, here are some small examples.

```julia
using DelaunayTriangulation, CairoMakie, StableRNGs

fig = Figure(fontsize=24)

## Unconstrained example: Just some random points 
rng = StableRNG(2)
pts = randn(rng, 2, 500)
tri = triangulate(pts; rng)
ax = Axis(fig[1, 1], title = "(a): Unconstrained", titlealign = :left, width=400,height=400)
triplot!(ax, tri) 

## Constrained example: Generate some points, convert to indices, triangulate 
points = [(0.0, 0.0), (1.0, 0.0), (1.0, 0.3), (0.8, 0.3), 
(0.8, 0.5), (1.0, 0.5), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
boundary_nodes, points = convert_boundary_points_to_indices(points)
edges = Set(((2, 8), (5, 7)))
tri = triangulate(points; boundary_nodes, edges, rng)
ax = Axis(fig[1, 2], title = "(b): Constrained, pre-refinement", titlealign = :left, width=400,height=400)
triplot!(ax, tri)

## Dynamic updating (Warning: Less supported on non-convex geometries)
n = num_points(tri)
add_point!(tri, (0.2, 0.2); rng)
add_point!(tri, (0.8, 0.8); rng)
add_point!(tri, (0.3, 0.2); rng)
add_point!(tri, (0.6, 0.2); rng)
add_point!(tri, (0.3, 0.3); rng)
add_edge!(tri, (n+1, n+4); rng)
delete_point!(tri, n+2; rng)
ax = Axis(fig[1, 3], title = "(c): Updated constrained, pre-refinement", titlealign = :left, width=400,height=400)
triplot!(ax, tri)

## Refinement example
A = get_total_area(tri)
refine!(tri; max_area = 1e-2A, min_angle = 28.7, rng)
ax = Axis(fig[1, 4], title = "(d): Constrained, post-refinement", titlealign = :left, width=400,height=400)
triplot!(ax, tri, show_convex_hull = false)

## Constrained example with holes
outer_boundary = [[(0.0, 0.0), (1.0, 0.0), (1.0, 0.3), (0.8, 0.3), 
(0.8, 0.5), (1.0, 0.5), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]]
inner_boundary = [[(0.2, 0.2), (0.2, 0.6), (0.4, 0.6), (0.4, 0.2), (0.2, 0.2)]]
boundary_nodes, points = convert_boundary_points_to_indices([outer_boundary, inner_boundary])
edges = Set(((2, 8), (5, 7)))
tri = triangulate(points; boundary_nodes, edges, rng)
A = get_total_area(tri)
refine!(tri; max_area = 1e-3A, min_angle = 31.5, rng)
ax = Axis(fig[2, 1], title = "(e): Multiply-connected", titlealign = :left, width=400,height=400)
triplot!(ax, tri, show_convex_hull = false)

## Voronoi tessellation: Make tessellations from their dual triangulation
pts = 25randn(rng, 2, 500)
tri = triangulate(pts; rng)
vorn = voronoi(tri)
ax = Axis(fig[2, 2], title = "(f): Voronoi tessellation", titlealign = :left, width=400,height=400)
voronoiplot!(ax, vorn, show_generators = false) 
xlims!(ax, -120, 120); ylims!(ax, -120, 120)

## Clipped Voronoi tessellation 
vorn = voronoi(tri, true)
cmap = cgrad(:jet)
colors = get_polygon_colors(vorn, cmap)
ax = Axis(fig[2, 3], title = "(g): Clipped Voronoi tessellation", titlealign = :left, width=400,height=400)
voronoiplot!(ax, vorn, show_generators = false, strokecolor = :red, strokewidth = 0.2, polygon_color = colors)

## Centroidal Voronoi tessellation (CVT)
smooth_vorn = centroidal_smooth(vorn; maxiters = 2500)
ax = Axis(fig[2, 4], title = "(h): Centroidal Voronoi tessellation", titlealign = :left, width=400,height=400)
voronoiplot!(ax, smooth_vorn, show_generators = false, strokecolor = :red, strokewidth = 0.2, polygon_color = colors)

resize_to_layout!(fig)
fig
```

![example_triangulations](https://user-images.githubusercontent.com/95613936/235439717-6e0258e4-c071-4578-b7bc-7f33b9ad401d.png)

## Animations

Just as a nice demonstration of the incremental behaviour of the algorithms in this package, here's an example of how we build a triangulation of the Julia logo.

https://user-images.githubusercontent.com/95613936/232210266-615e08bd-d4f2-4a40-849e-e832778a4b71.mp4

Here is also a nice animation showing the computation of a centroidal Voronoi tessellation.

https://user-images.githubusercontent.com/95613936/234076094-41489505-8fe2-431e-a1a2-e8b25ecfdb8c.mp4

## Similar Packages

This is not the only Delaunay triangulation package available. Some others are:

- [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl): A pure Julia library that constructs planar triangulations and tessellations like in this package, although no support for constrained triangulations / mesh refinement or clipped / centroid tessellations. Restricts points to $[1, 2] \times [1, 2]$.
- [VoronoiCells.jl](https://github.com/JuliaGeometry/VoronoiCells.jl): A pure Julia library that extends VoronoiDelaunay.jl. This package provides useful tools for constructing and working with Voronoi tessellations. Supports clipping Voronoi cells to a specified rectangle. Like VoronoiDelaunay.jl, restricts points to $[1, 2] \times [1, 2]$.
- [Delaunay.jl](https://github.com/eschnett/Delaunay.jl): Wraps Python's main Delaunay triangulation library, [`scipy.spatial.Delaunay`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html), for computing Delaunay triangulations in $\mathbb R^N$. I don't believe constrained triangulations or mesh refinement is available here.
- [MiniQhull.jl](https://github.com/gridap/MiniQhull.jl): Wraps [Qhull](http://www.qhull.org/) for computing unconstrained Delaunay triangulations in $\mathbb R^N$. No support is provided for mesh refinement.
- [DirectQhull.jl](https://github.com/JuhaHeiskala/DirectQhull.jl/): Similar to MiniQhull.jl, although also provides support for convex hulls and Voronoi tessellations from Qhull.
- [Delaunator.jl](https://github.com/JuliaGeometry/Delaunator.jl): A pure Julia library modelled after the [JavaScript Delaunator library](https://github.com/mapbox/delaunator). This package can construct unconstrained triangulations of planar point sets. No support is available for constrained triangulations or mesh refinement, although support exists for computing the dual Voronoi tessellation. Centroidal tessellations are not implemented, although the Voronoi cells can be clipped to a bounding box. 
- [TriangleMesh.jl](https://github.com/konsim83/TriangleMesh.jl), [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl), [Triangle.jl](https://github.com/cvdlab/Triangle.jl): Interfaces to [Shewchuk's Triangle library](https://www.cs.cmu.edu/~quake/triangle.html).
- [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl): This is for Delaunay tetrahedralisation, wrapping [TetGen](https://wias-berlin.de/software/index.jsp?id=TetGen).


