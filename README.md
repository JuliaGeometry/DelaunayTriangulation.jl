# DelaunayTriangulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/dev/)
[![Coverage](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl)
[![DOI](https://zenodo.org/badge/540660309.svg)](https://zenodo.org/badge/latestdoi/540660309)

This is a package for constructing Delaunay triangulations and Voronoi tessellations of planar point sets, with support for both unconstrained and constrained triangulations, and for mesh refinement. All geometric predicates are computed via [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl). Many features are available, some of these being:

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
- [Computation of Voronoi tessellations](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/tessellations/voronoi.md), including [clipping of polygons to the convex hull](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/tessellations/clipped.md). I hope to get this working for constrained triangulations, but it's difficult.
- Computation of [centroidal Voronoi tessellations](https://danielvandh.github.io/DelaunayTriangulation.jl/dev/tessellations/lloyd.md) using Lloyd's algorithm.

Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013). 

See the docs for plenty of examples. Feel free to use the issues tab for any suggestions, feedback, or if you have any questions about using the package, internals, etc.

Just as a nice demonstration of the incremental behaviour of the algorithms in this package, here's an example of how we build a triangulation of the Julia logo.

https://user-images.githubusercontent.com/95613936/232210266-615e08bd-d4f2-4a40-849e-e832778a4b71.mp4

Here is also a nice animation showing the computation of a centroidal Voronoi tessellation.

https://user-images.githubusercontent.com/95613936/234076094-41489505-8fe2-431e-a1a2-e8b25ecfdb8c.mp4

