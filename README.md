# DelaunayTriangulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGeometry.github.io/DelaunayTriangulation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeometry.github.io/DelaunayTriangulation.jl/dev/)
[![Coverage](https://codecov.io/gh/JuliaGeometry/DelaunayTriangulation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGeometryH/DelaunayTriangulation.jl)
[![DOI](https://zenodo.org/badge/540660309.svg)](https://zenodo.org/badge/latestdoi/540660309)

This is a package for constructing Delaunay triangulations and Voronoi tessellations of planar point sets. Supports unconstrained and constrained triangulations, mesh refinement, triangulation of curve bounded domains, Voronoi tessellations, and clipped and centroidal Voronoi tessellations. All geometric predicates are computed via [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl). To install the package, do

```julia
julia>] add DelaunayTriangulation
```

Many features are available, some of these being:

- Unconstrained and constrained Delaunay triangulations, supporting many types of domains.
- Computation of Voronoi tessellations, clipped Voronoi tessellations where the Voronoi tiles get clipped to the convex hull, and centroidal Voronoi tessellations where each Voronoi tile's generator is the tile's centroid.
- Mesh refinement, with support for custom angle and area constraints, as well as refinement of curve-bounded domains.
- Dynamic point insertion, point deletion, and segment insertion, amongst many other operations.
- Computation of convex hulls.
- Triangulation of convex polygons.
- Point location.
- Computation of the pole of inaccessibility.
- The interface for defining geometric primitives is fully customisable.

To ensure that the algorithms are robust, we use [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl) to define all geometric predicates in this package. Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013). Please see the documentation for much more information.

Some examples are below (and in the documentation), but if you would also like to see how DelaunayTriangulation.jl is used in other packages, see [FiniteVolumeMethod.jl](https://github.com/SciML/FiniteVolumeMethod.jl) (solving 2D PDEs) and [NaturalNeighbours.jl](https://github.com/DanielVandH/NaturalNeighbours.jl) (scattered data interpolation).

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
- [GMT.jl](https://github.com/GenericMappingTools/GMT.jl): A wrapper of [GMT](https://github.com/GenericMappingTools/gmt), allowing for [unconstrained Delaunay triangulations in two dimensions](https://www.generic-mapping-tools.org/GMTjl_doc/documentation/modules/triangulate/index.html#triangulate), and for [spherical triangulation, i.e. triangulation of points lying on a sphere](https://www.generic-mapping-tools.org/GMTjl_doc/documentation/modules/sphtriangulate/index.html#sphtriangulate).

