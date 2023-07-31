```@meta
CurrentModule = DelaunayTriangulation
```

# Introduction 

This is the documentation for DelaunayTriangulation.jl. [Click here to go back to the GitHub repository](DanielVandH/DelaunayTriangulation.jl).

This is a package for computing Delaunay triangulations and Voronoi tessellations of points in two dimensions, amongst many other features:

- Unconstrained and constrained Delaunay triangulations, supporting many types of domains.
- Computation of Voronoi tessellations, clipped Voronoi tessellations where the Voronoi tiles get clipped to the convex hull, and centroidal Voronoi tessellations where each Voronoi tile's generator is the tile's centroid.
- Mesh refinement, with support for custom angle and area constraints.
- Dynamic point insertion, point deletion, and segment insertion, amongst many other operations.
- Computation of convex hulls.
- Triangulation of convex polygons.
- Point location.
- Computation of the pole of inaccessibility.
- The interface for defining geometric primitives is fully customisable.

To ensure that the algorithms are robust, we use [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl) to define all geometric predicates in this package. Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013).

# Documentation Structure 

The documentation for this package is broken into several sections.

- **Tutorials**: These are examples that introduce the functionality of the package, demonstrating how to perform many operations. 
- **Manual**: Some particular descriptions of how certain things work in the package, such as representing boundaries and working with the triangulation and tessellation data structures. 
- **API Reference**: This section lists what is available from the public API, along with all relevant docstrings for these functions. There is some overlap between this section and the **Manual** section, with the difference being that this section is for providing the code itself and **Manual** is for describing how certain things work.
- **Extended Manual**: This section is for providing more details about the internals of the package, meaning not in the public API. The material in this section is subject to change across any version, but may be useful for you to understand how certain functions work. For example, all utility functions are documented here.
- **Mathematical and Implementation Details**: This section is for describing in detail the mathematics that underpins the algorithms in this package, for example walking you through the theory behind the algorithm for computing constrained Delaunay triangulations. In addition, how some of these algorithms are actually implemented is described, giving you a sense of how we translate theory into the actual code.
- **Example Applications**: In case you want to use this package for some of your applications, it might be useful to see how it has been applied in certain situations. In this section, we provide many examples of how the package can be applied, ranging from path-finding with obstacles to cellular biology simulations to finding roots of complex functions.

# Citing DelaunayTriangulation.jl

If you use DelaunayTriangulation.jl, please cite it. For now, we have a Zenodo record available at [10.5281/zenodo.7456989](https://doi.org/10.5281/zenodo.7456989).

# Installation

You can install DelaunayTriangulation.jl using the package manager:

```julia
import Pkg 
Pkg.add("DelaunayTriangulation")
```

Alternatively, the Pkg REPL can be used (accessed by typing `]` at the `julia>` prompt):

```julia
(@v1.9) pkg> add DelaunayTriangulation
```

# License 

DelaunayTriangulation.jl is provided under a [MIT license](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/LICENSE).

# Similar Packages 

This is not the only package available in Julia for working with Delaunay triangulations and Voronoi tessellations, although DelaunayTriangulation.jl is currently the most developed native Julia package available. Some other packages are:

- [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl): A pure Julia library that constructs planar triangulations and tessellations like in this package, although no support for constrained triangulations / mesh refinement or clipped / centroid tessellations. Restricts points to $[1, 2] \times [1, 2]$.
- [VoronoiCells.jl](https://github.com/JuliaGeometry/VoronoiCells.jl): A pure Julia library that extends VoronoiDelaunay.jl. This package provides useful tools for constructing and working with Voronoi tessellations. Supports clipping Voronoi cells to a specified rectangle. Like VoronoiDelaunay.jl, restricts points to $[1, 2] \times [1, 2]$.
- [Delaunay.jl](https://github.com/eschnett/Delaunay.jl): Wraps Python's main Delaunay triangulation library, [`scipy.spatial.Delaunay`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html), for computing Delaunay triangulations in $\mathbb R^N$. Constrained triangulations or mesh refinement are not available here.
- [MiniQhull.jl](https://github.com/gridap/MiniQhull.jl): Wraps [Qhull](http://www.qhull.org/) for computing unconstrained Delaunay triangulations in $\mathbb R^N$. No support is provided for mesh refinement.
- [DirectQhull.jl](https://github.com/JuhaHeiskala/DirectQhull.jl/): Similar to MiniQhull.jl, although also provides support for convex hulls and Voronoi tessellations from Qhull. 
- [Delaunator.jl](https://github.com/JuliaGeometry/Delaunator.jl): A pure Julia library modelled after the [JavaScript Delaunator library](https://github.com/mapbox/delaunator). This package can construct unconstrained triangulations of planar point sets. No support is available for constrained triangulations or mesh refinement, although support exists for computing the dual Voronoi tessellation. Centroidal tessellations are not implemented, although the Voronoi cells can be clipped to a bounding box. 
- [TriangleMesh.jl](https://github.com/konsim83/TriangleMesh.jl), [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl), [Triangle.jl](https://github.com/cvdlab/Triangle.jl): Interfaces to [Shewchuk's Triangle library](https://www.cs.cmu.edu/~quake/triangle.html).
- [TetGen.jl](https://github.com/JuliaGeometry/TetGen.jl): This is for Delaunay tetrahedralisation, wrapping [TetGen](https://wias-berlin.de/software/index.jsp?id=TetGen).

