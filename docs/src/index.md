```@meta
CurrentModule = DelaunayTriangulation
```

# Introduction 

This is the documentation for DelaunayTriangulation.jl. [Click here to go back to the GitHub repository](https://github.com/JuliaGeometry/DelaunayTriangulation.jl).

This is a package for computing Delaunay triangulations and Voronoi tessellations of points in two dimensions, amongst many other features:

- Unconstrained and constrained Delaunay triangulations, supporting many types of domains.
- Computation of Voronoi tessellations, clipped Voronoi tessellations where the Voronoi tiles get clipped to the convex hull, and centroidal Voronoi tessellations where each Voronoi tile's generator is the tile's centroid.
- Mesh refinement, with support for custom angle and area constraints, as well as refinement of curve-bounded domains.
- Dynamic point insertion, point deletion, and segment insertion, amongst many other operations.
- Computation of convex hulls.
- Triangulation of convex polygons.
- Point location.
- Computation of the pole of inaccessibility.
- The interface for defining geometric primitives is fully customisable.

To ensure that the algorithms are robust, we use [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl) to define all geometric predicates in this package. (ExactPredicates.jl can be disabled, as described [here](manual/disabling_ea.md).) Much of the work in this package is derived from the book [*Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013)](https://people.eecs.berkeley.edu/~jrs/meshbook.html).

# Documentation Structure 

The documentation for this package is broken into several sections.

- **Tutorials**: These are examples that introduce the functionality of the package, demonstrating how to perform many operations. 
- **Manual**: Some particular descriptions of how certain things work in the package, such as representing boundaries and working with the triangulation and tessellation data structures. 
- **API Reference**: This section lists docstrings for all functions in the public API.
- **Extended Reference**: This section is for providing more details about the internals of the package, meaning not in the public API.The material in this section is subject to change across any version, but may be useful for you to understand how certain functions work. For example, all utility functions are documented here.
- **Mathematical Details**: This section is for describing in detail the mathematics that underpins the algorithms in this package, for example walking you through the theory behind the algorithm for computing constrained Delaunay triangulations. 
- **Example Applications**: In case you want to use this package for some of your applications, it might be useful to see how it has been applied in certain situations. Only a few applications are considered here, but more could be added in the future.
- **Terminology**: This section provides a short description of some of the terminology used throughout the package.

If you see anything missing in any of these sections, or anything you think could be improved, feel free to file an issue.

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
pkg> add DelaunayTriangulation
```

# License 

DelaunayTriangulation.jl is provided under a [MIT license](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/blob/main/LICENSE).

# Contributing

Contributions are very welcome. If you encounter any issues or would like to make any feature requests and suggestions, you are encouraged to open an issue - any bug reports should be accompanied by a minimal example that reproduces your bug. You are also highly encouraged to make simple pull requests fixing any grammar or spelling issues you see in the documentation, or fixing any unclear explanations (or make an issue raising your concern instead if you prefer).

The issues tab also lists features that would be nice to have in the package. If you would like to contribute towards any of those features, or towards any other significant enhancements, you are recommended to first post on that issue about your ideas before committing towards a complete implementation. 

## Writing a Pull Request 

When contributing in the form of a pull request, there are a few important features that should be present, listed below. The point of these requirements is not to make you concerned about the amount of work needed to fulfill them, but to ensure that your contribution can be accepted more readily without the reviewer and yourself needing to go back-and-forth to meet the package's standards. If you do not meet them, you may still be fine depending on what you are contributing. You can always ask for help.

1. **You should outline what your pull request does in the description**. You should also link back to an associated issue if applicable.
2. **Document any new functions**. If your contribution involves any new functions, make sure they are well documented even if they are only for internal use.
3. **You need to include appropriate tests for your contribution**. For small changes, simple tests are fine, but for larger changes that implement or change an algorithm, understanding what tests you need to write is crucial when you are implementing a new feature. For example, if you wanted to implement a new algorithm for inserting a curve into a triangulation, you need to make sure that you test things like (1) point sets in general position, (2) point sets with many cocircular points, (3) collinear edges, and so on. In particular, degeneracies are some of the most important things to test for. This is not simple to do if you are inexperienced, so do feel free to ask for guidance.
4. **Your contribution should only do one thing**. To make sure that it is easy to track your changes, and to make it easier to track regressions across versions, please try to only make your contribution do one thing. For example, do not move around a bunch of files while at the same time implementing a function somewhere else - make two pull requests.
5. **Add to NEWS.md**. Describe your change in the NEWS.md file.

If you are just contributing something minimal, for example a typo, even a blank PR is fine so long as it is obvious what you have done.

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
