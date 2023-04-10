```@meta
CurrentModule = DelaunayTriangulation
```

# DelaunayTriangulation

Documentation for [DelaunayTriangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl).

This is a package for computing Delaunay triangulations of planar point sets. We support both unconstrained and constrained Delaunay triangulations. An interface for computing  constrained Delaunay triangulations with Gmsh is also available if needed; see the Gmsh discussion in the sidebar. Unconstrained Delaunay triangulations are computed with the Bowyer-Watson algorithm, and constrained Delaunay triangulations are computed with the incremental algorithm given by [https://doi.org/10.1016/j.comgeo.2015.04.006](https://doi.org/10.1016/j.comgeo.2015.04.006).

To ensure that the triangulations are robust to degeneracies, we use ExactPredicates.jl for all geometrical predicates. The results from these predicates are handled through a Certificate module, as outlined in the predicates section in the sidebar.

Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013).

For the future, the top priorities would be:

1. Spatial sorting.
2. (Clipped) Voronoi tessellations.

Of course, another priority is performance, as I am certain there are places where certain steps could be optimised further, e.g. for constrained triangulations where we might be able to interleave the process of walking across a triangulation when increasing a segment and triangulating the polygons on each side. Performance has of course been considered during development, but the priority has been on correctness and just getting a working algorithm.