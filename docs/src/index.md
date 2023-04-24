```@meta
CurrentModule = DelaunayTriangulation
```

# DelaunayTriangulation

Documentation for [DelaunayTriangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl).

This is a package for computing Delaunay triangulations of planar point sets. We support both unconstrained and constrained Delaunay triangulations, as well as mesh refinement with Rupper's algorithm. We also support Voronoi tessellations, clipped Voronoi tessellations, and central Voronoi tessellations; for these latter two cases, the triangulation is treated as constrained with the convex hull edges, but we do not support general boundary constraints for tessellations - if you do know about this, get in touch!

An interface for computing  constrained Delaunay triangulations with Gmsh is also available if needed; see the Gmsh discussion in the sidebar. Unconstrained Delaunay triangulations are computed with the Bowyer-Watson algorithm, and constrained Delaunay triangulations are computed with the incremental algorithm given by [https://doi.org/10.1016/j.comgeo.2015.04.006](https://doi.org/10.1016/j.comgeo.2015.04.006).

To ensure that the triangulations are robust to degeneracies, we use ExactPredicates.jl for all geometrical predicates. The results from these predicates are handled through a Certificate module, as outlined in the predicates section in the sidebar.

Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013).
