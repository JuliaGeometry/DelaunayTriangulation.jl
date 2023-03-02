```@meta
CurrentModule = DelaunayTriangulation
```

# DelaunayTriangulation

Documentation for [DelaunayTriangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl).

This is a package for computing Delaunay triangulations of planar point sets. We currently only support unconstrained Delaunay triangulations, but we do have an interface that supports constrained Delaunay triangulations. Constrained Delaunay triangulations can be computed using Gmsh if needed; see the Gmsh discussion in the sidebar. Unconstrained Delaunay triangulations are computed with the Bowyer-Watson algorithm.

To ensure that the triangulations are robust to degeneracies, we use ExactPredicates.jl for all geometrical predicates. The results from these predicates are handled through a Certificate module, as outlined in the predicates section in the sidebar.

Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013).

For the future, the top priorities would be:

1. Triangulation of convex polygons, and hence support for vertex deletion.
2. Constrained Delaunay triangulation.
3. Spatial sorting.
4. (Clipped) Voronoi tessellations.