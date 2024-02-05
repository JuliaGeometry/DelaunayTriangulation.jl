```@meta
CurrentModule = DelaunayTriangulation
```

# DelaunayTriangulation

Documentation for [DelaunayTriangulation](https://github.com/JuliaGeometry/DelaunayTriangulation.jl).

This is a package for computing Delaunay triangulations of planar point sets. We support both unconstrained and constrained Delaunay triangulations, as well as mesh refinement with Rupper's algorithm. We also support Voronoi tessellations, clipped Voronoi tessellations, and central Voronoi tessellations; for these latter two cases, the triangulation is treated as constrained with the convex hull edges, but we do not support general boundary constraints for tessellations - if you do know about this, get in touch!

An interface for computing  constrained Delaunay triangulations with Gmsh is also available if needed; see the Gmsh discussion in the sidebar. Unconstrained Delaunay triangulations are computed with the Bowyer-Watson algorithm, and constrained Delaunay triangulations are computed with the incremental algorithm given by [https://doi.org/10.1016/j.comgeo.2015.04.006](https://doi.org/10.1016/j.comgeo.2015.04.006).

To ensure that the triangulations are robust to degeneracies, we use ExactPredicates.jl for all geometrical predicates. The results from these predicates are handled through a Certificate module, as outlined in the predicates section in the sidebar.

Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013).

Some of the short-term plans for future features are outlined in the issues page, two important ones being:

- Clipped Voronoi tessellations on general domains.
- Spatial sorting.

But we also have some long-term plans, like

- Supporting curved boundaries. This is a bit complex since the boundary node interface must be extended to support not only segments, but general forms such as splines, while not slowing down the code in general for straight lines.
- Efficient reconstruction of a perturbed unconstrained or constrained triangulations, using e.g. star-splaying (that currently only works on unconstrained).
- Weighted Delaunay triangulations. These are not actually that complicated, but I just don't have plans to do it yet. Solutions to Exercise 2 of the *Delaunay Mesh Generation* book above in Chapter 4 would also be useful here.
- Anisotropic mesh generation.

One other feature would be three-dimensional Delaunay tetrahedralisation as a replacement for `TetGen`, just as this package (hopefully) can be used as a replacement for the `Triangle` software. This is definitely the furthest away, but it is something I'm interested in eventually giving my (very limited) time to.