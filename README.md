# DelaunayTriangulation

![Logo](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/delaunay_triangulation_logo.png?raw=true)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/dev/)
[![Coverage](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl)
[![DOI](https://zenodo.org/badge/540660309.svg)](https://zenodo.org/badge/latestdoi/540660309)

This is a package for constructing Delaunay triangulations of planar point sets, with partial support for constrained Delaunay triangulations. The support is partial in that we support them when creating via Gmsh (see the `generate_mesh` function), but we do not yet have a method for constructing constrained Delaunay triangulations with a given piecewise linear complex; such a method is in development. All geometric predicates are computed via ExactPredicates.jl.

Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013). 

See the docs for plenty of examples.