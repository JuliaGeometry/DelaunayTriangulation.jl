# DelaunayTriangulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/DelaunayTriangulation.jl/dev/)
[![Build Status](https://github.com/DanielVandH/DelaunayTriangulation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DanielVandH/DelaunayTriangulation.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DanielVandH/DelaunayTriangulation.jl)
[![DOI](https://zenodo.org/badge/540660309.svg)](https://zenodo.org/badge/latestdoi/540660309)

This is a package for constructing Delaunay triangulations of planar point sets, with partial support for constrained Delaunay triangulations. The support is partial in that we support them when creating via Gmsh (see the `generate_mesh` function), but we do not yet have a method for constructing constrained Delaunay triangulations with a given piecewise linear complex; such a method is in development. All geometric predicates are computed via ExactPredicates.jl.

Much of the work in this package is derived from the book *Delaunay Mesh Generation* by Cheng, Dey, and Shewchuk (2013).

## Simple example 

Here we give a simple example, but you should refer to the docs or tests for more details.

```julia
using DelaunayTriangulation, CairoMakie
pts = rand(2, 500) # This is a matrix, but we could just as well do a vector of vectors or a vector of tuples
tri = triangulate(pts)
fig, ax, sc = triplot(tri; show_convex_hull = true, convex_hull_linewidth=4)
```

![A triangulation](https://github.com/DanielVandH/DelaunayTriangulation.jl/blob/main/test/figures/custom_interface_testing.png?raw=true)
