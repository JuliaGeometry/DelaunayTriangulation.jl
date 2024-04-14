```@meta
CurrentModule = DelaunayTriangulation
```

# Overview

This section presents a series of basic tutorials. The point of these tutorials is to give a quick overview of how many applications can be performed, rather than providing a detailed introduction into how the package works and, for example, what inputs each important function takes. See the later sections in the sidebar for more detail than what we give here, and the introduction for a description of each other section. More complicated applications of how this package can be used are given in the [Example Applications](../applications/overview.md) section.

In the tutorials, we make consistent use of the package [`StableRNGs.jl`](https://github.com/JuliaRandom/StableRNGs.jl) to make the random results reproducible across Julia versions. In practice, you could either

1. Just ignore worrying about the random number generation entirely,
2. Do something like `Random.seed!(seed)` (after doing `using Random`).

The first option is fine for most applications where the points sets are in general position (meaning no four points are collinear), and the second result is fine and easier to apply more generally.

Lastly, in all but one of the tutorials we consider only a basic interface for defining the geometric primitives, e.g. matrices or vectors for points and sets for edges. This interface is customisable, as described in [Representing Primitives](../manual/primitives.md) and [Representing Boundaries](../manual/boundaries.md) (with the API defined in [Primitive Interfaces](../api/primitives.md)). The last tutorial will demonstrate how to use custom structs for 
representing these primitives.

The tutorials that we consider are as follows:

- [Unconstrained Triangulations](unconstrained.md): How Delaunay triangulations can be computed and accessed. We also give examples here of how we iterate over the triangulation's vertices, edges, and triangles, some warnings about having to be careful about ghost objects defined in the [manual](../manual/ghost_triangles.md), and how to query neighbour information.
- Constrained Triangulations: How constrained Delaunay triangulations can be computed, and how different types of boundaries can be represented for triangulating. This tutorial is broken into multiple sub-tutorials, starting with considering [constrained segments](constrained_edges.md) and ending with considering [multipolygons](constrained_multipolygon.md).
- Triangulation Operations: How certain operations such as vertex insertion and deletion can be applied to an existing triangulation. This tutorial is broken into multiple sub-tutorials, demonstrating the multiple operations available for use, starting with [vertex insertion and deletion](operations_vertex_insertion_deletion.md).
- [Mesh Refinement](refinement.md): How triangulations can be refined to meet certain quality constraints.
- [Triangulating Rectangular Regions](lattice.md): A simple example of how rectangular regions in the plane can be triangulated quickly.
- [Triangulating Convex Polygons](convex.md): How triangulations of convex polygons can be computed.
- [Triangulation Curve-Bounded Domains](curve_bounded.md): How triangulations of curve-bounded domains can be computed.
- [Voronoi Tessellations](voronoi.md): How Voronoi tessellations can be computed. We also give examples of how you can iterate over the edges and generators in the tessellation, and other features.
- [Clipped Voronoi Tessellations](clipped.md): How to compute a Voronoi tessellation of a point set such that the polygons are [clipped to the point set's convex hull](clipped.md), or [to an arbitrary rectangle](clipped_rectangle.md).
- [Centroidal Voronoi Tessellation](centroidal.md): How to compute centroidal Voronoi tessellations, in particular how to shift a given set of points so that each point's associated Voronoi tile is that tile's centroid.
- [Point Location](point_location.md): How to use a triangulation to perform point location.
- [Nearest Neighbour Queries](nearest.md): How to use a Voronoi tessellation to find a point's nearest neighbour.
- [Convex Hulls](convex_hull.md): How to compute a convex hull of a point set, using either an existing triangulation or from scratch.
- [Pole of Inaccessibility](pole_of_inaccessibility.md): How to compute the [pole of inaccessibility](https://blog.mapbox.com/a-new-algorithm-for-finding-a-visual-center-of-a-polygon-7c77e6492fbc) of a polygon, also considering cases of multipolygons (polygons composed of multiple disjoint polygons) and multiply-connected polygons.
- [Point-in-Polygon Testing](point_in_polygon.md): How to perform point-in-polygon testing. The method we give is obviously not going to be the fastest out of all available methods, but the algorithms in this package naturally provide a nice reasonably fast way to do this testing.
- [Using Custom Structs for Primitives and Boundaries](custom_primitive.md): How to use custom structs for defining the geometric primitives.