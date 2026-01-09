# Overview 

This section provides docstrings for all functions in the public API, split into sections. If you want a written description of some of the main features of this package, see the [manual](../manual/overview.md). We will describe the following:

- [Data Structures](data_structures.md): Some of the main data structures used in this package.
- [Triangulations](triangulation.md): Functions for creating and working with triangulations.
- [Triangulation Operations](operations.md): Functions for modifying triangulations.
- [Voronoi Tessellations](voronoi.md): Functions for creating and working with Voronoi tessellations.
- [Convex Hull](convex_hull.md): Functions for creating convex hulls.
- [Curves](curves.md): Functions for defining curves for use with curve-bounded domains.
- [Iterators](iterators.md): Functions for iterating over triangulation and tessellation features.
- [Point Location](point_location.md): Functions for performing point location.
- [Predicates](predicates.md): Functions for computing geometric predicates.
- [Triangulation Statistics](statistics.md): Functions for computing statistics about triangulations.
- [Primitive Interfaces](primitives.md): Functions for defining primitive interfaces.
- [Other](other.md): Other functions that don't fit into the above categories.

We emphasise that a function having a docstring does **NOT** guarantee it being in the public API - only it being listed here will guarantee this.

We also note that, for more structs such as `Triangulation` and `VoronoiTessellations`, their fields are not intended to be in the public API. Rather, their 
associated accessor functions are intended to be used, e.g. prefer `get_X(tri)` instead of `tri.X` for `tri::Triangulation`.

Each section will first start with the list of all functions to be listed, and then the docstrings of those functions will be given. There will be some docstrings that fit into multiple categories, in which case one is chosen. Here is an index of all the functions listed in the above pages.

```@index 
Pages = ["convex_hull.md", "curves.md", "data_structures.md", "iterators.md", "operations.md", "other.md", "point_location.md", "predicates.md", "primitives.md", "statistics.md", "triangulation.md", "voronoi.md"]
```