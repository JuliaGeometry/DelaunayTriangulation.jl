# Overview 

This section shows the docstrings for all functions not in the public API; for those in the public API, see the [API](../api/overview.md). We will describe the following:

- [Data Structures](data_structures.md): All data structures used inside the package, along with functions related to working with them. Both public and private data structures are listed here.
- [Algorithm Internals](algorithms.md): All internal functions used for algorithms inside the package.
- [Utility Functions](utils.md): All internal utility or other miscellaneous functions used inside the package.

Some of the functions shown here will also appear in the public API, although only very few; for example, `Triangulation` and its accessors are repeated here. Functions listed here are **NOT** necessarily in the public API (similarly, a function having a docstring does **NOT** guarantee that it is in the public API) - only those listed in the [API](../api/overview.md) have this guarantee. Here is an index of all the functions listed in the above pages.

```@index
Pages = ["data_structures.md", "algorithm_internals.md", "utils.md"]
```