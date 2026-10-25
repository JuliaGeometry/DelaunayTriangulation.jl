```@meta
CurrentModule = DelaunayTriangulation
```

# Package Structure

This page describes how the repository and the source code of the package are organised, so that you can find where something is implemented, tested, and documented. This is mostly useful if you want to read the source code or contribute to the package; please also see the [CONTRIBUTING.md](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/blob/main/CONTRIBUTING.md) file in the repository. As for the rest of this section, none of this is part of the public API, and it may change across any version.

If you are looking for a specific function rather than a feature, the quickest way to find it is to use the _source_ link on its docstring (in the [API](../api/overview.md) or in this section), or `@which` or `methods` in the REPL, e.g. `@which triangulate(rand(2, 10))`.

The layout of the repository, at a glance:

```
DelaunayTriangulation.jl/
├── assets/                     reference images for the tests and the tutorials, and the README's images
├── .github/                    the GitHub workflows
└── DelaunayTriangulation/      the package
    ├── src/
    │   ├── DelaunayTriangulation.jl   the module, which includes all the other files
    │   ├── primitives/         the interface for points, edges, triangles, and boundary nodes
    │   ├── predicates/         predicate kernels, orientation and incircle tests, and predicates built on them
    │   ├── utils/              queues and other utilities
    │   ├── geometry/           convex hulls, poles of inaccessibility, and polygon clipping
    │   ├── trees/              the R-tree and the polygon hierarchy
    │   ├── curves/             the curves for curve-bounded domains, and the BoundaryEnricher
    │   ├── triangulation/      the Triangulation type and its data structures
    │   │   ├── construction/   triangulate and the algorithms it uses
    │   │   ├── methods/        functions for working with the fields of a Triangulation
    │   │   └── operations/     add_point!, delete_point!, add_segment!, flip_edge!, ...
    │   ├── point_location/     find_triangle, find_polygon, and nearest neighbours
    │   ├── refinement/         refine!
    │   ├── statistics/         statistics of triangulations and triangles
    │   └── voronoi/            Voronoi tessellations, clipping, and centroidal tessellations
    ├── test/                   the tests, with the same folders as src/
    ├── docs/                   this documentation
    └── benchmark/              the benchmark suite
```

## The Repository

The package does not live at the top level of the repository, but in the `DelaunayTriangulation` folder. This is so that the reference images used by the tests and the documentation, which live in the top-level `assets` folder, are not downloaded with the package. The top level of the repository has:

- `DelaunayTriangulation/`: The package itself.
- `assets/`: The reference images for the tests (`assets/test/`) and for the tutorials (`assets/docs/figures/`), and the images in the README.
- `.github/`: The GitHub workflows, described in [Continuous Integration](@ref package-structure-ci) below, and the templates for issues.
- `NEWS.md`: The changes made in each version.
- `README.md`, `CONTRIBUTING.md`, `CITATION.cff`, and `LICENSE`.

The `DelaunayTriangulation` folder has the usual layout of a Julia package:

- `src/`: The source code, described in [The Source Code](@ref package-structure-src).
- `test/`: The tests, described in [The Tests](@ref package-structure-test).
- `docs/`: The documentation, described in [The Documentation](@ref package-structure-docs).
- `benchmark/`: The benchmarks, described in [The Benchmarks](@ref package-structure-benchmark).

## [The Source Code](@id package-structure-src)

The module is defined in `src/DelaunayTriangulation.jl`, which includes all the other files. The files are included in groups, in this order:

1. The data structures, e.g. the [`Triangulation`](@ref) and [`VoronoiTessellation`](@ref) types. These come first since the rest of the package uses their types.
2. The geometric predicates.
3. The interface for geometric primitives.
4. The utility functions.
5. The algorithms.
6. [`validate_triangulation`](@ref).

A new file should be included in the group that it belongs to. The other files at the top level of `src` are:

- `setup.jl`: Some constants used throughout the package, e.g. `∅` for an edge with no adjacent vertex and `𝒢` for the ghost vertex, and the switches [`toggle_inf_warn!`](@ref) and [`toggle_warn_on_dupes!`](@ref) for some of the package's warnings.
- `exports.jl`: The exported names.
- `public.jl`: The names that are in the public API but not exported, which are marked `public` on Julia 1.11 and later.

The rest of the source code is in folders, one for each feature of the package:

- `primitives/`: The interface for representing points, edges, triangles, and boundary nodes, e.g. [`getxy`](@ref) and [`construct_edge`](@ref). See [Representing Primitives](../manual/primitives.md).
- `predicates/`: The [predicate kernels](../manual/predicate_kernels.md) (`kernels.jl`), the [`Certificate`](@ref)s returned by the predicates (`certificate.jl`), the orientation and incircle predicates themselves (`predicate_definitions.jl`), and the predicates built on them, for example for testing triangles and edges of a triangulation (`predicates.jl`) or for testing for ghost vertices and boundary edges (`boundaries_and_ghosts.jl`).
- `utils/`: The queues used by the package and some miscellaneous utility functions.
- `geometry/`: Geometry that does not need a triangulation: [`convex_hull`](@ref) (`convex_hull/`), [`pole_of_inaccessibility`](@ref) (`pole_of_inaccessibility/`), polygon clipping with the Sutherland-Hodgman and Liang-Barsky algorithms (`polygon_clipping/`), and other utility functions for working with polygons (`geometry_utils.jl`).
- `trees/`: The R-tree used for finding intersections with the boundary of curve-bounded domains (`rtree.jl` and `rtree_intersections.jl`), and the [`PolygonHierarchy`](@ref) that records which boundaries are inside which (`polygon_hierarchy.jl`).
- `curves/`: The curves for [curve-bounded domains](../manual/curve_bounded.md), one file for each curve along with the interface for all curves in `abstract.jl`, and the [`BoundaryEnricher`](@ref) used for discretising the boundary (`boundary_enricher.jl`).
- `triangulation/`: Triangulations, described below.
- `point_location/`: Point location with [`find_triangle`](@ref) (`jump_and_march.jl`) and [`brute_force_search`](@ref) (`brute_force.jl`), [`find_polygon`](@ref) for finding the polygon of the boundary that contains a point, and nearest neighbour queries with [`get_nearest_neighbour`](@ref).
- `refinement/`: Mesh refinement with [`refine!`](@ref) (`mesh_refinement.jl`), along with the data structures it uses.
- `statistics/`: Statistics of a triangulation with [`statistics`](@ref), and statistics of individual triangles.
- `voronoi/`: Voronoi tessellations and power diagrams: the [`VoronoiTessellation`](@ref) type (`tessellation.jl`), its construction with [`voronoi`](@ref) (`construction.jl` and `unbounded.jl`), its clipping (`clipped.jl` and `clipped_coordinates.jl`), and [`centroidal_smooth`](@ref) (`centroidal.jl`).

The `triangulation` folder is the largest. In it are:

- `triangulation.jl`: The [`Triangulation`](@ref) type. The data structures that it stores have their own files, e.g. `adjacent.jl` for [`Adjacent`](@ref), `adjacent2vertex.jl` for [`Adjacent2Vertex`](@ref), `graph.jl` for [`Graph`](@ref), and `triangulation_cache.jl` for the [`TriangulationCache`](@ref).
- `construction/`: [`triangulate`](@ref) (`triangulate.jl`) and the checks of its arguments (`check_args.jl`), together with the algorithms that it uses: the Bowyer-Watson algorithm (`unconstrained.jl`) along with the spatial sorting of its insertion order (`spatial_sorting.jl`), the insertion of segments and boundaries (`constrained.jl`), and the discretisation of curve-bounded domains (`curve_bounded.jl`). The other ways of constructing a triangulation are also here: [`triangulate_convex`](@ref) (`convex.jl` and `shuffled_polygon_linked_list.jl`) and [`triangulate_rectangle`](@ref) (`rectangle.jl`).
- `methods/`: The functions for working with the fields of a `Triangulation`, with one file for each field. For example, `methods/adjacent.jl` has [`get_adjacent`](@ref), and `methods/boundary_nodes.jl` has the functions for working with the boundary nodes. `methods/iterators.jl` has the iterators like [`each_solid_triangle`](@ref), and `methods/checks.jl` has queries like [`is_constrained`](@ref).
- `operations/`: The operations that modify a triangulation, with one file for each operation, e.g. [`add_point!`](@ref) is in `add_point.jl` and [`flip_edge!`](@ref) is in `flip_edge.jl`. See the [operations tutorials](../tutorials/operations_vertex_insertion_deletion.md).
- `validation.jl`: [`validate_triangulation`](@ref), for checking that a triangulation is valid.

### What `triangulate` Does

The main steps of [`triangulate`](@ref), with the files in `src/triangulation/` that they are in, are:

1. For a curve-bounded domain, the curves are discretised and the boundary is enriched (`construction/curve_bounded.jl`), and `triangulate` is then called again with the new points and boundary nodes.
2. The input is checked (`construction/check_args.jl`).
3. The points are inserted with the Bowyer-Watson algorithm (`construction/unconstrained.jl`), in the order from `construction/spatial_sorting.jl`, finding the triangle containing each new point with [`find_triangle`](@ref).
4. The segments and boundary edges are inserted with [`add_segment!`](@ref) (`operations/add_segment.jl` and `construction/constrained.jl`), and the triangles outside of the domain are deleted (`operations/delete_holes.jl`).
5. The ghost triangles, empty features, and representative points are updated (`construction/triangulate.jl`).

## [The Tests](@id package-structure-test)

The folders in `test` mirror those in `src`, e.g. the tests for `src/voronoi/` are in `test/voronoi/`. All the tests are run by `test/runtests.jl`, in groups, e.g. `"Voronoi"` and `"Refinement"`. The full test suite takes over 20 minutes, so while working on a change you might want to only run some of the groups, which you can do by passing their names as test arguments, e.g.

```julia
import Pkg
Pkg.test("DelaunayTriangulation"; test_args = ["Voronoi", "Point Location"])
```

The names of all the groups are listed at the start of `test/runtests.jl`. Some other files in `test` are:

- `helper_functions.jl`: The `HelperFunctions` module of functions used throughout the tests. They are split by purpose into the files in `test/helpers/`, e.g. `examples.jl` has example geometries, `validation.jl` has validators like `validate_statistics`, `comparisons.jl` has comparisons that do not depend on the order of iteration over a `Set` or `Dict`, and `reference_implementations.jl` has slow but simple implementations of some of the package's algorithms that the tests compare against.
- `docstring_examples.jl`: The examples shown in the docstrings.
- `readme_example.jl`: The example in the README.
- `makie/makie.jl`: Tests of the plotting recipes, like `triplot` and `voronoiplot`, that Makie.jl defines for this package.

The tutorials and applications in the documentation are also run by the tests, in the `"Run the documentation examples"` group. The lines in those files that end with `#src` do not appear in the documentation, and they are used for testing the tutorials, e.g. by comparing their figures with the reference images in `assets/docs/figures/`.

## [The Documentation](@id package-structure-docs)

The documentation is built by `docs/make.jl`. Every page is listed in its `_PAGES`, and the build fails if a Markdown file in `docs/src` is not listed there, if a docstring is not on any page, or if a cross-reference such as `[triangulate](@ref)` is broken. The pages are in `docs/src`, in a folder for each section, e.g. `manual/`, `api/`, `extended/`, and `math/`. The exceptions are the tutorials and applications, which are written as [Literate.jl](https://github.com/fredrikekre/Literate.jl) scripts in `docs/src/literate_tutorials/` and `docs/src/literate_applications/`. These are converted by `docs/make.jl` into the pages in `docs/src/tutorials/` and `docs/src/applications/`, so it is the scripts that should be edited.

The docstrings of the public API are on the [API](../api/overview.md) pages, listed by name in `@docs` blocks. All other docstrings are in this section, mostly in `@autodocs` blocks that list the source files whose docstrings they show. A new source file that has docstrings therefore needs to be added to one of these lists, or its docstrings to an `@docs` block.

The timings on the [Performance Measurements](../manual/performance_measurements.md) page are not measured when the documentation is built, but by the script `docs/performance_measurements.jl`, which saves them to `docs/src/manual/performance_measurements.toml`. Run it again to update them, e.g. for a new version of the package.

To view the documentation locally while you edit it, run `docs/liveserver.jl`. This uses [LiveServer.jl](https://github.com/JuliaDocs/LiveServer.jl) to build the documentation and update it in your browser whenever a file changes.

## [The Benchmarks](@id package-structure-benchmark)

The `benchmark` folder has a benchmark suite, `benchmarks.jl`, that covers most of the package, as well as the scripts `run.jl`, for running the benchmarks for any version of the package, and `compare.jl`, for comparing the results of two runs. See the CONTRIBUTING.md file for how to use them.

## [Continuous Integration](@id package-structure-ci)

The GitHub workflows are in `.github/workflows/`:

- `CI.yml`: Runs the tests on several Julia versions on Linux, and on the latest release on Windows and macOS, and builds the documentation. For each Julia version and operating system, the test groups are split between several jobs that run in parallel, and a final job, using `.github/test_summary.jl`, puts their results together and checks that every test file was run. A new test group needs to be added to one of these jobs.
- `Benchmark.yml`: For a pull request that changes the source code or the benchmarks, compares the benchmarks of the pull request with those of the branch it is merging into. The results are in the summary of the workflow's run.
- `Downstream.yml`: Runs the tests of some packages that depend on this package, such as [NaturalNeighbours.jl](https://github.com/DanielVandH/NaturalNeighbours.jl) and [FiniteVolumeMethod.jl](https://github.com/SciML/FiniteVolumeMethod.jl), against the new version.
- `enforce_news.yml`: Checks that a pull request adds an entry to `NEWS.md`.
- `DocCleanup.yml`, `CompatHelper.yml`, and `TagBot.yml`: Remove the documentation previews of closed pull requests, update the compat bounds of the dependencies, and tag new releases, respectively.
