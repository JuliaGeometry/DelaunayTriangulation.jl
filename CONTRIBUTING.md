# Contributing to DelaunayTriangulation.jl 

Contributions are more than welcome. For the advice below, please take it more as loose guidelines rather than a strict set of rules you must follow. For most cases, please remember to update the NEWS.md file if you have made a pull request.

## Documentation Improvements 

An easy way to contribute to the package is to make pull requests fixing any issues you see with documentation, for example correcting any typos or rewriting an unclear sentence. In these cases, a simple pull request correcting the issue is sufficient. If the improvement is a bit more substantial, make sure to explain the issue clearly in the description of your pull request.

## Writing a Pull Request

For changes beyond simple documentation improvements, your pull request description should have a bit more detail. Some loose requirements are listed below. The point of these requirements is not to make you concerned about the amount of work needed to fulfill them, but to ensure that your contribution can be accepted easily without the reviewer and yourself needing to go back-and-forth to meet the package's standards. If you do not meet them, you may still be fine depending on what you are contributing. You can always ask for help.

### 1. Outline what your pull request does in the description

For example, if you are fixing a bug, you should link to an associated issue (if there is one) and explain your fix. It may not be easy to understand your contribution if you have an empty description and a complicated edit.

### 2. Document any new functions

Unlike most packages, this package tries to document all internal functions. If your contribution involves any new functions, please make sure they are well documented - regardless of whether they are only for internal use. The form of these docstrings should be similar to other docstrings in the package. A good general format to use is along the lines of 
```julia
"""
    fnc(args...; kwargs...) -> ReturnType 

<description of function>

# Arguments 
- `arg1`: <description>
- ...

# Keyword Arguments 
- `kwarg1`: <description>
- ...

# Output 
- `out1`: <description>
- ...
"""
fnc(args..; kwargs...) = ...
```

For shorter functions, a short-form docstring like

```julia
"""
    add_nums(a, b) -> Number 

Given `a` and `b`, returns `a + b`.
"""
add_nums(a, b) = a + b 
```

is fine.

### 3. Include appropriate tests 

Tests are crucial, regardless of whether your change is big or small. For small tests, simple tests are fine, but for larger changes that implement, for example, a new feature, understanding what tests you need to write is extremely important.

For example, suppose you wanted to implement a new algorithm for inserting a curve into a triangulation. The tests involved with this change need to test things like:

1. Inserting a curve into an unconstrained triangulation.
2. Inserting a curve into a constrained triangulation.
3. Inserting a curve into a constrained triangulation with multiple boundaries and pre-existing segments.
4. Inserting multiple curves.
5. Points sets with problems such as cocircular/collinear points, collinear segments, collinear points on the inserted curve, etc.

This final point often takes the most time. Understanding how to test for these things may not be straightforward if you are inexperienced, so please feel free to ask for guidance. 

Please also ensure to make your tests easy to understand.

Some practical notes on the test suite:

- The folders in `DelaunayTriangulation/test` mirror those in `DelaunayTriangulation/src`, e.g. the tests for `src/voronoi/` are in `test/voronoi/`. New test files need to be added to `test/runtests.jl`.
- The full test suite takes over 20 minutes. To only run some of its groups, e.g. while working on a change, pass their names as test arguments: `Pkg.test("DelaunayTriangulation"; test_args = ["Voronoi", "Point Location"])`. The group names are listed in `test/runtests.jl`. On CI, the groups are split between parallel jobs in `.github/workflows/CI.yml`, so a new group needs to be added there as well; the "Test summary" job, which puts the results of the jobs together, fails if any test script was not run.
- `test/helper_functions.jl` defines a module of helpers for the tests, e.g. example geometries, validators like `validate_statistics`, comparisons that do not depend on hash order like `same_polygon`, and slow reference implementations. The helpers are split by purpose over the files in `test/helpers/`. To use them outside of the test suite, e.g. in the REPL with the test environment active, run `include(joinpath(pkgdir(DelaunayTriangulation), "test", "helper_functions.jl")); using .HelperFunctions`.
- Pass a seeded RNG to everything random in a test, both to calls like `rand(rng, 2, 50)` and to functions with an `rng` keyword like `triangulate(points; rng)`, so that a failure can be reproduced by running that testset on its own. The convention is an `rng = StableRNG(k)` at the start of each top-level testset (or before a loop at the top level), which loops inside it then draw from. Do not define a new `rng` with a fixed seed inside a loop that repeats a test for different insertion orders, since each iteration would then be the same.
- Reference images for tests that compare figures are in the top-level `assets/` folder rather than in the package, so that they are not downloaded with the package. See `reference_path` in `test/helpers/testing.jl` for how images that differ between Julia versions are handled.

### 4. Try to only do one thing per pull request

Ideally, a single pull request should only do one thing. This makes it easy to track your changes and track regressions between versions. For example, do not move around a bunch of files while at the same time implementing a function somewhere else - make two pull requests.

### 5. Add to NEWS.md 

As stated at the start, please describe your changes in the NEWS.md file.

## Performance Improvements

In addition to the guidance above for making a pull request, if your contribution is aimed at improving performance, benchmarks comparing the results with and without your change are extremely useful.

The `DelaunayTriangulation/benchmark` directory has a benchmark suite, `benchmarks.jl`, covering most of the package: triangulating (unconstrained, weighted, constrained, multiply connected, disjoint and curve-bounded domains), refinement, operations such as `add_point!` and `delete_point!`, point location, Voronoi tessellations and power diagrams, and some geometry utilities. To run it, from the repository's root:

```sh
julia DelaunayTriangulation/benchmark/run.jl --output new.json
```

To compare your changes with another version, e.g. the `main` branch, check that version out somewhere else, e.g. with `git worktree add ../baseline main`, benchmark it with the same script, and compare the two results:

```sh
julia DelaunayTriangulation/benchmark/run.jl --package ../baseline/DelaunayTriangulation --output baseline.json
julia DelaunayTriangulation/benchmark/compare.jl baseline.json new.json
```

The comparison is printed as a Markdown table. You can run only some of the benchmarks with e.g. `--filter voronoi`, and see the comments at the start of `run.jl` and `compare.jl` for their other options. For pull requests that change the source code, the Benchmarks workflow does this comparison automatically, and its results are in the summary of its run in the Actions tab. If you add a feature, please consider adding a benchmark for it to `benchmarks.jl`.

## Feature Contributions

If you want to contribute a new feature then, depending on how substantial the feature is, please make an issue or a discussion thread (in the discussions tab) to discuss your ideas first. If there is already an existing feature request in the discussions tab, comment there with some ideas. This will help make sure we get the design of your feature right and that we can help each other.