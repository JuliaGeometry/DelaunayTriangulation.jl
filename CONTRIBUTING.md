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

### 4. Try to only do one thing per pull request

Ideally, a single pull request should only do one thing. This makes it easy to track your changes and track regressions between versions. For example, do not move around a bunch of files while at the same time implementing a function somewhere else - make two pull requests.

### 5. Add to NEWS.md 

As stated at the start, please describe your changes in the NEWS.md file.

## Performance Improvements

In addition to the guidance above for making a pull request, if your contribution is aimed at improving performance, benchmarks comparing the results with and without your change are extremely useful.

## Feature Contributions

If you want to contribute a new feature then, depending on how substantial the feature is, please make an issue or a discussion thread (in the discussions tab) to discuss your ideas first. If there is already an existing feature request in the discussions tab, comment there with some ideas. This will help make sure we get the design of your feature right and that we can help each other.