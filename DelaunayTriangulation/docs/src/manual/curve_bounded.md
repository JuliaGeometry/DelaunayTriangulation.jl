```@meta 
CurrentModule = DelaunayTriangulation
```

# Defining Curve-Bounded Domains 

This section discusses how curve-bounded domains, and curves, are defined in this package. A good demonstration of how these domains are worked with is in the [curve-bounded tutorial](../tutorials/curve_bounded.md).

## Curves

To start, let us discuss curves. All functions that work with curves in this package treat them as being subtypes of the [`AbstractParametricCurve`](@ref) type. These curves must:

1. Be defined as parametric curves, parametrised over $0 \leq t \leq 1$.
2. Not be self-intersecting, with the exception of allowing for closed curves.
3. Implement [`differentiate`](@ref), [`twice_differentiate`](@ref), and [`thrice_differentiate`](@ref).
4. Be defined as a callable struct.
5. Either implement [`point_position_relative_to_curve`](@ref) or [`get_closest_point`](@ref). Alternatively, the struct should have a `lookup_table` field so that `lookup_table[i]` is the value of the curve at `t = (i - 1) / (length(lookup_table) - 1)`.

With these specifications, a curve can be fully compatible with the functions in this package, as other functions such as [`arc_length`](@ref), [`curvature`](@ref), and [`total_variation`](@ref) are automatically defined for such a curve. The curves that are defined in this package can be found using `subtypes`:

```@example 
using DelaunayTriangulation, InteractiveUtils
subtypes(DelaunayTriangulation.AbstractParametricCurve)
```

## Curve-Bounded Domains

For defining curve-bounded domains, defining them is similar to defining piecewise linear boundaries as described [here](boundaries.md). You still need to be careful about orientation, and you of course need to make sure that curves belonging to the same boundary connect appropriately. For defining curve-bounded domains that also have a piecewise linear section, the piecewise linear section should be defined as a vector of vertices just like [here](boundaries.md).