```@meta 
CurrentModule = DelaunayTriangulation
```

# Disabling Exact Predicates

For performance reasons, you may find it useful to want to disable exact predicates using [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl). This can be easily done using a setup with [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl), but before you consider disabling exact predicates, there are a few things to be aware of. If you just want to disable them without reading a lot of information warning you about the consequences, please skip to the end.

## Why use exact predicates?

Three great resources for understanding why we need exact predicates are

1. Jonathan Shewchuk's paper on adaptive precision floating-point arithmetic [here](https://doi.org/10.1007/PL00009321).
2. Jonathan Shewchuk's lecture notes on geometric robustness [here](https://people.eecs.berkeley.edu/~jrs/meshpapers/robnotes.pdf).
3. [This paper](https://doi.org/10.1016/j.comgeo.2007.06.003) by Kettner et al. (2008) on some examples of issues with inexact arithmetic.

We give a simple summary here. A big component of the algorithms used in this package are what are known as _geometric predicates_, some of these being:

- `orient(p, q, r)`: Is `r` left, right, or on the line through `pq`?
- `incircle(p, q, r, s)`: Is `s` inside, outside, or on the circle through `p`, `q`, and `r`?

These predicates can be computed using determinants:

```math 
\begin{align*}
O_{pqr} &:= \textrm{orient}(p, q, r) = \begin{vmatrix} p_x - r_x & p_y - r_y \\ q_x - r_x & q_y - r_y \end{vmatrix}, \\
C_{pqrs} &:= \textrm{incircle}(p, q, r, s) = \begin{vmatrix} p_x - s_x & p_y - s_y & (p_x - s_x)^2 + (p_y - s_y)^2 \\ q_x - s_x & q_y - s_y & (q_x - s_x)^2 + (q_y - s_y)^2 \\ r_x - s_x & r_y - s_y & (r_x - s_x)^2 + (r_y - s_y)^2 \end{vmatrix}.
\end{align*}
```

The signs of these determinants $O_{pqr}$ and $C_{pqrs}$ are used to determinant the answers to the above questions. In inexact arithmetic, it is common that the sign picked is wrong when the determinants are close to zero. The consequences of this can be catastrophic:

1. The algorithms may hang or crash.
2. The final triangulation may be completely invalid. For example, if a point is being added into a triangulation right onto an existing edge, then in exact arithmetic we would know to split the edge to the left and to the right. In inexact arithmetic, the point may be to the left of the edge but detected as being to the right of it, thus adding a triangle that crosses an edge.
3. You may encounter `BoundsError`s from bad `Adjacent` queries where a triangle is expected to exist but doesn't.

Another issue is due to the fact that floating point arithmetic is not associative. In exact arithmetic, we would expect for example that 

```math 
O_{pqr} = O_{qrp} = O_{rpq},
```

but this is not true in floating point arithmetic. This causes issues with consistency - a point may be found to be both left and right of an edge depending on the order of the points given to the `orient` predicate, inevitably leading to an invalid triangulation. With the use of exact predicates, this property is guaranteed to hold, ensuring that all the predicate results are consistent with each other. This has the following consequence: **Even if you think exact predicates are not necessary for you because none of your inputs are exact (for example), you still want them to guarantee consistency with predicates regardless of the input order**.

## Will disabling exact predicates give me better performance?

It is also not even the case that using inexact predicates will give you better performance than if you were to us exact predicates. ExactPredicates.jl uses clever filters to ensure that the fast path, meaning an inexact computation, is performed first. This means that, in most cases, computing the `orient` predicate would have the same speed as if you were to compute it using the determinant definition. 

One case where you could see improvements would be if there were many collinear points, where the slow computations from exact predicates would actually be needed - but this is exactly the case where you expect inexact predicates to cause you problems! 

You should always benchmark your problems to see if disabling exact predicates, if you choose to do, will actually give you better performance.

## How do I disable exact predicates?

If you still want to disable exact predicates, here is how you can do so. Before doing `using DelaunayTriangulation`, you can use 

```julia
using Preferences: set_preferences! 
set_preferences!("DelaunayTriangulation", "USE_EXACTPREDICATES" => false)
using DelaunayTriangulation # only load after setting the preference 
```

The `set_preferences!` call will make a `LocalPreferences.toml` file in your directory that sets this preference. Once this file exists you are free to delete `Preferences.jl`. If you want to skip using Preferences.jl entirely, you can also just create `LocalPreferences.toml` in your working directory manually and put 

```
[DelaunayTriangulation]
USE_EXACTPREDICATES = false
```

into it. If you later want to re-enable exact predicates, either delete the file or write `USE_EXACTPREDICATES = true` instead. This setup ensures that there is no slowdown in the package from checking if ExactPredicates.jl is being used at runtime, as it is all done at compile time instead.

## Can I check if my computed triangulation is valid?

When you are not using exact predicates, you may want to check if your computed triangulation is actually a valid Delaunay triangulation. We provide the function `DelaunayTriangulation.validate_triangulation` for this purpose. This functionality is quite slow to use and is not currently optimised or well-documented (contributions towards addressing these issues are welcome), but it will work. One important note is that this check does actually use predicates in certain areas, so this check is still not guaranteed to be 100% accurate. Here is an example of its use.

```julia
using DelaunayTriangulation
tri = triangulate(rand(2, 50))
DelaunayTriangulation.validate_triangulation(tri)
```
```julia 
true
```
```julia
T = first(each_solid_triangle(tri)) 
DelaunayTriangulation.delete_triangle!(tri, T) # break the triangulation for this example
DelaunayTriangulation.validate_triangulation(tri)
```
```julia
The edge (12, 40) does not have two incident triangles.
The edge (12, 40) appears as an edge in the graph but it and its reverse are not both a key of the adjacent map.

false
```
```julia
DelaunayTriangulation.validate_triangulation(tri; print_result = false)
```
```julia 
false
```
