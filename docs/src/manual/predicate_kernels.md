```@meta 
CurrentModule = DelaunayTriangulation
```

# Predicate Kernels

By default, this package uses adaptive arithmetic via [AdaptivePredicates.jl](https://github.com/JuliaGeometry/AdaptivePredicates.jl) for computing predicates.
In total, there are three different kernels offered for computing predicates:
- `Fast()`: Predicates will be computed without any adaptive or exact arithmetic. 
- `Adaptive()`: Predicates will be computed using adaptive arithmetic via [AdaptivePredicates.jl](https://github.com/JuliaGeometry/AdaptivePredicates.jl).
- `Exact()`: Predicates will be computed using exact arithmetic via [ExactPredicates.jl](https://github.com/lairez/ExactPredicates.jl). 

There are clear strengths and weaknesses to each of these choices. To summarise them, here is when each kernel should be considered:
- `Fast()`: This kernel can be safely used when you know that there are no collinear points or cocircular points in your data set.
It may still work even in those cases, but it can not be safely relied upon. If you trust that there are no issues, this should 
be the kernel you use as it is the fastest. If you run into issues while using this kernel, please use try `Adaptive()`.
- `Adaptive()`: This is the kernel we use by default. It has performance that is reasonably close to what is offered by `Fast()`,
except it also guarantees that predicates will return the correct result even with collinear points or cocircular points, or in other 
degenerate cases where one typically expects predicates to be problematic. If you are using `Fast()` and run into issues, this should be the next 
kernel you try.
- `Exact()`: This is the slowest kernel, but it is the safest. This kernel works on a much wider range of numbers than `Adaptive()`, and is guaranteed
to satisfy certain combinatorial properties such as `orient(a, b, c) == orient(b, c, a) == orient(c, a, b)`. I have not seen any examples 
where `Adaptive()` fails but `Exact()` works, though, so you should only consider using this kernel if you do actually encounter such a case, i.e.
treat this kernel as a fallback for `Adaptive()`.

We give a discussion below about why robust arithmetic is actually important, to help you understand these choices. A key point is that it is highly advised that you do not use `Fast()`.

## Why use robust predicates?

Three great resources for understanding why we need robust predicates are

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

but this is not true in floating point arithmetic. This causes issues with consistency - a point may be found to be both left and right of an edge depending on the order of the points given to the `orient` predicate, inevitably leading to an invalid triangulation. With the use of robust predicates, this property is guaranteed to hold, ensuring that all the predicate results are consistent with each other. (This identity does not always hold with adaptive arithmetic,
although this is less problematic due to the design of the adaptive predicates; read Shewchuk's paper for more information.) This has the following consequence: **Even if you think robust predicates are not necessary for you because none of your inputs are exact (for example), you still want them to guarantee consistency with predicates regardless of the input order**.

## Will disabling exact predicates give me better performance?

It is also not even the case that using inexact predicates will give you better performance than if you were to use robust predicates. `Adaptive()`'s performance is typically similar to `Fast()`, with the exception of queries on collinear points. This exception is irrelevant, though, as `Fast()` is not even 
reliable when used on collinear points. `Exact()` is a bit slower, but its performance is still not terrible compared to `Fast()` since ExactPredicates.jl
uses clever filters that typically do as much work as `Fast()` or `Adaptive()` would. Thus, the only cases where performance is improved significantly using `Fast()` is exactly in the cases where you do not want to be using `Fast()`.

You should always benchmark your problems to see if using `Fast()` over the robust kernels `Adaptive()` or `Exact()`, if you choose to do, will actually give you better performance.

## Can I check if my computed triangulation is valid?

When you are not using robust predicates, you may want to check if your computed triangulation is actually a valid Delaunay triangulation. We provide the function `DelaunayTriangulation.validate_triangulation` for this purpose. This functionality is quite slow to use and is not currently optimised or well-documented (contributions towards addressing these issues are welcome), but it will work. One important note is that this check does actually use predicates in certain areas, so this check is still not guaranteed to be 100% accurate without robust predicates (by default, `validate_triangulation` will use the `Exact()` kernel). Here is an example of its use.

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
