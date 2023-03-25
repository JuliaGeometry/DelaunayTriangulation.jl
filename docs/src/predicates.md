```@meta
CurrentModule = DelaunayTriangulation
```

# Predicates 

The predicates that we use in this package are all built from ExactPredicates.jl, avoiding degeneracies from predicates owing to floating point arithmetic. The results from predicates are based on certificates, coming from a `Certificate` type defined with EnumX.jl. The definition of this is below.

```@docs 
Certificate
```

## General 

Below we list some general predicates. The core ones that all other predicates are based on are:

```@docs 
orient_predicate
incircle_predicate 
```

The mathematical definitions for these predicates are:

- `orient_predicate`: Let $O(p, q, r)$ denote `orient_predicate(p, q, r)`. The definition is 

```math
O(p, q, r) = \text{sgn}\left(\begin{vmatrix} p_x - r_x & p_y - r_y \\ q_x - r_x & q_y - r_y \end{vmatrix}\right).
```

With this definition, $O(p, q, r) > 0$ means $r$ is left of the line $\overrightarrow{pq}$ or that the triangle $pqr$ is positively oriented; $O(p, q, r) = 0$ means $r$ is collinear with $\overrightarrow{pq}$ or that the triangle $pqr$ is degenerate; $O(p, q, r) < 0$ means $r$ is to the right of $\overrightarrow{pq}$ or that the triangle $pqr$ is negatively oriented.

- `incircle_predicate`: Let $O(p, q, r, s)$ denote `incircle_predicate(p, q, r, s)`. The definition is

```math
O(p, q, r, s) = \text{sgn}\left(\begin{vmatrix} 
p_x - s_x & p_y - s_y & (p_x - s_x)^2 + (p_y - s_y)^2 \\
q_x - s_x & q_y - s_y & (q_x - s_x)^2 + (q_y - s_y)^2 \\
r_x - s_x & r_y - s_y & (r_x - s_x)^2 + (r_y - s_y)^2
\end{vmatrix}\right).
```

With this definition, $O(p, q, r, s) > 0$ means $s$ is inside the circle through $p$, $q$, and $r$; $O(p,q r)=0$ means $s$ is cocircular with $p$, $q$, and $r$; $O(p, q, r) < 0$ means $s$ is outside the circle through $p$, $q$, and $r$.

In code, these two predicates could be defined by (the actual definition with ExactPredicates.jl is much more involved):

```julia
_det(a, b, c, d) = a * d - b * c
_det(a, b, c, d, e, f, g, h, i) = a * _det(e, f, h, i) - d * _det(b, c, h, i) + g * _det(b, c, e, f) # cofactor expansion 
function orient_predicate(a, b, c)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    o = _det(ax - cx, ay - cy, bx - cx, by - cy)
    return Int(sign(o)) # need Int for xor
end
function incircle_predicate(a, b, c, d)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    dx, dy = getxy(d)
    o = _det(ax - dx, ay - dy, (ax - dx)^2 + (ay - dy)^2,
        bx - dx, by - dy, (bx - dx)^2 + (by - dy)^2,
        cx - dx, cy - dy, (cx - dx)^2 + (cy - dy)^2)
    return Int(sign(o)) # need Int for xor
end
```

You could use this as a reference if you want to disconnect from using ExactPredicates.jl (or e.g. use the predicates also defined in GeometricalPredicates.jl). This could be useful if you are not too worried about robustness (although you should typically care about this, so be careful) and just want fast code. Let's see what happens if we randomly triangulate some set of $100,000$ points using ExactPredicates.jl versus the definitions above.

```julia
using DelaunayTriangulation
using BenchmarkTools 
n = 100_000
pts = 20randn(2, n)

## Benchmark the original definition 
b1 = @benchmark triangulate($pts)

## Now change the definitions 
_det(a, b, c, d) = a * d - b * c
_det(a, b, c, d, e, f, g, h, i) = a * _det(e, f, h, i) - d * _det(b, c, h, i) + g * _det(b, c, e, f) # cofactor expansion 
function DelaunayTriangulation.orient_predicate(a, b, c)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    o = _det(ax - cx, ay - cy, bx - cx, by - cy)
    return Int(sign(o)) 
end
function DelaunayTriangulation.incircle_predicate(a, b, c, d)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    dx, dy = getxy(d)
    o = _det(ax - dx, ay - dy, (ax - dx)^2 + (ay - dy)^2,
        bx - dx, by - dy, (bx - dx)^2 + (by - dy)^2,
        cx - dx, cy - dy, (cx - dx)^2 + (cy - dy)^2)
    return Int(sign(o))
end

## Benchmark these new definitions 
b2 = @benchmark triangulate($pts)
```
```julia-repl
julia> b1
BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  2.218 s …   2.386 s  ┊ GC (min … max): 5.61% … 7.91%
 Time  (median):     2.360 s              ┊ GC (median):    5.28%
 Time  (mean ± σ):   2.321 s ± 90.182 ms  ┊ GC (mean ± σ):  5.81% ± 2.02%

  █                                               █       █
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█ ▁
  2.22 s         Histogram: frequency by time        2.39 s <

 Memory estimate: 418.00 MiB, allocs estimate: 4079091.

julia> b2
BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  2.343 s …   2.447 s  ┊ GC (min … max): 3.93% … 7.49%
 Time  (median):     2.370 s              ┊ GC (median):    6.31%
 Time  (mean ± σ):   2.387 s ± 53.921 ms  ┊ GC (mean ± σ):  5.93% ± 1.82%

  █             █                                         █
  █▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  2.34 s         Histogram: frequency by time        2.45 s <

 Memory estimate: 417.99 MiB, allocs estimate: 4079731.
```

Not much difference - ExactPredicates.jl probably never has to run the slow definition in this case. What if the numbers are all very small?

```julia
using ExactPredicates
## Go back to the original definitions 
DelaunayTriangulation.orient_predicate(p, q, r) = orient(getxy(p), getxy(q), getxy(r))
DelaunayTriangulation.incircle_predicate(a, b, c, p) = incircle(getxy(a), getxy(b), getxy(c), getxy(p))

## Get another set of points 
pts = 1e-8rand(2, 100_000)

## Do the benchmarks again 
b1 = @benchmark triangulate($pts)
function DelaunayTriangulation.orient_predicate(a, b, c)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    o = _det(ax - cx, ay - cy, bx - cx, by - cy)
    return Int(sign(o)) 
end
function DelaunayTriangulation.incircle_predicate(a, b, c, d)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    cx, cy = getxy(c)
    dx, dy = getxy(d)
    o = _det(ax - dx, ay - dy, (ax - dx)^2 + (ay - dy)^2,
        bx - dx, by - dy, (bx - dx)^2 + (by - dy)^2,
        cx - dx, cy - dy, (cx - dx)^2 + (cy - dy)^2)
    return Int(sign(o))
end
b2 = @benchmark triangulate($pts)
```
```julia-repl 
julia> b1
BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  2.301 s …   2.423 s  ┊ GC (min … max): 2.34% … 3.96%
 Time  (median):     2.305 s              ┊ GC (median):    3.64%
 Time  (mean ± σ):   2.343 s ± 69.581 ms  ┊ GC (mean ± σ):  3.32% ± 0.85%

  ██                                                      █
  ██▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  2.3 s          Histogram: frequency by time        2.42 s <

 Memory estimate: 430.78 MiB, allocs estimate: 4079516.

julia> b2
BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  2.218 s …   2.269 s  ┊ GC (min … max): 4.03% … 2.51%
 Time  (median):     2.239 s              ┊ GC (median):    4.00%
 Time  (mean ± σ):   2.242 s ± 25.344 ms  ┊ GC (mean ± σ):  3.75% ± 1.14%

  █                      █                                █
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  2.22 s         Histogram: frequency by time        2.27 s <

 Memory estimate: 417.98 MiB, allocs estimate: 4079100.
```
Still not much of a difference, so I would not really recommend bothering changing these definitions -- but the option is there if your application calls for it. (With a million points in the above example, the changed definition is about two seconds faster at a total of 38 seconds).

The other predicates are:

```@docs
sameside_predicate 
meet_predicate
triangle_orientation(::Any, ::Any, ::Any)
point_position_relative_to_circle 
point_position_relative_to_line(::Any, ::Any, ::Any) 
point_position_on_line_segment(::Any, ::Any, ::Any) 
line_segment_intersection_type(::Any, ::Any, ::Any, ::Any) 
point_position_relative_to_triangle(::Any, ::Any, ::Any, ::Any) 
point_position_relative_to_oriented_outer_halfplane
is_legal(::Any, ::Any, ::Any, ::Any)
```

## Boundaries and Ghosts 

Below we list some predicates for working with boundaries and ghost triangles. 

```@docs 
is_boundary_index 
is_boundary_edge(::Any, ::Adjacent) 
is_boundary_triangle(::Any, ::Any, ::Any, ::Any) 
is_ghost_edge 
is_ghost_triangle 
is_interior_curve
is_outer_boundary_index(::Any, ::Any) 
is_outer_ghost_triangle 
is_outer_ghost_edge
is_outer_boundary_node(::Any, ::Graph{I}, ::Any) where {I} 
edge_exists(::I) where {I}
edge_exists(::Any, ::Adjacent{I,E}) where {I,E}
has_ghost_triangles(::Adjacent{I,E}, ::Any) where {I,E} 
```

## Index and Ghost Handling

Below we list methods for working with predicates that are used when we provide indices for points rather than points directly.

```@docs 
triangle_orientation(::Any, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_circumcircle(::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_line(::Any, ::Any, ::Any, ::Any, ::Any)
point_position_on_line_segment(::Any, ::Any, ::Any, ::Any)
line_segment_intersection_type(::Any, ::Any, ::Any, ::Any, ::Any)
point_position_relative_to_triangle(::Any, ::Any, ::Any, ::Any, ::Any, ::AbstractDict)
point_position_relative_to_triangle(::Any, ::Any, ::Any, ::AbstractDict)
```