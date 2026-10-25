```@meta
CurrentModule = DelaunayTriangulation
```

# [Performance Tips](@id performance-tips)

The default settings of this package are chosen to be robust and reasonably fast for most problems, but there are some things you can do to make your code faster, and some things to avoid. The timings below are only to give an idea of the differences. They were measured on a laptop with Julia 1.13, and will depend on your machine and your problem. As for any Julia code, the general [performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/) in the Julia manual also apply, e.g. you should put loops that work with a triangulation inside of a function rather than at the top level of a script. The [Performance Measurements](performance_measurements.md) page shows how the differences discussed here grow with the number of points, together with the code used to measure them.

## Constructing Triangulations

### Keep the Insertion Order Randomised

By default, [`triangulate`](@ref) inserts the points in a random order that is also spatially sorted, as described in [Insertion Order](@ref delaunay-insertion-order). With `randomise = false`, the points are instead inserted in the order that they are given in. This is rarely faster, and it can be much slower, even when the points look sorted. For example:

- For 100,000 random points, `triangulate` takes 1.1 s by default and 2.2 s with `randomise = false`.
- For a $300 \times 300$ lattice of points listed row by row, `triangulate` takes 0.84 s by default and 21.5 s with `randomise = false`, since inserting a lattice row by row takes $\mathcal O(n^3)$ time for an $n \times n$ lattice.

So, you should only use `randomise = false` if you need the points to be inserted in exactly the order they are given in. To get the same triangulation each time, pass an `rng` instead, e.g. `triangulate(points; rng = StableRNG(123))` with [StableRNGs.jl](https://github.com/JuliaRandom/StableRNGs.jl).

### Use `triangulate_rectangle` for Rectangular Grids

For a regular grid of points in a rectangle, [`triangulate_rectangle`](@ref) is faster than constructing the points yourself and using `triangulate`, since it builds the triangles directly from the grid rather than inserting the points one at a time: for a $300 \times 300$ grid, it takes 0.46 s rather than 0.84 s. It also gives you the boundary of the rectangle as boundary nodes. See [this tutorial](../tutorials/lattice.md).

### Use a Smaller Integer Type

By default, the vertices of a triangulation are `Int`s, i.e. `Int64` on a 64-bit system. If your triangulation has fewer than two billion or so vertices, you can store them as `Int32`s instead with the `IntegerType` keyword argument:

```julia
tri = triangulate(points; IntegerType = Int32)
```

This reduces the memory used by the triangulation by more than a third. For 100,000 random points, the triangulation then uses 85 MiB rather than 136 MiB, and `triangulate` is around 10% faster. [`triangulate_rectangle`](@ref) also has an `IntegerType` keyword argument.

### The Predicate Kernel Matters Less than You Might Expect

The robust predicate kernels are not much slower than `FastKernel()`. For 100,000 random points, `triangulate` takes about 1 s with each of `AdaptiveKernel()` (the default), `ExactKernel()`, and `FastKernel()`, with the times within about 10% of each other. Since `FastKernel()` can give invalid triangulations for points that are collinear or cocircular, or nearly so, it is very rarely worth using. See [Predicate Kernels](predicate_kernels.md) for more details.

### Insert All the Points at Once

If you have all the points from the start, it is faster to triangulate them all at once with `triangulate` than to add them one at a time with [`add_point!`](@ref). For 20,000 random points, `triangulate` takes 0.12 s, while adding them one at a time in a random order takes 0.20 s. Most of the difference is in finding the triangle that contains each new point, which `add_point!` starts from the closest of a few randomly sampled vertices, whereas `triangulate` starts from the last inserted point, which its insertion order keeps close to the next point.

If you do need to add points one at a time, and each new point is close to the previous one, you can similarly start the search from the previous point with the `initial_search_point` keyword argument:

```julia
function add_points!(tri, new_points)
    k = first(each_solid_vertex(tri)) # any vertex to start the first search from
    for p in new_points
        add_point!(tri, p; initial_search_point = k)
        k = DelaunayTriangulation.num_points(tri) # the vertex of the point that was just added
    end
    return tri
end
```

For the same 20,000 points sorted into strips with the `strip_sort` function below, this takes 0.17 s rather than 0.20 s. Only do this if consecutive points are close together, though: for points in a random order, it takes 0.29 s instead.

## Point Location and Nearest Neighbours

[`find_triangle`](@ref) finds the triangle containing a point by walking through the triangulation, starting from the closest of a few randomly sampled vertices (about the cube root of the number of vertices). If you have many points to locate, you can make this much faster by sorting the points so that consecutive points are close together, and then starting each search from a vertex of the triangle found for the previous point using the `k` keyword argument. For 100,000 random points in a triangulation of 100,000 random points, this takes 0.14 s, of which 0.02 s is for the sorting, rather than 1.2 s:

```julia
# Sort the points into horizontal strips, alternating the direction along each strip, so that 
# consecutive points are close together
function strip_sort(points)
    ymin, ymax = extrema(last, points)
    m = round(Int, sqrt(length(points))) # the number of strips
    function key(p)
        x, y = p
        s = min(m - 1, floor(Int, m * (y - ymin) / (ymax - ymin))) # the strip that p is in
        return (s, isodd(s) ? -x : x)
    end
    return sort(points; by = key)
end

# Returns the sorted queries and the triangle containing each of them
function locate_all(tri, queries)
    sorted_queries = strip_sort(queries)
    triangles = Vector{DelaunayTriangulation.triangle_type(tri)}(undef, length(sorted_queries))
    k = first(each_solid_vertex(tri))
    for (n, q) in enumerate(sorted_queries)
        V = find_triangle(tri, q; k)
        triangles[n] = V
        # If q is outside of the triangulation, V is a ghost triangle, so one of its vertices is a ghost vertex,
        # which can not be used as the start of the next search
        i, j, _ = triangle_vertices(V)
        k = DelaunayTriangulation.is_ghost_vertex(i) ? j : i
    end
    return sorted_queries, triangles
end
```

The sorting is what makes this work: without it, starting each search from the previous triangle is about four times slower than the default, since consecutive points are then far apart. The same approach works for [`get_nearest_neighbour`](@ref), which passes its keyword arguments to `find_triangle` and returns a vertex that can be used for the next `k` directly. For the same points, this takes 0.15 s rather than 1.2 s:

```julia
# Returns the sorted queries and the nearest neighbour of each of them
function nearest_neighbours(tri, queries)
    sorted_queries = strip_sort(queries)
    neighbours = Vector{Int}(undef, length(sorted_queries))
    k = first(each_solid_vertex(tri))
    for (n, q) in enumerate(sorted_queries)
        k = get_nearest_neighbour(tri, q; k)
        neighbours[n] = k
    end
    return sorted_queries, neighbours
end
```

## Other Tips

- `copy(tri)` makes a deep copy of a triangulation, so that the copy is independent of `tri`. This takes 0.14 s for a triangulation of 100,000 points, about an eighth of the time it takes to construct it, so it is best avoided inside of loops.
- [`validate_triangulation`](@ref) checks that a triangulation is valid, e.g. that it is Delaunay and that its fields are consistent with each other, and takes more than twice as long as constructing the triangulation does (0.13 s rather than 0.06 s for 10,000 random points). It is useful for testing and for debugging, but you should not need it otherwise.
