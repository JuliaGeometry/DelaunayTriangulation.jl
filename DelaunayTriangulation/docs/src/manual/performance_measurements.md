```@meta
CurrentModule = DelaunayTriangulation
```

# Performance Measurements

This page shows measurements of the effects described in the [Performance Tips](@ref performance-tips), to show how they grow with the number of points.

```@eval
using DelaunayTriangulation, TOML, Markdown
metadata = TOML.parsefile(joinpath(pkgdir(DelaunayTriangulation), "docs", "src", "manual", "performance_measurements.toml"))["metadata"]
Markdown.parse("These measurements were made on $(metadata["date"]) with Julia $(metadata["julia"]) and DelaunayTriangulation.jl v$(metadata["DelaunayTriangulation"]), on $(metadata["os"]) with the processor $(metadata["cpu"]).")
```

The timings are not measured when the documentation is built. Instead, they are measured by the script [`docs/performance_measurements.jl`](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/blob/main/DelaunayTriangulation/docs/performance_measurements.jl), which saves them to `docs/src/manual/performance_measurements.toml` for this page to plot. Each timing is the minimum over as many runs as fit into ten seconds, measured with [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl). To update them, e.g. for a new version of the package, run the script from the `DelaunayTriangulation` directory of the repository:

```sh
julia --project=docs docs/performance_measurements.jl
```

The code below each heading shows what is compared in each figure, and the full code is in the script. The points are generated with

```julia
random_points(n; seed = n) = (rng = StableRNG(seed); [(rand(rng), rand(rng)) for _ in 1:n])
lattice_points(m) = [(i / m, j / m) for j in 0:(m - 1) for i in 0:(m - 1)]
```

```@setup perf
using DelaunayTriangulation
using CairoMakie
using TOML

data = TOML.parsefile(joinpath(pkgdir(DelaunayTriangulation), "docs", "src", "manual", "performance_measurements.toml"))

# A log-log plot of the times in results[key] for each key => label in series, against the numbers of points
function plot_times(results, series; title, ylabel = "Time (s)")
    fig = Figure(size = (650, 420))
    ax = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "Number of points", ylabel, title)
    for (key, label) in series
        scatterlines!(ax, results["n"], results[key]; label)
    end
    axislegend(ax; position = :lt)
    return fig
end

# Plots the times in results[key] for each key => label in series, relative to the times in results[reference]
function plot_relative!(ax, results, reference, series)
    for (i, (key, label)) in enumerate(series)
        scatterlines!(ax, results["n"], results[key] ./ results[reference]; label, color = Makie.wong_colors()[i])
    end
    hlines!(ax, 1; color = :black, linestyle = :dash)
    ylims!(ax, 0, 1.5)
    axislegend(ax; position = :lb)
    return ax
end

function plot_integer_type(results)
    fig = Figure(size = (900, 380))
    ax = Axis(fig[1, 1]; xscale = log10, xlabel = "Number of points", ylabel = "Time relative to IntegerType = Int", title = "Time")
    plot_relative!(ax, results, "time_Int", ["time_Int" => "IntegerType = Int (default)", "time_Int32" => "IntegerType = Int32"])
    ax = Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "Number of points", ylabel = "Size of the triangulation (MiB)", title = "Memory")
    scatterlines!(ax, results["n"], results["memory_Int"]; label = "IntegerType = Int (default)")
    scatterlines!(ax, results["n"], results["memory_Int32"]; label = "IntegerType = Int32")
    axislegend(ax; position = :lt)
    return fig
end

function plot_kernels(results)
    fig = Figure(size = (650, 420))
    ax = Axis(fig[1, 1]; xscale = log10, xlabel = "Number of points", ylabel = "Time relative to AdaptiveKernel()", title = "Predicate kernels")
    plot_relative!(ax, results, "AdaptiveKernel", ["AdaptiveKernel" => "AdaptiveKernel() (default)", "ExactKernel" => "ExactKernel()", "FastKernel" => "FastKernel()"])
    return fig
end
```

## Insertion Order

By default, `triangulate` inserts the points in a random order that is spatially sorted. With `randomise = false`, it inserts them in the order they are given in instead. For random points, this is slower, increasingly so as the number of points grows:

```julia
points = random_points(n)
triangulate(points; rng = StableRNG(1)) # randomise = true, the default
triangulate(points; randomise = false)
```

```@example perf
plot_times(data["insertion_order"], ["randomise_true" => "randomise = true (default)", "randomise_false" => "randomise = false"]; title = "Random points") # hide
```

For the points of an $m \times m$ lattice listed row by row, which might look like a sorted order, `randomise = false` is much worse: the time grows like $m^3$, i.e. like $n^{3/2}$ for $n = m^2$ points, rather than like $n$. For a lattice, [`triangulate_rectangle`](@ref) is faster still, since it builds the triangles directly from the lattice.

```julia
points = lattice_points(m)
triangulate(points; rng = StableRNG(1)) # randomise = true, the default
triangulate(points; randomise = false)
triangulate_rectangle(0.0, 1.0, 0.0, 1.0, m, m)
```

```@example perf
plot_times(data["lattice"], ["randomise_true" => "triangulate, randomise = true (default)", "randomise_false" => "triangulate, randomise = false", "triangulate_rectangle" => "triangulate_rectangle"]; title = "Points of a lattice, listed row by row") # hide
```

## Integer Type

Storing the vertices as `Int32`s rather than `Int64`s reduces the memory used by the triangulation by more than a third, and the time by around 10%:

```julia
tri = triangulate(points; rng = StableRNG(1)) # IntegerType = Int, the default
tri = triangulate(points; IntegerType = Int32, rng = StableRNG(1))
Base.summarysize(tri) # the memory used by tri
```

```@example perf
plot_integer_type(data["integer_type"]) # hide
```

## Predicate Kernels

The robust predicate kernels are not much slower than `FastKernel()`, which is not robust:

```julia
triangulate(points; predicates = AdaptiveKernel(), rng = StableRNG(1)) # the default
triangulate(points; predicates = ExactKernel(), rng = StableRNG(1))
triangulate(points; predicates = FastKernel(), rng = StableRNG(1))
```

```@example perf
plot_kernels(data["kernels"]) # hide
```

## Inserting Points One at a Time

Adding the points one at a time with [`add_point!`](@ref) is slower than triangulating them all at once, mostly because each search for the triangle containing a new point starts from a few randomly sampled vertices. If consecutive points are close together, starting each search from the previous point with `initial_search_point` makes up most of the difference. Here, the points are sorted with the `strip_sort` function from the [Performance Tips](@ref performance-tips), and the time for the sorting is not included. In a random order, starting from the previous point is slower than the default.

```julia
function add_points!(tri, new_points; start_from_previous = false)
    k = first(each_solid_vertex(tri))
    for p in new_points
        if start_from_previous
            add_point!(tri, p; initial_search_point = k)
            k = DelaunayTriangulation.num_points(tri)
        else
            add_point!(tri, p)
        end
    end
    return tri
end

first_points, other_points = points[1:3], points[4:end]
triangulate(points; rng = StableRNG(1))
add_points!(triangulate(first_points), other_points)
add_points!(triangulate(first_points), strip_sort(other_points); start_from_previous = true)
add_points!(triangulate(first_points), other_points; start_from_previous = true)
```

```@example perf
plot_times(data["one_at_a_time"], ["triangulate" => "triangulate", "random_order" => "add_point!, random order", "sorted_previous" => "add_point!, sorted, initial_search_point = previous", "random_previous" => "add_point!, random order, initial_search_point = previous"]; title = "Inserting points one at a time") # hide
```

## Locating Many Points

The same applies to [`find_triangle`](@ref): for many points, sorting them and starting each search from a vertex of the triangle found for the previous point, using the `k` keyword argument, is much faster than the default, while without the sorting it is slower. Here, `n` random points are located in a triangulation of `n` random points, and the time for the sorting is included.

```julia
function locate_all(tri, queries; sorted = true)
    k = first(each_solid_vertex(tri))
    for q in (sorted ? strip_sort(queries) : queries)
        V = find_triangle(tri, q; k)
        i, j, _ = triangle_vertices(V)
        k = DelaunayTriangulation.is_ghost_vertex(i) ? j : i
    end
    return k
end

foreach(q -> find_triangle(tri, q), queries) # the default
locate_all(tri, queries)
locate_all(tri, queries; sorted = false)
```

```@example perf
plot_times(data["point_location"], ["default" => "default", "sorted_previous" => "sorted, k = previous", "random_previous" => "random order, k = previous"]; title = "Locating n random points in a triangulation of n random points") # hide
```

## Copying and Validating

Finally, the costs of `copy` and [`validate_triangulation`](@ref), compared with constructing the triangulation:

```julia
tri = triangulate(points; rng = StableRNG(1))
copy(tri)
DelaunayTriangulation.validate_triangulation(tri; print_result = false)
```

```@example perf
plot_times(data["copy_validate"], ["triangulate" => "triangulate", "copy" => "copy", "validate_triangulation" => "validate_triangulation"]; title = "Copying and validating") # hide
```
