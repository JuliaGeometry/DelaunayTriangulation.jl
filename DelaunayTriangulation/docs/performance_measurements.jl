# Measures the timings shown on the Performance Measurements page of the documentation (src/manual/performance_measurements.md)
# and saves them to src/manual/performance_measurements.toml, which that page plots. They are not measured when the documentation
# is built, since that would make the build much slower and the shared machines that build it give noisy timings. To update them,
# e.g. for a new version of the package, run this script from the DelaunayTriangulation directory with
#     julia --project=docs docs/performance_measurements.jl
# It takes about half an hour. Close anything else that is running, so that it does not affect the timings.

using DelaunayTriangulation
using BenchmarkTools
using StableRNGs
using Dates
using TOML

# The minimum time, in seconds, of f(setup()) over as many runs as fit in `seconds` seconds (at least one). The setup is not
# timed, so it can e.g. create a triangulation for f to modify
function measure(f, setup = () -> nothing; seconds = 10)
    GC.gc() # so that the garbage from the previous measurements does not affect this one
    trial = @benchmark $f(x) setup = (x = $setup()) evals = 1 seconds = seconds
    return minimum(trial).time / 1.0e9
end

random_points(n; seed = n) = (rng = StableRNG(seed); [(rand(rng), rand(rng)) for _ in 1:n])
lattice_points(m) = [(i / m, j / m) for j in 0:(m - 1) for i in 0:(m - 1)]

# Sorts the points into horizontal strips, alternating the direction along each strip, so that consecutive points are close
# together. This is the same as in the Performance Tips
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

function add_points!(tri, new_points; start_from_previous = false)
    rng = StableRNG(1)
    k = first(each_solid_vertex(tri))
    for p in new_points
        if start_from_previous
            add_point!(tri, p; initial_search_point = k, rng)
            k = DelaunayTriangulation.num_points(tri)
        else
            add_point!(tri, p; rng)
        end
    end
    return tri
end

function locate_all(tri, queries; sorted = true)
    k = first(each_solid_vertex(tri))
    for q in (sorted ? strip_sort(queries) : queries)
        V = find_triangle(tri, q; k)
        i, j, _ = triangle_vertices(V)
        k = DelaunayTriangulation.is_ghost_vertex(i) ? j : i
    end
    return k
end

function locate_all_default(tri, queries)
    rng = StableRNG(1)
    for q in queries
        find_triangle(tri, q; rng)
    end
    return nothing
end

results = Dict{String, Any}()
results["metadata"] = Dict(
    "date" => string(today()),
    "julia" => string(VERSION),
    "DelaunayTriangulation" => string(pkgversion(DelaunayTriangulation)),
    "cpu" => String(strip(Sys.cpu_info()[1].model)),
    "os" => Sys.iswindows() ? "Windows" : Sys.isapple() ? "macOS" : Sys.islinux() ? "Linux" : string(Sys.KERNEL),
)

println("Insertion order: random points")
ns = [1_000, 3_000, 10_000, 30_000, 100_000]
r = results["insertion_order"] = Dict("n" => ns, "randomise_true" => Float64[], "randomise_false" => Float64[])
for n in ns
    points = random_points(n)
    push!(r["randomise_true"], measure(_ -> triangulate(points; rng = StableRNG(1))))
    push!(r["randomise_false"], measure(_ -> triangulate(points; randomise = false)))
end

println("Insertion order: lattices")
ms = [20, 40, 80, 120, 160]
r = results["lattice"] = Dict("n" => ms .^ 2, "randomise_true" => Float64[], "randomise_false" => Float64[], "triangulate_rectangle" => Float64[])
for m in ms
    points = lattice_points(m)
    push!(r["randomise_true"], measure(_ -> triangulate(points; rng = StableRNG(1))))
    push!(r["randomise_false"], measure(_ -> triangulate(points; randomise = false)))
    push!(r["triangulate_rectangle"], measure(_ -> triangulate_rectangle(0.0, 1.0, 0.0, 1.0, m, m)))
end

println("Integer type")
ns = [1_000, 3_000, 10_000, 30_000, 100_000]
r = results["integer_type"] = Dict("n" => ns, "time_Int" => Float64[], "time_Int32" => Float64[], "memory_Int" => Float64[], "memory_Int32" => Float64[])
for n in ns
    points = random_points(n)
    push!(r["time_Int"], measure(_ -> triangulate(points; rng = StableRNG(1))))
    push!(r["time_Int32"], measure(_ -> triangulate(points; IntegerType = Int32, rng = StableRNG(1))))
    push!(r["memory_Int"], Base.summarysize(triangulate(points; rng = StableRNG(1))) / 2^20)
    push!(r["memory_Int32"], Base.summarysize(triangulate(points; IntegerType = Int32, rng = StableRNG(1))) / 2^20)
end

println("Predicate kernels")
ns = [1_000, 3_000, 10_000, 30_000, 100_000]
r = results["kernels"] = Dict("n" => ns, "AdaptiveKernel" => Float64[], "ExactKernel" => Float64[], "FastKernel" => Float64[])
for n in ns
    points = random_points(n)
    for kernel in (AdaptiveKernel(), ExactKernel(), FastKernel())
        push!(r[string(nameof(typeof(kernel)))], measure(_ -> triangulate(points; predicates = kernel, rng = StableRNG(1))))
    end
end

println("Inserting points one at a time")
ns = [1_000, 3_000, 10_000, 30_000]
r = results["one_at_a_time"] = Dict("n" => ns, "triangulate" => Float64[], "random_order" => Float64[], "sorted_previous" => Float64[], "random_previous" => Float64[])
for n in ns
    points = random_points(n)
    first_points, other_points = points[1:3], points[4:end]
    sorted_points = strip_sort(other_points)
    push!(r["triangulate"], measure(_ -> triangulate(points; rng = StableRNG(1))))
    push!(r["random_order"], measure(tri -> add_points!(tri, other_points), () -> triangulate(copy(first_points))))
    push!(r["sorted_previous"], measure(tri -> add_points!(tri, sorted_points; start_from_previous = true), () -> triangulate(copy(first_points))))
    push!(r["random_previous"], measure(tri -> add_points!(tri, other_points; start_from_previous = true), () -> triangulate(copy(first_points))))
end

println("Locating many points")
ns = [1_000, 3_000, 10_000, 30_000, 100_000]
r = results["point_location"] = Dict("n" => ns, "default" => Float64[], "sorted_previous" => Float64[], "random_previous" => Float64[])
for n in ns
    tri = triangulate(random_points(n); rng = StableRNG(1))
    queries = random_points(n; seed = n + 1)
    push!(r["default"], measure(_ -> locate_all_default(tri, queries)))
    push!(r["sorted_previous"], measure(_ -> locate_all(tri, queries)))
    push!(r["random_previous"], measure(_ -> locate_all(tri, queries; sorted = false)))
end

println("Copying and validating")
ns = [1_000, 3_000, 10_000, 30_000, 100_000]
r = results["copy_validate"] = Dict("n" => ns, "triangulate" => Float64[], "copy" => Float64[], "validate_triangulation" => Float64[])
for n in ns
    points = random_points(n)
    tri = triangulate(points; rng = StableRNG(1))
    push!(r["triangulate"], measure(_ -> triangulate(points; rng = StableRNG(1))))
    push!(r["copy"], measure(_ -> copy(tri)))
    push!(r["validate_triangulation"], measure(_ -> DelaunayTriangulation.validate_triangulation(tri; print_result = false)))
end

path = joinpath(@__DIR__, "src", "manual", "performance_measurements.toml")
open(path, "w") do io
    println(io, "# Written by docs/performance_measurements.jl, which describes how these were measured. Times are in seconds, and memory in MiB.")
    TOML.print(io, results; sorted = true)
end
println("Saved the results to ", path)
