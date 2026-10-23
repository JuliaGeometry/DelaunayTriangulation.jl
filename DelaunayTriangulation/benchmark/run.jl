# Runs the benchmarks in benchmarks.jl for a version of DelaunayTriangulation.jl and saves the results. Usage:
#
#     julia benchmark/run.jl [--package PATH] [--output FILE] [--filter TEXT] [--seconds SECONDS]
#
# - `--package`: the directory of the version of the package to benchmark, i.e. the directory containing its Project.toml. By
#   default this is the package that contains this script. To compare against another version, e.g. the main branch, check
#   that version out somewhere else and benchmark it with this script, so that both versions run the same benchmarks.
# - `--output`: the JSON file to save the results to, default `benchmark_results.json`.
# - `--filter`: only run the benchmarks whose name, including its groups, contains this text, e.g. `--filter voronoi`.
# - `--seconds`: the time budget for each benchmark, default `2`. Each benchmark is run at least twice, once to compile it.
#
# The benchmarks are run in a temporary environment made from benchmark/Project.toml, with DelaunayTriangulation
# taken from `--package`. Use benchmark/compare.jl to compare two results files.

using Pkg, TOML

function parse_arguments(args)
    options = Dict("package" => dirname(@__DIR__), "output" => "benchmark_results.json", "filter" => "", "seconds" => "2")
    for i in 1:2:length(args)
        key = args[i]
        (startswith(key, "--") && haskey(options, key[3:end]) && i < length(args)) || error("Invalid argument $key. See the comments at the start of $(@__FILE__).")
        options[key[3:end]] = args[i + 1]
    end
    return options
end
const OPTIONS = parse_arguments(ARGS)

# A temporary copy of the benchmark environment, with DelaunayTriangulation taken from --package
let dir = mktempdir()
    project = TOML.parsefile(joinpath(@__DIR__, "Project.toml"))
    project["sources"] = Dict("DelaunayTriangulation" => Dict("path" => abspath(OPTIONS["package"])))
    open(io -> TOML.print(io, project), joinpath(dir, "Project.toml"), "w")
    Pkg.activate(dir)
    Pkg.instantiate()
end

using BenchmarkTools
include(joinpath(@__DIR__, "benchmarks.jl"))

# Adds trial to results under the groups and name given by path, creating the groups as needed
function add_result!(results, path, trial)
    group = results
    for key in path[1:(end - 1)]
        group = haskey(group, key) ? group[key] : (group[key] = BenchmarkGroup())
    end
    group[path[end]] = trial
    return results
end

function run_benchmarks(suite, options)
    results = BenchmarkGroup()
    seconds = parse(Float64, options["seconds"])
    @info "Benchmarking DelaunayTriangulation v$(pkgversion(DelaunayTriangulation)) from $(pkgdir(DelaunayTriangulation))"
    for (path, benchmark) in BenchmarkTools.leaves(suite)
        name = join(path, " / ")
        occursin(options["filter"], name) || continue
        print(name, ": ")
        try
            warmup(benchmark; verbose = false)
            trial = run(benchmark; seconds)
            add_result!(results, path, trial)
            println(BenchmarkTools.prettytime(minimum(trial).time), ", ", BenchmarkTools.prettymemory(memory(trial)), " (", length(trial), " samples)")
            flush(stdout)
        catch err
            # e.g. a benchmark that uses something that the benchmarked version of the package does not have
            println("failed: ", sprint(showerror, err))
            flush(stdout)
        end
    end
    BenchmarkTools.save(options["output"], results)
    @info "Saved the results to $(abspath(options["output"]))"
    return results
end

run_benchmarks(SUITE, OPTIONS)
