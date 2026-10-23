# Compares two results files saved by benchmark/run.jl and prints a Markdown table of the results. Usage:
#
#     julia benchmark/compare.jl BASELINE.json NEW.json [TOLERANCE]
#
# The times are the minimum times over the samples of each benchmark. A ratio above 1 means that NEW is slower than BASELINE.
# Changes in time larger than TOLERANCE, default 0.1 (i.e. 10%), are marked as slower or faster, since timings on shared
# machines, such as those used by GitHub Actions, can easily vary by a few percent.

using Pkg
Pkg.activate(; temp = true, io = devnull)
Pkg.add(name = "BenchmarkTools", version = "1"; io = devnull)
using BenchmarkTools

function compare_benchmarks(baseline_file, new_file, tolerance)
    baseline = BenchmarkTools.load(baseline_file)[1]
    new = BenchmarkTools.load(new_file)[1]
    baseline_trials = Dict(path => trial for (path, trial) in BenchmarkTools.leaves(baseline))
    new_trials = Dict(path => trial for (path, trial) in BenchmarkTools.leaves(new))
    names = sort!(collect(union(keys(baseline_trials), keys(new_trials))))
    num_slower = num_faster = 0
    rows = String[]
    for name in names
        old_trial, new_trial = get(baseline_trials, name, nothing), get(new_trials, name, nothing)
        old_time = isnothing(old_trial) ? "–" : BenchmarkTools.prettytime(minimum(old_trial).time)
        new_time = isnothing(new_trial) ? "–" : BenchmarkTools.prettytime(minimum(new_trial).time)
        old_memory = isnothing(old_trial) ? "–" : BenchmarkTools.prettymemory(memory(old_trial))
        new_memory = isnothing(new_trial) ? "–" : BenchmarkTools.prettymemory(memory(new_trial))
        ratio, change = "–", ""
        if !isnothing(old_trial) && !isnothing(new_trial)
            r = minimum(new_trial).time / minimum(old_trial).time
            ratio = string(round(r; digits = 2))
            if r > 1 + tolerance
                change, num_slower = "slower", num_slower + 1
            elseif r < 1 - tolerance
                change, num_faster = "faster", num_faster + 1
            end
        end
        push!(rows, "| $(join(name, " / ")) | $old_time | $new_time | $ratio | $change | $old_memory | $new_memory |")
    end
    println("## Benchmarks\n")
    println("$num_slower slower and $num_faster faster, out of $(length(names)) benchmarks (minimum times, with a tolerance of $(round(Int, 100tolerance))%).\n")
    println("| Benchmark | Baseline | New | Ratio | | Baseline memory | New memory |")
    println("|---|---|---|---|---|---|---|")
    foreach(println, rows)
    return nothing
end

length(ARGS) ∈ (2, 3) || error("Usage: julia benchmark/compare.jl BASELINE.json NEW.json [TOLERANCE]")
compare_benchmarks(ARGS[1], ARGS[2], length(ARGS) == 3 ? parse(Float64, ARGS[3]) : 0.1)
