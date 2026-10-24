# Puts together the test reports written by the CI test jobs (see write_test_report in DelaunayTriangulation/test/runtests.jl).
# Each job only runs some of the test groups, so this is also where we check that every test script was run on each Julia version
# and operating system.
#
# Usage:
#   julia .github/test_summary.jl test_reports
# where test_reports has a directory test-report-<version>-<os>-<job index> for each downloaded artifact, containing its test_report.tsv.
# The summary is written to the job summary page on CI (and printed). The script fails if any test failed or errored, if any
# test script was not run, or if TESTS_RESULT (the result of the test jobs) is given and is not "success".

const REPORT_DIR = only(ARGS)

# The results are kept for each configuration, which is the version in the matrix, e.g. "lts", followed by the operating system
# if it is not Linux, e.g. "1 on windows"
const JULIA_VERSIONS = Dict{String, String}() # the Julia version that each configuration resolved to
const GROUPS = String[]
const GROUP_RESULTS = Dict{Tuple{String, String}, Vector{String}}() # (version, group) => [passes, fails, errors, broken, duration]
const FAILED_TESTSETS = Tuple{String, String, String, Int, Int}[] # (version, group, testset, fails, errors)
const INCLUDED_SCRIPTS = Dict{String, Set{String}}()
const EXPECTED_SCRIPTS = Set{String}()

for dir in readdir(REPORT_DIR)
    m = match(r"^test-report-(.+)-([a-z]+)-latest-\d+$", dir)
    report = joinpath(REPORT_DIR, dir, "test_report.tsv")
    (isnothing(m) || !isfile(report)) && continue
    version = m[2] == "ubuntu" ? m[1] : "$(m[1]) on $(m[2])"
    included = get!(INCLUDED_SCRIPTS, version, Set{String}())
    for line in eachline(report)
        fields = split(line, '\t')
        kind = fields[1]
        if kind == "julia"
            JULIA_VERSIONS[version] = fields[2]
        elseif kind == "groups"
            isempty(GROUPS) && append!(GROUPS, fields[2:end])
        elseif kind == "group"
            GROUP_RESULTS[(version, fields[2])] = fields[3:end]
        elseif kind == "file"
            fails, errors = parse(Int, fields[5]), parse(Int, fields[6])
            fails + errors > 0 && push!(FAILED_TESTSETS, (version, fields[2], fields[3], fails, errors))
        elseif kind == "included"
            push!(included, fields[2])
        elseif kind == "expected"
            push!(EXPECTED_SCRIPTS, fields[2])
        end
    end
end

function count_summary(fails, errors, broken)
    counts = String[]
    fails > 0 && push!(counts, "$fails failed")
    errors > 0 && push!(counts, "$errors errored")
    broken > 0 && push!(counts, "$broken broken")
    return join(counts, ", ")
end

function result_cell(result)
    isnothing(result) && return "not run"
    fails, errors, broken = parse.(Int, result[2:4])
    duration = result[5]
    if fails + errors == 0
        return broken == 0 ? "✅ $duration" : "✅ $duration, $(count_summary(fails, errors, broken))"
    else
        return "❌ $(count_summary(fails, errors, broken))"
    end
end

const VERSIONS = sort!(collect(keys(INCLUDED_SCRIPTS)))
const MISSING_SCRIPTS = Dict(version => sort!(collect(setdiff(EXPECTED_SCRIPTS, INCLUDED_SCRIPTS[version]))) for version in VERSIONS)
const TESTS_RESULT = get(ENV, "TESTS_RESULT", "success")

io = IOBuffer()
println(io, "## Tests\n")
if isempty(VERSIONS)
    println(io, "No test reports were found.")
else
    println(io, "Julia versions: ", join(("`$version` is $(JULIA_VERSIONS[version])" for version in VERSIONS), ", "), ".\n")
    println(io, "| Group | ", join(("`$version`" for version in VERSIONS), " | "), " |")
    println(io, "|---", "|---"^length(VERSIONS), "|")
    for group in GROUPS
        println(io, "| $group | ", join((result_cell(get(GROUP_RESULTS, (version, group), nothing)) for version in VERSIONS), " | "), " |")
    end
    totals = map(VERSIONS) do version
        results = [parse.(Int, result[1:4]) for ((v, _), result) in GROUP_RESULTS if v == version]
        passes, fails, errors, broken = isempty(results) ? (0, 0, 0, 0) : sum(results)
        return join(filter(!isempty, ["$passes passed", count_summary(fails, errors, broken)]), ", ")
    end
    println(io, "| **Total** | ", join(totals, " | "), " |")

    if !isempty(FAILED_TESTSETS)
        println(io, "\n### Failed testsets\n")
        for (version, group, testset, fails, errors) in sort!(FAILED_TESTSETS)
            println(io, "- `$version`, $group › $testset: ", count_summary(fails, errors, 0))
        end
    end

    println(io, "\n### Test scripts\n")
    if all(isempty, values(MISSING_SCRIPTS))
        println(io, "All $(length(EXPECTED_SCRIPTS)) test scripts were run on every Julia version and operating system.")
    else
        for version in VERSIONS
            isempty(MISSING_SCRIPTS[version]) && continue
            println(io, "These test scripts were not run on `$version`:")
            foreach(script -> println(io, "- `$script`"), MISSING_SCRIPTS[version])
        end
        println(io, "\nEach test group needs to be listed in the test matrix in .github/workflows/CI.yml. If a job failed before it ",
            "finished, the scripts of its groups are missing too.")
    end
end
TESTS_RESULT == "success" || println(io, "\nThe result of the test jobs was `$TESTS_RESULT`; see their logs.")

summary = String(take!(io))
print(summary)
haskey(ENV, "GITHUB_STEP_SUMMARY") && open(file -> write(file, summary), ENV["GITHUB_STEP_SUMMARY"], "a")

passed = !isempty(VERSIONS) && isempty(FAILED_TESTSETS) && all(isempty, values(MISSING_SCRIPTS)) && TESTS_RESULT == "success" &&
    all(parse(Int, result[2]) + parse(Int, result[3]) == 0 for result in values(GROUP_RESULTS))
exit(passed ? 0 : 1)
