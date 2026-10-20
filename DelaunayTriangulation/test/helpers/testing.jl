# Test utilities: logging of failing tests, reference image paths, random kernels, and comparisons of types and approximate values.

export @test_with_log, @check, ⊢, ⪧, rt, reference_path

# Like @test, but also prints msg if the test does not pass. This used to wrap each test in a `@testset let ctx = msg`,
# but creating a testset for every test made e.g. validate_statistics about 14 times slower, since it runs millions of tests.
macro test_with_log(msg, expr, kws...)
    if msg == ""
        return esc(:(@test $expr $(kws...)))
    else
        return esc(quote
            let result = @test $expr $(kws...)
                result isa Test.Pass || println("  Context: ", $msg)
                result
            end
        end)
    end
end

"""
    @check failures key expr [kws...]

Evaluates `expr` like `@test` does, so that e.g. `a ≈ b rtol = 1.0e-6` works, and if it is `false`, records `key` in the
`Dict` `failures` under the text of `expr` (keeping only the first `key` for each check). This is for checking every
element of a large collection, e.g. every triangle of a triangulation, where a `@test` for each element would be slow:
use `@check` for each element, and then `@test isempty(failures)`.
"""
macro check(failures, key, expr, kws...)
    if isempty(kws)
        cond = expr
    else
        expr isa Expr && expr.head === :call && expr.args[1] in (:≈, :isapprox) || error("keyword arguments are only supported for ≈")
        kwargs = [Expr(:kw, kw.args[1], kw.args[2]) for kw in kws]
        cond = Expr(:call, :isapprox, Expr(:parameters, kwargs...), expr.args[2], expr.args[3])
    end
    label = join((string(expr), string.(kws)...), " ")
    return esc(:(($cond) || get!($failures, $label, $key); nothing))
end

"""
    reference_path(dir, name)

Returns the path to the reference image `name` in `dir`. Some figures depend on hash iteration
order (e.g. which triangles mesh refinement splits first), which changed in Julia 1.13. For those
figures, the Julia 1.13 version of the image is stored in `dir/julia-1.13/` and used on Julia ≥ 1.13.

To regenerate the Julia 1.13 images, e.g. after a change that affects these figures:

1. Run the tests on Julia ≥ 1.13 with the environment variable `JULIA_REFERENCETESTS_UPDATE=true`. This overwrites
   only the reference images whose tests fail.
2. Move each changed image into the `julia-1.13/` folder next to it (replacing any existing copy there).
3. Restore the original images with `git checkout -- <file>`, since those are still the references for Julia < 1.13.

If a figure changes on Julia < 1.13 instead, regenerate it on that version in the same way, but keep the updated image in
place rather than moving it.
"""
function reference_path(dir, name)
    if VERSION ≥ v"1.13-"
        path = joinpath(dir, "julia-1.13", name)
        isfile(path) && return path
    end
    return joinpath(dir, name)
end

⊢(x, y) = x == y && typeof(x) == typeof(y)

⪧(a::Vector{Tuple}, b::Vector{Tuple}; kwargs...) = ⪧(collect.(a), collect.(b); kwargs...)

⪧(a::Vector{<:Union{Vector,Number}}, b::Vector; kwargs...) = isapprox(a, b; kwargs...)

⪧(a, b; kwargs...) = isapprox(collect(collect.(a)), collect(collect.(b)); kwargs...)

function ⪧(x::DT.CatmullRomSplineSegment, y::DT.CatmullRomSplineSegment; kwargs...)
    ⪧(x.a, y.a; kwargs...) || return false
    ⪧(x.b, y.b; kwargs...) || return false
    ⪧(x.c, y.c; kwargs...) || return false
    ⪧(x.d, y.d; kwargs...) || return false
    ⪧(x.p₁, y.p₁; kwargs...) || return false
    ⪧(x.p₂, y.p₂; kwargs...) || return false
    return true
end

rt() = rand((DT.FastKernel(), DT.ExactKernel(), DT.AdaptiveKernel()))
