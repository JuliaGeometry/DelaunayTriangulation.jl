# Test utilities: logging of failing tests, reference image paths, random kernels, and comparisons of types and approximate values.

export @test_with_log, ⊢, ⪧, rt, reference_path

macro test_with_log(msg, expr, kws...)
    if msg == ""
        return esc(:(@test $expr $(kws...)))
    else
        return esc(quote
            @testset let ctx = $msg
                @test $expr $(kws...)
            end
        end)
    end
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
