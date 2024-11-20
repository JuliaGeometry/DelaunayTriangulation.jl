#=
Much of the code in the file is derived or taken from ExactPredicates.jl (https://github.com/lairez/ExactPredicates.jl). 
The license for ExactPredicates.jl is printed below, taken from https://github.com/lairez/ExactPredicates.jl/blob/a9dce0334ee62104a14930fe206764bf802c23db/LICENSE.

Copyright (c) 2019 Pierre Lairez

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#
struct FastKernel <: AbstractPredicateKernel end
struct ExactKernel <: AbstractPredicateKernel end
struct AdaptiveKernel <: AbstractPredicateKernel end

const DEFAULT_KERNEL = AdaptiveKernel()

@inline ext(ux, uy, vx, vy) = ux * vy - uy * vx

@inline inp(ux, uy, vx, vy) = ux * vx + uy * vy

@inline det(a, b, c, d) = a * d - b * c

@inline det(a, b, c, d, e, f, g, h, i) = a * det(e, f, h, i) - b * det(d, f, g, i) + c * det(d, e, g, h)

@inline opposite_signs(x, y) == xor(x, y) == -2

@inline sgn(x) = Int(sign(x))

@inline @optarg1 DEFAULT_KERNEL function orient_predicate(kernel::AbstractPredicateKernel, p, q, r; ctr)
    add_orient2!(ctr)
    return orient(kernel, getxy(p), getxy(q), getxy(r))
end

@inline orient(::FastKernel, p, q, r) = sgn(AdaptivePredicates.orient2fast(p, q, r))
@inline orient(::ExactKernel, p, q, r) = ExactPredicates.orient(f64_getxy(p), f64_getxy(q), f64_getxy(r))
@inline orient(::AdaptiveKernel, p, q, r) = AdaptivePredicates.orient2p(p, q, r)

@inline @optarg1 DEFAULT_KERNEL function orient_predicate(kernel::AbstractPredicateKernel, p, q, r, s; cache::PredicateCacheType=nothing, ctr)
    add_orient3!(ctr)
    return orient(kernel, getxyz(p), getxyz(q), getxyz(r), getxyz(s), cache)
end

@inline orient(::FastKernel, p, q, r, s, ::PredicateCacheType=nothing) = sgn(AdaptivePredicates.orient3fast(p, q, r, s))
@inline orient(::ExactKernel, p, q, r, s, ::PredicateCacheType=nothing) = ExactPredicates.orient(f64_getxyz(p), f64_getxyz(q), f64_getxyz(r), f64_getxyz(s))
@inline orient(::AdaptiveKernel, p, q, r, s, cache::PredicateCacheType=nothing) = AdaptivePredicates.orient3p(p, q, r, s, cache)

@inline @optarg1 DEFAULT_KERNEL function incircle_predicate(kernel::AbstractPredicateKernel, a, b, c, p; cache::PredicateCacheType=nothing, ctr)
    add_incircle!(ctr)
    return incircle(kernel, getxy(a), getxy(b), getxy(c), getxy(p), cache)
end

@inline incircle(::FastKernel, a, b, c, p, _::PredicateCacheType=nothing) = sgn(AdaptivePredicates.incirclefast(a, b, c, p))
@inline incircle(::ExactKernel, a, b, c, p, _::PredicateCacheType=nothing) = ExactPredicates.incircle(f64_getxy(a), f64_getxy(b), f64_getxy(c), f64_getxy(p))
@inline incircle(::AdaptiveKernel, a, b, c, p, cache::PredicateCacheType=nothing) = AdaptivePredicates.incirclep(a, b, c, p, cache)

@inline @optarg1 DEFAULT_KERNEL function parallelorder_predicate(kernel::AbstractPredicateKernel, a, b, p, q; ctr)
    add_parallelorder!(ctr)
    return parallelorder(kernel, getxy(a), getxy(b), getxy(p), getxy(q))
end

@inline function parallelorder(::FastKernel, a, b, p, q)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    px, py = getxy(p)
    qx, qy = getxy(q)
    bax, bay = bx - ax, by - ay
    qpx, qpy = qx - px, qy - py
    return sgn(ext(bax, bay, qpx, qpy))
end
@inline parallelorder(::ExactKernel, a, b, p, q) = ExactPredicates.parallelorder(f64_getxy(a), f64_getxy(b), f64_getxy(p), f64_getxy(q))
@inline parallelorder(::AdaptiveKernel, a, b, p, q) = parallelorder(ExactKernel(), a, b, p, q) # not implemented yet 

@inline @optarg1 DEFAULT_KERNEL function angle_is_acute_predicate(kernel::AbstractPredicateKernel, p, q, r; ctr)
    add_angle_is_acute!(ctr)
    return angle_is_acute(kernel, getxy(p), getxy(q), getxy(r))
end

ExactPredicates.Codegen.@genpredicate function _angle_is_acute(p::2, q::2, r::2)
    pr = p - r
    qr = q - r
    EP.Codegen.group!(pr...)
    EP.Codegen.group!(qr...)
    return pr[1] * qr[1] + pr[2] * qr[2]
end
@inline function angle_is_acute(::FastKernel, p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
    prx, pry = px - rx, py - ry
    qrx, qry = qx - rx, qy - ry
    return sgn(inp(prx, pry, qrx, qry))
end
@inline angle_is_acute(::ExactKernel, p, q, r) = _angle_is_acute(f64_getxy(p), f64_getxy(q), f64_getxy(r))
@inline angle_is_acute(::AdaptiveKernel, p, q, r) = angle_is_acute(ExactKernel(), p, q, r) # not implemented yet 

@inline function sameside_predicate(a, b, p; ctr)
    add_sameside!(ctr)
    p, a, b = getxy(p), getxy(a), getxy(b)
    if a < p && b < p || a > p && b > p
        return 1
    elseif a < p && b > p || a > p && b < p
        return -1
    else
        return 0
    end
end

@inline @optarg1 DEFAULT_KERNEL function meet_predicate(kernel::AbstractPredicateKernel, p, q, a, b; ctr)
    add_meet!(ctr)
    pqa = orient_predicate(kernel, p, q, a; ctr)
    pqb = orient_predicate(kernel, p, q, b; ctr)
    abp = orient_predicate(kernel, a, b, p; ctr)
    abq = orient_predicate(kernel, a, b, q; ctr)
    if opposite_signs(pqa, pqb) && opposite_signs(abp, abq)
        return 1
    elseif (!iszero(pqa) && pqa == pqb) || (!iszero(abq) && abp == abq)
        return -1
    elseif iszero(pqa) && iszero(pqb)
        if isone(sameside_predicate(a, b, p; ctr)) &&
           isone(sameside_predicate(a, b, q; ctr)) &&
           isone(sameside_predicate(p, q, a; ctr)) &&
           isone(sameside_predicate(p, q, b; ctr))
            return -1
        else
            return 0
        end
    else
        return 0
    end
end
