#=
Much of the code in the file is derived or taken from ExactPredicates.jl (https://github.com/lairez/ExactPredicates.jl). 
The license for ExactPredicates.jl is printed below, taken from https://github.com/lairez/ExactPredicates.jl/blob/master/LICENSE.

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

"""
    ext(ux, uy, vx, vy)

Computes `ExactPredicates.ext((ux, uy), (vx, vy))`, i.e. 
returns `ux * vy - uy * vx`.
"""
@inline ext(ux, uy, vx, vy) = ux * vy - uy * vx

"""
    inp(ux, uy, vx, vy)

Computes `ExactPredicates.inp((ux, uy), (vx, vy))`, i.e.
returns `ux * vx + uy * vy`.
"""
@inline inp(ux, uy, vx, vy) = ux * vx + uy * vy

"""
    det(a, b, c, d)

Computes `ExactPredicates.det(a, b, c, d)`,
i.e. returns `a*d - b*c`.
"""
@inline det(a, b, c, d) = a * d - b * c

"""
    det(a, b, c, d, e, f, g, h, i)

Computes `ExactPredicates.det(a, b, c, d, e, f, g, h, i)`,
i.e. returns the determinant of 
```math
\\det \\begin{bmatrix} a & b & c \\\\ d & e & f \\\\ g & h & i \\end{bmatrix}
```
"""
@inline det(a, b, c, d, e, f, g, h, i) = a * det(e, f, h, i) - d * det(b, c, h, i) + g * det(b, c, e, f)

"""
    opposite_signs(x, y) -> Bool

From ExactPredicates.jl, returns `true` if `x` and `y` have opposite signs, and `false` otherwise.
Assumes that `x` and `y` are in `[-1, 0, 1]`.
"""
@inline opposite_signs(x, y) = xor(x, y) == -2

"""
    sgn(x) -> Int 

Returns `Int(sign(x))`.
"""
@inline sgn(x) = Int(sign(x))

"""
    orient_predicate([kernel::AbstractPredicateKernel,] p, q, r) -> Integer 

Returns 

- `1`: `(p, q, r)` is positively oriented.
- `0`: `(p, q, r)` is collinear / degenerate.
- `-1`: `(p, q, r)` is negatively oriented.

The `kernel` argument determines how this result is computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.
"""
@inline function orient_predicate(kernel::AbstractPredicateKernel, p, q, r)
    return orient(kernel, getxy(p), getxy(q), getxy(r))
end
@inline orient_predicate(p, q, r) = orient_predicate(AdaptiveKernel(), p, q, r)

@inline orient(::FastKernel, p, q, r) = sgn(AP.orient2fast(p, q, r))
@inline orient(::ExactKernel, p, q, r) = EP.orient(_getxy(p), _getxy(q), _getxy(r))
@inline orient(::AdaptiveKernel, p, q, r) = AP.orient2p(p, q, r)

"""
    orient_predicate([kernel::AbstractPredicateKernel,] p, q, r; cache = nothing) -> Integer

Returns 

- `1`: `(p, q, r, s)` is positively oriented.
- `0`: `(p, q, r, s)` is collinear / degenerate.
- `-1`: `(p, q, r, s)` is negatively oriented.

Here, a positively oriented tetrahedron `(p, q, r, s)` takes the form

                                   z.
                                 .
                               ,/
                             s
                           ,/|'\\
                         ,/  |  '\\
                       ,/    '.   '\\
                     ,/       |     '\\                 
                   ,/         |       '\\              
                  p-----------'.--------q --> x
                   '\\.         |      ,/              
                      '\\.      |    ,/                 
                         '\\.   '. ,/    
                            '\\. |/      
                               'r       
                                 '\\.    

The `kernel` argument determines how this result is computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.

The optional `cache` keyword argument can be used for preallocating memory for intermediate results, passing the argument from `AdaptivePredicates.orient3adapt_cache(T)`,
where `T` is the number type of the input points. If `nothing` is passed, no cache is used. This is only needed if an `AdaptiveKernel()` is used.
"""
@inline function orient_predicate(kernel::AbstractPredicateKernel, p, q, r, s; cache::PredicateCacheType = nothing)
    return orient(kernel, getxyz(p), getxyz(q), getxyz(r), getxyz(s), cache)
end
@inline orient_predicate(p, q, r, s; cache::PredicateCacheType = nothing) = orient_predicate(AdaptiveKernel(), p, q, r, s; cache)

@inline orient(::FastKernel, p, q, r, s, _::PredicateCacheType = nothing) = sgn(AP.orient3fast(p, q, r, s))
@inline orient(::ExactKernel, p, q, r, s, _::PredicateCacheType = nothing) = EP.orient(_getxyz(p), _getxyz(q), _getxyz(r), _getxyz(s))
@inline orient(::AdaptiveKernel, p, q, r, s, cache::PredicateCacheType = nothing) = AP.orient3p(p, q, r, s, cache)

"""
    incircle_predicate([kernel::AbstractPredicateKernel,] a, b, c, p; cache = nothing) -> Integer

Assuming that `(a, b, c)` is a positively oriented triangle, returns

- `1`: If `p` is inside the circle defined by `(a, b, c)`.
- `0`: If `p` is on the circle defined by `(a, b, c)`.
- `-1`: If `p` is outside the circle defined by `(a, b, c)`.

The `kernel` argument determines how this result is computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.

The optional `cache` keyword argument can be used for preallocating memory for intermediate results, passing the argument from `AdaptivePredicates.incircleadapt_cache(T)`,
where `T` is the number type of the input points. If `nothing` is passed, no cache is used. This is only needed if an `AdaptiveKernel()` is used.
"""
@inline function incircle_predicate(kernel::AbstractPredicateKernel, a, b, c, p; cache::PredicateCacheType = nothing)
    return incircle(kernel, getxy(a), getxy(b), getxy(c), getxy(p), cache)
end
@inline incircle_predicate(a, b, c, p; cache::PredicateCacheType = nothing) = incircle_predicate(AdaptiveKernel(), a, b, c, p; cache)

@inline incircle(::FastKernel, a, b, c, p, _::PredicateCacheType = nothing) = sgn(AP.incirclefast(a, b, c, p))
@inline incircle(::ExactKernel, a, b, c, p, _::PredicateCacheType = nothing) = EP.incircle(_getxy(a), _getxy(b), _getxy(c), _getxy(p))
@inline incircle(::AdaptiveKernel, a, b, c, p, cache::PredicateCacheType = nothing) = AP.incirclep(a, b, c, p, cache)

"""
    parallelorder_predicate([kernel::AbstractPredicateKernel,] a, b, p, q) -> Integer

Returns

- `1`: `q` is closer to the line `(a, b)` than `p`.
- `0`: `p` and `q` are equidistant from the line `(a, b)`.
- `-1`: `p` is closer to the line `(a, b)` than `q`.

The `kernel` argument determines how this result is computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.
"""
@inline function parallelorder_predicate(kernel::AbstractPredicateKernel, a, b, p, q)
    return parallelorder(kernel, getxy(a), getxy(b), getxy(p), getxy(q))
end
@inline parallelorder_predicate(a, b, p, q) = parallelorder(AdaptiveKernel(), a, b, p, q)

@inline function parallelorder(::FastKernel, a, b, p, q)
    ax, ay = getxy(a)
    bx, by = getxy(b)
    px, py = getxy(p)
    qx, qy = getxy(q)
    bax, bay = bx - ax, by - ay
    qpx, qpy = qx - px, qy - py
    return sgn(ext(bax, bay, qpx, qpy))
end
@inline parallelorder(::ExactKernel, a, b, p, q) = EP.parallelorder(_getxy(a), _getxy(b), _getxy(p), _getxy(q))
@inline parallelorder(::AdaptiveKernel, a, b, p, q) = parallelorder(ExactKernel(), a, b, p, q) # not implemented yet 

"""
    angle_is_acute([kernel::AbstractPredicateKernel,] p, q, r)

Tests if the angle opposite `(p, q)` in the triangle `(p, q, r)`, 
meaning `∠prq`, is acute, returning:

- `1`: `∠prq` is acute.
- `0`: `∠prq` is a right angle.
- `-1`: `∠prq` is obtuse.

The `kernel` argument determines how this result is computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.
"""
@inline function angle_is_acute_predicate(kernel::AbstractPredicateKernel, p, q, r)
    return angle_is_acute(kernel, getxy(p), getxy(q), getxy(r))
end
@inline angle_is_acute_predicate(p, q, r) = angle_is_acute_predicate(AdaptiveKernel(), p, q, r)

EP.Codegen.@genpredicate function _angle_is_acute(p::2, q::2, r::2)
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
@inline angle_is_acute(::ExactKernel, p, q, r) = _angle_is_acute(_getxy(p), _getxy(q), _getxy(r))
@inline angle_is_acute(::AdaptiveKernel, p, q, r) = angle_is_acute(ExactKernel(), p, q, r) # not implemented yet 

"""
    sameside_predicate(a, b, p) -> Integer

Assuming all of `a, b, p` are collinear, returns

- `1`: `a` and `b` are on the same side of `p` on the line.
- `0`: `a == p` or `b == p`.
- `-1`: `a` and `b` are on different sides of `p` on the line.

!!! note 

    The difference in the argument order to ExactPredicates.jl is to match the convention that the 
    main point being tested is the last argument.
"""
@inline function sameside_predicate(a, b, p)
    _p = getxy(p)
    _a = getxy(a)
    _b = getxy(b)
    if _a < _p && _b < _p || _a > _p && _b > _p
        return 1
    elseif _a < _p && _b > _p || _a > _p && _b < _p
        return -1
    else
        return 0
    end
end

"""
    meet_predicate([kernel::AbstractPredicateKernel], p, q, a, b) -> Integer

Returns

- `1`: The open line segments `(p, q)` and `(a, b)` meet in a single point.
- `0`: The closed line segments `[p, q]` and `[a, b]` meet in one or several points. 
- `-1`: Otherwise.

The `kernel` argument determines how this result is computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.
"""
@inline function meet_predicate(kernel::AbstractPredicateKernel, p, q, a, b)
    pqa = orient_predicate(kernel, p, q, a)
    pqb = orient_predicate(kernel, p, q, b)
    abp = orient_predicate(kernel, a, b, p)
    abq = orient_predicate(kernel, a, b, q)
    if opposite_signs(pqa, pqb) && opposite_signs(abp, abq)
        return 1
    elseif (pqa ≠ 0 && pqa == pqb) || (abq ≠ 0 && abp == abq)
        return -1
    elseif pqa == 0 && pqb == 0
        if sameside_predicate(a, b, p) == 1 &&
                sameside_predicate(a, b, q) == 1 &&
                sameside_predicate(p, q, a) == 1 &&
                sameside_predicate(p, q, b) == 1
            return -1
        else
            return 0
        end
    else
        return 0
    end
end
@inline meet_predicate(p, q, a, b) = meet_predicate(AdaptiveKernel(), p, q, a, b)
