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
ext(ux, uy, vx, vy) = ux * vy - uy * vx

"""
    inp(ux, uy, vx, vy)

Computes `ExactPredicates.inp((ux, uy), (vx, vy))`, i.e.
returns `ux * vx + uy * vy`.
"""
inp(ux, uy, vx, vy) = ux * vx + uy * vy

"""
    det(a, b, c, d)

Computes `ExactPredicates.det(a, b, c, d)`,
i.e. returns `a*d - b*c`.
"""
det(a, b, c, d) = a * d - b * c

"""
    det(a, b, c, d, e, f, g, h, i)

Computes `ExactPredicates.det(a, b, c, d, e, f, g, h, i)`,
i.e. returns the determinant of 
```math
\\det \\begin{bmatrix} a & b & c \\\\ d & e & f \\\\ g & h & i \\end{bmatrix}
```
"""
det(a, b, c, d, e, f, g, h, i) = a * det(e, f, h, i) - d * det(b, c, h, i) + g * det(b, c, e, f)

"""
    opposite_signs(x, y) -> Bool

From ExactPredicates.jl, returns `true` if `x` and `y` have opposite signs, and `false` otherwise.
Assumes that `x` and `y` are in `[-1, 0, 1]`.
"""
opposite_signs(x, y) = xor(x, y) == -2

"""
    sgn(x) -> Int 

Returns `Int(sign(x))`.
""" 
sgn(x) = Int(sign(x))

"""
    orient_predicate([method::AbstractPredicateType,] p, q, r) -> Integer 

Returns 

- `1`: `(p, q, r)` is positively oriented.
- `0`: `(p, q, r)` is collinear / degenerate.
- `-1`: `(p, q, r)` is negatively oriented.

The `method` argument determines how this result is computed, and should be 
one of [`Exact`](@ref), [`Fast`](@ref), and [`Adaptive`](@ref) (the default).
See the documentation for more information about these choices.
"""
function orient_predicate(method::AbstractPredicateType, p, q, r) 
    return orient(method, getxy(p), getxy(q), getxy(r))
end
orient_predicate(p, q, r) = orient_predicate(def_alg222(), p, q, r)

orient(::Fast, p, q, r) = sgn(AP.orient2fast(p, q, r))
orient(::Exact, p, q, r) = EP.orient(_getxy(p), _getxy(q), _getxy(r))
orient(::Adaptive, p, q, r) = AP.orient2p(p, q, r)

"""
    orient_predicate([method::AbstractPredicateType,] p, q, r) -> Integer

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

The `method` argument determines how this result is computed, and should be 
one of [`Exact`](@ref), [`Fast`](@ref), and [`Adaptive`](@ref) (the default).
See the documentation for more information about these choices.
"""
function orient_predicate(method::AbstractPredicateType, p, q, r, s) 
    return orient(method, getxy(p), getxy(q), getxy(r), getxy(s))
end
orient_predicate(p, q, r, s) = orient_predicate(def_alg222(), p, q, r, s)

orient(::Fast, p, q, r, s) = sgn(AP.orient3fast(p, q, r, s))
orient(::Exact, p, q, r, s) = EP.orient(_getxy(p), _getxy(q), _getxy(r), _getxy(s))
orient(::Adaptive, p, q, r, s) = AP.orient3(p, q, r, s)

"""
    incircle_predicate([method::AbstractPredicateType,] a, b, c, p) -> Integer

Assuming that `(a, b, c)` is a positively oriented triangle, returns

- `1`: If `p` is inside the circle defined by `(a, b, c)`.
- `0`: If `p` is on the circle defined by `(a, b, c)`.
- `-1`: If `p` is outside the circle defined by `(a, b, c)`.

The `method` argument determines how this result is computed, and should be 
one of [`Exact`](@ref), [`Fast`](@ref), and [`Adaptive`](@ref) (the default).
See the documentation for more information about these choices.
"""
function incircle_predicate(method::AbstractPredicateType, a, b, c, p) 
    return incircle(method, getxy(p), getxy(b), getxy(c), getxy(p))
end
incircle_predicate(a, b, c, p) = incircle_predicate(def_alg222(), a, b, c, p)

incircle(::Fast, a, b, c, p) = sgn(AP.incirclefast(a, b, c, p))
incircle(::Exact, a, b, c, p) = EP.incircle(_getxy(a), _getxy(b), _getxy(c), _getxy(p))
incircle(::Adaptive, a, b, c, p) = AP.incirclep(a, b, c, p)

"""
    parallelorder_predicate([method::AbstractPredicateType,] a, b, p, q) -> Integer

Returns

- `1`: `q` is closer to the line `(a, b)` than `p`.
- `0`: `p` and `q` are equidistant from the line `(a, b)`.
- `-1`: `p` is closer to the line `(a, b)` than `q`.

The `method` argument determines how this result is computed, and should be 
one of [`Exact`](@ref), [`Fast`](@ref), and [`Adaptive`](@ref) (the default).
See the documentation for more information about these choices.
"""
function parallelorder_predicate(method::AbstractPredicateType, a, b, p, q)
    return parallelorder(method, getxy(a), getxy(b), getxy(p), getxy(q))
end
parallelorder_predicate(a, b, p, q) = parallelorder(def_alg222(), a, b, p, q)

function parallelorder(::Fast, a, b, p, q)
    ax, ay = getxy(a) 
    bx, by = getxy(b)
    px, py = getxy(p) 
    qx, qy = getxy(q) 
    bax, bay = bx - ax, by - ay 
    qpx, qpy = qx - px, qy - py 
    return sgn(ext(bax, bay, qpx, qpy))
end
parallelorder(::Exact, a, b, p, q) =  EP.parallelorder(_getxy(a), _getxy(b), _getxy(p), _getxy(q))
parallelorder(::Adaptive, a, b, p, q) = parallelorder(def_alg222(), a, b, p, q) # not implemented yet 

"""
    angle_is_acute([method::AbstractPredicateType,] p, q, r)

Tests if the angle opposite `(p, q)` in the triangle `(p, q, r)`, 
meaning `∠prq`, is acute, returning:

- `1`: `∠prq` is acute.
- `0`: `∠prq` is a right angle.
- `-1`: `∠prq` is obtuse.

The `method` argument determines how this result is computed, and should be 
one of [`Exact`](@ref), [`Fast`](@ref), and [`Adaptive`](@ref) (the default).
See the documentation for more information about these choices.
"""
function angle_is_acute_predicate(method::AbstractPredicateType, p, q, r) 
    return angle_is_acute(method, getxy(p), getxy(q), getxy(r))
end
angle_is_acute_predicate(p, q, r) = angle_is_acute_predicate(def_alg222(), p, q, r)

EP.Codegen.@genpredicate function _angle_is_acute(p::2, q::2, r::2)
    pr = p - r
    qr = q - r
    EP.Codegen.group!(pr...)
    EP.Codegen.group!(qr...)
    return pr[1] * qr[1] + pr[2] * qr[2]
end
function angle_is_acute(::Fast, p, q, r)
    px, py = getxy(p)
    qx, qy = getxy(q) 
    rx, ry = getxy(r) 
    prx, pry = px - rx, py - ry 
    qrx, qry = qx - rx, qy - ry 
    return sgn(inp(prx, pry, qrx, qry))
end
angle_is_acute(::Exact, p, q, r) = _angle_is_acute(_getxy(p), _getxy(q), _getxy(r))
angle_is_acute(::Adaptive, p, q, r) = angle_is_acute(def_alg222(), p, q, r) # not implemented yet 

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
function sameside_predicate(a, b, p)
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
    meet_predicate([method::AbstractPredicateType], p, q, a, b) -> Integer

Returns

- `1`: The open line segments `(p, q)` and `(a, b)` meet in a single point.
- `0`: The closed line segments `[p, q]` and `[a, b]` meet in one or several points. 
- `-1`: Otherwise.

The `method` argument determines how this result is computed, and should be 
one of [`Exact`](@ref), [`Fast`](@ref), and [`Adaptive`](@ref) (the default).
See the documentation for more information about these choices.
"""
function meet_predicate(method::AbstractPredicateType, p, q, a, b)
    pqa = orient_predicate(method, p, q, a)
    pqb = orient_predicate(method, p, q, b)
    abp = orient_predicate(method, a, b, p)
    abq = orient_predicate(method, a, b, q)
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
meet_predicate(p, q, a, b) = meet_predicate(def_alg222(), p, q, a, b)