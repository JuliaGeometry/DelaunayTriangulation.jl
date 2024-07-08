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
    orient_predicate(p, q, r) -> Integer 

Returns `ExactPredicates.orient(p, q, r)`, in particular we return:

- `1`: `(p, q, r)` is positively oriented.
- `0`: `(p, q, r)` is collinear / degenerate.
- `-1`: `(p, q, r)` is negatively oriented.

If ExactPredicates.jl has been disabled using Preferences.jl, the 
determinant defining `orient` is evaluated numerically 
without exact arithmetic.

# Extended help 
The orient predicate is defined by the determinant 

```math 
\\text{orient}(p, q, r) = \\text{sgn} \\det \\begin{bmatrix} p_x & p_y & 1 \\\\ q_x & q_y & 1 \\\\ r_x & r_y & 1 \\end{bmatrix} = \\text{sgn} \\det \\begin{bmatrix} p_x-r_x & p_y-r_y \\\\ q_x-r_x & q_y-r_y \\end{bmatrix}.
```
"""
orient_predicate(::Any, ::Any, ::Any)

@static if USE_EXACTPREDICATES
    orient_predicate(p, q, r) = orient(_getxy(p), _getxy(q), _getxy(r))
elseif USE_INEXACTPREDICATES
    orient_predicate(p, q, r) = _orient_predicate(_getxy(p), _getxy(q), _getxy(r))
end
function _orient_predicate(p, q, r)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    rx, ry = _getxy(r)
    ux, uy = px - rx, py - ry
    vx, vy = qx - rx, qy - ry
    return sgn(ext(ux, uy, vx, vy))
end

"""
    orient_predicate(p, q, r, s) -> Integer

Returns `ExactPredicates.orient(p, q, r, s)`, in particular we return:

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

If ExactPredicates.jl has been disabled using Preferences.jl, the 
determinant defining `orient` is evaluated numerically 
without exact arithmetic.

# Extended help 
The orient predicate is defined by the determinant 

```math 
\\text{orient}(p, q, r, s) = \\text{sgn} \\det \\begin{bmatrix} p_x & p_y & p_z & 1 \\\\ q_x & q_y & q_z & 1 \\\\ r_x & r_y & r_z & 1 \\\\ s_x & s_y & s_z & 1 \\end{bmatrix} = \\text{sgn} \\det \\begin{bmatrix} p_x - s_x & p_y - s_y & p_z - s_y \\\\ q_x - s_x & q_y - s_y & q_z - s_z \\\\ r_x - s_x & r_y - s_y & r_z - s_z \\end{bmatrix}.
```
"""
orient_predicate(::Any, ::Any, ::Any, ::Any)

@static if USE_EXACTPREDICATES
    orient_predicate(p, q, r, s) = orient(_getxyz(p), _getxyz(q), _getxyz(r), _getxyz(s))
elseif USE_INEXACTPREDICATES
    orient_predicate(p, q, r, s) = _orient_predicate(_getxyz(p), _getxyz(q), _getxyz(r), _getxyz(s))
end
function _orient_predicate(p, q, r, s)
    px, py, pz = _getxyz(p)
    qx, qy, qz = _getxyz(q)
    rx, ry, rz = _getxyz(r)
    sx, sy, sz = _getxyz(s)
    a, b, c = px - sx, py - sy, pz - sz
    d, e, f = qx - sx, qy - sy, qz - sz
    g, h, i = rx - sx, ry - sy, rz - sz
    return sgn(det(a, b, c, d, e, f, g, h, i))
end

"""
    incircle_predicate(a, b, c, p) -> Integer

Returns `ExactPredicates.incircle(a, b, c, p)`, in particular we return:

- `1`: If `p` is inside the circle defined by `(a, b, c)`.
- `0`: If `p` is on the circle defined by `(a, b, c)`.
- `-1`: If `p` is outside the circle defined by `(a, b, c)`.

If ExactPredicates.jl has been disabled using Preferences.jl, the 
determinant defining `incircle` is evaluated numerically 
without exact arithmetic.

# Extended help 
The incircle predicate is defined by the determinant

```math
\\text{incircle}(a, b, c, d) = \\text{sgn} \\det \\begin{bmatrix} a_x & a_y & a_x^2 + a_y^2 & 1 \\\\ b_x & b_y & b_x^2 + b_y^2 & 1 \\\\ c_x & c_y & c_x^2 + c_y^2 & 1 \\\\ d_x & d_y & d_x^2 + d_y^2 & 1 \\end{bmatrix} = \\text{sgn} \\det \\begin{bmatrix} a_x - d_x & a_y - d_y & (a_x - d_x)^2 + (a_y - d_y)^2 \\\\ b_x - d_x & b_y - d_y & (b_x - d_x)^2 + (b_y - d_y)^2 \\\\ c_x - d_x & c_y - d_y & (c_x - d_x)^2 + (c_y - d_y)^2 \\end{bmatrix}.
```
"""
incircle_predicate

@static if USE_EXACTPREDICATES
    incircle_predicate(a, b, c, p) = incircle(_getxy(a), _getxy(b), _getxy(c), _getxy(p))
elseif USE_INEXACTPREDICATES
    incircle_predicate(a, b, c, p) = _incircle_predicate(_getxy(a), _getxy(b), _getxy(c), _getxy(p))
end
function _incircle_predicate(p, q, r, a)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    rx, ry = _getxy(r)
    ax, ay = _getxy(a)
    qpx, qpy = qx - px, qy - py
    rpx, rpy = rx - px, ry - py
    apx, apy = ax - px, ay - py
    aqx, aqy = ax - qx, ay - qy
    rqx, rqy = rx - qx, ry - qy
    val1 = ext(qpx, qpy, apx, apy) * inp(rpx, rpy, rqx, rqy)
    val2 = ext(qpx, qpy, rpx, rpy) * inp(apx, apy, aqx, aqy)
    return sgn(val1 - val2)
end

"""
    parallelorder_predicate(a, b, p, q) -> Integer

Returns `ExactPredicates.parallelorder(a, b, p, q)`, in particular we return:

- `1`: `q` is closer to the line `(a, b)` than `p`.
- `0`: `p` and `q` are equidistant from the line `(a, b)`.
- `-1`: `p` is closer to the line `(a, b)` than `q`.

If ExactPredicates.jl has been disabled using Preferences.jl, the 
the predicates defining `parallelorder` are all evaluated numerically 
without exact arithmetic.
"""
parallelorder_predicate

@static if USE_EXACTPREDICATES
    parallelorder_predicate(a, b, p, q) = parallelorder(_getxy(a), _getxy(b), _getxy(p), _getxy(q))
elseif USE_INEXACTPREDICATES
    parallelorder_predicate(a, b, p, q) = _parallelorder_predicate(_getxy(a), _getxy(b), _getxy(p), _getxy(q))
end
function _parallelorder_predicate(a, b, p, q)
    ax, ay = _getxy(a)
    bx, by = _getxy(b)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    δx, δy = bx - ax, by - ay
    qpx, qpy = qx - px, qy - py
    return sgn(ext(δx, δy, qpx, qpy))
end

"""
    sameside_predicate(a, b, p) -> Integer

Returns `ExactPredicates.sameside(p, a, b)` where all three points are collinear, in particular we return:

- `1`: `a` and `b` are on the same side of `p` on the line.
- `0`: `a == p` or `b == p`.
- `-1`: `a` and `b` are on different sides of `p` on the line.

!!! note 

    The difference in the argument order to ExactPredicates.jl is to match the convention that the 
    main point being tested is the last argument.
"""
function sameside_predicate(a, b, p)
    _p = _getxy(p)
    _a = _getxy(a)
    _b = _getxy(b)
    if _a < _p && _b < _p || _a > _p && _b > _p
        return 1
    elseif _a < _p && _b > _p || _a > _p && _b < _p
        return -1
    else
        return 0
    end
end

"""
    meet_predicate(p, q, a, b) 

Returns `ExactPredicates.meet(p, q, a, b)`, in particular we return:

- `1`: The open line segments `(p, q)` and `(a, b)` meet in a single point.
- `0`: The closed line segments `[p, q]` and `[a, b]` meet in one or several points. 
- `-1`: Otherwise.

If ExactPredicates.jl has been disabled using Preferences.jl, then all 
the `orient` predicates computed for this result are evaluated numerically 
without exact arithmetic.
"""
function meet_predicate(p, q, a, b)
    pqa = orient_predicate(p, q, a)
    pqb = orient_predicate(p, q, b)
    abp = orient_predicate(a, b, p)
    abq = orient_predicate(a, b, q)
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

@static if USE_EXACTPREDICATES
    ExactPredicates.Codegen.@genpredicate function angle_is_acute(p::2, q::2, r::2)
        pr = p - r
        qr = q - r
        ExactPredicates.Codegen.group!(pr...)
        ExactPredicates.Codegen.group!(qr...)
        return pr[1] * qr[1] + pr[2] * qr[2]
    end
elseif USE_INEXACTPREDICATES
    function angle_is_acute(p, q, r)
        px, py = _getxy(p)
        qx, qy = _getxy(q)
        rx, ry = _getxy(r)
        ux, uy = px - rx, py - ry
        vx, vy = qx - rx, qy - ry
        return sgn(inp(ux, uy, vx, vy))
    end
end

@doc """
    angle_is_acute(p, q, r)

Tests if the angle opposite `(p, q)` in the triangle `(p, q, r)`, 
meaning `∠prq`, is acute, returning:

- `1`: `∠prq` is acute.
- `0`: `∠prq` is a right angle.
- `-1`: `∠prq` is obtuse.
"""
angle_is_acute