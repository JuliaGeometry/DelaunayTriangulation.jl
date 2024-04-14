"""
    orient_predicate(p, q, r) -> Integer 

Returns `ExactPredicates.orient(p, q, r)`, in particular we return:

- `1`: `(p, q, r)` is positively oriented.
- `0`: `(p, q, r)` is collinear / degenerate.
- `-1`: `(p, q, r)` is negatively oriented.

# Extended help 
The orient predicate is defined by the determinant 

```math 
\\text{orient}(p, q, r) = \\text{sgn} \\det \\begin{bmatrix} p_x & p_y & 1 \\\\ q_x & q_y & 1 \\\\ r_x & r_y & 1 \\end{bmatrix} = \\text{sgn} \\det \\begin{bmatrix} p_x-r_x & p_y-r_y \\\\ q_x-r_x & q_y-r_y \\end{bmatrix}.
```
"""
orient_predicate(p, q, r) = orient(_getxy(p), _getxy(q), _getxy(r))

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

# Extended help 
The orient predicate is defined by the determinant 

```math 
\\text{orient}(p, q, r, s) = \\text{sgn} \\det \\begin{bmatrix} p_x & p_y & p_y & 1 \\\\ q_x & q_y & q_r & 1 \\\\ r_x & r_y & r_z & 1 \\\\ s_x & s_y & s_z & 1 \\end{bmatrix} = \\text{sgn} \\det \\begin{bmatrix} p_x - s_x & p_y - s_y & p_z - s_y \\\\ q_x - s_x & q_y - s_y & q_z - s_z \\\\ r_x - s_x & r_y - s_y & r_z - s_z \\end{bmatrix}.
```
"""
orient_predicate(p, q, r, s) = orient(p, q, r, s)

"""
    incircle_predicate(a, b, c, p) -> Integer

Returns `ExactPredicates.incircle(a, b, c, p)`, in particular we return:

- `1`: If `p` is inside the circle defined by `(a, b, c)`.
- `0`: If `p` is on the circle defined by `(a, b, c)`.
- `-1`: If `p` is outside the circle defined by `(a, b, c)`.

# Extended help 
The incircle predicate is defined by the determinant

```math
\\text{incircle}(a, b, c, d) = \\text{sgn} \\det \\begin{bmatrix} a_x & a_y & a_x^2 + a_y^2 & 1 \\\\ b_x & b_y & b_x62 + b_y^2 & 1 \\\\ c_x & c_y & c_x^2 + c_y^2 & 1 \\\\ d_x & d_y & d_x^2 + d_y^2 & 1 \\end{bmatrix} = \\text{sgn} \\det \\begin{bmatrix} a_x - d_x & a_y - d_y & (a_x - d_x)^2 + (a_y - d_y)^2 \\\\ b_x - d_x & b_y - d_y & (b_x - d_x)^2 + (b_y - d_y)^2 \\\\ c_x - d_x & c_y - d_y & (c_x - d_x)^2 + (c_y - d_y)^2 \\end{bmatrix}.
```
"""
incircle_predicate(a, b, c, p) = incircle(_getxy(a), _getxy(b), _getxy(c), _getxy(p))

"""
    parallelorder_predicate(a, b, p, q) -> Integer

Returns `ExactPredicates.parallelorder(a, b, p, q)`, in particular we return:

- `1`: `q` is closer to the line `(a, b)` than `p`.
- `0`: `p` and `q` are equidistant from the line `(a, b)`.
- `-1`: `p` is closer to the line `(a, b)` than `q`.
"""
parallelorder_predicate(a, b, p, q) = parallelorder(_getxy(a), _getxy(b), _getxy(p), _getxy(q))

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
    opposite_signs(x, y) -> Bool

From ExactPredicates.jl, returns `true` if `x` and `y` have opposite signs, and `false` otherwise.
"""
opposite_signs(x, y) = xor(x, y) == -2

"""
    meet_predicate(p, q, a, b) 

Returns `ExactPredicates.meet(p, q, a, b)`, in particular we return:

- `1`: The open line segments `(p, q)` and `(a, b)` meet in a single point.
- `0`: The closed line segments `[p, q]` and `[a, b]` meet in one or several points. 
- `-1`: Otherwise.
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