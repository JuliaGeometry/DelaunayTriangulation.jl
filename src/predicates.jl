function ExactPredicates.orient(𝒯::TriangleType, pts)
    u, v, w = 𝒯
    return orient(pts[u], pts[v], pts[w])
end
function ExactPredicates.incircle(pts, i, j, k, ℓ)
    return incircle(pts[i], pts[j], pts[k], pts[ℓ])
end

"""
    is_point_higher(p, q)
    is_point_lower(p, q)

Tests if `p` is lexicographically higher than `q`, where we say that `p = (xp, yp)`
is lexicographically higher than `q = (xq, yp)` if `yp > yq` or 
`yp = yq` and `xq > xp`.
"""
is_point_higher(p, q) = (p[2] > q[2]) || (p[2] == q[2] && q[1] > p[1])
is_point_lower(p, q) = is_point_higher(q, p)

"""
    partial_highest_point_sort!(v, k)

Partially sorts the first `v` so that the first `k` entries are the highest points, with the first 
being the highest (acccording to `is_point_higher`).
"""
partial_highest_point_sort!(v, k) = partialsort!(v, k, lt = is_point_higher)

"""
    leftofline(p, pᵢ, pⱼ)

Tests if `p` is to the left of the line from `pᵢ` to `pⱼ`.
"""
leftofline(p, pᵢ, pⱼ) = orient(pᵢ, pⱼ, p)
"""
    leftofline(pts, p, i, j)

Tests if the point `p` is to the left of the oriented line through 
`pts[i]` to `pts[j]`. Checks are made for non-positive indices.
"""
function leftofline(pts, p, i, j)
    if j == LargeRightIdx && i > LargeRightIdx      # p₀ → pᵢ
        return is_point_higher(p, pts[i]) ? 1 : -1
    elseif j == LargeLeftIdx && i > LargeRightIdx   # p₋₁ → pᵢ
        return is_point_lower(p, pts[i]) ? 1 : -1
    elseif i == LargeRightIdx && j > LargeRightIdx  # p₀ → pᵢ
        return is_point_lower(p, pts[j]) ? 1 : -1
    elseif i == LargeLeftIdx && j > LargeRightIdx   # p₋₁ → pᵢ
        return is_point_higher(p, pts[j]) ? 1 : -1
    elseif i == LargeRightIdx && j == LargeLeftIdx  # p₀ → p₋₁
        return -1
    elseif i == LargeLeftIdx && j == LargeRightIdx  # p₋₁ → p₀
        return 1
    end
    return leftofline(pts[i], pts[j], p) # If a line is from i → j, then k is left of i → j is (i, j, k) is a counter-clockwise triangle
end

"""
    intriangle(𝒯::TriangleType, pts, p)

Tests if the point `p` is in the triangle `𝒯`, where the vertices of 
`𝒯 = (i, j, k)` are `(pts[i], pts[j], pts[k])`. It is assumed that `𝒯` is 
positively oriented.
"""
function intriangle(e1, e2, e3) # https://stackoverflow.com/a/2049593
    if e1 == 0 || e2 == 0 || e3 == 0
        return 0
    end
    ininterior = e1 > 0 && e2 > 0 && e3 > 0
    if ininterior
        return 1
    else
        return -1
    end
end
function intriangle(𝒯::TriangleType, pts, p)
    i, j, k = 𝒯
    e1 = leftofline(pts, p, i, j)
    e2 = leftofline(pts, p, j, k)
    e3 = leftofline(pts, p, k, i)
    return intriangle(e1, e2, e3)
end

"""
    edge_on_large_triangle(i, j)

Returns true if `(i, j)` is an edge of the triangle `(1, -1, 0)`.
"""
function edge_on_large_triangle(i, j)
    if i > (LargeRightIdx + 1) || j > (LargeRightIdx + 1) # = 1 case can be p₀ 
        return false
    elseif (i, j) == (LargeRightIdx + 1, LargeRightIdx) ||
           (i, j) == (LargeRightIdx, LargeLeftIdx) ||
           (i, j) == (LargeLeftIdx, LargeRightIdx + 1) ||
           (i, j) == (LargeRightIdx, LargeRightIdx + 1) ||
           (i, j) == (LargeLeftIdx, LargeRightIdx) ||
           (i, j) == (LargeRightIdx + 1, LargeLeftIdx)
        return true
    else
        return false
    end
end

"""
    is_legal(i, j, 𝒜, pts)

Tests if the edge `(i, j)` is a legal edge. `𝒜` is the adjacency list of the triangulation, and `pts` is the point set.
Returns `true` if the edge is legal.
"""
function is_legal(i, j, k, ℓ, pts)
    if i > LargeRightIdx && j > LargeRightIdx && k > LargeRightIdx && ℓ > LargeRightIdx
        return incircle(pts, i, j, k, ℓ) ≤ 0
    else
        return min(k, ℓ) < min(i, j)
    end
end
function is_legal(i, j, 𝒜::Adjacent, pts)
    edge_on_large_triangle(i, j) && return true
    k, ℓ = 𝒜(i, j), 𝒜(j, i)
    return is_legal(i, j, k, ℓ, pts)
end

