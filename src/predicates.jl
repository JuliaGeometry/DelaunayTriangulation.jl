############################################
##
## PREDICATES 
##
############################################
"""
    ExactPredicates.orient(T::AbstractTriangle, pts::Points)

Tests if the triangle `T` is positively oriented.
"""
function ExactPredicates.orient(T::AbstractTriangle, pts::Points)
    u, v, w = T
    pu, pv, pw = get_point(pts, u), get_point(pts, v), get_point(pts, w)
    return orient(pu, pv, pw)
end

"""
    ExactPredicates.incircle(pts::Points, i, j, k, ℓ)

Tests if `get_point(pts, ℓ)` is in the circle through `(get_point(pts, i), get_point(pts, j), get_point(pts, k))`.
"""
function ExactPredicates.incircle(pts::Points, i, j, k, ℓ)
    pti, ptj, ptk, ptℓ = get_point(pts, i), get_point(pts, j), get_point(pts, k), get_point(pts, ℓ)
    return incircle(pti, ptj, ptk, ptℓ)
end

"""
    leftofline(p, pᵢ, pⱼ)

Tests if `p` is to the left of the line from `pᵢ` to `pⱼ`.
"""
leftofline(p::AbstractPoint, pᵢ::AbstractPoint, pⱼ::AbstractPoint) = orient(pᵢ, pⱼ, p)
"""
    leftofline(pts, p, i, j)

Tests if the point `p` is to the left of the oriented line through 
`pts[i]` to `pts[j]`. Checks are made for non-positive indices.
"""
function leftofline(pts::Points, p::AbstractPoint, i, j)
    if i == LowerRightBoundingIndex && j == LowerLeftBoundingIndex
        return -1
    elseif i == LowerRightBoundingIndex && j == UpperBoundingIndex
        return 1
    elseif i == LowerLeftBoundingIndex && j == LowerRightBoundingIndex
        return 1
    elseif i == LowerLeftBoundingIndex && j == UpperBoundingIndex
        return -1
    elseif i == UpperBoundingIndex && j == LowerLeftBoundingIndex
        return 1
    elseif i == UpperBoundingIndex && j == LowerRightBoundingIndex
        return -1
    end
    return leftofline(get_point(pts, i), get_point(pts, j), p) # If a line is from i → j, then k is left of i → j is (i, j, k) is a counter-clockwise triangle
end

"""
    intriangle(T::AbstractTriangle, pts::Points, p::AbstractPoint)

Tests if the point `p` is in the triangle `T`, where the vertices of 
`T = (i, j, k)` are `(pts[i], pts[j], pts[k])`. It is assumed that `T` is 
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
function intriangle(T::AbstractTriangle, pts::Points, p::AbstractPoint)
    i, j, k = T
    if i < FirstPointIndex && j < FirstPointIndex && k < FirstPointIndex
        return 1 # all poiints are in the bounding triangle
    end
    e1 = leftofline(pts, p, i, j)
    e2 = leftofline(pts, p, j, k)
    e3 = leftofline(pts, p, k, i)
    return intriangle(e1, e2, e3)
end

"""
    edge_on_bounding_triangle(i, j)

Returns true if `(i, j)` is an edge of bounding triangle $(BoundingTriangle).
"""
function edge_on_bounding_triangle(i, j)
    return i < FirstPointIndex && j < FirstPointIndex
end

"""
    islegal(i, j, k, ℓ, pts::Points)
    islegal(i, j, adj::Adjacent, pts::Points)

Returns `true` if the edge `(i, j)` is legal. It is assumed that 
`(i, j, k)` and `(j, i, ℓ)` are positively oriented triangles.
"""
function islegal(i, j, k, ℓ, pts::Points)
    #=
    if i ≥ FirstPointIndex && j ≥ FirstPointIndex && k ≥ FirstPointIndex && ℓ ≥ FirstPointIndex
        return incircle(pts, i, j, k, ℓ) ≤ 0
    else
        num_neg = num_less(FirstPointIndex, (i, j, k, ℓ))
        if num_neg == 1
            return !(i < FirstPointIndex || j < FirstPointIndex)
        elseif num_neg == 2
            return min(i, j) < min(k, ℓ)
        else
            throw("Error occured.")
        end
    end
    throw("Error occured.")
    =#
    return incircle(pts, i, j, k, ℓ) ≤ 0
end
@doc (@doc islegal(::Any, ::Any, ::Any, ::Any, ::Points))
function islegal(i, j, adj::Adjacent, pts::Points)
    edge_on_bounding_triangle(i, j) && return true
    k = get_edge(adj, i, j)
    ℓ = get_edge(adj, j, i)
    return islegal(i, j, k, ℓ, pts)
end