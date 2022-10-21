"""
    isoriented(T, pts)

Tests if the triangle `T` is positively oriented, returning
`1` if the triangle is positively oriented, `0` 
if the triangle is degenerate (i.e. the points are all 
collinear), and `-1` if the triangle is negatively 
oriented.

The argument `pts` should be a collection of points 
such that the `i`th vertex of `T` corresponds to 
`get_point(pts, i)`.

The check is exact, making use of `ExactPredicates.orient` (unless 
you define a new method for it for your point type, obviously).
"""
function isoriented(T, pts)
    i, j, k = indices(T)
    pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
    return orient(pᵢ, pⱼ, pₖ)
end

"""
    isincircle(pts, i, j, k, ℓ)

Tests if the `ℓ`th point is in the circle through the 
`i`th, `j`th, and `k` points (as obtained via `get_point(pts, i)`, etc.).
Returns `1` if the point is inside, `0` if the point is on the circle,
and `-1` if the point is outside. It is assumed that the 
`i`th, `j`th`, and `k`th points are provided counter-clockwise. 

The check is exact, making use of `ExactPredicates.incircle` (unless 
you define a new method for it for your point type, obviously).

Checks are made for ghost triangles, which are triangles of the form `(i, j, $BoundaryIndex)`. We say 
that a point is in the circumcircle of `(i, j, $BoundaryIndex)` if it is to the left of `(i, j)`.
"""
function isincircle(pts, i::I, j::I, k::I, ℓ::I) where {I}
    if i == I(BoundaryIndex)
        return I(isleftofline(pts, j, k, ℓ))
    elseif j == I(BoundaryIndex)
        return I(isleftofline(pts, k, i, ℓ))
    elseif k == I(BoundaryIndex)
        return I(isleftofline(pts, i, j, ℓ))
    end
    pᵢ, pⱼ, pₖ, pₗ = get_point(pts, i, j, k, ℓ)
    return I(incircle(pᵢ, pⱼ, pₖ, pₗ))
end
"""
    isincircle(T, pts, ℓ) 

Tests if the `ℓ`th point in `pts` is inside the circumdisk of the 
triangle `T`. 
"""
function isincircle(T, pts, ℓ)
    i, j, k = indices(T)
    return isincircle(pts, i, j, k, ℓ)
end

"""
    isleftofline(pts, i, j, k)

Tests if the `k`th point if left of the line from the `i`th point 
to the `j`th point (as obtained via `get_point(pts, i)`, etc.). 
Returns `1` if the point is to the left, `0` if the points are 
collinear, and `-1` if the point is to the right.

The check is exact, making use of `ExactPredicates.orient` (unless 
you define a new method for it for your point type, obviously).
"""
function isleftofline(pts, i::I, j::I, k::I) where {I}
    if i == I(LowerRightBoundingIndex) && j == I(LowerLeftBoundingIndex)
        return I(-1)
    elseif i == I(LowerRightBoundingIndex) && j == I(UpperBoundingIndex)
        return I(1)
    elseif i == I(LowerLeftBoundingIndex) && j == I(LowerRightBoundingIndex)
        return I(1)
    elseif i == I(LowerLeftBoundingIndex) && j == I(UpperBoundingIndex)
        return I(-1)
    elseif i == I(UpperBoundingIndex) && j == I(LowerLeftBoundingIndex)
        return I(1)
    elseif i == I(UpperBoundingIndex) && j == I(LowerRightBoundingIndex)
        return I(-1)
    end
    pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
    return I(orient(pᵢ, pⱼ, pₖ))
end

"""
    isintriangle(T, pts, ℓ)

Tests if the `ℓ`th point of `pts` is in the triangle `T`. Returns `1`
if the point is inside, `0` if a point is on an edge of `T`, and `-1` 
if the point is outside. It is assumed that `T` is positively oriented. 

The check is exact, making use of `isleftofline`. 
"""
function isintriangle(T, pts, ℓ::I) where {I}
    i, j, k = indices(T)
    if i < I(FirstPointIndex) && j < I(FirstPointIndex) && k < I(FirstPointIndex)
        return I(1)
    end
    is_left_of_edge_ij = isleftofline(pts, i, j, ℓ)
    is_left_of_edge_jk = isleftofline(pts, j, k, ℓ)
    is_left_of_edge_ki = isleftofline(pts, k, i, ℓ)
    return isintriangle(
        is_left_of_edge_ij,
        is_left_of_edge_jk,
        is_left_of_edge_ki
    )
end
function isintriangle(e1::I, e2::I, e3::I) where {I}
    if e1 == I(0) || e2 == I(0) || e3 == I(0)
        return I(0)
    end
    isininterior = e1 == I(1) && e2 == I(1) && e3 == I(1)
    if isininterior
        return I(1)
    else
        return I(-1)
    end
end

"""
    find_edge(T, pts, ℓ)

Given a triangle `T` and a set of points `pts`, with the `ℓ`th point of `pts` 
on an edge of `T`, find the edge of `T` that the point is on.
"""
function find_edge(T, pts, ℓ)
    i, j, k = indices(T)
    isleftofline(pts, i, j, ℓ) == 0 && return (i, j)
    isleftofline(pts, j, k, ℓ) == 0 && return (j, k)
    isleftofline(pts, k, i, ℓ) == 0 && return (k, i)
    throw("The point $p is not on an edge of $T.")
end

"""
    edge_on_bounding_triangle(i, j)

Returns true if `(i, j)` is an edge of bounding triangle $(BoundingTriangle).
"""
edge_on_bounding_triangle(i, j) = i < FirstPointIndex && j < FirstPointIndex

"""
    islegal(i, j, adj, pts)

Tests if the edge `(i, j)` is legal. `adj` is the adjacent map of the 
triangulation, and `pts` is the point set.
"""
function islegal(i, j, adj, pts)
    edge_on_bounding_triangle(i, j) && return true
    is_boundary_edge(i, j, adj) && return true
    is_boundary_edge(j, i, adj) && return true
    !edge_exists(i, j, adj) && return true
    !edge_exists(j, i, adj) && return true
    k = get_edge(adj, i, j)
    ℓ = get_edge(adj, j, i)
    return islegal(i, j, k, ℓ, pts)
end
function islegal(i::I, j::I, k::I, ℓ::I, pts) where {I}
    incirc = isincircle(pts, i, j, k, ℓ)
    return incirc == I(-1) || incirc == I(0)
end