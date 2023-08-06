"""
    opposite_angle_is_obtuse(p, q, r)

Returns `Certificate.Obtuse` if the angle opposite to the edge `pq` is obtuse, `Certificate.Right` if 
it is a right angle, and `Certificate.Acute` otherwise, letting the vertex opposite to `pq` be `r`.

The test is done by checking if `dot(p - r, q - r) ≤ 0`, where `< 0` means obtuse, `> 0` means acute, and 
`== 0` means right.

!!! note 

    While it would be easy to define a predicate for this that is robust via ExactPredicates.jl's codegen, 
    we do not employ this here.
"""
function opposite_angle_is_obtuse(p, q, r) # https://math.stackexchange.com/a/701682/861404 
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    rx, ry = _getxy(r)
    d = (px - rx) * (qx - rx) + (py - ry) * (qy - ry)
    if d < zero(d)
        return Certificate.Obtuse
    elseif d > zero(d)
        return Certificate.Acute
    else
        return Certificate.Right
    end
end

"""
    point_position_relative_to_diametral_circle(p, q, r)

Given an edge `pq` and a point `r`, returns the position of `r` relative to the diametral circle of `pq`. This test 
is done by noting that a point is in the diametral circle if, and only if, the angle at `r` is obtuse. Returns 

- `Certificate.Inside`: `r` is inside the diametral circle.
- `Certificate.On`: `r` is on the diametral circle.
- `Certificate.Outside`: `r` is outside the diametral circle.

!!! warning 

    This predicate does not use ExactPredicates.jl.
"""
function point_position_relative_to_diametral_circle(p, q, r)
    d = opposite_angle_is_obtuse(p, q, r)
    if is_obtuse(d)
        return Certificate.Inside
    elseif is_right(d)
        return Certificate.On
    else
        return Certificate.Outside
    end
end

#=
Need to eventually get around to implementing the diametral lens definition 
of encroachment. Just not yet, sorry. I think the check is something like 
dot(p-r, q-r)² ≥ (2cos²(θₘᵢₙ) - 1)² * norm(p-r)^2 * norm(q-r)^2, with 
θₘᵢₙ the user-specified angle, but I've not got around to confirming.
function point_position_relative_to_diametral_lens(p, q, r)
end
=#

"""
    encroaches_upon(p, q, r)

Returns `true` if the point `r` encroaches upon the edge `pq`, and `false` otherwise. Here, encroachment means 
`r` is in the diametral circle of `pq`, or, equivalently, the angle opposite to `pq` is obtuse. 

This will eventually be extended to allow for a user to use the diametral lens definition of encroachment,
but for now, this is the only definition.

!!! warning 

    This predicate does not use ExactPredicates.jl. 
"""
function encroaches_upon(p, q, r)
    d = point_position_relative_to_diametral_circle(p, q, r)
    return !is_outside(d)
end

"""
    is_encroached(tri::Triangulation, e)

Tests if the edge `e` in the triangulation `tri` is encroached. See also [`encroaches_upon`](@ref).
"""
function is_encroached(tri::Triangulation, e)
    e′ = reverse_edge(e)
    if !edge_exists(tri, e) && !edge_exists(tri, e′)
        return true
    elseif is_ghost_edge(e)
        return false
    else
        u, v = edge_indices(e)
        p, q = get_point(tri, u, v)
        e_opp = get_adjacent(tri, e)
        if !is_boundary_index(e_opp)
            r = get_point(tri, e_opp)
            encroaches_upon(p, q, r) && return true
        end
        e′_opp = get_adjacent(tri, e′)
        if !is_boundary_index(e′_opp)
            r = get_point(tri, e′_opp)
            encroaches_upon(p, q, r) && return true
        end
        return false
    end
end

"""
compute_concentric_shell_ternary_split_position(p, q)
compute_concentric_shell_quarternary_split_position(p, q)

Computes the position of the split point for a concentric shell split of the edge `pq`. Following 
Ruppert's original paper [here](https://www.nas.nasa.gov/assets/pdf/techreports/1994/rnr-94-002.pdf)
and the book by Cheng, Dey, and Shewchuk.
"""
function compute_concentric_shell_ternary_split_position(p, q)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    ℓ² = (px - qx)^2 + (py - qy)^2
    ℓ = sqrt(ℓ²)
    #=
    The idea is to split a subsegment not at its midpoint, but at one of the circular shells 
    surrounding the shared vertex, assuming the encroached segment is adjacent to another 
    segment. As recommended by Cheng, Dey, and Shewchuk:
        Choose the shell that gives the best-balanced split, so the two new subsegments 
        produced by the split are between one third and two thirds the length of the 
        split subsegment.
    So, the power of 2 for the shell's radius must somehow define a split between one third 
    and two thirds of ℓ. Thus, the relevant quantity is not log₂(ℓ) as it is in Ruppert's 
    original paper, but log₂(ℓ/3) and log₂(ℓ/1.5). We could do e.g. 
        x = floor(log2(ℓ/3))
        y = floor(log2(ℓ/1.5))
        if abs(x - log2(ℓ)) < abs(y - log2(ℓ))
            return exp2(x)
        else 
            return exp2(y)
        end
    but this is slow. Instead, note that this is the same as finding the first power of 2  
    exceeding ℓ/3, and then adjusting it so that we get the nearest power of 2 not exceeding 
    two-thirds of ℓ, i.e. ℓ/1.5. So, let's just do a loop instead.
    =#
    s = balanced_power_of_two_ternary_split(ℓ)
    t = s / ℓ # We want the split to a distance s from p, but our line parametrisation uses p + t(q-p), t ∈ [0, 1]. Thus, normalise s.
    return t
end
function compute_concentric_shell_quarternary_split_position(p, q)
    px, py = _getxy(p)
    qx, qy = _getxy(q)
    ℓ² = (px - qx)^2 + (py - qy)^2
    ℓ = sqrt(ℓ²)
    s = balanced_power_of_two_ternary_split(ℓ)
    t = s / ℓ
    return t
end