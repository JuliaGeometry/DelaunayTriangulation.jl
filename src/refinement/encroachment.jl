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
    px, py = getxy(p)
    qx, qy = getxy(q)
    rx, ry = getxy(r)
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