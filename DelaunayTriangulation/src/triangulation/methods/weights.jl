"""
    add_weight!(weights, w) 

Pushes the weight `w` into `weights`. The default definition for this is `push!(weights, w)`.
"""
add_weight!(weights, w) = push!(weights, w)

"""
    get_weight(weights, i) -> Number

Gets the `i`th weight from `weights`. The default definition for this is `weights[i]`,
but this can be extended - e.g., [`ZeroWeight`](@ref) uses `get_weight(weights, i) = 0.0`.

If `i` is not an integer, then the default definition is `is_point3(i) ? i[3] : zero(number_type(weights))`.
"""
get_weight(weights, i::Integer) = weights[i]
get_weight(weights, i) = is_point3(i) ? i[3] : zero(number_type(weights))

"""
    ZeroWeight{T}()

Struct used for indicating that a triangulation has zero weights. The weights have type `T`,
which is `Float64` by default.
"""
struct ZeroWeight{T} end
ZeroWeight() = ZeroWeight{Float64}()
get_weight(weights::ZeroWeight{T}, i) where {T} = is_point3(i) ? i[3] : zero(T)
get_weight(weights::ZeroWeight{T}, i::Integer) where {T} = zero(T)

Base.copy(w::ZeroWeight) = w

"""
    add_weight!(tri::Triangulation, w)

Pushes the weight `w` into the weights of `tri`.
"""
function add_weight!(tri::Triangulation, w)
    weights = get_weights(tri)
    add_weight!(weights, w)
    return tri
end

"""
    get_weight(tri::Triangulation, i) -> Number

Gets the `i`th weight from `tri`.
"""
function get_weight(tri::Triangulation, i)
    weights = get_weights(tri)
    return get_weight(weights, i)
end
function get_weight(tri::Triangulation, i::Integer) # ambiguity
    weights = get_weights(tri)
    return get_weight(weights, i)
end

"""
    is_weighted(weights) -> Bool

Returns `true` if `weights` represents a set of `weights` that are not all zero, and `false` otherwise.
Note that even for vectors like `zeros(n)`, this will return `true`; by default, `false` is returned only for 
`weights = ZeroWeight()`.
"""
is_weighted(weights::ZeroWeight) = false
is_weighted(weights) = true

"""
    get_lifted_point(p, w) -> NTuple{3, Number}

Returns the lifted companion of the point `p`, in particular `(x, y, x^2 + y^2 - w)`, where `(x, y)` is `p`.
"""
function get_lifted_point(p, w)
    x, y = getxy(p)
    # The coordinates are promoted to the same type, since the predicates need it, e.g. for Float32 points with Float64 weights
    return promote(x, y, norm_sqr((x, y)) - w)
end

"""
    get_lifted_point(tri::Triangulation, i) -> NTuple{3, Number}

Returns the lifted companion of the `i`th vertex of `tri`, in particular `(x, y, x^2 + y^2 - w)`, where `w` is the
`i`th weight of `tri` and `(x, y)` is the `i`th point of `tri`.
"""
function get_lifted_point(tri::Triangulation, i)
    p = get_point(tri, i)
    w = get_weight(tri, i)
    return get_lifted_point(p, w)
end

"""
    get_lifted_point(p, w, o) -> NTuple{3, Number}
    get_lifted_point(tri::Triangulation, i, o) -> NTuple{3, Number}

Returns the lifted companion of the point `p` with weight `w`, or of the `i`th vertex of `tri`, lifted relative to the point `o`,
in particular `(x, y, (x - ox)^2 + (y - oy)^2 - w)`, where `(x, y)` is `p` and `(ox, oy)` is `o`. See also [`get_lifting_origin`](@ref).
"""
function get_lifted_point(p, w, o)
    x, y = getxy(p)
    ox, oy = getxy(o)
    # Only the height is relative to o. Translating x and y too would round them, which could change the orientation of nearly collinear points
    return promote(x, y, norm_sqr((x - ox, y - oy)) - w)
end
function get_lifted_point(tri::Triangulation, i, o)
    p = get_point(tri, i)
    w = get_weight(tri, i)
    return get_lifted_point(p, w, o)
end

"""
    get_lifting_origin(tri::Triangulation) -> NTuple{2, Number}

Returns the point that the vertices of `tri` are lifted relative to in [`point_position_relative_to_witness_plane`](@ref). This
is the first point of `tri`, or the origin if that point is not finite, e.g. if it is a skipped point.

# Extended help
Lifting each point relative to a fixed point `o`, i.e. to the height `(x - ox)^2 + (y - oy)^2 - w` rather than `x^2 + y^2 - w`, changes
the heights by an affine function of `x` and `y`, which does not change the orientation of any four lifted points. Far from the origin,
though, it keeps the heights small, so that rounding them does not lose the weights or the distances between nearby points: near
`(1e5, 1e5)`, for example, the heights relative to the origin are about `2e10`, and are rounded by about `1e-6`. The coordinates `x` and
`y` themselves are not changed, so that the lifted points project exactly onto the points. Since `o` is the same for every test, the
tests are all exact for the same lifted points, and so they are consistent with each other.
"""
function get_lifting_origin(tri::Triangulation)
    x, y = getxy(get_point(tri, first(each_point_index(tri))))
    return isfinite(x) && isfinite(y) ? (x, y) : (zero(x), zero(y))
end

"""
    get_power_distance(tri::Triangulation, i, j) -> Number

Returns the power distance between vertices `i` and `j`, defined by 
`||pᵢ - pⱼ||^2 - wᵢ - wⱼ`, where `wᵢ` and `wⱼ` are the respective weights.
"""
function get_power_distance(tri::Triangulation, i, j)
    dij² = edge_length_sqr(tri, i, j)
    wᵢ = get_weight(tri, i)
    wⱼ = get_weight(tri, j)
    return dij² - wᵢ - wⱼ
end

"""
    get_distance_to_witness_plane([kernel::AbstractPredicateKernel = AdaptiveKernel(), ] tri::Triangulation, i, V; cache = nothing) -> Number

Computes the distance between the lifted companion of the vertex `i` and the witness plane to the triangle `V`. If `V` is a ghost triangle 
and `i` is not on its solid edge, then the distance is `-Inf` if it is below the ghost triangle's witness plane and `Inf` if it is above. If `V` is a ghost triangle and `i` 
is on its solid edge, then the distance returned is the distance associated with the solid triangle adjoining `V`.

In general, the distance is positive if the lifted vertex is above the witness plane, negative if it is below, 
and zero if it is on the plane.
    
The `kernel` argument determines how this result is computed, and should be 
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.

The `cache` keyword argument is passed to [`point_position_relative_to_circumcircle`]. Please see the documentation for that function for more information.

See also [`point_position_relative_to_witness_plane`](@ref) and [`get_distance_to_plane`](@ref).
"""
function get_distance_to_witness_plane(kernel::AbstractPredicateKernel, tri::Triangulation, i, V; cache::PredicateCacheType = nothing)
    if !is_ghost_triangle(V)
        u, v, w = triangle_vertices(V)
        p⁺ = get_lifted_point(tri, u)
        q⁺ = get_lifted_point(tri, v)
        r⁺ = get_lifted_point(tri, w)
        s⁺ = get_lifted_point(tri, i)
        return get_distance_to_plane(p⁺, q⁺, r⁺, s⁺)
    else
        cert = point_position_relative_to_circumcircle(kernel, tri, V, i; cache)
        if is_inside(cert)
            return -Inf
        elseif is_outside(cert)
            return Inf
        else # is_on(cert) 
            V′ = replace_ghost_triangle_with_boundary_triangle(tri, V)
            return get_distance_to_witness_plane(kernel, tri, i, V′; cache)
        end
    end
end
get_distance_to_witness_plane(tri::Triangulation, i, V; cache::PredicateCacheType = nothing) = get_distance_to_witness_plane(AdaptiveKernel(), tri, i, V; cache)


"""
    get_weighted_nearest_neighbour(tri::Triangulation, i, j = rand(each_solid_vertex(tri))) -> Vertex 

Using a greedy search, finds the closest vertex in `tri` to the vertex `i` (which might not already be in `tri`), 
measuring distance in lifted space (i.e., using the power distance - see [`get_power_distance`](@ref)). 
The search starts from the vertex `j` which should be in `tri`. 
"""
function get_weighted_nearest_neighbour(tri::Triangulation, i, j = rand(each_solid_vertex(tri)))
    if has_vertex(tri, i)
        return i
    else
        return _get_weighted_nearest_neighbour(tri, i, j)
    end
end
function _get_weighted_nearest_neighbour(tri, i, j)
    # This is a loop rather than recursion, since the walk can take as many steps as there are vertices
    while true
        min_δ² = get_power_distance(tri, i, j)
        min_j = j
        for k in get_neighbours(tri, j)
            if !is_ghost_vertex(k)
                δ² = get_power_distance(tri, i, k)
                if δ² < min_δ²
                    min_δ² = δ²
                    min_j = k
                end
            end
        end
        min_j == j && return j
        j = min_j
    end
end

@doc """
    is_submerged([kernel::AbstractPredicateKernel = AdaptiveKernel(), ] tri::Triangulation, i; cache = nothing) -> Bool 
    is_submerged([kernel::AbstractPredicateKernel = AdaptiveKernel(), ] tri::Triangulation, i, V; cache = nothing) -> Bool

Returns `true` if the vertex `i` is submerged in `tri` and `false` otherwise. In the 
second method, `V` is a triangle containing `tri`.

The `kernel` argument determines how this result is computed, and should be
one of [`ExactKernel`](@ref), [`FastKernel`](@ref), and [`AdaptiveKernel`](@ref) (the default).
See the documentation for more information about these choices.

The `cache` keyword argument is passed to [`point_position_relative_to_circumcircle`]. Please see the documentation for that function for more information.
"""
is_submerged
function is_submerged(kernel::AbstractPredicateKernel, tri::Triangulation, i; cache::PredicateCacheType = nothing)
    # A source that mentions that testing if `i` is submerged only needs to consider the triangle that contains it 
    # is given in https://otik.uk.zcu.cz/bitstream/11025/21574/1/Zemek.pdf on p.17. 
    # (If the link dies, it is the PhD thesis of Michal Zemek, "Regular Triangulation in 3D and Its Applications".)
    is_ghost_vertex(i) && return false
    q = get_point(tri, i)
    V = find_triangle(tri, q)
    return is_submerged(kernel, tri, i, V; cache)
end
function is_submerged(kernel::AbstractPredicateKernel, tri::Triangulation, i, V; cache::PredicateCacheType = nothing)
    is_ghost_vertex(i) && return false
    cert = point_position_relative_to_circumcircle(kernel, tri, V, i; cache)
    return is_outside(cert)
end
is_submerged(tri::Triangulation, i; cache::PredicateCacheType = nothing) = is_submerged(AdaptiveKernel(), tri, i; cache)
is_submerged(tri::Triangulation, i, V; cache::PredicateCacheType = nothing) = is_submerged(AdaptiveKernel(), tri, i, V; cache)