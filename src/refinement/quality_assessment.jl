"""
    RefinementTargets{A,R,P}

A struct containing the user-specified refinement targets. 

# Fields

- `max_area::A`

The maximum area of a triangle. This can also be a function of the form `f(T, p, q, r, A)`, where `T` is the triangle` with area and coordinates `p`, `q`, `r`, returning `true` if the triangle should be refined.

- `max_radius_edge_ratio::R`

The maximum permitted radius-edge ratio. This can also be a function of the form `f(T, p, q, r, ρ)`, where `T` is the triangle` with ratio `ρ` and coordinates `p`, `q`, `r`, returning `true` if the triangle should be refined. Defaults to 1.0, corresponding to a minimum angle of 30°.
- `max_points::P`

The maximum number of points in the mesh. Defaults to Inf.
# Constructors 

The constructor is 

    RefinementTargets(;
            max_area=typemax(Float64),
            max_radius_edge_ratio=nothing,
            max_points=typemax(Int64),
            min_angle = nothing
            )

which allows for a user to specify either the maximum radius-edge ratio `ρ` or the minimum angle `min_angle` in degrees, using the relationship `min_angle = asin[1/(2ρ)]`. If both are provided, the minimum angle is ignored. If neither are provided, the default value of `ρ = 1.0` is used, corresponding to a minimum angle of 30°.

Note that you cannot use `min_angle` as a function, unlike `max_radius_edge_ratio`. If you want to use a function, use `max_radius_edge_ratio` instead.
"""
struct RefinementTargets{A,R,P}
    max_area::A
    max_radius_edge_ratio::R
    max_points::P
    function RefinementTargets(;
        max_area=typemax(Float64),
        max_radius_edge_ratio=nothing,
        max_points=typemax(Int64),
        min_angle = nothing
        )
        if min_angle isa Function 
            throw(ArgumentError("Cannot provide min_angle as a function."))
        end
        if max_radius_edge_ratio isa Number && max_radius_edge_ratio < sqrt(3)/3
            @warn "The provided maximum radius-edge ratio, $max_radius_edge_ratio, is below the theoretical limit of 1/√3. Replacing it with 1/sqrt(3)."
            max_radius_edge_ratio = sqrt(3)/3
        end
        if max_radius_edge_ratio === nothing && min_angle === nothing
            max_radius_edge_ratio = 1.0
        elseif max_radius_edge_ratio === nothing && min_angle !== nothing
            max_radius_edge_ratio = cscd(min_angle) / 2 
        elseif max_radius_edge_ratio !== nothing && min_angle !== nothing
            @warn "Both max_radius_edge_ratio and min_angle are provided. Ignoring min_angle."
        end
        min_angle = max_radius_edge_ratio isa Number ? asind(1/(2max_radius_edge_ratio)) : nothing 
        if !isnothing(min_angle) && (33.9 < min_angle < 60)
            @warn "The provided max_radius_edge_ratio, ρ = $max_radius_edge_ratio, corresponds to a minimum angle of $(min_angle)°. The algorithm may fail to halt with a minimum angle constraint exceeding 33.9° (corresponding to ρ <
             0.9) or higher (meaning lower ρ). Consider using a smaller minimum angle constraint or a larger maximum radius-edge ratio. Will proceed with the provided maximum radius-edge ratio."
        end
        if max_area isa Number && max_area < 0.0
            @warn "The provided maximum area constraint, $max_area, is negative. Replacing with Inf."
            max_area = typemax(Float64)
        end
        if !isnothing(min_angle) && (min_angle < 0.0 || min_angle > 60.00005) # 60.00005 for the =sqrt(3)/3 case
            @warn "The provided max_radius_edge_ratio, ρ = $max_radius_edge_ratio, corresponds to a minimum angle of $(min_angle)°, which is outside the range [0, 60]. Replacing ρ with 1.0, corresponding to a minimum angle of 30°."
            max_radius_edge_ratio = 1.0
        end
        if max_points isa Number && max_points < 0
            @warn "The provided maximum number of points, $max_points, is negative. Replacing with Inf."
            max_points = typemax(Int64)
        end
        return new{typeof(max_area), typeof(max_radius_edge_ratio), typeof(max_points)}(max_area, max_radius_edge_ratio, max_points)
    end
end
function compare_area(targets::RefinementTargets, T, A, p, q, r)
    flag = if targets.max_area isa Function
        targets.max_area(T, p, q, r, A)
    else
        A > targets.max_area
    end
    return flag
end
function compare_ratio(targets::RefinementTargets, T, ρ, p, q, r)
    flag = if targets.max_radius_edge_ratio isa Function
        targets.max_radius_edge_ratio(T, p, q, r, ρ)
    else
        ρ > targets.max_radius_edge_ratio
    end
    return flag
end
function compare_points(targets::RefinementTargets, n)
    return n > targets.max_points
end

"""
    assess_triangle_quality(tri::Triangulation, T, targets::RefinementTargets)

Assess the quality of a triangle `T` in a `Triangulation` `tri` with respect to the provided `RefinementTargets`. 

Returns `(ρ, flag)`, where `ρ` is the radius-edge ratio of the triangle and `flag` is `true` if the triangle should be refined.
"""
function assess_triangle_quality(tri::Triangulation, T, targets::RefinementTargets)
    u, v, w = indices(T)
    p, q, r = get_point(tri, u, v, w)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    r = triangle_circumradius(A, ℓmin², ℓmed², ℓmax²)
    ρ = r / sqrt(ℓmin²)
    area_flag = compare_area(targets, T, A, p, q, r)
    ratio_flag = compare_ratio(targets, T, ρ, p, q, r)
    return ρ, (area_flag || ratio_flag)
end

"""
    RefinementQueue{T,E,F}

A queue for storing encroachment and triangle refinement priority queues. 

# Fields

- `encroachment_queue::PriorityQueue{T,E,F,Base.Order.ReverseOrdering}`

A priority queue for storing encroached segments to be split. The keys are the squared edge lengths.
- `triangle_queue::PriorityQueue{T,E,F,Base.Order.ReverseOrdering}`

A priority queue for storing triangles to be refined. The keys are radius-edge ratio.
"""
struct RefinementQueue{T,E,F}
    encroachment_queue::PriorityQueue{E,F,Base.Order.ReverseOrdering}
    triangle_queue::PriorityQueue{T,F,Base.Order.ReverseOrdering}
    function RefinementQueue{T,E,F}() where {T,E,F}
        return new{T,E,F}(
            PriorityQueue{E,F,Base.Order.ReverseOrdering}(Base.Order.Reverse),
            PriorityQueue{T,F,Base.Order.ReverseOrdering}(Base.Order.Reverse))
    end
end

function encroachment_enqueue!(queue::RefinementQueue, e, e_length²)
    encroachment_queue = queue.encroachment_queue
    existing_segments = keys(encroachment_queue)
    if e ∈ existing_segments
        existing_length² = encroachment_queue[e]
        if e_length² > existing_length²
            encroachment_queue[e] = e_length²
        end
    elseif reverse_edge(e) ∈ existing_segments
        existing_length² = encroachment_queue[reverse_edge(e)]
        if e_length² > existing_length²
            encroachment_queue[reverse_edge(e)] = e_length²
        end
    else
        enqueue!(encroachment_queue, e, e_length²)
    end
    return nothing
end

function triangle_enqueue!(queue::RefinementQueue, T, ρ)
    triangle_queue = queue.triangle_queue
    existing_triangles = keys(triangle_queue)
    T, flag = contains_triangle(T, existing_triangles)
    if !flag
        enqueue!(triangle_queue, T, ρ)
    else
        existing_ρ = triangle_queue[T]
        if ρ > existing_ρ
            triangle_queue[T] = ρ
        end
    end
    return nothing
end

function encroachment_dequeue!(queue::RefinementQueue)
    encroachment_queue = queue.encroachment_queue
    return dequeue!(encroachment_queue)
end

function triangle_dequeue!(queue::RefinementQueue)
    triangle_queue = queue.triangle_queue
    return dequeue!(triangle_queue)
end

function encroachment_queue_is_empty(queue::RefinementQueue)
    encroachment_queue = queue.encroachment_queue
    return isempty(encroachment_queue)
end

function triangle_queue_is_empty(queue::RefinementQueue)
    triangle_queue = queue.triangle_queue
    return isempty(triangle_queue)
end

"""
    initialise_refinement_queue(tri::Triangulation, targets::RefinementTargets)

Initialise a `RefinementQueue` for a `Triangulation` `tri` with respect to the provided `RefinementTargets`. 
"""
function initialise_refinement_queue(tri::Triangulation, targets::RefinementTargets)
    T = triangle_type(tri)
    F = number_type(tri)
    E = edge_type(tri)
    queue = RefinementQueue{T,E,F}()
    for T in each_solid_triangle(tri)
        ρ, refine_flag = assess_triangle_quality(tri, T, targets)
        refine_flag && triangle_enqueue!(queue, T, ρ)
        for e in triangle_edges(T)
            encroachment_flag = is_encroached(tri, e)
            if encroachment_flag
                u, v = edge_indices(e)
                p, q = get_point(tri, u, v)
                px, py = getxy(p)
                qx, qy = getxy(q)
                e_length² = (px - qx)^2 + (py - qy)^2
                encroachment_enqueue!(queue, e, e_length²)
            end
        end
    end
    return queue
end

"""
    assess_added_triangles!(tri::Triangulation, queue, events, targets)

Assess the quality of the triangles added to the triangulation `tri` after an insertion event,
as stored in `events` and add them to the `queue` if they should be refined according to the mesh 
targets defined in `targets`.
"""
function assess_added_triangles!(tri::Triangulation, queue, events, targets)
    for T in each_added_triangle(events)
        ρ, flag = assess_triangle_quality(tri, T, targets)
        if flag
            triangle_enqueue!(queue, T, ρ)
        end
        for e in triangle_edges(T)
            if is_encroached(tri, e)
                u, v = edge_indices(e)
                p, q = get_point(tri, u, v)
                px, py = getxy(p)
                qx, qy = getxy(q)
                ℓ² = (qx - px)^2 + (qy - py)^2
                encroachment_enqueue!(queue, e, ℓ²)
            end
        end
    end 
    return nothing
end
