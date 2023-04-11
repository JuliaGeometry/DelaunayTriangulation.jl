"""
    convert_minimum_angle(θ)

Converts the user-specified minimum angle `θ` to `sin(ψ)^2`, `ψ = rad2deg(θ)`.
"""
convert_angle(θ) = sind(θ)^2

"""
    RefinementTargets{A,B,C,D,E}

A struct containing the user-specified refinement targets. 

# Fields

- `min_area=0.0` (not currently used)

The minimum area of a triangle. This can also be a function of the form `f(p, q, r, A)`, where `A` is a triangle's area with coordinates `p`, `q`, `r`, returning `true` if the triangle should be refined.
- `max_area=Inf`

The maximum area of a triangle. This can also be a function of the form `f(p, q, r, A)`, where `A` is a triangle's area with coordinates `p`, `q`, `r`, returning `true` if the triangle should be refined.
- `min_angle=20.0`

The minimum angle of a triangle. While the user should provide this as an angle in degrees between 0 and 60,
the internal representation is `sin(ψ)^2`, where `ψ = rad2deg(min_angle)`. This can also be a function of the form `f(p, q, r, θ)`, where `θ` is a triangle's minimum angle with coordinates `p`, `q`, `r`, returning `true` if the triangle should be refined.
- `max_angle=180.0` (not currently used)

The maximum angle of a triangle. While the user should provide this as an angle in degrees between 60 and 180,
the internal representation is `sin(ψ)^2`, where `ψ = rad2deg(max_angle)`. This can also be a function of the form `f(p, q, r, θ)`, where `θ` is a triangle's maximum angle with coordinates `p`, `q`, `r`, returning `true` if the triangle should be refined.
- `max_points=Inf`

The maximum number of points in the mesh.
"""
struct RefinementTargets{A,B,C,D,E}
    min_area::A
    max_area::B
    min_angle::C
    max_angle::D
    max_points::E
    function MeshQualities(;
        min_area::A=zero(Float64),
        max_area::B=typemax(Float64),
        min_angle::C=20.0,
        max_angle::D=180.0,
        max_points::E=typemax(Int64)) where {A,B,C,D,E}
        if min_area isa Number && min_area < 0.0
            @warn "The provided minimum area constraint, $min_area, is negative. Replacing with zero."
            min_area = 0.0
        end
        if max_area isa Number && (max_area < 60.0 || max_area > 180.0)
            @warn "The provided maximum area constraint, $max_area, is outside the range [60, 180]. Replacing with 180. Note that areas must be provided in degrees."
            max_area = 180.0
        end
        if min_angle isa Number && (min_angle < 0.0 || min_angle > 60.0)
            @warn "The provided min_angle, $min_angle, is outside the range [0, 60]. Replacing with 20. Note that angles must be provided in degrees."
            min_angle = 20.0
            if min_angle > 33.9
                @warn "The algorithm may fail to halt with a minimum angle constraint of 33.9° or higher. Consider using a smaller minimum angle constraint. Will proceed with the provided minimum angle constraint."
            end
        end
        if max_angle isa Number && (max_angle < 60.0 || max_angle > 180.0)
            @warn "The provided max_angle, $max_angle, is outside the range [60, 180]. Replacing with 180. Note that angles must be provided in degrees."
            max_angle = 180.0
        end
        if max_points isa Number && max_points < 0
            throw(ArgumentError("The provided maximum number of points, $max_points, is negative."))
        end
        min_angle = convert_angle(min_angle)
        max_angle = convert_angle(max_angle)
        return new{A,B,C,D,E}(min_area, max_area, min_angle, max_angle, max_points)
    end
end
function compare_area(targets::RefinementTargets, A, p, q, r)
    min_flag = if targets.min_area isa Function
        targets.min_area(p, q, r, A)
    else
        A < targets.min_area
    end
    max_flag = if targets.max_area isa Function
        targets.max_area(p, q, r, A)
    else
        A > targets.max_area
    end
    return (min_flag, max_flag)
end
function compare_angle(targets::RefinementTargets, sinθₘᵢₙ², sinθₘₐₓ², p, q, r)
    min_flag = if targets.min_angle isa Function
        targets.min_angle(p, q, r, sinθₘᵢₙ²)
    else
        sinθₘᵢₙ² < targets.min_angle
    end
    max_flag = if targets.max_angle isa Function
        targets.max_angle(p, q, r, sinθₘₐₓ²)
    else
        sinθₘₐₓ² > targets.max_angle
    end
    return (min_flag, max_flag)
end
function compare_points(targets::RefinementTargets, n)
    return n > targets.max_points
end

"""
    assess_triangle_quality(tri::Triangulation, T, targets::RefinementTargets)

Assess the quality of a triangle `T` in a `Triangulation` `tri` with respect to the provided `RefinementTargets`. 

Returns `(sinθₘᵢₙ², flag)`, where `sinθₘᵢₙ²` is the squared sine of the minimum angle of the triangle and `flag` is `true` if the triangle should be refined.
"""
function assess_triangle_quality(tri::Triangulation, T, targets::RefinementTargets)
    u, v, w = indices(T)
    p, q, r = get_point(tri, u, v, w)
    ℓmin², ℓmed², ℓmax² = squared_triangle_lengths(p, q, r)
    A = triangle_area(ℓmin², ℓmed², ℓmax²)
    sinθₘᵢₙ² = triangle_sine_minimum_angle_squared(A, ℓmin², ℓmed²)
    sinθₘₐₓ² = triangle_sine_maximum_angle_squared(A, ℓmed², ℓmax²)
    min_area_flag, max_area_flag = compare_area(targets, A, p, q, r)
    min_angle_flag, max_angle_flag = compare_angle(targets, sinθₘᵢₙ², sinθₘₐₓ², p, q, r)
    return sinθₘᵢₙ², (min_area_flag || max_area_flag || min_angle_flag || max_angle_flag)
end

"""
    RefinementQueue{T,E,F}

A queue for storing encroachment and triangle refinement priority queues. 

# Fields

- `encroachment_queue::PriorityQueue{T,E,F,Base.Order.ForwardOrdering}`

A priority queue for storing encroached segments to be split. The keys are the squared edge lengths.
- `triangle_queue::PriorityQueue{T,E,F,Base.Order.ForwardOrdering}`

A priority queue for storing triangles to be refined. The keys are the squared sine of the minimum angle of the triangle.
"""
struct RefinementQueue{T,E,F}
    encroachment_queue::PriorityQueue{T,E,F,Base.Order.ForwardOrdering}
    triangle_queue::PriorityQueue{T,E,F,Base.Order.ForwardOrdering}
    function RefinementQueue{T,E,F}() where {T,F}
        return new{T,E,F}(
            PriorityQueue{T,E,Base.Order.ForwardOrdering}(Base.Order.Forward),
            PriorityQueue{T,F,Base.Order.ForwardOrdering}(Base.Order.Forward))
    end
end

function encroachment_enqueue!(queue::RefinementQueue, e, e_length²)
    encroachment_queue = queue.encroachment_queue
    existing_segments = keys(encroachment_queue)
    if e ∈ existing_segments
        existing_length² = encroachment_queue[e]
        if e_length² < existing_length²
            encroachment_queue[e] = e_length²
        end
    elseif reverse_edge(e) ∈ existing_segments
        existing_length² = encroachment_queue[reverse_edge(e)]
        if e_length² < existing_length²
            encroachment_queue[reverse_edge(e)] = e_length²
        end
    else
        enqueue!(encroachment_queue, e, e_length²)
    end
    return nothing

end

function triangle_enqueue!(queue::RefinementQueue, T, sinθₘᵢₙ²)
    triangle_queue = queue.triangle_queue
    existing_triangles = keys(triangle_queue)
    T, flag = contains_triangle(T, existing_triangles)
    if !flag
        enqueue!(triangle_queue, T, sinθₘᵢₙ²)
    else
        existing_sinθₘᵢₙ² = triangle_queue[T]
        if sinθₘᵢₙ² < existing_sinθₘᵢₙ²
            triangle_queue[T] = sinθₘᵢₙ²
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
        sinθₘᵢₙ², refine_flag = assess_triangle_quality(tri, T, targets)
        refine_flag && triangle_enqueue!(queue, T, sinθₘᵢₙ²)
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
