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

function Base.isempty(queue::RefinementQueue)
    return encroachment_queue_is_empty(queue) && triangle_queue_is_empty(queue)
end

peek_triangle_ρ(queue::RefinementQueue) = peek(queue.triangle_queue)[2]

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
    end
    for e in each_constrained_edge(tri)
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
    return queue
end