"""
    RefinementQueue{T,E,F}

Struct defining a pair of priority queues for encroached segments and poor-quality triangles. 

# Fields 
- `segments::MaxPriorityQueue{E,F}`: A priority queue of encroached segments, where the priorities are the squared edge lengths.
- `triangles::MaxPriorityQueue{T,F}`: A priority queue of poor-quality triangles, where the priorities are the radius-edge ratios.

# Constructor 
The default constructor is available, but we also provide 

    RefinementQueue(tri::Triangulation)

which will initialise this struct with empty queues with the appropriate types.
"""
struct RefinementQueue{T, E, F}
    segments::MaxPriorityQueue{E, F}
    triangles::MaxPriorityQueue{T, F}
    function RefinementQueue{T, E, F}() where {T, E, F}
        return new{T, E, F}(MaxPriorityQueue{E, F}(), MaxPriorityQueue{T, F}())
    end
end
function Base.show(io::IO, ::MIME"text/plain", queue::RefinementQueue)
    println(io, "RefinementQueue")
    println(io, "   $(length(queue.segments)) segments: ", queue.segments)
    print(io, "   $(length(queue.triangles)) triangles: ", queue.triangles)
end

function RefinementQueue(tri::Triangulation)
    T = triangle_type(tri)
    E = edge_type(tri)
    F = number_type(tri)
    return RefinementQueue{T, E, F}()
end

#=
    haskey(queue::RefinementQueue{T,E,F}, segment::E) -> Bool

Return `true` if `queue` has `segment` or its reverse, and `false` otherwise.
=#
function Base.haskey(queue::RefinementQueue{T, E, F}, segment::E) where {T, E, F}
    return haskey(queue.segments, segment) || haskey(queue.segments, reverse_edge(segment))
end

#=
    haskey(queue::RefinementQueue{T,E,F}, triangle::T) -> Bool

Return `true` if `queue` has `triangle` or any of its counter-clockwise rotations, and `false` otherwise.
=#
function Base.haskey(queue::RefinementQueue{T, E, F}, triangle::T) where {T, E, F}
    return haskey(queue.triangles, triangle) ||
        haskey(queue.triangles, rotate_triangle(triangle, Val(1))) ||
        haskey(queue.triangles, rotate_triangle(triangle, Val(2)))
end

#=
    getindex(queue::RefinementQueue{T,E,F}, triangle::T) -> F
    queue[triangle] -> F 

Return the radius-edge ratio of `triangle` in `queue`.
=#
function Base.getindex(queue::RefinementQueue{T, E, F}, triangle::T) where {T, E, F}
    if haskey(queue.triangles, triangle)
        return queue.triangles[triangle]
    elseif haskey(queue.triangles, rotate_triangle(triangle, Val(1)))
        return queue.triangles[rotate_triangle(triangle, Val(1))]
    else #if haskey(queue.triangles, rotate_triangle(triangle, Val(2)))
        return queue.triangles[rotate_triangle(triangle, Val(2))]
    end
end


#=
    setindex!(queue::RefinementQueue{T,E,F}, ℓ²::F, segment::E) where {T,E,F}
    queue[segment] = ℓ²

Add a `segment` to `queue` whose squared length is `ℓ²`. If the `segment` is already in the `queue`, its priority is updated to `ℓ`.
=#
function Base.setindex!(queue::RefinementQueue{T, E, F}, ℓ²::F, segment::E) where {T, E, F}
    segments = queue.segments
    if haskey(segments, reverse_edge(segment))
        segments[reverse_edge(segment)] = ℓ²
    else
        segments[segment] = ℓ²
    end
    return segments
end

#=
    setindex!(queue::RefinementQueue{T,E,F}, ρ::F, triangle::T) where {T,E,F}
    queue[triangle] = ρ

Add a `triangle` to `queue` whose radius-edge ratio is `ρ`. If the `triangle` is already in the `queue`, its priority is updated to `ρ`.
=#
function Base.setindex!(queue::RefinementQueue{T, E, F}, ρ, triangle::T) where {T, E, F}
    triangles = queue.triangles
    if haskey(triangles, triangle)
        triangles[triangle] = ρ
    elseif haskey(triangles, rotate_triangle(triangle, Val(1)))
        triangles[rotate_triangle(triangle, Val(1))] = ρ
    elseif haskey(triangles, rotate_triangle(triangle, Val(2)))
        triangles[rotate_triangle(triangle, Val(2))] = ρ
    else
        triangles[triangle] = ρ
    end
    return triangles
end

"""
    popfirst_segment!(queue::RefinementQueue) -> Edge 

Dequeue the next segment from `queue`, returning the segment and its squared length.
"""
popfirst_segment!(queue::RefinementQueue) = popfirst!(queue.segments)

"""
    popfirst_triangle!(queue::RefinementQueue) -> (Triangle, Number)

Dequeue the next triangle from `queue`, returning the triangle and its radius-edge ratio.
"""
popfirst_triangle!(queue::RefinementQueue) = popfirst!(queue.triangles)

"""
    has_segments(queue::RefinementQueue) -> Bool

Return `true` if `queue` has any segments, `false` otherwise.
"""
has_segments(queue::RefinementQueue) = !isempty(queue.segments)

"""
    has_triangles(queue::RefinementQueue) -> Bool

Return `true` if `queue` has any triangles, `false` otherwise.
"""
has_triangles(queue::RefinementQueue) = !isempty(queue.triangles)

#=
    isempty(queue::RefinementQueue) -> Bool

Return `true` if `queue` has no segments or triangles, `false` otherwise.
=#
Base.isempty(queue::RefinementQueue) = !has_segments(queue) && !has_triangles(queue)
