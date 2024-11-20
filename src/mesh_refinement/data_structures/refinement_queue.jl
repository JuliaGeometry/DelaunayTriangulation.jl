struct RefinementQueue{I,E,F}
    segments::MaxPriorityQueue{E,F}
    triangles::MaxPriorityQueue{NTuple{3,I},F}
    function RefinementQueue{I,E,F}() where {I,E,F}
        return new{I,E,F}(MaxPriorityQueue{E,F}(), MaxPriorityQueue{NTuple{3,I},F}())
    end
end
function Base.show(io::IO, ::MIME"text/plain", queue::RefinementQueue)
    println(io, "RefinementQueue")
    println(io, "   $(length(queue.segments)) segments: ", queue.segments)
    print(io, "   $(length(queue.triangles)) triangles: ", queue.triangles)
end

@inline function RefinementQueue(tri::Triangulation)
    I = integer_type(tri)
    E = edge_type(tri)
    F = number_type(tri)
    return RefinementQueue{I,E,F}()
end

@inline function Base.haskey(queue::RefinementQueue{T,E,F}, segment::E) where {T,E,F}
    return haskey(queue.segments, segment) || haskey(queue.segments, reverse_edge(segment))
end

@inline function Base.haskey(queue::RefinementQueue{T,E,F}, triangle::T) where {T,E,F}
    return haskey(queue.triangles, triangle) ||
           haskey(queue.triangles, rotate_triangle(triangle, Val(1))) ||
           haskey(queue.triangles, rotate_triangle(triangle, Val(2)))
end

@inline function Base.getindex(queue::RefinementQueue{T,E,F}, triangle::T) where {T,E,F}
    if haskey(queue.triangles, triangle)
        return queue.triangles[triangle]
    elseif haskey(queue.triangles, rotate_triangle(triangle, Val(1)))
        return queue.triangles[rotate_triangle(triangle, Val(1))]
    else #if haskey(queue.triangles, rotate_triangle(triangle, Val(2)))
        return queue.triangles[rotate_triangle(triangle, Val(2))]
    end
end

function Base.setindex!(queue::RefinementQueue{T,E,F}, ℓ²::F, segment::E) where {T,E,F}
    segments = queue.segments
    if haskey(segments, reverse_edge(segment))
        segments[reverse_edge(segment)] = ℓ²
    else
        segments[segment] = ℓ²
    end
    return segments
end

function Base.setindex!(queue::RefinementQueue{T,E,F}, ρ, triangle::T) where {T,E,F}
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

@inline popfirst_segment!(queue::RefinementQueue) = popfirst!(queue.segments)
@inline popfirst_triangle!(queue::RefinementQueue) = popfirst!(queue.triangles)
@inline has_segments(queue::RefinementQueue) = !isempty(queue.segments)
@inline has_triangles(queue::RefinementQueue) = !isempty(queue.triangles)
@inline Base.isempty(queue::RefinementQueue) = !has_segments(queue) && !has_triangles(queue)
