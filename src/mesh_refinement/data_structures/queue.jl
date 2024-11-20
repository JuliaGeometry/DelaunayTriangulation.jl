struct Queue{T}
    data::Vector{T}
end
Queue{T}() where {T} = Queue{T}(T[])
Base.:(==)(q1::Queue, q2::Queue) = q1.data == q2.data

@inline Base.copy(queue::Queue) = Queue(copy(queue.data))
@inline Base.isempty(queue::Queue) = isempty(queue.data)
@inline Base.length(queue::Queue) = length(queue.data)
@inline Base.eltype(::Type{Queue{T}}) where {T} = T
@inline Base.push!(queue::Queue, item) = (push!(queue.data, item); return queue)
@inline Base.popfirst!(queue::Queue) = popfirst!(queue.data)
@inline enqueue_all!(queue::Queue, data) = append!(queue.data, data)
