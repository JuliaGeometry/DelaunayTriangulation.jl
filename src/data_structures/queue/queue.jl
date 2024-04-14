"""
    Queue{T}

Struct for a first-in first-out queue. 

!!! note "Performance"

    Under the hood, `Queue` simply uses a `Vector`. This may not be as optimised compared to other implementations,
    e.g. DataStructure.jl's block-based approach with a `Dequeue`.
""" 
struct Queue{T}
    data::Vector{T}
end
Queue{T}() where T = Queue{T}(T[])
Base.:(==)(q1::Queue, q2::Queue) = q1.data == q2.data

"""
    isempty(queue::Queue) -> Bool 

Returns `true` if the `queue` is empty, `false` otherwise.
"""
Base.isempty(queue::Queue) = isempty(queue.data)

"""
    length(queue::Queue) -> Int

Returns the number of elements in the `queue`.
"""
Base.length(queue::Queue) = length(queue.data)

"""
    eltype(queue::Queue{T}) -> Type{T}

Returns the type of elements stored in `q`.
"""
Base.eltype(queue::Type{Queue{T}}) where {T} = T

"""
    push!(queue::Queue, item) 

Adds `item` to the end of the `queue`.
"""
Base.push!(queue::Queue, item) = (push!(queue.data, item); return queue)

"""
    popfirst!(queue::Queue) 

Removes the element from the front of the `queue` and returns it.
"""
Base.popfirst!(queue::Queue) = popfirst!(queue.data)

"""
    enqueue_all!(queue::Queue, data) 

Adds all `data` to the end of the `queue`.
"""
enqueue_all!(queue::Queue, data) = append!(queue.data, data)

